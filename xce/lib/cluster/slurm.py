# -*- coding: utf-8 -*-
import os
import ssl
import json
import httplib
import paramiko
import time
import traceback
import gtk
import getpass
from PyQt4 import QtGui
from datetime import datetime
from xce.lib.XChemLog import updateLog
from uuid import uuid4

CLUSTER_USER = (os.environ.get("CLUSTER_USER")
                or os.environ.get("USER")
                or getpass.getuser())

CLUSTER_BASTION = os.environ.get("CLUSTER_BASTION", "")
CLUSTER_HOST = os.environ.get("CLUSTER_HOST", "")
CLUSTER_PORT = int(os.environ.get("CLUSTER_PORT", "6820"))
CLUSTER_PARTITION = os.environ.get("CLUSTER_PARTITION", "regular")

TOKEN = None
TOKEN_EXPIRY = None

POPUP_TITLE = "SLURM Authentication"


def fetch_password_qt(password_prompt):
    password, ok = QtGui.QInputDialog.getText(
        None, POPUP_TITLE, password_prompt, mode=QtGui.QLineEdit.Password
    )
    return password if ok else None


def fetch_password_gtk(password_prompt):
    dialog = gtk.MessageDialog(
        None,
        gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,
        gtk.MESSAGE_QUESTION,
        gtk.BUTTONS_OK_CANCEL,
        None,
    )
    dialog.set_title(POPUP_TITLE)
    dialog.set_markup(password_prompt)

    entry = gtk.Entry()
    entry.set_visibility(False)
    dialog.vbox.pack_end(entry)
    dialog.show_all()

    dialog.run()
    password = entry.get_text()
    dialog.destroy()
    return password


def get_token(fetch_password, error=None):
    global TOKEN
    global TOKEN_EXPIRY

    if TOKEN is None or TOKEN_EXPIRY is None or TOKEN_EXPIRY < time.time() + 60:
        password_prompt = error + "\n" + "Password:" if error else "Password:"
        password = fetch_password(password_prompt)
        if password is None:
            return None

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.load_system_host_keys()
        try:
            ssh.connect(CLUSTER_BASTION, username=CLUSTER_USER, password=str(password))
        except paramiko.AuthenticationException:
            print(traceback.format_exc())
            return get_token(fetch_password, error="SSH Authentication Failed")
        stdin, stdout, stderr = ssh.exec_command("scontrol token lifespan=3600")
        if stdout.channel.recv_exit_status() != 0:
            return get_token(fetch_password, error="Token Acquisition Failed")
        final_line = stdout.next()
        for final_line in stdout:
            continue
        TOKEN = final_line.split("=")[1].strip()
        TOKEN_EXPIRY = time.time() + 3600
    return TOKEN


def construct_headers(token):
    return {
        "Content-Type": "application/json",
        "X-SLURM-USER-NAME": CLUSTER_USER,
        "X-SLURM-USER-TOKEN": token,
    }


def submit_cluster_job(
    name, file, xce_logfile, token, array=None, exclusive=False, memory=None, tasks=None
):
    with open(file) as script_file:
        script = "\n".join(script_file.readlines())
    payload = dict(
        script=script,
        job=dict(
            partition=CLUSTER_PARTITION,
            name=str(name),
            account=CLUSTER_USER,
            current_working_directory=os.getcwd(),
            environment=[
                "USER={}".format(os.environ["USER"]),
                "XChemExplorer_DIR={}".format(os.environ["XChemExplorer_DIR"]),
            ],
            standard_output=os.path.join(os.getcwd(), "{}.stdout".format(name)),
            standard_error=os.path.join(os.getcwd(), "{}.stderr".format(name)),
        ),
    )
    if array is not None:
        payload["job"]["array"] = array
    if exclusive is True:
        payload["job"]["exclusive"] = "mcs"
        payload["job"]["mcs_label"] = str(uuid4())
    if memory is not None:
        payload["job"]["memory_per_node"] = dict()
        payload["job"]["memory_per_node"]["set"] = True
        payload["job"]["memory_per_node"]["number"] = memory
    if tasks is not None:
        payload["job"]["tasks_per_node"] = tasks
    body = json.dumps(payload)
    logfile = updateLog(xce_logfile)
    logfile.insert("Submitting job, '{}', to Slurm with body: {}".format(name, body))

    if not CLUSTER_HOST:
        # No REST API configured - use plain sbatch (e.g. WEHI Milton).
        # subprocess pipes fail inside Apptainer on restricted kernels (ENOTCONN),
        # so we use os.system() which forks without creating pipes.
        sbatch_bin = "sbatch"
        for _candidate in [
            "/usr/local/slurm/bin/sbatch",
            "/usr/local/bin/sbatch",
            "/usr/bin/sbatch",
            "/opt/slurm/bin/sbatch",
        ]:
            if os.path.isfile(_candidate) and os.access(_candidate, os.X_OK):
                sbatch_bin = _candidate
                break

        sbatch_cmd = (
            "{sbatch} --partition {partition} --job-name {name}"
            " --output {out} --error {err}"
        ).format(
            sbatch=sbatch_bin,
            partition=CLUSTER_PARTITION,
            name=str(name),
            out=os.path.join(os.getcwd(), "{}.stdout".format(name)),
            err=os.path.join(os.getcwd(), "{}.stderr".format(name)),
        )
        if array is not None:
            sbatch_cmd += " --array {!s}".format(array)
        if exclusive:
            sbatch_cmd += " --exclusive"
        if memory is not None:
            sbatch_cmd += " --mem {!s}".format(memory)
        if tasks is not None:
            sbatch_cmd += " --ntasks-per-node {!s}".format(tasks)
        sbatch_cmd += " " + file
        logfile.insert("No CLUSTER_HOST set; submitting via sbatch: " + sbatch_cmd)
        ret = os.system(sbatch_cmd)
        if ret != 0:
            logfile.insert("sbatch exited with code {!s}".format(ret))
            raise OSError("sbatch failed with exit code {!s}".format(ret))
        return

    connection = httplib.HTTPSConnection(
        CLUSTER_HOST, CLUSTER_PORT, context=ssl._create_unverified_context()
    )
    connection.request(
        "POST", "/slurm/v0.0.40/job/submit", body=body, headers=construct_headers(token)
    )
    response = connection.getresponse().read()
    logfile.insert("Got response: {}".format(response))


def _parse_squeue_time(time_str):
    """Parse squeue elapsed time string (D-HH:MM:SS, HH:MM:SS or MM:SS) to timedelta."""
    from datetime import timedelta
    try:
        days = 0
        if "-" in time_str:
            d, time_str = time_str.split("-", 1)
            days = int(d)
        parts = time_str.split(":")
        if len(parts) == 3:
            h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
        elif len(parts) == 2:
            h, m, s = 0, int(parts[0]), int(parts[1])
        else:
            h, m, s = 0, 0, int(parts[0])
        return timedelta(days=days, hours=h, minutes=m, seconds=s)
    except Exception:
        from datetime import timedelta
        return timedelta(0)


def query_running_jobs(xce_logfile, token):
    """Query running jobs via squeue (local sbatch, no REST API needed at WEHI)."""
    import tempfile
    logfile = updateLog(xce_logfile)

    # Locate squeue on the host path (bind-mounted into the container)
    squeue_bin = "squeue"
    for candidate in [
        "/usr/local/slurm/bin/squeue",
        "/usr/local/bin/squeue",
        "/usr/bin/squeue",
        "/opt/slurm/bin/squeue",
    ]:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            squeue_bin = candidate
            break

    # Write to a temp file - avoids subprocess pipe issues inside Apptainer
    tmp = tempfile.mktemp(suffix="_xce_squeue.txt")
    ret = os.system(
        "{squeue} --user={user} --noheader --format='%i %j %T %M' > {out} 2>/dev/null".format(
            squeue=squeue_bin, user=CLUSTER_USER, out=tmp
        )
    )

    jobs = []
    if ret != 0:
        logfile.insert("squeue exited with code {!s}".format(ret))
        return jobs

    try:
        with open(tmp) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(None, 3)
                if len(parts) < 4:
                    continue
                job_id, job_name, job_status, elapsed = parts
                run_time = _parse_squeue_time(elapsed)
                jobs.append((job_id, job_name, job_status, run_time))
    except Exception as exc:
        logfile.insert("Error reading squeue output: {!s}".format(exc))
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass

    return jobs

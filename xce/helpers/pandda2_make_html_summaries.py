"""
pandda2_make_html_summaries.py  --  generate XCE-compatible HTML summary files
from pandda2 output CSVs.

Reads:
  <pandda_dir>/analyses/pandda_analyse_events.csv
  <pandda_dir>/analyses/pandda_analyse_sites.csv  (optional)

Writes:
  <pandda_dir>/analyses/html_summaries/pandda_initial.html   (Dataset Summary tab)
  <pandda_dir>/analyses/html_summaries/pandda_analyse.html   (Processing Output tab)

Compatible with Python 2.7 (ccp4-python) and Python 3.
Usage:
  ccp4-python pandda2_make_html_summaries.py <pandda_dir>
"""
from __future__ import print_function
import csv
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

_STYLE = """
<style>
  body { font-family: Arial, sans-serif; font-size: 13px; margin: 20px; }
  h2   { color: #2c5f8a; }
  p    { color: #444; }
  table { border-collapse: collapse; width: 100%; }
  th { background: #2c5f8a; color: white; padding: 6px 10px; text-align: left; }
  td { padding: 5px 10px; border-bottom: 1px solid #ddd; }
  tr:nth-child(even) { background: #f4f8fc; }
  tr:hover { background: #ddeeff; }
  .hit { background: #d4edda !important; }
</style>
"""


def _html_page(title, body):
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        "<title>{title}</title>{style}</head><body>"
        "<h2>{title}</h2>{body}</body></html>"
    ).format(title=title, style=_STYLE, body=body)


def _table(headers, rows, row_class_fn=None):
    head_cells = "".join("<th>{0}</th>".format(h) for h in headers)
    head = "<thead><tr>{0}</tr></thead>".format(head_cells)
    body_rows = []
    for row in rows:
        cls = row_class_fn(row) if row_class_fn else ""
        cells = "".join("<td>{0}</td>".format(c) for c in row)
        body_rows.append("<tr{0}>{1}</tr>".format(
            " class='{0}'".format(cls) if cls else "", cells))
    return "<table>{0}<tbody>{1}</tbody></table>".format(head, "".join(body_rows))


# ---------------------------------------------------------------------------
# Read CSVs
# ---------------------------------------------------------------------------

def _read_csv(path):
    """Return list of dicts; empty list if file missing."""
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Build pandda_initial.html  (Dataset Summary tab)
# ---------------------------------------------------------------------------

def make_initial_html(events, out_path):
    """One row per dataset: resolution, R-free, R-work, #events, max z_peak."""
    datasets = defaultdict(lambda: {
        "resolution": "", "r_free": "", "r_work": "",
        "n_events": 0, "max_z_peak": 0.0,
    })

    for row in events:
        dtag = row.get("dtag", "")
        d = datasets[dtag]
        if not d["resolution"]:
            d["resolution"] = row.get("high_resolution", "")
            d["r_free"]     = row.get("r_free", "")
            d["r_work"]     = row.get("r_work", "")
        d["n_events"] += 1
        try:
            zp = float(row.get("z_peak", 0) or 0)
            if zp > d["max_z_peak"]:
                d["max_z_peak"] = zp
        except (ValueError, TypeError):
            pass

    headers = ["Dataset", "Resolution (A)", "R-free", "R-work",
               "Events found", "Max z-peak"]

    rows = []
    for dtag in sorted(datasets.keys()):
        d = datasets[dtag]
        rows.append([
            dtag,
            d["resolution"],
            d["r_free"],
            d["r_work"],
            str(d["n_events"]),
            "{0:.2f}".format(d["max_z_peak"]) if d["max_z_peak"] else "0.00",
        ])

    def row_class(row):
        try:
            return "hit" if int(row[4]) > 0 else ""
        except (ValueError, IndexError):
            return ""

    n_with_events = sum(1 for d in datasets.values() if d["n_events"] > 0)
    summary = "<p>{0} datasets processed; {1} with events.</p>".format(
        len(datasets), n_with_events)
    body = summary + _table(headers, rows, row_class)

    _ensure_dir(os.path.dirname(out_path))
    with open(out_path, "w") as fh:
        fh.write(_html_page("PanDDA 2 - Dataset Summary", body))
    print("Written: {0}".format(out_path))


# ---------------------------------------------------------------------------
# Build pandda_analyse.html  (Processing Output tab)
# ---------------------------------------------------------------------------

def make_analyse_html(events, sites, out_path):
    """One row per event sorted by z_peak descending."""
    def sort_key(row):
        try:
            return -float(row.get("z_peak", 0) or 0)
        except (ValueError, TypeError):
            return 0

    sorted_events = sorted(events, key=sort_key)

    headers = ["Dataset", "Event", "Site", "z-peak", "z-mean",
               "Cluster size", "BDC", "1-BDC", "x", "y", "z"]

    rows = []
    for row in sorted_events:
        bdc_val = row.get("bdc", "")
        try:
            one_bdc = "{0:.2f}".format(1.0 - float(bdc_val))
        except (ValueError, TypeError):
            one_bdc = row.get("1-BDC", "")
        rows.append([
            row.get("dtag", ""),
            row.get("event_idx", ""),
            row.get("site_idx", ""),
            row.get("z_peak", ""),
            row.get("z_mean", ""),
            row.get("cluster_size", ""),
            bdc_val,
            one_bdc,
            row.get("x", ""),
            row.get("y", ""),
            row.get("z", ""),
        ])

    def row_class(row):
        try:
            return "hit" if float(row[3]) >= 3.0 else ""
        except (ValueError, IndexError):
            return ""

    n_datasets = len(set(r.get("dtag", "") for r in events))
    summary = "<p>{0} events across {1} datasets.</p>".format(
        len(events), n_datasets)
    body = summary + _table(headers, rows, row_class)

    _ensure_dir(os.path.dirname(out_path))
    with open(out_path, "w") as fh:
        fh.write(_html_page("PanDDA 2 - Processing Output", body))
    print("Written: {0}".format(out_path))


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _ensure_dir(path):
    if path and not os.path.isdir(path):
        os.makedirs(path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(pandda_dir):
    analyses   = os.path.join(pandda_dir, "analyses")
    html_dir   = os.path.join(analyses, "html_summaries")
    events_csv = os.path.join(analyses, "pandda_analyse_events.csv")
    sites_csv  = os.path.join(analyses, "pandda_analyse_sites.csv")

    if not os.path.isfile(events_csv):
        print("Events CSV not found: {0}".format(events_csv))
        sys.exit(1)

    events = _read_csv(events_csv)
    sites  = _read_csv(sites_csv)

    make_initial_html(
        events,
        os.path.join(html_dir, "pandda_initial.html"),
    )
    make_analyse_html(
        events,
        sites,
        os.path.join(html_dir, "pandda_analyse.html"),
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: pandda2_make_html_summaries.py <pandda_dir>")
        sys.exit(1)
    main(sys.argv[1])

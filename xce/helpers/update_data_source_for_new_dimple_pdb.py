import os
import sys

from iotbx import mtz

if __name__ == "__main__":
    sys.path.insert(0, os.environ["XChemExplorer_DIR"])
    from xce.lib import XChemDB
    from xce.lib.XChemUtils import parse

    db_file = sys.argv[1]
    xtal = sys.argv[2]
    inital_model_directory = sys.argv[3]

    db = XChemDB.data_source(db_file)
    if os.path.isfile(os.path.join(inital_model_directory, xtal, "dimple.pdb")):
        db_dict = {
            "DimplePathToPDB": os.path.join(inital_model_directory, xtal, "dimple.pdb")
        }
        dimple_ran_successfully = False
        if os.path.isfile(os.path.join(inital_model_directory, xtal, "dimple.mtz")):
            db_dict["DimplePathToMTZ"] = os.path.join(
                inital_model_directory, xtal, "dimple.mtz"
            )
            dimple_ran_successfully = True
            db_dict["DataProcessingDimpleSuccessful"] = "True"
            db_dict["DimpleStatus"] = "finished"
        if not dimple_ran_successfully:
            db_dict["DataProcessingDimpleSuccessful"] = "False"
            db_dict["DimpleStatus"] = "failed"
        pdb = parse().PDBheader(
            os.path.join(inital_model_directory, xtal, "dimple.pdb")
        )
        db_dict["DimpleRcryst"] = pdb["Rcryst"]
        db_dict["DimpleRfree"] = pdb["Rfree"]
        db_dict["RefinementOutcome"] = "1 - Analysis Pending"
        db_dict["RefinementSpaceGroup"] = pdb["SpaceGroup"]

        # setting free.mtz file
        os.chdir(os.path.join(inital_model_directory, xtal))
        os.system("/bin/rm -f %s.free.mtz" % xtal)
        mtzFree = None
        db_dict["RefinementMTZfree"] = ""
        if os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple",
                "prepared2.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple",
                "prepared2.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple",
                "prepared.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple",
                "prepared.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory, xtal, "dimple", "dimple", "prepared.mtz"
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory, xtal, "dimple", "dimple", "prepared.mtz"
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory, xtal, "dimple", "dimple", "prepared2.mtz"
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory, xtal, "dimple", "dimple", "prepared2.mtz"
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple",
                "free.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple",
                "free.mtz",
            )
        elif os.path.isfile(
            os.path.join(inital_model_directory, xtal, "dimple", "dimple", "free.mtz")
        ):
            mtzFree = os.path.join(
                inital_model_directory, xtal, "dimple", "dimple", "free.mtz"
            )

        if mtzFree is not None:
            if "F_unique" in mtz.object(mtzFree).column_labels():
                cmd = (
                    "cad hklin1 %s hklout %s.free.mtz << eof\n" % (mtzFree, xtal)
                    + " monitor BRIEF\n"
                    " labin file 1 E1=F E2=SIGF E3=FreeR_flag\n"
                    " labout file 1 E1=F E2=SIGF E3=FreeR_flag\n"
                    "eof\n"
                )

                os.system(cmd)
            else:
                os.symlink(mtzFree, xtal + ".free.mtz")

            db_dict["RefinementMTZfree"] = xtal + ".free.mtz"

        print("==> xce: updating data source after DIMPLE run")
        db.update_data_source(xtal, db_dict)

    else:
        # the actual dimple script creates symbolic links regardless if dimple was
        # successful or not python os.path.isfile is False if symbolic link points to
        # non existing file so we remove all of them!
        os.chdir(os.path.join(inital_model_directory, xtal))
        os.system("/bin/rm dimple.pdb")
        os.system("/bin/rm dimple.mtz")
        os.system("/bin/rm 2fofc.map")
        os.system("/bin/rm fofc.map")

        # For phenix.ligand_pipeline: refine.pdb is the pipeline output.
        # Populate R-factors, space group, and point group so Maps tab columns
        # (Refinement Rcryst/Rfree, Dimple Rcryst/Rfree, DataProcessing SpaceGroup)
        # are all kept up to date with the Phaser-solved space group.
        refine_pdb = os.path.join(inital_model_directory, xtal, "refine.pdb")
        if os.path.isfile(refine_pdb):
            pdb = parse().PDBheader(refine_pdb)
            ref_db_dict = {
                # Refinement R-factors — shown in the new Maps tab columns
                "RefinementRcryst": pdb["Rcryst"],
                "RefinementRcrystTraficLight": pdb["RcrystTL"],
                "RefinementRfree": pdb["Rfree"],
                "RefinementRfreeTraficLight": pdb["RfreeTL"],
                "RefinementResolution": pdb["ResolutionHigh"],
                "RefinementSpaceGroup": pdb["SpaceGroup"],
                "RefinementStatus": "finished",
                "RefinementOutcome": "1 - Analysis Pending",
                # Dimple R-factor columns — already shown in Maps tab; fill from
                # refine.pdb so they are not left blank for pipeline runs
                "DimpleRcryst": pdb["Rcryst"],
                "DimpleRfree": pdb["Rfree"],
                "DimpleStatus": "finished",
                "DimplePathToPDB": os.path.realpath(refine_pdb),
                # Keep DataProcessing space group / point group in sync with the
                # Phaser-solved space group stored in the CRYST1 record of refine.pdb
                "DataProcessingSpaceGroup": pdb["SpaceGroup"],
                "DataProcessingPointGroup": pdb["PointGroup"],
            }
            print("==> xce: updating data source with phenix.ligand_pipeline results from refine.pdb")
            db.update_data_source(xtal, ref_db_dict)

# -*- coding: utf-8 -*-
# plugin for COOT to inspect and model pandda.analyse results
# Python 3 / GTK3 / Coot 0.9.x (CCP4 7.1) compatible version
#
# Based on the original by Tobias Krojer, MAX IV Laboratory
# MIT License - see inspect_pandda_analyse.py for full licence text

import os
import glob
import sys
import shutil
import csv
import logging

# ---------------------------------------------------------------------------
# GTK2 / GTK3 import -- CCP4 7.1 ships a GTK2-based Coot environment.
# gi.require_version raises ValueError (not ImportError) when GTK3 is absent,
# so we must catch both.
# ---------------------------------------------------------------------------
try:
    import gi
    gi.require_version('Gtk', '3.0')
    from gi.repository import Gtk
    _GTK2 = False
except (ImportError, ValueError):
    import gtk as Gtk
    _GTK2 = True


# ---------------------------------------------------------------------------
# GTK2 / GTK3 widget factory helpers
# ---------------------------------------------------------------------------

def _VBox(spacing=0):
    if _GTK2:
        return Gtk.VBox(False, spacing)
    return Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=spacing)


def _HBox(spacing=0):
    if _GTK2:
        return Gtk.HBox(False, spacing)
    return Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=spacing)


def _Label(text=''):
    # GTK2 gtk.Label uses 'str' as param name, not 'label'
    if _GTK2:
        return Gtk.Label(str(text))
    return Gtk.Label(label=text)


def _ComboBoxText():
    if _GTK2:
        return Gtk.combo_box_new_text()
    return Gtk.ComboBoxText()


def _make_info_grid(nrows, ncols):
    if _GTK2:
        t = Gtk.Table(nrows, ncols, False)
        t.set_row_spacings(2)
        t.set_col_spacings(2)
        return t
    g = Gtk.Grid()
    g.set_row_homogeneous(False)
    g.set_column_homogeneous(True)
    g.set_row_spacing(2)
    g.set_column_spacing(2)
    return g


def _grid_attach(container, widget, col, row):
    if _GTK2:
        container.attach(widget, col, col + 1, row, row + 1)
    else:
        container.attach(widget, col, row, 1, 1)


def _RadioButton(label, group_widget=None):
    if _GTK2:
        # GTK2: first arg is a RadioButton instance (or None), not a group list
        return Gtk.RadioButton(group_widget, label)
    if group_widget is None:
        return Gtk.RadioButton(label=label)
    return Gtk.RadioButton.new_with_label_from_widget(group_widget, label)


def _CheckButton(label=''):
    if _GTK2:
        return Gtk.CheckButton(label)
    return Gtk.CheckButton(label=label)


if _GTK2:
    _FOLDER_ACTION = Gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER
    _RESP_CANCEL   = Gtk.RESPONSE_CANCEL
    _RESP_OK       = Gtk.RESPONSE_OK
else:
    _FOLDER_ACTION = Gtk.FileChooserAction.SELECT_FOLDER
    _RESP_CANCEL   = Gtk.ResponseType.CANCEL
    _RESP_OK       = Gtk.ResponseType.OK


import coot
import sys as _sys

# Python 2/3 compatibility helpers for csv file open
if _sys.version_info[0] >= 3:
    def _csv_open_r(path):
        return open(path, newline='')
    def _csv_open_w(path):
        return open(path, 'w', newline='')
else:
    def _csv_open_r(path):
        return open(path, 'rb')
    def _csv_open_w(path):
        return open(path, 'wb')


# ---------------------------------------------------------------------------
# Coot API compatibility wrappers
# ---------------------------------------------------------------------------
# In Coot 0.9.x (Python 3) several helpers that lived in __main__ moved into
# the coot module itself.  These thin wrappers insulate the rest of the code.

def _molecule_number_list():
    try:
        return coot.molecule_number_list()
    except AttributeError:
        import __main__
        return __main__.molecule_number_list()


def _set_map_displayed(imol, show):
    """Show (1) or hide (0) a map molecule."""
    try:
        coot.set_map_displayed(imol, show)
    except AttributeError:
        import __main__
        __main__.toggle_display_map(imol, show)


def _move_molecule_here(imol):
    """Move a molecule to the current screen centre / pointer."""
    try:
        coot.move_molecule_to_screen_centre_py(imol)
    except AttributeError:
        try:
            import __main__
            __main__.move_molecule_here(imol)
        except AttributeError:
            pass


def _merge_molecules(imol_list, target_imol):
    try:
        coot.merge_molecules(imol_list, target_imol)     # Coot 0.9.x
    except AttributeError:
        coot.merge_molecules_py(imol_list, target_imol)  # Coot 0.8.x


# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

def init_logger(logfile):
    logger = logging.getLogger('pandda_inspect')
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter(
        '%(asctime)s | %(levelname)s - INSPECT | %(message)s',
        '%m-%d-%Y %H:%M:%S',
    )
    if not logger.handlers:
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(fmt)
        logger.addHandler(sh)
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


# ---------------------------------------------------------------------------
# Main inspect GUI class
# ---------------------------------------------------------------------------

class inspect_gui(object):

    def __init__(self):
        self.logger = init_logger('inspect.log')
        self.logger.info('starting new session of pandda event map inspection')

        self.index = -1
        self.Todo = []
        self.cb_list = []
        self.mol_dict = {'protein': None, 'emap': None, 'ligand': None}

        self.panddaDir = None
        self.eventCSV  = None
        self.reset_params()
        self.merged = False

        # CSV column indices; set in parsepanddaDir
        self.comment_index            = None
        self.interesting_index        = None
        self.ligand_placed_index      = None
        self.ligand_confidence_index  = None
        self.slist_site_num_index     = None
        self.slist_site_name_index    = None
        self.slist_site_comment_index = None

        # Suppress widget callbacks during programmatic updates
        self._updating_widgets = False

        self.selection_criteria = [
            'show all events',
            'show all events - sort by cluster size',
            'show all events - sort alphabetically',
            'show not viewed events',
            'show unassigned',
            'show no ligands bound',
            'show low confidence ligands',
            'show medium confidence ligands',
            'show high confidence ligands',
        ]
        self.selected_selection_criterion = None

    # ------------------------------------------------------------------
    # GUI construction
    # ------------------------------------------------------------------

    def startGUI(self):
        self.window = Gtk.Window()
        self.window.connect("delete-event", lambda w, e: Gtk.main_quit())
        self.window.set_border_width(6)
        self.window.set_default_size(750, 900)
        self.window.set_title("PANDDA inspect")

        self.vbox = _VBox(4)

        # ---- Row 1: Quit / Summary / Update HTML / Progress / Go to Dataset ----
        row1 = _HBox(6)
        for lbl, cb in [("Quit",        lambda w: Gtk.main_quit()),
                        ("Summary",     self.show_summary),
                        ("Update HTML", self.update_html)]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            row1.pack_start(b, False, False, 0)
        prog_vbox = _VBox(1)
        prog_vbox.pack_start(_Label("Overall Inspection Event/Site Progress:"),
                             False, False, 0)
        self.crystal_progressbar = Gtk.ProgressBar()
        prog_vbox.pack_start(self.crystal_progressbar, False, False, 0)
        row1.pack_start(prog_vbox, True, True, 4)
        goto_hbox = _HBox(2)
        goto_hbox.pack_start(_Label("Go to Dataset:"), False, False, 0)
        self.goto_entry = Gtk.Entry()
        self.goto_entry.set_width_chars(18)
        self.goto_entry.connect("activate", self.go_to_dataset)
        goto_hbox.pack_start(self.goto_entry, False, False, 0)
        go_btn = Gtk.Button(label="Go")
        go_btn.connect("clicked", self.go_to_dataset)
        goto_hbox.pack_start(go_btn, False, False, 0)
        row1.pack_start(goto_hbox, False, False, 0)
        self.vbox.pack_start(row1, False, False, 0)

        # ---- Row 2: Event counter + site navigation ----
        row2 = _HBox(6)
        self.event_counter_label = _Label("Event  -  of  -")
        row2.pack_start(self.event_counter_label, False, False, 4)
        for lbl, cb in [("<<< Go to Prev Site <<<", self.previous_site),
                        (">>> Go to Next Site >>>", self.next_site)]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            row2.pack_start(b, True, True, 0)
        self.vbox.pack_start(row2, False, False, 0)

        # ---- Row 3: Site counter + unviewed / modelled ----
        row3 = _HBox(6)
        self.site_counter_label = _Label("Site  -  of  -")
        row3.pack_start(self.site_counter_label, False, False, 4)
        for lbl, cb in [(">>> Go to Next Unviewed >>>", self.next_unviewed_event),
                        (">>> Go to Next Modelled >>>", self.next_modelled_event)]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            row3.pack_start(b, True, True, 0)
        self.vbox.pack_start(row3, False, False, 0)

        # ---- Row 4: Prev / Next navigation ----
        row4 = _HBox(4)
        for lbl, cb in [
            ("<<< Prev <<<\n(Don't Save Model)", self.previous_event),
            (">>> Next >>>\n(Don't Save Model)", self.next_event_no_save),
            (">>> Next >>>\n(Save Model)",        self.next_event_save),
        ]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            row4.pack_start(b, True, True, 0)
        self.vbox.pack_start(row4, False, False, 0)

        # ---- Info grid + Modelling buttons ----
        info_outer = _HBox(6)
        _INFO_ROWS = [
            ('Dataset ID',      'dataset_id'),
            ('Event #',         'event_num'),
            ('Site #',          'site_num'),
            ('1 - BDC',         'bdc'),
            ('Resolution',      'resolution'),
            ('R-Free / R-Work', 'rfree_rwork'),
        ]
        info_grid = _make_info_grid(len(_INFO_ROWS), 2)
        self.info_labels = {}
        for row, (name, key) in enumerate(_INFO_ROWS):
            lf = Gtk.Frame()
            lf.add(_Label(name))
            _grid_attach(info_grid, lf, 0, row)
            val = _Label('')
            vf = Gtk.Frame()
            vf.add(val)
            _grid_attach(info_grid, vf, 1, row)
            self.info_labels[key] = val
        self.xtal_label       = self.info_labels['dataset_id']
        self.event_label      = self.info_labels['event_num']
        self.site_label       = self.info_labels['site_num']
        self.bdc_label        = self.info_labels['bdc']
        self.resolution_label = self.info_labels['resolution']
        self.r_free_label     = self.info_labels['rfree_rwork']
        self.r_work_label     = self.info_labels['rfree_rwork']
        info_outer.pack_start(info_grid, True, True, 0)
        mod_vbox = _VBox(2)
        for lbl, cb in [
            ("Merge Ligand\nWith Model",    self.merge_ligand_into_protein),
            ("Save Model",                  self.save_model),
            ("Move New\nLigand Here",       self.place_ligand_here),
            ("Reload Last\nSaved Model",    self.reload_last_saved_model),
            ("Open Next Ligand",            self.open_next_ligand),
            ("Reset to\nUnfitted Model",    self.reset_to_unfitted),
        ]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            mod_vbox.pack_start(b, False, False, 0)
        info_outer.pack_start(mod_vbox, False, False, 0)
        self.vbox.pack_start(info_outer, False, False, 0)

        # ---- Record Event Information ----
        ev_frame = Gtk.Frame(
            label='Record Event Information (this event only)')
        ev_vbox = _VBox(4)
        comment_hbox = _HBox(4)
        comment_hbox.pack_start(_Label("Event Comment:"), False, False, 0)
        self.comment_entry = Gtk.Entry()
        self.comment_entry.connect("activate", self.save_comment)
        comment_hbox.pack_start(self.comment_entry, True, True, 0)
        ev_vbox.pack_start(comment_hbox, False, False, 0)
        ann_hbox = _HBox(8)
        # Interesting
        int_vbox = _VBox(2)
        self.interesting_radio_yes = _RadioButton("Mark Event as Interesting")
        self.interesting_radio_no  = _RadioButton(
            "Mark Event as Not Interesting", self.interesting_radio_yes)
        self.interesting_radio_no.set_active(True)
        self.interesting_radio_yes.connect("toggled", self.set_interesting_radio)
        self.interesting_radio_no.connect("toggled",  self.set_interesting_radio)
        int_vbox.pack_start(self.interesting_radio_yes, False, False, 0)
        int_vbox.pack_start(self.interesting_radio_no,  False, False, 0)
        ann_hbox.pack_start(int_vbox, True, True, 0)
        # Ligand placed
        lig_vbox = _VBox(2)
        self.ligand_placed_radio_yes = _RadioButton("Ligand Placed")
        self.ligand_placed_radio_no  = _RadioButton(
            "No Ligand Placed", self.ligand_placed_radio_yes)
        self.ligand_placed_radio_no.set_active(True)
        self.ligand_placed_radio_yes.connect("toggled", self.set_ligand_placed_radio)
        self.ligand_placed_radio_no.connect("toggled",  self.set_ligand_placed_radio)
        lig_vbox.pack_start(self.ligand_placed_radio_yes, False, False, 0)
        lig_vbox.pack_start(self.ligand_placed_radio_no,  False, False, 0)
        ann_hbox.pack_start(lig_vbox, True, True, 0)
        # Confidence
        conf_vbox = _VBox(2)
        self.confidence_radio_high   = _RadioButton("Model: High Confidence")
        self.confidence_radio_medium = _RadioButton(
            "Model: Medium Confidence", self.confidence_radio_high)
        self.confidence_radio_low    = _RadioButton(
            "Model: Low Confidence",    self.confidence_radio_high)
        self.confidence_radio_low.set_active(True)
        for r in (self.confidence_radio_high,
                  self.confidence_radio_medium,
                  self.confidence_radio_low):
            r.connect("toggled", self.set_confidence_radio)
            conf_vbox.pack_start(r, False, False, 0)
        ann_hbox.pack_start(conf_vbox, True, True, 0)
        ev_vbox.pack_start(ann_hbox, False, False, 0)
        ev_frame.add(ev_vbox)
        self.vbox.pack_start(ev_frame, False, False, 0)

        # ---- Record Site Information ----
        site_frame = Gtk.Frame(
            label='Record Site Information (for all events with this site)')
        site_grid = _make_info_grid(2, 2)
        lf = Gtk.Frame()
        lf.add(_Label("Name:"))
        _grid_attach(site_grid, lf, 0, 0)
        self.site_name_entry = Gtk.Entry()
        self.site_name_entry.connect("activate", self.save_site_info)
        _grid_attach(site_grid, self.site_name_entry, 1, 0)
        lf2 = Gtk.Frame()
        lf2.add(_Label("Comment:"))
        _grid_attach(site_grid, lf2, 0, 1)
        self.site_comment_entry = Gtk.Entry()
        self.site_comment_entry.connect("activate", self.save_site_info)
        _grid_attach(site_grid, self.site_comment_entry, 1, 1)
        site_frame.add(site_grid)
        self.vbox.pack_start(site_frame, False, False, 0)

        # ---- Dataset selector + PanDDA folder ----
        sel_frame = Gtk.Frame(label='Dataset selector')
        sel_hbox = _HBox(4)
        self.cb = _ComboBoxText()
        self.cb.connect("changed", self.select_crystal)
        sel_hbox.pack_start(self.cb, True, True, 0)
        folder_btn = Gtk.Button(label="Select pandda directory")
        folder_btn.connect("clicked", self.select_pandda_folder)
        sel_hbox.pack_start(folder_btn, False, False, 0)
        sel_frame.add(sel_hbox)
        self.vbox.pack_start(sel_frame, False, False, 0)

        # ---- Event filter ----
        filter_frame = Gtk.Frame(label='Event filter')
        filter_hbox = _HBox(4)
        self.select_events_combobox = _ComboBoxText()
        for c in self.selection_criteria:
            self.select_events_combobox.append_text(c)
        filter_hbox.pack_start(self.select_events_combobox, True, True, 0)
        filter_btn = Gtk.Button(label="Apply")
        filter_btn.connect("clicked", self.select_events)
        filter_hbox.pack_start(filter_btn, False, False, 0)
        filter_frame.add(filter_hbox)
        self.vbox.pack_start(filter_frame, False, False, 0)

        # ---- Toggle maps ----
        maps_frame = Gtk.Frame(label='Toggle Maps')
        maps_hbox = _HBox(4)
        for lbl, cb in [("event map",    self.toggle_emap),
                        ("Z-map",        self.toggle_zmap),
                        ("(2)fofc maps", self.toggle_x_ray_maps)]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            maps_hbox.pack_start(b, True, True, 0)
        self.toggle_average_map_button = Gtk.Button(label="average map")
        self.toggle_average_map_button.connect("clicked", self.toggle_average_map)
        maps_hbox.pack_start(self.toggle_average_map_button, True, True, 0)
        maps_frame.add(maps_hbox)
        self.vbox.pack_start(maps_frame, False, False, 0)

        # ---- Miscellaneous buttons ----
        misc_frame = Gtk.Frame(label='Miscellaneous buttons')
        misc_hbox = _HBox(4)
        for lbl, cb in [
            ("Load input\nmtz file",                    self.load_input_mtz_clicked),
            ("Load average map",                        self.load_average_map_clicked),
            ("Load unfitted model\n(for comparison only)", self.load_unfitted_model_clicked),
            ("Create new ligand",                       self.create_new_ligand),
        ]:
            b = Gtk.Button(label=lbl)
            b.connect("clicked", cb)
            misc_hbox.pack_start(b, True, True, 0)
        misc_frame.add(misc_hbox)
        self.vbox.pack_start(misc_frame, False, False, 0)

        self.window.add(self.vbox)
        self.window.show_all()

    # ------------------------------------------------------------------
    # CSV persistence
    # ------------------------------------------------------------------

    def save_pandda_inspect_events_csv_file(self):
        path = os.path.join(self.analysis_folder, 'pandda_inspect_events.csv')
        self.logger.info('updating {0!s}'.format(path))
        with _csv_open_w(path) as f:
            csv.writer(f).writerows(self.elist)

    # ------------------------------------------------------------------
    # Annotation callbacks
    # ------------------------------------------------------------------

    def set_confidence_radio(self, widget):
        """Save confidence radio selection to CSV."""
        if not widget.get_active():
            return
        if getattr(self, '_updating_widgets', False) or self.index < 1:
            return
        if self.ligand_confidence_index is None:
            return
        if widget is self.confidence_radio_high:
            val = 'High Confidence'
        elif widget is self.confidence_radio_medium:
            val = 'Medium Confidence'
        else:
            val = 'Low Confidence'
        self.elist[self.index][self.ligand_confidence_index] = val
        self.save_pandda_inspect_events_csv_file()
        self.logger.info('confidence: {0!s}'.format(val))

    def set_interesting_radio(self, widget):
        """Save interesting flag to CSV."""
        if not widget.get_active():
            return
        if getattr(self, '_updating_widgets', False) or self.index < 1:
            return
        if self.interesting_index is None:
            return
        val = 'True' if widget is self.interesting_radio_yes else 'False'
        self.elist[self.index][self.interesting_index] = val
        self.save_pandda_inspect_events_csv_file()
        self.logger.info('interesting: {0!s}'.format(val))

    def set_ligand_placed_radio(self, widget):
        """Save ligand-placed flag to CSV."""
        if not widget.get_active():
            return
        if getattr(self, '_updating_widgets', False) or self.index < 1:
            return
        if self.ligand_placed_index is None:
            return
        val = 'True' if widget is self.ligand_placed_radio_yes else 'False'
        self.elist[self.index][self.ligand_placed_index] = val
        self.save_pandda_inspect_events_csv_file()
        self.logger.info('ligand placed: {0!s}'.format(val))

    def _set_annotation_buttons(self):
        """Update annotation widgets to reflect current event data."""
        self._updating_widgets = True
        lc = self.ligand_confidence or ''
        # Normalise old-style values
        if lc == '2fofc map':
            lc = 'High Confidence'
        elif lc == 'event map only':
            lc = 'Medium Confidence'
        if lc == 'High Confidence':
            self.confidence_radio_high.set_active(True)
        elif lc == 'Medium Confidence':
            self.confidence_radio_medium.set_active(True)
        else:
            self.confidence_radio_low.set_active(True)
        self.interesting_radio_yes.set_active(self.interesting)
        self.interesting_radio_no.set_active(not self.interesting)
        self.ligand_placed_radio_yes.set_active(self.ligand_placed)
        self.ligand_placed_radio_no.set_active(not self.ligand_placed)
        self._updating_widgets = False

    def save_event_as_viewed(self):
        self.elist[self.index][self.viewed_index] = 'True'
        self.save_pandda_inspect_events_csv_file()

    def save_comment(self, widget):
        """Save the comment entry text to the current event's CSV row."""
        if getattr(self, '_updating_widgets', False) or self.index < 1:
            return
        if self.comment_index is None:
            return
        comment = self.comment_entry.get_text()
        self.elist[self.index][self.comment_index] = comment
        self.save_pandda_inspect_events_csv_file()
        self.logger.info("saved comment: '{0!s}'".format(comment))

    # ------------------------------------------------------------------
    # Folder / CSV initialisation
    # ------------------------------------------------------------------

    def select_pandda_folder(self, widget):
        dlg = Gtk.FileChooserDialog(
            title="Select PanDDA directory",
            parent=None,
            action=_FOLDER_ACTION,
        )
        dlg.add_button("_Cancel", _RESP_CANCEL)
        dlg.add_button("_Open",   _RESP_OK)

        response = dlg.run()
        if response != _RESP_OK:
            dlg.destroy()
            return

        self.panddaDir = dlg.get_filename()
        dlg.destroy()

        self.analysis_folder = ''
        for candidate in ('results', 'analyses', 'analysis'):
            p = os.path.join(self.panddaDir, candidate)
            if os.path.isdir(p):
                self.analysis_folder = p
                break

        self.eventCSV = os.path.realpath(
            os.path.join(self.analysis_folder, 'pandda_inspect_events.csv'))
        self.siteCSV = os.path.realpath(
            os.path.join(self.analysis_folder, 'pandda_inspect_sites.csv'))

        if not os.path.isfile(self.eventCSV):
            analyse_csv = self.eventCSV.replace(
                'pandda_inspect_events.csv', 'pandda_analyse_events.csv')
            if not os.path.isfile(analyse_csv):
                self.logger.error('cannot find {0!s}'.format(analyse_csv))
                return
            self.initialize_inspect_events_csv_file(analyse_csv)

        if not os.path.isfile(self.eventCSV):
            self.logger.error('cannot find {0!s}'.format(self.eventCSV))
            return

        if not os.path.isfile(self.siteCSV):
            analyse_csv = self.siteCSV.replace(
                'pandda_inspect_sites.csv', 'pandda_analyse_sites.csv')
            if not os.path.isfile(analyse_csv):
                self.logger.error('cannot find {0!s}'.format(analyse_csv))
                return
            self.initialize_inspect_sites_csv_file(analyse_csv)

        if not os.path.isfile(self.siteCSV):
            self.logger.error('cannot find {0!s}'.format(self.siteCSV))
            return

        self.parsepanddaDir()

    def make_secure_copy_of_original_csv(self, csv_file):
        backup = csv_file + '.original'
        if not os.path.isfile(backup):
            self.logger.info('backing up {0!s}'.format(csv_file))
            shutil.copy(csv_file, backup)

    def initialize_inspect_events_csv_file(self, analyse_csv):
        self.make_secure_copy_of_original_csv(analyse_csv)
        with _csv_open_r(analyse_csv) as f:
            rows = list(csv.reader(f))
        for i, row in enumerate(rows):
            if i == 0:
                rows[i].extend(['Interesting', 'Ligand Placed',
                                 'Ligand Confidence', 'Comment', 'Viewed'])
            else:
                rows[i].extend(['False', 'False', 'Low', 'None', 'False'])
        out = os.path.join(self.analysis_folder, 'pandda_inspect_events.csv')
        with _csv_open_w(out) as f:
            csv.writer(f).writerows(rows)

    def initialize_inspect_sites_csv_file(self, analyse_csv):
        self.make_secure_copy_of_original_csv(analyse_csv)
        with _csv_open_r(analyse_csv) as f:
            rows = list(csv.reader(f))
        for i, row in enumerate(rows):
            if i == 0:
                rows[i].extend(['Name', 'Comment'])
            else:
                rows[i].extend(['None', 'None'])
        out = os.path.join(self.analysis_folder, 'pandda_inspect_sites.csv')
        with _csv_open_w(out) as f:
            csv.writer(f).writerows(rows)

    def parsepanddaDir(self):
        self.logger.info("reading {0!s}".format(self.eventCSV))
        with _csv_open_r(self.eventCSV) as f:
            self.elist = list(csv.reader(f))

        self.logger.info("reading {0!s}".format(self.siteCSV))
        with _csv_open_r(self.siteCSV) as f:
            self.slist = list(csv.reader(f))

        for n, item in enumerate(self.elist[0]):
            if item == 'dtag':
                self.xtal_index = n
            elif item == 'Ligand Confidence':
                self.ligand_confidence_index = n
            elif item in ('event_num', 'event_idx'):
                self.event_index = n
            elif item in ('site_num', 'site_idx'):
                self.site_index = n
            elif item in ('bdc', '1-BDC'):
                self.bdc_index = n
            elif item == 'x':
                self.x_index = n
            elif item == 'y':
                self.y_index = n
            elif item == 'z':
                self.z_index = n
            elif item in ('analysed_resolution', 'high_resolution'):
                self.resolution_index = n
            elif item == 'r_work':
                self.r_work_index = n
            elif item == 'r_free':
                self.r_free_index = n
            elif item == 'Viewed':
                self.viewed_index = n
            elif item == 'cluster_size':
                self.cluster_size_index = n
            elif item == 'Ligand Placed':
                self.ligand_placed_index = n
            elif item == 'Comment':
                self.comment_index = n
            elif item == 'Interesting':
                self.interesting_index = n

        # Parse slist (sites) column indices
        if self.slist:
            for n, item in enumerate(self.slist[0]):
                if item in ('site_num', 'site_idx'):
                    self.slist_site_num_index = n
                elif item == 'Name':
                    self.slist_site_name_index = n
                elif item == 'Comment':
                    self.slist_site_comment_index = n

        self.show_content_of_event_csv_file()

    def show_content_of_event_csv_file(self):
        self.logger.info("contents of {0!s}:".format(self.eventCSV))
        for n in range(1, len(self.elist)):
            x = round(float(self.elist[n][self.x_index]), 1)
            y = round(float(self.elist[n][self.y_index]), 1)
            z = round(float(self.elist[n][self.z_index]), 1)
            self.logger.info(
                ' xtal: {0!s} - event/site: {1!s}/{2!s}'
                ' - BDC: {3!s} - x,y,z: {4!s},{5!s},{6!s}'
                ' - Res: {7!s} - Rwork/Rfree: {8!s}/{9!s}'
                ' - viewed: {10!s} - confidence: {11!s}'.format(
                    self.elist[n][self.xtal_index],
                    self.elist[n][self.event_index],
                    self.elist[n][self.site_index],
                    self.elist[n][self.bdc_index],
                    x, y, z,
                    self.elist[n][self.resolution_index],
                    self.elist[n][self.r_work_index],
                    self.elist[n][self.r_free_index],
                    self.elist[n][self.viewed_index],
                    self.elist[n][self.ligand_confidence_index],
                ))
        self.init_crystal_selection_combobox()

    # ------------------------------------------------------------------
    # Crystal / event combobox
    # ------------------------------------------------------------------

    def init_crystal_selection_combobox(self):
        self.logger.info('rebuilding crystal selection combobox')
        model = self.cb.get_model()
        if model is not None:
            while len(model) > 0:
                self.cb.remove(0)
        self.cb_list = []
        for n in range(1, len(self.elist)):
            text = '{0!s} - event: {1!s} - site: {2!s}'.format(
                self.elist[n][self.xtal_index],
                self.elist[n][self.event_index],
                self.elist[n][self.site_index],
            )
            self.cb_list.append(text)
            self.cb.append_text(text)

    def update_crystal_selection_combobox(self):
        x = self.elist[self.index][self.xtal_index]
        e = self.elist[self.index][self.event_index]
        s = self.elist[self.index][self.site_index]
        text = '{0!s} - event: {1!s} - site: {2!s}'.format(x, e, s)
        for n, i in enumerate(self.cb_list):
            if i == text:
                self.cb.set_active(n)
                break

    def select_crystal(self, widget):
        tmp = str(widget.get_active_text())
        if not tmp:
            return
        self.logger.info('selected: {0!s}'.format(tmp))
        tmpx = tmp.replace(' - event: ', ' ').replace(' - site: ', ' ')
        parts = tmpx.split()
        xtal, event, site = parts[0], parts[1], parts[2]
        index_increment = 0
        for n in range(len(self.elist)):
            if (self.elist[n][self.xtal_index] == xtal and
                    self.elist[n][self.event_index] == event and
                    self.elist[n][self.site_index] == site):
                index_increment = n - self.index
                break
        self.change_event(index_increment)

    # ------------------------------------------------------------------
    # File location helpers
    # ------------------------------------------------------------------

    def get_pdb(self, missing_files):
        pdb = ''
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        modelled = os.path.join(
            ds, 'modelled_structures',
            '{0!s}-pandda-model.pdb'.format(self.xtal))
        input_pdb = os.path.join(
            ds, '{0!s}-pandda-input.pdb'.format(self.xtal))
        if os.path.isfile(modelled):
            pdb = modelled
            self.logger.info('found pdb (modelled): {0!s}'.format(pdb))
        elif os.path.isfile(input_pdb):
            pdb = input_pdb
            self.logger.info('found pdb: {0!s}'.format(pdb))
        else:
            self.logger.error('did not find pdb file')
            missing_files = True
        return pdb, missing_files

    def load_pdb(self):
        coot.set_nomenclature_errors_on_read("ignore")
        imol = coot.handle_read_draw_molecule_with_recentre(self.pdb, 0)
        self.mol_dict['protein'] = imol
        coot.set_show_symmetry_master(1)
        coot.set_show_symmetry_molecule(imol, 1)

    def get_emap(self, missing_files):
        emap = ''
        new_pandda_output = False
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        event_number = (3 - len(str(self.event))) * '0' + str(self.event)
        candidates = [
            (os.path.join(ds, '{0!s}-event_{1!s}_1-BDC_{2!s}_map.native.mtz'.format(
                self.xtal, self.event, self.bdc)), False),
            (os.path.join(ds, '{0!s}-event_{1!s}_1-BDC_{2!s}_map.native.ccp4'.format(
                self.xtal, self.event, self.bdc)), False),
            (os.path.join(ds, '{0!s}-pandda-output-event-{1!s}.mtz'.format(
                self.xtal, event_number)), True),
        ]
        for path, is_new in candidates:
            if os.path.isfile(path):
                emap = path
                new_pandda_output = is_new
                self.logger.info('found event map: {0!s}'.format(emap))
                break
        if not emap:
            self.logger.error(
                'cannot find event map for {0!s} event {1!s}'.format(
                    self.xtal, self.event))
            missing_files = True
        return emap, new_pandda_output, missing_files

    def load_emap(self):
        if self.new_pandda_output:
            imol = coot.make_and_draw_map(
                self.emap, "FEVENT", "PHEVENT", "1", 0, 0)
        elif self.emap.endswith(".ccp4"):
            imol = coot.read_ccp4_map(self.emap, 0)
        else:
            imol = coot.make_and_draw_map(
                self.emap, "FWT", "PHWT", "1", 0, 0)
        self.mol_dict['emap'] = imol
        coot.set_colour_map_rotation_on_read_pdb(0)
        coot.set_last_map_colour(0, 0, 1)
        self.show_emap = 1
        coot.set_contour_level_in_sigma(
            self.mol_dict['emap'], 2.0 * (1.0 - float(self.bdc)))

    def get_zmap(self, missing_files):
        zmap = ''
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        for fname in (
            '{0!s}-z_map.native.mtz'.format(self.xtal),
            '{0!s}-z_map.native.ccp4'.format(self.xtal),
            '{0!s}-pandda-output.mtz'.format(self.xtal),
        ):
            p = os.path.join(ds, fname)
            if os.path.isfile(p):
                zmap = p
                self.logger.info('found z-map: {0!s}'.format(zmap))
                break
        if not zmap:
            self.logger.error('cannot find z-map')
            missing_files = True
        return zmap, missing_files

    def load_zmap(self):
        coot.set_default_initial_contour_level_for_difference_map(3)
        if self.new_pandda_output:
            imol = coot.make_and_draw_map(
                self.zmap, "FZVALUES", "PHZVALUES", "1", 0, 1)
            self.mol_dict['zmap'] = imol
            coot.set_map_is_difference_map(imol, True)
        elif self.zmap.endswith(".ccp4"):
            self.mol_dict['zmap'] = coot.read_ccp4_map(self.zmap, 1)
        else:
            imol = coot.auto_read_make_and_draw_maps(self.zmap)
            self.mol_dict['zmap'] = (
                imol[0] if isinstance(imol, (list, tuple)) else imol)
            coot.set_contour_level_in_sigma(self.mol_dict['zmap'], 3)
        self.show_zmap = 1

    def get_xraymap(self, missing_files):
        p = os.path.join(
            self.panddaDir, 'processed_datasets', self.xtal,
            '{0!s}-pandda-input.mtz'.format(self.xtal))
        if os.path.isfile(p):
            self.logger.info('found xray map: {0!s}'.format(p))
            return p, missing_files
        self.logger.error('did not find xray map')
        return '', True

    def load_xraymap(self):
        imol = coot.auto_read_make_and_draw_maps(self.xraymap)
        self.mol_dict['xraymap'] = imol
        coot.set_colour_map_rotation_on_read_pdb(0)
        if isinstance(imol, (list, tuple)):
            _set_map_displayed(imol[0], self.show_xraymap)
            _set_map_displayed(imol[1], self.show_xraymap)
        else:
            _set_map_displayed(imol, self.show_xraymap)

    def get_averagemap(self):
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        for fname in (
            '{0!s}-ground-state-average-map.native.mtz'.format(self.xtal),
            '{0!s}-ground-state-average-map.native.ccp4'.format(self.xtal),
            '{0!s}-pandda-output.mtz'.format(self.xtal),
        ):
            p = os.path.join(ds, fname)
            if os.path.isfile(p):
                self.logger.info('found average map: {0!s}'.format(p))
                self.toggle_average_map_button.set_sensitive(True)
                return p
        self.logger.warning('did not find average map; disabling button')
        self.toggle_average_map_button.set_sensitive(False)
        return ''

    def load_averagemap(self):
        if self.new_pandda_output:
            imol = coot.make_and_draw_map(
                self.zmap, "FGROUND", "PHGROUND", "1", 0, 0)
            self.mol_dict['averagemap'] = imol
        elif self.averagemap.endswith(".ccp4"):
            self.mol_dict['averagemap'] = coot.read_ccp4_map(
                self.averagemap, 0)
        else:
            imol = coot.auto_read_make_and_draw_maps(self.averagemap)
            self.mol_dict['averagemap'] = (
                imol[0] if isinstance(imol, (list, tuple)) else imol)
        coot.set_colour_map_rotation_on_read_pdb(0)
        _set_map_displayed(self.mol_dict['averagemap'], self.show_averagemap)
        coot.set_last_map_colour(0, 0, 1)

    def get_ligcif(self):
        ds = os.path.join(self.panddaDir, 'processed_datasets', self.xtal)
        if self.event:
            rhofit = os.path.join(
                ds, str(self.event), 'rhofit', 'best.cif')
            if os.path.isfile(rhofit):
                self.logger.info(
                    'found ligand cif (rhofit): {0!s}'.format(rhofit))
                return rhofit
        for ligcif in glob.glob(os.path.join(ds, 'ligand_files', '*.cif')):
            self.logger.info('found ligand cif: {0!s}'.format(ligcif))
            return ligcif
        self.logger.warning('no ligand cif in: {0!s}'.format(
            os.path.join(ds, 'ligand_files')))
        return ''

    def load_ligcif(self):
        if self.ligcif and os.path.isfile(self.ligcif):
            coot.read_cif_dictionary(self.ligcif)
            pdb = self.ligcif.replace('.cif', '.pdb')
            imol = coot.handle_read_draw_molecule_with_recentre(pdb, 0)
            self.mol_dict['ligand'] = imol
            coot.set_b_factor_residue_range(imol, "X", 1, 1, 20.0)
            coot.set_occupancy_residue_range(imol, "X", 1, 1, float(self.bdc))

    # ------------------------------------------------------------------
    # Event state
    # ------------------------------------------------------------------

    def reset_params(self):
        self.xtal = self.event = self.bdc = self.site = None
        self.pdb = self.emap = self.zmap = self.xraymap = self.averagemap = None
        self.new_pandda_output = False
        self.ligcif = ''
        self.x = self.y = self.z = None
        self.resolution = self.r_free = self.r_work = None
        self.ligand_confidence = None
        self.comment = ''
        self.interesting = False
        self.ligand_placed = False
        self.site_name = ''
        self.site_comment = ''
        self.merged = False

    def update_params(self):
        missing_files = False
        self.xtal  = self.elist[self.index][self.xtal_index]
        self.event = self.elist[self.index][self.event_index]
        self.bdc   = self.elist[self.index][self.bdc_index]
        self.site  = self.elist[self.index][self.site_index]
        self.logger.info(
            'checking files for {0!s} event {1!s} site {2!s}'.format(
                self.xtal, self.event, self.site))
        self.pdb, missing_files = self.get_pdb(missing_files)
        self.emap, self.new_pandda_output, missing_files = \
            self.get_emap(missing_files)
        self.zmap, missing_files = self.get_zmap(missing_files)
        self.xraymap, missing_files = self.get_xraymap(missing_files)
        self.averagemap = self.get_averagemap()
        self.ligcif = self.get_ligcif()
        self.x = float(self.elist[self.index][self.x_index])
        self.y = float(self.elist[self.index][self.y_index])
        self.z = float(self.elist[self.index][self.z_index])
        self.logger.info('event coords: {0!s}, {1!s}, {2!s}'.format(
            self.x, self.y, self.z))
        self.resolution = self.elist[self.index][self.resolution_index]
        self.r_free     = self.elist[self.index][self.r_free_index]
        self.r_work     = self.elist[self.index][self.r_work_index]
        self.ligand_confidence = \
            self.elist[self.index][self.ligand_confidence_index]
        if self.comment_index is not None:
            raw = self.elist[self.index][self.comment_index]
            self.comment = '' if raw in ('None', 'none', '') else raw
        if self.interesting_index is not None:
            self.interesting = (
                self.elist[self.index][self.interesting_index].lower() == 'true')
        if self.ligand_placed_index is not None:
            self.ligand_placed = (
                self.elist[self.index][self.ligand_placed_index].lower() == 'true')
        self._load_site_info()
        return missing_files

    def update_labels(self):
        self.info_labels['dataset_id'].set_label(str(self.xtal))
        self.info_labels['event_num'].set_label(str(self.event))
        self.info_labels['site_num'].set_label(str(self.site))
        self.info_labels['bdc'].set_label(str(self.bdc))
        self.info_labels['resolution'].set_label(str(self.resolution))
        self.info_labels['rfree_rwork'].set_label(
            '{0!s} / {1!s}'.format(self.r_free, self.r_work))
        self._updating_widgets = True
        self.comment_entry.set_text(self.comment)
        self.site_name_entry.set_text(self.site_name)
        self.site_comment_entry.set_text(self.site_comment)
        self._updating_widgets = False
        self._set_annotation_buttons()
        self._update_counters()

    def recentre_on_event(self):
        coot.set_rotation_centre(self.x, self.y, self.z)

    def current_sample_matches_selection_criteria(self):
        sc = self.selected_selection_criterion
        if sc is None or sc.startswith("show all events"):
            return True
        lc = self.ligand_confidence or ''
        if sc == 'show unassigned':
            return lc in ('unassigned', '', 'Low Confidence')
        if sc == 'show no ligands bound':
            return 'no ligand bound' in lc or lc == 'Low Confidence'
        if sc == 'show low confidence ligands':
            return lc in ('Low Confidence', 'low confidence', 'unassigned',
                          'no ligand bound', 'unknown ligand', 'ambiguous density')
        if sc == 'show medium confidence ligands':
            return lc in ('Medium Confidence', 'event map only')
        if sc == 'show high confidence ligands':
            return lc in ('High Confidence', '2fofc map')
        if sc == 'show not viewed events':
            return self.elist[self.index][self.viewed_index].lower() != 'true'
        return False

    # ------------------------------------------------------------------
    # Refresh (loads all files for the current event index)
    # ------------------------------------------------------------------

    def RefreshData(self):
        self.reset_params()

        for imol in _molecule_number_list():
            coot.close_molecule(imol)

        self.mol_dict = {
            'pdb': None, 'emap': None, 'zmap': None,
            'ligand': None, 'xraymap': None, 'averagemap': None,
        }
        self.show_emap = self.show_zmap = \
            self.show_xraymap = self.show_averagemap = 0

        if self.index < 1:
            self.index = 1
        if self.index > len(self.elist) - 1:
            self.index = len(self.elist) - 1
            self.logger.warning('reached end of events!')
            return None

        missing_files = self.update_params()

        if self.current_sample_matches_selection_criteria() and not missing_files:
            self.logger.info('loading {0!s} event {1!s}'.format(
                self.xtal, self.event))
            self._set_annotation_buttons()
            self.update_labels()
            self.recentre_on_event()
            self.load_ligcif()
            self.load_pdb()
            self.load_emap()
            self.load_zmap()
            self.load_xraymap()
            if self.averagemap:
                self.load_averagemap()
            self.logger.info('setting event map as RSR map')
            coot.set_imol_refinement_map(self.mol_dict['emap'])

        elif self.current_sample_matches_selection_criteria() and missing_files:
            self.logger.error('essential files missing; skipping...')
            self.change_event(1)
        else:
            self.logger.warning(
                '{0!s} event {1!s} does not match selection; skipping'.format(
                    self.xtal, self.event))
            self.change_event(1)

    # ------------------------------------------------------------------
    # Ligand modelling callbacks
    # ------------------------------------------------------------------

    def place_ligand_here(self, widget):
        self.logger.info('moving ligand to screen centre')
        _move_molecule_here(self.mol_dict['ligand'])

    def merge_ligand_into_protein(self, widget):
        self.logger.info('merging ligand into protein')
        _merge_molecules([self.mol_dict['ligand']], self.mol_dict['protein'])
        coot.close_molecule(self.mol_dict['ligand'])
        self.merged = True

    def reset_to_unfitted(self, widget):
        for imol in _molecule_number_list():
            if 'pandda-model.pdb' in coot.molecule_name(imol):
                self.pdb = os.path.join(
                    self.panddaDir, 'processed_datasets', self.xtal,
                    '{0!s}-pandda-input.pdb'.format(self.xtal))
                coot.close_molecule(imol)
                self.load_pdb()
                break

    def check_if_modelled_structures_folder_exists(self):
        p = os.path.join(self.panddaDir, 'processed_datasets',
                         self.xtal, 'modelled_structures')
        if not os.path.isdir(p):
            self.logger.info('creating {0!s}'.format(p))
            os.mkdir(p)

    def save_model(self, widget):
        """Save current protein model to modelled_structures (without advancing)."""
        if self.mol_dict.get('protein') is None or self.panddaDir is None:
            self.logger.warning('no protein loaded; cannot save model')
            return
        self.check_if_modelled_structures_folder_exists()
        base = os.path.join(self.panddaDir, 'processed_datasets',
                            self.xtal, 'modelled_structures')
        existing = sorted(glob.glob(os.path.join(base, 'fitted-v*.pdb')))
        if existing:
            nums = [int(os.path.basename(p)[8:12]) for p in existing]
            new_n = max(nums) + 1
        else:
            new_n = 1
        new = 'fitted-v{0:04d}.pdb'.format(new_n)
        coot.write_pdb_file(self.mol_dict['protein'],
                            os.path.join(base, new))
        model_link = os.path.join(
            base, '{0!s}-pandda-model.pdb'.format(self.xtal))
        if os.path.isfile(model_link):
            os.remove(model_link)
        os.chdir(base)
        os.system('/bin/cp {0!s} {1!s}-pandda-model.pdb'.format(new, self.xtal))
        self.logger.info('saved model: {0!s}'.format(new))
        if self.merged:
            if self.ligand_placed_index is not None:
                self.elist[self.index][self.ligand_placed_index] = 'True'
                self._updating_widgets = True
                self.ligand_placed_radio_yes.set_active(True)
                self._updating_widgets = False
            self.save_pandda_inspect_events_csv_file()

    def save_next(self, widget):
        """Backward-compat alias for save_model."""
        self.save_model(widget)

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------

    def select_events(self, widget):
        self.selected_selection_criterion = \
            self.select_events_combobox.get_active_text()
        self.crystal_progressbar.set_fraction(0)
        header = self.elist[0]
        data   = self.elist[1:]
        sc = self.selected_selection_criterion
        if sc == "show all events - sort by cluster size":
            data = sorted(data, key=lambda x: x[self.cluster_size_index])
        elif sc == "show all events - sort alphabetically":
            self.logger.info("sorting alphabetically")
            data = sorted(data, key=lambda x: x[self.xtal_index])
        self.elist = [header] + data
        self.init_crystal_selection_combobox()
        self.logger.info("selection: {0!s}".format(sc))
        self.index = -1

    def previous_event(self, widget):
        self.save_comment(None)
        self.change_event(-1)

    def next_event_no_save(self, widget):
        """Advance to next event, saving annotation but not the model."""
        self.save_comment(None)
        self.save_event_as_viewed()
        self.change_event(1)

    def next_event_save(self, widget):
        """Save model then advance to next event."""
        self.save_comment(None)
        self.save_model(None)
        self.save_event_as_viewed()
        self.change_event(1)

    # kept for backward compatibility
    def next_event(self, widget):
        self.next_event_no_save(widget)

    def previous_site(self, widget):
        self.logger.info('moving to previous site')
        self.change_site(-1)

    def next_site(self, widget):
        self.logger.info('moving to next site')
        self.change_site(1)

    def change_site(self, n):
        current_site = int(self.site)
        new_site = current_site + n
        index_increment = 0
        for i, item in enumerate(self.elist):
            if item[self.site_index] == str(new_site):
                index_increment = i - self.index
                break
        self.change_event(index_increment)

    def change_event(self, n):
        self.index += n
        total = max(len(self.elist) - 1, 1)
        self.crystal_progressbar.set_fraction(
            min(1.0, float(self.index) / float(total)))
        self.update_crystal_selection_combobox()
        self.RefreshData()

    # ------------------------------------------------------------------
    # Map toggle callbacks
    # ------------------------------------------------------------------

    def toggle_emap(self, widget):
        if self.mol_dict['emap'] is not None:
            self.show_emap = 1 - self.show_emap
            _set_map_displayed(self.mol_dict['emap'], self.show_emap)

    def toggle_zmap(self, widget):
        if self.mol_dict['zmap'] is not None:
            self.show_zmap = 1 - self.show_zmap
            _set_map_displayed(self.mol_dict['zmap'], self.show_zmap)

    def toggle_x_ray_maps(self, widget):
        if self.mol_dict['xraymap'] is not None:
            self.show_xraymap = 1 - self.show_xraymap
            xm = self.mol_dict['xraymap']
            if isinstance(xm, (list, tuple)):
                _set_map_displayed(xm[0], self.show_xraymap)
                _set_map_displayed(xm[1], self.show_xraymap)
            else:
                _set_map_displayed(xm, self.show_xraymap)

    def toggle_average_map(self, widget):
        if self.mol_dict['averagemap'] is not None:
            self.show_averagemap = 1 - self.show_averagemap
            _set_map_displayed(
                self.mol_dict['averagemap'], self.show_averagemap)

    # ------------------------------------------------------------------
    # New helper methods
    # ------------------------------------------------------------------

    def _load_site_info(self):
        """Populate self.site_name / self.site_comment from slist."""
        self.site_name = ''
        self.site_comment = ''
        if not hasattr(self, 'slist') or not self.slist:
            return
        if self.slist_site_num_index is None or self.site is None:
            return
        for row in self.slist[1:]:
            if str(row[self.slist_site_num_index]) == str(self.site):
                if self.slist_site_name_index is not None:
                    raw = row[self.slist_site_name_index]
                    self.site_name = '' if raw in ('None', 'none', '') else raw
                if self.slist_site_comment_index is not None:
                    raw = row[self.slist_site_comment_index]
                    self.site_comment = '' if raw in ('None', 'none', '') else raw
                break

    def _update_counters(self):
        """Refresh event/site counter labels and progress bar."""
        if not hasattr(self, 'elist') or not self.elist:
            return
        total_events = max(len(self.elist) - 1, 0)
        self.event_counter_label.set_label(
            'Event  {0!s}  of  {1!s}'.format(max(self.index, 0), total_events))
        total_sites = max(len(self.slist) - 1, 0) if hasattr(self, 'slist') and self.slist else 0
        self.site_counter_label.set_label(
            'Site  {0!s}  of  {1!s}'.format(
                self.site if self.site is not None else '-', total_sites))
        if total_events > 0:
            self.crystal_progressbar.set_fraction(
                min(1.0, float(max(self.index, 0)) / float(total_events)))

    def go_to_dataset(self, widget):
        """Jump to the first event for the dataset name in goto_entry."""
        if not hasattr(self, 'elist') or not self.elist:
            return
        target = self.goto_entry.get_text().strip()
        if not target:
            return
        for n in range(1, len(self.elist)):
            if self.elist[n][self.xtal_index] == target:
                self.change_event(n - self.index)
                return
        self.logger.warning('dataset not found: {0!s}'.format(target))

    def show_summary(self, widget):
        """Log a count of events in each annotation state."""
        if not hasattr(self, 'elist') or not self.elist:
            return
        total  = len(self.elist) - 1
        viewed = sum(1 for r in self.elist[1:]
                     if r[self.viewed_index].lower() == 'true')
        high = medium = low = 0
        if self.ligand_confidence_index is not None:
            for r in self.elist[1:]:
                lc = r[self.ligand_confidence_index]
                if lc in ('High Confidence', '2fofc map'):
                    high += 1
                elif lc in ('Medium Confidence', 'event map only'):
                    medium += 1
                else:
                    low += 1
        self.logger.info(
            'Summary: {0} total, {1} viewed | '
            'High: {2}  Medium: {3}  Low/unassigned: {4}'.format(
                total, viewed, high, medium, low))

    def update_html(self, widget):
        """Hook for XCE HTML update (no-op outside XCE context)."""
        self.logger.info('update_html: no-op outside XCE context')

    def next_unviewed_event(self, widget):
        """Jump to the next event with Viewed=False."""
        if not hasattr(self, 'elist') or not self.elist:
            return
        start = self.index + 1 if self.index >= 1 else 1
        for n in range(start, len(self.elist)):
            if self.elist[n][self.viewed_index].lower() != 'true':
                self.change_event(n - self.index)
                return
        self.logger.info('no more unviewed events')

    def next_modelled_event(self, widget):
        """Jump to the next event that has a saved pandda-model.pdb."""
        if not hasattr(self, 'elist') or not self.elist or self.panddaDir is None:
            return
        start = self.index + 1 if self.index >= 1 else 1
        for n in range(start, len(self.elist)):
            xtal = self.elist[n][self.xtal_index]
            model = os.path.join(
                self.panddaDir, 'processed_datasets', xtal,
                'modelled_structures', '{0!s}-pandda-model.pdb'.format(xtal))
            if os.path.isfile(model):
                self.change_event(n - self.index)
                return
        self.logger.info('no more modelled events')

    def reload_last_saved_model(self, widget):
        """Close current protein and reload the latest fitted-v*.pdb."""
        if self.panddaDir is None or self.xtal is None:
            return
        base = os.path.join(self.panddaDir, 'processed_datasets',
                            self.xtal, 'modelled_structures')
        existing = sorted(glob.glob(os.path.join(base, 'fitted-v*.pdb')))
        if not existing:
            self.logger.warning('no saved model for {0!s}'.format(self.xtal))
            return
        latest = existing[-1]
        if self.mol_dict.get('protein') is not None:
            coot.close_molecule(self.mol_dict['protein'])
        self.pdb = latest
        self.load_pdb()
        self.logger.info('reloaded: {0!s}'.format(latest))

    def open_next_ligand(self, widget):
        """Advance to the next event without saving the model."""
        self.save_comment(None)
        self.save_event_as_viewed()
        self.change_event(1)

    def save_site_info(self, widget):
        """Save site Name and Comment entries to pandda_inspect_sites.csv."""
        if getattr(self, '_updating_widgets', False) or self.site is None:
            return
        if not hasattr(self, 'slist') or not self.slist:
            return
        if self.slist_site_num_index is None:
            return
        name    = self.site_name_entry.get_text()
        comment = self.site_comment_entry.get_text()
        for row in self.slist[1:]:
            if str(row[self.slist_site_num_index]) == str(self.site):
                if self.slist_site_name_index is not None:
                    row[self.slist_site_name_index] = name
                if self.slist_site_comment_index is not None:
                    row[self.slist_site_comment_index] = comment
                break
        path = os.path.join(self.analysis_folder, 'pandda_inspect_sites.csv')
        with _csv_open_w(path) as f:
            csv.writer(f).writerows(self.slist)
        self.logger.info('saved site info for site {0!s}'.format(self.site))

    def load_input_mtz_clicked(self, widget):
        """Reload the input MTZ for the current crystal."""
        if self.xraymap and os.path.isfile(self.xraymap):
            imol = coot.auto_read_make_and_draw_maps(self.xraymap)
            self.mol_dict['xraymap'] = imol
            self.logger.info('loaded input mtz: {0!s}'.format(self.xraymap))
        else:
            self.logger.warning('input mtz not found')

    def load_average_map_clicked(self, widget):
        """Load the ground-state average map."""
        if self.averagemap and os.path.isfile(self.averagemap):
            self.load_averagemap()
        else:
            self.logger.warning('average map not found')

    def load_unfitted_model_clicked(self, widget):
        """Load the input (unfitted) PDB alongside the current model."""
        if self.panddaDir is None or self.xtal is None:
            return
        input_pdb = os.path.join(
            self.panddaDir, 'processed_datasets', self.xtal,
            '{0!s}-pandda-input.pdb'.format(self.xtal))
        if os.path.isfile(input_pdb):
            coot.handle_read_draw_molecule_with_recentre(input_pdb, 0)
            self.logger.info('loaded unfitted model: {0!s}'.format(input_pdb))
        else:
            self.logger.warning('unfitted model not found: {0!s}'.format(input_pdb))

    def create_new_ligand(self, widget):
        """Load a fresh copy of the ligand from its CIF file."""
        if self.ligcif and os.path.isfile(self.ligcif):
            self.load_ligcif()
            self.logger.info('created new ligand from {0!s}'.format(self.ligcif))
        else:
            self.logger.warning('no ligand CIF for {0!s}'.format(self.xtal))

    def CANCEL(self, widget):
        self.window.destroy()


# ---------------------------------------------------------------------------
# Entry point  (load via Coot scripting window or .coot.py)
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    inspect_gui().startGUI()

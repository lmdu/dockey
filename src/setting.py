from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from param import *
from utils import *
from widget import *

__all__ = ['DockeyGlobalSettingDialog']

class AutodockSettingTabWidget(QTabWidget):
	params = AutodockParameter()

	def __init__(self, parent):
		super().__init__(parent)
		self.settings = QSettings()
		self.input_widgets = []

		self.algorithms = {
			'global': "General",
			'LGA': "Lamarckian GA",
			'GA': "Genetic Algorithm",
			'SA': "Simulated Annealing",
			'LS': "Local Search"
		}

		self.create_algorithm_widgets()
		self.read_settings()

	def register_widget(self, widget, wgtype, option, default, convert, index=-1):
		self.input_widgets.append(AttrDict(
			widget = widget,
			wgtype = wgtype,
			option = option,
			default = default,
			convert = convert,
			index = index
		))

	def create_algorithm_widgets(self):
		for scope in ['global', 'LGA', 'GA', 'SA', 'LS']:
			widget = QWidget(self)
			layout = QFormLayout()
			#if scope != 'global':
			#	layout.addRow(QLabel("Set parameters for {}".format(self.algorithms[scope])))
			#layout.setContentsMargins(0, 0, 0, 0)
			widget.setLayout(layout)
			self.create_parameter_widgets(layout, scope)
			self.addTab(widget, self.algorithms[scope])

	def create_parameter_widgets(self, layout, scope='global'):
		group = "ADT{}".format(scope)
		#self.settings.beginGroup(group)

		for p, m in self.params.items():
			if not m.user:
				continue

			if scope not in m.scope:
				continue

			option = '{}/{}'.format(group, p)

			if m.type is int:
				editor = QSpinBox(self)
				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)
				editor.setValue(self.settings.value(option, m.default, int))
				self.register_widget(editor, 'spin', option, m.default, m.type)
				layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)
				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)
				editor.setValue(self.settings.value(option, m.default, float))
				self.register_widget(editor, 'spin', option, m.default, m.type)
				layout.addRow(m.comment, editor)

			elif p in ('geometric_schedule', 'set_sw1'):
				pass

			elif p == 'linear_schedule':
				btn_layout = QHBoxLayout()
				btn_group = QButtonGroup()
				line_btn = QRadioButton("Linear", self)
				line_btn.setChecked(self.settings.value(option, m.default, m.type))
				self.register_widget(line_btn, 'radio', option, m.default, m.type)
				geome_btn = QRadioButton("Geometric", self)
				option1 = '{}/{}'.format(group, 'geometric_schedule')
				geome_btn.setChecked(self.settings.value(option1, not m.default, m.type))
				self.register_widget(geome_btn, 'radio', option1, not m.default, m.type)
				btn_layout.addWidget(line_btn)
				btn_group.addButton(line_btn)
				btn_layout.addWidget(geome_btn)
				btn_group.addButton(geome_btn)
				btn_group.setExclusive(True)
				layout.addRow(m.comment, btn_layout)

			elif p == 'set_psw1':
				btn_layout = QHBoxLayout()
				btn_group = QButtonGroup()
				sw_btn = QRadioButton("classic", self)
				option1 = '{}/{}'.format(group, 'set_sw1')
				sw_btn.setChecked(self.settings.value(option1, not m.default, m.type))
				self.register_widget(sw_btn, 'radio', option1, not m.default, m.type)
				btn_group.addButton(sw_btn)
				btn_layout.addWidget(sw_btn)

				psw_btn = QRadioButton("pseudo", self)
				psw_btn.setChecked(self.settings.value(option, m.default, m.type))
				self.register_widget(psw_btn, 'radio', option, m.default, m.type)
				btn_group.addButton(psw_btn)
				btn_layout.addWidget(psw_btn)
				btn_group.setExclusive(True)
				layout.addRow(m.comment, btn_layout)

			elif p == 'parameter_file':
				editor = BrowseLineEdit(self)
				editor.set_text(self.settings.value(option, m.default, m.type))
				self.register_widget(editor, 'browse', option, m.default, m.type)
				layout.addRow(m.comment, editor)

			elif m.type is bool:
				editor = QCheckBox(self)
				editor.setChecked(self.settings.value(option, m.default, m.type))
				self.register_widget(editor, 'check', option, m.default, m.type)
				layout.addRow(m.comment, editor)

			elif m.type is str:
				editor = QLineEdit(self)
				editor.setText(self.settings.value(option, m.default, m.type))
				self.register_widget(editor, 'line', option, m.default, m.type)
				layout.addRow(m.comment, editor)

			elif 'choices' in m:
				editor = QComboBox(self)
				editor.addItems(m.choices)
				idx = editor.findText(self.settings.value(option, m.default, m.type))
				editor.setCurrentIndex(idx)
				self.register_widget(editor, 'select', option, m.default, m.type)
				layout.addRow(m.comment, editor)

			elif isinstance(m.type, list):
				vals = self.settings.value(option, m.default)
				defs = m.default

				for i, t in enumerate(m.type):
					if t is int:
						editor = QSpinBox(self)
						self.register_widget(editor, 'spin', option, defs[i], m.type, i)
						editor.setValue(int(vals[i]))
					elif t is float:
						editor = QDoubleSpinBox(self)
						self.register_widget(editor, 'spin', option, defs[i], m.type, i)
						editor.setValue(float(vals[i]))
					else:
						editor = QLineEdit(self)
						self.register_widget(editor, 'spin', option, defs[i], m.type, i)
						editor.setText(vals[i])

					if 'range' in m:
						_min, _max = m.range
						editor.setRange(_min, _max)

					if i == 0:
						layout.addRow(m.comment, editor)
					else:
						layout.addRow('', editor)

			else:
				layout.addRow(QLabel(m.comment))

		#self.settings.endGroup()

	def read_settings(self):
		for w in self.input_widgets:
			if w.convert is int:
				w.widget.setValue(self.settings.value(w.option, w.default, int))

			elif w.convert is float:
				w.widget.setValue(self.settings.value(w.option, w.default, float))

			elif w.convert is bool:
				w.widget.setChecked(self.settings.value(w.option, w.default, bool))

			elif 'parameter_file' in w.option:
				w.widget.set_text(self.settings.value(w.option, w.default, str))

			elif w.convert is str:
				w.widget.setText(self.settings.value(w.option, w.default, str))

			elif w.wgtype == 'select':
				idx = w.widget.findText(self.settings.value(w.option, w.default, str))
				w.widget.setCurrentIndex(idx)

			elif isinstance(w.convert, list):
				vals = self.settings.value(w.option)

				if vals:
					val = vals[w.index]
				else:
					val = w.default

				if w.convert[w.index] is str:
					w.widget.setText(val)
				elif w.convert[w.index] is float:
					w.widget.setValue(float(val))
				else:
					w.widget.setValue(int(val))

	def write_settings(self):
		for w in self.input_widgets:
			if w.convert is int:
				self.settings.setValue(w.option, w.widget.value())

			elif w.convert is float:
				self.settings.setValue(w.option, w.widget.value())

			elif w.convert is bool:
				self.settings.setValue(w.option, w.widget.isChecked())

			elif 'parameter_file' in w.option:
				self.settings.setValue(w.option, w.widget.get_text())

			elif w.convert is str:
				self.settings.setValue(w.option, w.widget.text())

			elif w.wgtype == 'select':
				self.settings.setValue(w.option, w.widget.currentText())

			elif isinstance(w.convert, list):
				vals = self.settings.value(w.option)

				if vals is None:
					vals = [''] * len(w.convert)

				if w.convert[w.index] is str:
					val = w.widget.text()
				else:
					val = w.widget.value()

				vals[w.index] = val

				self.settings.setValue(w.option, vals)

	def reset_settings(self):
		for scope in ['global', 'LGA', 'GA', 'SA', 'LS']:
			group = "ADT{}".format(scope)
			self.settings.remove(group)

		self.read_settings()

class AutodockVinaSettingWidget(QWidget):
	params = AutodockVinaParameter()

	def __init__(self, parent=None):
		super().__init__(parent)
		self.settings = QSettings()
		self.input_widgets = []

		self.create_algorithm_widgets()
		self.read_settings()

	def register_widget(self, widget, wgtype, option, default, convert, index=-1):
		self.input_widgets.append(AttrDict(
			widget = widget,
			wgtype = wgtype,
			option = option,
			default = default,
			convert = convert,
			index = index
		))

	def create_algorithm_widgets(self):
		param_layout = QFormLayout()
		self.setLayout(param_layout)
		self.create_parameter_widgets(param_layout)

	def create_parameter_widgets(self, param_layout):
		for p, m in self.params.items():
			if not m.user:
				continue

			option = "VINA/{}".format(p)

			if m.type is int:
				editor = QSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				self.register_widget(editor, 'spin', option, m.default, m.type)
				param_layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				self.register_widget(editor, 'spin', option, m.default, m.type)
				param_layout.addRow(m.comment, editor)

			elif p == 'scoring':
				editor = QComboBox(self)
				editor.addItems(m.choices)
				text = self.settings.value(option, m.default, m.type)

				if text == 'ad4':
					text = 'autodock4'

				idx = editor.findText(text)
				editor.setCurrentIndex(idx)
				self.register_widget(editor, 'select', option, m.default, m.type)
				param_layout.addRow(m.comment, editor)

	def read_settings(self):
		for w in self.input_widgets:
			if w.convert is int:
				w.widget.setValue(self.settings.value(w.option, w.default, int))

			elif w.convert is float:
				w.widget.setValue(self.settings.value(w.option, w.default, float))

			elif w.wgtype == 'select':
				text = self.settings.value(w.option, w.default, str)

				if text == 'ad4':
					text = 'autodock4'

				idx = w.widget.findText(text)
				w.widget.setCurrentIndex(idx)

	def write_settings(self):
		for w in self.input_widgets:
			if w.convert is int:
				self.settings.setValue(w.option, w.widget.value())

			elif w.convert is float:
				self.settings.setValue(w.option, w.widget.value())

			elif w.wgtype == 'select':
				text = w.widget.currentText()

				if text == 'autodock4':
					text = 'ad4'

				self.settings.setValue(w.option, text)

	def reset_settings(self):
		self.settings.remove('VINA')
		self.read_settings()

class QuickVinaSettingWidget(QWidget):
	params = QuickVinaParameter()

	def __init__(self, parent=None):
		super().__init__(parent)
		self.settings = QSettings()
		self.input_widgets = []

		self.create_algorithm_widgets()
		self.read_settings()

	def register_widget(self, widget, wgtype, option, default, convert, index=-1):
		self.input_widgets.append(AttrDict(
			widget = widget,
			wgtype = wgtype,
			option = option,
			default = default,
			convert = convert,
			index = index
		))

	def create_algorithm_widgets(self):
		param_layout = QFormLayout()
		self.setLayout(param_layout)
		self.create_parameter_widgets(param_layout)

	def create_parameter_widgets(self, param_layout):
		for p, m in self.params.items():
			if not m.user:
				continue

			option = "QVINA/{}".format(p)

			if m.type is int:
				editor = QSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				self.register_widget(editor, 'spin', option, m.default, m.type)
				param_layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				self.register_widget(editor, 'spin', option, m.default, m.type)
				param_layout.addRow(m.comment, editor)

	def read_settings(self):
		for w in self.input_widgets:
			if w.convert is int:
				w.widget.setValue(self.settings.value(w.option, w.default, int))

			elif w.convert is float:
				w.widget.setValue(self.settings.value(w.option, w.default, float))

	def write_settings(self):
		for w in self.input_widgets:
			if w.convert is int:
				self.settings.setValue(w.option, w.widget.value())

			elif w.convert is float:
				self.settings.setValue(w.option, w.widget.value())

	def reset_settings(self):
		self.settings.remove('QVINA')
		self.read_settings()

class ADTLigandPreparationConfig(DockeyConfigPage):
	def create_config_inputs(self):
		#prepare_ligand4.py settings
		prelig_page = QWidget(self)
		prelig_layout = QVBoxLayout()
		prelig_page.setLayout(prelig_layout)
		#stack_widget.addWidget(prelig_page)
		self.main_layout.addWidget(prelig_page)

		b_h_radio = QRadioButton('Build bonds and add hydrogens', self)
		b_radio = QRadioButton('Build a single bond from each atom with no bonds to its closest neighbor', self)
		h_radio = QRadioButton('Add hydrogens', self)
		c_h_radio = QRadioButton('Add hydrogens only if there are none already', self)
		non_radio = QRadioButton('Do not make any repairs', self)

		self.register_widget(b_h_radio, 'radio', 'ADTLigand/repairs', 'bonds_hydrogens', str)
		self.register_widget(b_radio, 'radio', 'ADTLigand/repairs', 'bonds', str)
		self.register_widget(h_radio, 'radio', 'ADTLigand/repairs', 'hydrogens', str)
		self.register_widget(c_h_radio, 'radio', 'ADTLigand/repairs', 'checkhydrogens', str)
		self.register_widget(non_radio, 'radio', 'ADTLigand/repairs', '', str)

		repair_group = QButtonGroup(self)
		repair_group.addButton(b_h_radio)
		repair_group.addButton(b_radio)
		repair_group.addButton(h_radio)
		repair_group.addButton(c_h_radio)
		repair_group.addButton(non_radio)

		prelig_layout.addWidget(QLabel("<b>Types of repairs to make:</b>", self))
		repair_layout = QGridLayout()
		repair_layout.addWidget(non_radio, 0, 0)
		repair_layout.addWidget(b_h_radio, 0, 1)
		repair_layout.addWidget(h_radio, 1, 0)
		repair_layout.addWidget(c_h_radio, 1, 1)
		repair_layout.addWidget(b_radio, 2, 0, 1, 2)
		prelig_layout.addLayout(repair_layout)

		n_c_radio = QRadioButton("Preserve all input charges, do not add new charges", self)
		a_g_radio = QRadioButton("Add gasteiger charges", self)

		self.register_widget(n_c_radio, 'radio', 'ADTLigand/charges_to_add', None, str)
		self.register_widget(a_g_radio, 'radio', 'ADTLigand/charges_to_add', 'gasteiger', str)

		charge_group = QButtonGroup(self)
		charge_group.addButton(a_g_radio)
		charge_group.addButton(n_c_radio)

		charge_layout = QHBoxLayout()
		charge_layout.addWidget(a_g_radio)
		charge_layout.addWidget(n_c_radio)
		prelig_layout.addWidget(QLabel("<b>Charges:</b>", self))
		prelig_layout.addLayout(charge_layout)

		nphs_check = QCheckBox("Merge charges and remove non-polar hydrogens")
		lps_check = QCheckBox("Merge charges and remove lone pairs")

		self.register_widget(nphs_check, 'check', 'ADTLigand/cleanup', 'nphs', str)
		self.register_widget(lps_check, 'check', 'ADTLigand/cleanup', 'lps', str)
	
		prelig_layout.addWidget(QLabel("<b>Clean types:</b>", self))
		prelig_layout.addWidget(nphs_check)
		prelig_layout.addWidget(lps_check)

		bb_check = QCheckBox("backbone", self)
		am_check = QCheckBox("amide", self)
		gd_check = QCheckBox("guaninium", self)

		self.register_widget(bb_check, 'check', 'ADTLigand/allowed_bonds', 'backbone', str)
		self.register_widget(am_check, 'check', 'ADTLigand/allowed_bonds', 'amide', str)
		self.register_widget(gd_check, 'check', 'ADTLigand/allowed_bonds', 'guaninium', str)

		allow_layout = QHBoxLayout()
		allow_layout.addWidget(QLabel("Types of bonds to allow to rotate:", self))
		allow_layout.addWidget(bb_check)
		allow_layout.addWidget(am_check)
		allow_layout.addWidget(gd_check)
		prelig_layout.addWidget(QLabel("<b>Bounds:</b>", self))
		prelig_layout.addLayout(allow_layout)

		f_check = QCheckBox("Check for and use largest nonbonded fragment", self)
		z_check = QCheckBox("Inactivate all active torsions", self)
		g_check = QCheckBox("Attach all nonbonded fragments", self)
		s_check = QCheckBox("Attach all nonbonded singletons", self)

		self.register_widget(f_check, 'check', 'ADTLigand/check_for_fragments', False, bool)
		self.register_widget(z_check, 'check', 'ADTLigand/inactivate_all_torsions', False, bool)
		self.register_widget(g_check, 'check', 'ADTLigand/attach_nonbonded_fragments', False, bool)
		self.register_widget(s_check, 'check', 'ADTLigand/attach_singletons', False, bool)

		bound_layout = QGridLayout()
		bound_layout.addWidget(f_check, 0, 0)
		bound_layout.addWidget(z_check, 0, 1)
		bound_layout.addWidget(g_check, 1, 0)
		bound_layout.addWidget(s_check, 1, 1)
		prelig_layout.addLayout(bound_layout)

	def read_settings(self):
		for i in self.input_widgets:
			if i.wgtype == 'radio':
				if 'repairs' in i.option:
					val = self.settings.value(i.option, 'bonds_hydrogens')

				elif 'charges_to_add' in i.option:
					val = self.settings.value(i.option, 'gasteiger')

					if val == 'None':
						val = None

				else:
					val = self.settings.value(i.option, i.default, i.convert)

				i.widget.setChecked(val == i.default)

			elif i.wgtype == 'check':
				if 'cleanup' in i.option:
					val = self.settings.value(i.option, 'nphs_lps_waters_nonstdres')
					val = i.default in val

				elif i.option == 'Ligand/allowed_bonds':
					val = self.settings.value(i.option, 'backbone')
					val = i.default in val

				else:
					val = self.settings.value(i.option, i.default, i.convert)

				if val:
					i.widget.setCheckState(Qt.Checked)

				else:
					i.widget.setCheckState(Qt.Unchecked)

			elif i.wgtype == 'edit':
				val = self.settings.value(i.option, i.default)
				i.widget.setText(val)

			elif i.wgtype == 'spin':
				val = self.settings.value(i.option, i.default, i.convert)
				i.widget.setValue(val)

			elif i.wgtype == 'select':
				val = self.settings.value(i.option, i.default)
				index = i.widget.findText(val)
				i.widget.setCurrentIndex(index)

	def write_settings(self):
		for i in self.input_widgets:
			if i.wgtype == 'radio':
				if i.widget.isChecked():
					self.settings.setValue(i.option, i.default)

			elif i.wgtype == 'check':
				if 'cleanup' in i.option or 'allowed_bonds' in i.option:
					vals = self.settings.value(i.option, '')

					if vals:
						vals = vals.split('_')
					else:
						vals = []

					if i.widget.isChecked() and i.default not in vals:
						vals.append(i.default)
					elif not i.widget.isChecked() and i.default in vals:
						vals.remove(i.default)

					self.settings.setValue(i.option, '_'.join(vals))

				else:
					self.settings.setValue(i.option, i.widget.isChecked())

			elif i.wgtype == 'edit':
				val = i.widget.text()
				self.settings.setValue(i.option, val)

			elif i.wgtype == 'spin':
				val = i.widget.value()
				self.settings.setValue(i.option, val)

			elif i.wgtype == 'select':
				val = i.widget.currentText()
				self.settings.setValue(i.option, val)

	def reset_settings(self):
		self.settings.remove('ADTLigand')
		self.read_settings()

class MeekoLigandPreparationConfig(DockeyConfigPage):
	def create_config_inputs(self):
		#meeko settings
		meeko_page = QWidget(self)
		meeko_layout = QVBoxLayout()
		meeko_page.setLayout(meeko_layout)
		#stack_widget.addWidget(meeko_page)
		self.main_layout.addWidget(meeko_page)

		#ahc_check = QCheckBox("Add hydrogens and 3D coordinates", self)
		rmc_check = QCheckBox("Keep macrocycles rigid in input conformation", self)
		kcr_check = QCheckBox("Return all rings from exhaustive perception", self)
		ker_check = QCheckBox("Equivalent rings have the same size and neighbors", self)
		hyd_check = QCheckBox("Add water molecules for hydrated docking", self)
		#knh_check = QCheckBox("Keep non-polar hydrogens (default: merge onto heavy atom)", self)
		far_check = QCheckBox("Allow amide bonds to rotate and be non-planar, which is bad", self)
		aim_check = QCheckBox("Write map of atom indices from input to pdbqt", self)
		rms_check = QCheckBox("Do not write smiles as remark to pdbqt", self)

		#self.register_widget(ahc_check, 'check', 'Meeko/add_h_3d', True, bool)
		self.register_widget(rmc_check, 'check', 'MeekoLigand/rigid_macrocycles', False, bool)
		self.register_widget(kcr_check, 'check', 'MeekoLigand/keep_chorded_rings', False, bool)
		self.register_widget(ker_check, 'check', 'MeekoLigand/keep_equivalent_rings', False, bool)
		self.register_widget(hyd_check, 'check', 'MeekoLigand/hydrate', False, bool)
		#self.register_widget(knh_check, 'check', 'MeekoLigand/keep_nonpolar_hydrogens', False, bool)
		self.register_widget(far_check, 'check', 'MeekoLigand/flexible_amides', False, bool)
		self.register_widget(aim_check, 'check', 'MeekoLigand/add_index_map', False, bool)
		self.register_widget(rms_check, 'check', 'MeekoLigand/remove_smiles', False, bool)

		#meeko_layout.addWidget(ahc_check)

		meeko_layout.addWidget(rmc_check)
		meeko_layout.addWidget(kcr_check)
		meeko_layout.addWidget(ker_check)
		meeko_layout.addWidget(hyd_check)
		#meeko_layout.addWidget(knh_check)
		meeko_layout.addWidget(far_check)
		meeko_layout.addWidget(aim_check)
		meeko_layout.addWidget(rms_check)

		mat_input = QLineEdit(self)
		mat_layout = QHBoxLayout()
		mat_layout.addWidget(QLabel("Merge these atom types:", self))
		mat_layout.addWidget(mat_input)
		mat_layout.addWidget(QLabel("<font color='gray'>use space to separate multiple atoms</font>", self))
		meeko_layout.addLayout(mat_layout)
		rbs_input = QLineEdit(self)
		rbs_input.setPlaceholderText('Multiple values can be separated by comma')
		rbs_layout = QHBoxLayout()
		rbs_layout.addWidget(QLabel("SMARTS patterns to rigidify bonds:", self))
		rbs_layout.addWidget(rbs_input)
		meeko_layout.addLayout(rbs_layout)
		rbi_input = QLineEdit(self)
		rbi_input.setPlaceholderText("Multiple values can be separated by comma, e.g. 1 2, 3 4. Corresponding to SMARTS patterns")
		meeko_layout.addWidget(QLabel("Indices of two atoms (in the SMARTS) that define a bond (start at 1)", self))
		meeko_layout.addWidget(rbi_input)
		ats_input = QLineEdit(self)
		ats_layout = QHBoxLayout()
		ats_layout.addWidget(QLabel("SMARTS based atom typing (JSON format):", self))
		ats_layout.addWidget(ats_input)
		meeko_layout.addLayout(ats_layout)

		dbp_spin = QSpinBox()
		dbp_spin.setRange(0, 100000)
		dbp_layout = QHBoxLayout()
		dbp_layout.addWidget(QLabel("Double bond penalty:", self))
		dbp_layout.addWidget(dbp_spin)
		dbp_layout.addWidget(QLabel("<font color='gray'><small>penalty > 100 prevents breaking double bonds</small></font>"))
		meeko_layout.addLayout(dbp_layout)

		self.register_widget(mat_input, 'edit', 'MeekoLigand/merge_these_atom_types', 'H', str)
		self.register_widget(rbs_input, 'edit', 'MeekoLigand/rigidify_bonds_smarts', '', str)
		self.register_widget(rbi_input, 'edit', 'MeekoLigand/rigidify_bonds_indices', '', str)
		self.register_widget(ats_input, 'edit', 'MeekoLigand/atom_type_smarts', '', str)
		self.register_widget(dbp_spin, 'spin', 'MeekoLigand/double_bond_penalty', 50, int)

	def read_settings(self):
		for i in self.input_widgets:
			if i.wgtype == 'radio':
				if 'repairs' in i.option:
					val = self.settings.value(i.option, 'bonds_hydrogens')

				elif 'charges_to_add' in i.option:
					val = self.settings.value(i.option, 'gasteiger')

					if val == 'None':
						val = None

				else:
					val = self.settings.value(i.option, i.default, i.convert)

				i.widget.setChecked(val == i.default)

			elif i.wgtype == 'check':
				if 'cleanup' in i.option:
					val = self.settings.value(i.option, 'nphs_lps_waters_nonstdres')
					val = i.default in val

				elif i.option == 'Ligand/allowed_bonds':
					val = self.settings.value(i.option, 'backbone')
					val = i.default in val

				else:
					val = self.settings.value(i.option, i.default, i.convert)

				if val:
					i.widget.setCheckState(Qt.Checked)

				else:
					i.widget.setCheckState(Qt.Unchecked)

			elif i.wgtype == 'edit':
				val = self.settings.value(i.option, i.default)
				i.widget.setText(val)

			elif i.wgtype == 'spin':
				val = self.settings.value(i.option, i.default, i.convert)
				i.widget.setValue(val)

			elif i.wgtype == 'select':
				val = self.settings.value(i.option, i.default)
				index = i.widget.findText(val)
				i.widget.setCurrentIndex(index)

	def write_settings(self):
		for i in self.input_widgets:
			if i.wgtype == 'radio':
				if i.widget.isChecked():
					self.settings.setValue(i.option, i.default)

			elif i.wgtype == 'check':
				if 'cleanup' in i.option or 'allowed_bonds' in i.option:
					vals = self.settings.value(i.option, '')

					if vals:
						vals = vals.split('_')
					else:
						vals = []

					if i.widget.isChecked() and i.default not in vals:
						vals.append(i.default)
					elif not i.widget.isChecked() and i.default in vals:
						vals.remove(i.default)

					self.settings.setValue(i.option, '_'.join(vals))

				else:
					self.settings.setValue(i.option, i.widget.isChecked())

			elif i.wgtype == 'edit':
				val = i.widget.text()
				self.settings.setValue(i.option, val)

			elif i.wgtype == 'spin':
				val = i.widget.value()
				self.settings.setValue(i.option, val)

			elif i.wgtype == 'select':
				val = i.widget.currentText()
				self.settings.setValue(i.option, val)

	def reset_settings(self):
		self.settings.remove('MeekoLigand')
		self.read_settings()

class AutodockToolSettingTabWidget(QTabWidget):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.receptor_widget = ReceptorPreparationConfigPage(self)
		self.ligand_widget = ADTLigandPreparationConfig(self)
		self.addTab(self.receptor_widget , "Prepare Receptors")
		self.addTab(self.ligand_widget, "Prepare Ligands")

	def write_settings(self):
		self.receptor_widget.write_settings()
		self.ligand_widget.write_settings()

	def reset_settings(self):
		self.receptor_widget.reset_settings()
		self.ligand_widget.reset_settings()

class MeekoSettingTabWidget(QTabWidget):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.ligand_widget = MeekoLigandPreparationConfig(self)
		self.addTab(self.ligand_widget, "Prepare Ligands")

	def write_settings(self):
		self.ligand_widget.write_settings()

	def reset_settings(self):
		self.ligand_widget.reset_settings()

class PDBFixerSettingConfig(DockeyConfigPage):
	def create_config_inputs(self):
		self.create_pdbfix_widgets()
		self.main_layout.addSpacing(20)
		self.main_layout.addStretch(1)

	def create_pdbfix_widgets(self):
		pdbfix_grid = QGridLayout()
		pdbfix_grid.setColumnStretch(1, 1)
		pdbfix_grid.setColumnStretch(2, 1)
		self.main_layout.addLayout(pdbfix_grid)
		#pdbfix_label = QLabel("<b>Use PDBFixer to fix receptor PDB file</b>", self)
		#use_pdbfix = QCheckBox(self)
		#self.register_widget(use_pdbfix, 'check', 'PDBFixer/use_pdbfix', False, bool)

		#pdbfix_grid.addWidget(use_pdbfix, 0, 0)
		#pdbfix_grid.addWidget(pdbfix_label, 0, 1)
		
		replace_nonres = QCheckBox("Replace non-standard residues", self)
		#replace_nonres.setEnabled(False)
		self.register_widget(replace_nonres, 'check', 'PDBFixer/replace_nonres', True, bool)
		pdbfix_grid.addWidget(replace_nonres, 1, 1)
		#use_pdbfix.stateChanged.connect(replace_nonres.setEnabled)

		remove_heterogen = QCheckBox("Remove all heterogens including water", self)
		#remove_heterogen.setEnabled(False)
		self.register_widget(remove_heterogen, 'check', 'PDBFixer/remove_heterogen', True, bool)
		pdbfix_grid.addWidget(remove_heterogen, 2, 1)
		#use_pdbfix.stateChanged.connect(remove_heterogen.setEnabled)

		add_misheavy = QCheckBox("Add missing heavy atoms", self)
		#add_misheavy.setEnabled(False)
		self.register_widget(add_misheavy, 'check', 'PDBFixer/add_misheavy', True, bool)
		pdbfix_grid.addWidget(add_misheavy, 3, 1)
		#use_pdbfix.stateChanged.connect(add_misheavy.setEnabled)

		add_mishydrogen = QCheckBox("Add missing hydrogen atoms, the pH to use for adding hydrogens:", self)
		#add_mishydrogen.setEnabled(False)
		self.register_widget(add_mishydrogen, 'check', 'PDBFixer/add_mishydrogen', True, bool)
		pdbfix_grid.addWidget(add_mishydrogen, 4, 1)
		mishydrogen_ph = QDoubleSpinBox(self)
		#mishydrogen_ph.setEnabled(False)
		mishydrogen_ph.setRange(0, 14)
		mishydrogen_ph.setDecimals(1)
		self.register_widget(mishydrogen_ph, 'spin', 'PDBFixer/mishydrogen_ph', 7.0, float)
		pdbfix_grid.addWidget(mishydrogen_ph, 4, 2)
		#use_pdbfix.stateChanged.connect(add_mishydrogen.setEnabled)
		#use_pdbfix.stateChanged.connect(mishydrogen_ph.setEnabled)

	def read_settings(self):
		for i in self.input_widgets:
			val = self.settings.value(i.option, i.default, i.convert)

			if i.wgtype == 'check':
				if val:
					i.widget.setCheckState(Qt.Checked)

				else:
					i.widget.setCheckState(Qt.Unchecked)

			elif i.wgtype == 'spin':
				i.widget.setValue(val)

			elif i.wgtype == 'select':
				i.widget.setCurrentIndex(i.widget.findText(val))

	def write_settings(self):
		for i in self.input_widgets:
			if i.wgtype == 'check':
				val = i.widget.isChecked()

			elif i.wgtype == 'spin':
				val = i.widget.value()

			elif i.wgtype == 'select':
				val = i.widget.currentText()

			self.settings.setValue(i.option, val)

	def reset_settings(self):
		self.settings.remove('PDBFixer')
		self.read_settings()

class PDB2PQRSettingConfig(DockeyConfigPage):
	def create_config_inputs(self):
		self.create_pdbpqr_widgets()
		self.main_layout.addSpacing(20)
		self.main_layout.addStretch(1)

	def create_pdbpqr_widgets(self):
		pdbpqr_grid = QGridLayout()
		pdbpqr_grid.setColumnStretch(1, 1)
		pdbpqr_grid.setColumnStretch(2, 1)
		self.main_layout.addLayout(pdbpqr_grid)
		#pdbpqr_label = QLabel("<b>Use PDB2PQR to convert receptor PDB file to PQR file</b>", self)
		#use_pdbpqr = QCheckBox(self)
		#self.register_widget(use_pdbpqr, 'check', 'PDB2PQR/use_pdbpqr', False, bool)

		#pdbpqr_grid.addWidget(use_pdbpqr, 0, 0)
		#pdbpqr_grid.addWidget(pdbpqr_label, 0, 1)

		pdbpqr_grid.addWidget(QLabel("choose a forcefield to use:", self), 1, 1)
		force_field = QComboBox(self)
		force_field.addItems(['AMBER', 'CHARMM', 'PEOEPB', 'PARSE', 'SWANSON', 'TYL06'])
		force_field.setCurrentIndex(3)
		#force_field.setEnabled(False)
		pdbpqr_grid.addWidget(force_field, 1, 2)
		self.register_widget(force_field, 'select', 'PDB2PQR/force_field', 'PARSE', str)
		#use_pdbpqr.stateChanged.connect(force_field.setEnabled)

		use_propka = QCheckBox("Use PROPKA to assign protonation states at provided pH:", self)
		#use_propka.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(use_propka.setEnabled)
		pdbpqr_grid.addWidget(use_propka, 2, 1)
		self.register_widget(use_propka, 'check', 'PDB2PQR/use_propka', True, bool)

		propka_ph = QDoubleSpinBox(self)
		propka_ph.setRange(0, 14)
		propka_ph.setDecimals(1)
		#propka_ph.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(propka_ph.setEnabled)
		self.register_widget(propka_ph, 'spin', 'PDB2PQR/propka_ph', 7.0, float)
		pdbpqr_grid.addWidget(propka_ph, 2, 2)

		node_bump = QCheckBox("Ensure that new atoms are not rebuilt too close to existing atoms", self)
		#node_bump.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(node_bump.setEnabled)
		pdbpqr_grid.addWidget(node_bump, 3, 1, 1, 2)
		self.register_widget(node_bump, 'check', 'PDB2PQR/node_bump', False, bool)

		no_hopt = QCheckBox("Optimize the hydrogen bonding network", self)
		#no_hopt.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(no_hopt.setEnabled)
		pdbpqr_grid.addWidget(no_hopt, 4, 1)
		self.register_widget(no_hopt, 'check', 'PDB2PQR/no_hopt', False, bool)

		remove_water = QCheckBox("Remove the waters from the output file", self)
		#remove_water.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(remove_water.setEnabled)
		pdbpqr_grid.addWidget(remove_water, 5, 1)
		self.register_widget(remove_water, 'check', 'PDB2PQR/remove_water', True, bool)

		neutraln = QCheckBox("Make the N-terminus of a protein neutral (only for PARSE)", self)
		#neutraln.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(neutraln.setEnabled)
		pdbpqr_grid.addWidget(neutraln, 6, 1)
		self.register_widget(neutraln, 'check', 'PDB2PQR/neutraln', False, bool)

		neutralc = QCheckBox("Make the C-terminus of a protein neutral (only for PARSE)", self)
		#neutralc.setEnabled(False)
		#use_pdbpqr.stateChanged.connect(neutralc.setEnabled)
		pdbpqr_grid.addWidget(neutralc, 7, 1)
		self.register_widget(neutralc, 'check', 'PDB2PQR/neutralc', False, bool)

	def read_settings(self):
		for i in self.input_widgets:
			val = self.settings.value(i.option, i.default, i.convert)

			if i.wgtype == 'check':
				if val:
					i.widget.setCheckState(Qt.Checked)

				else:
					i.widget.setCheckState(Qt.Unchecked)

			elif i.wgtype == 'spin':
				i.widget.setValue(val)

			elif i.wgtype == 'select':
				i.widget.setCurrentIndex(i.widget.findText(val))

	def write_settings(self):
		for i in self.input_widgets:
			if i.wgtype == 'check':
				val = i.widget.isChecked()

			elif i.wgtype == 'spin':
				val = i.widget.value()

			elif i.wgtype == 'select':
				val = i.widget.currentText()

			self.settings.setValue(i.option, val)

	def reset_settings(self):
		self.settings.remove('PDB2PQR')
		self.read_settings()

class DockeyGlobalSettingDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.setWindowTitle("Global Settings")
		self.list_widget = QListWidget(self)
		self.stack_widget = QStackedWidget(self)
		self.list_widget.currentRowChanged.connect(self.stack_widget.setCurrentIndex)

		self.button_box = QDialogButtonBox(
			QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel
		)

		self.button_box.accepted.connect(self.save_settings)
		self.button_box.rejected.connect(self.reject)
		self.button_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.reset_settings)

		main_layout = QVBoxLayout()
		main_layout.setSpacing(20)
		widget_layout = QHBoxLayout()
		widget_layout.setSpacing(20)
		widget_layout.addWidget(self.list_widget)
		widget_layout.addWidget(self.stack_widget, 1)
		main_layout.addLayout(widget_layout)
		main_layout.addWidget(self.button_box)
		self.setLayout(main_layout)

		self.create_pages()

	def create_pages(self):
		pages = [
			('Docking tools', DockingToolConfigPage),
			('AutoDock4', AutodockSettingTabWidget),
			('AutoDock Vina', AutodockVinaSettingWidget),
			('QuickVina-W', QuickVinaSettingWidget),
			('AutoDockTools', AutodockToolSettingTabWidget),
			('Meeko', MeekoSettingTabWidget),
			#('OpenBabel', QWidget),
			('PDBFixer', PDBFixerSettingConfig),
			('PDB2PQR', PDB2PQRSettingConfig)
		]

		list_width = 0
		for text, pager in pages:
			item = QListWidgetItem(text)
			self.list_widget.addItem(item)
			page = pager(self)
			self.stack_widget.addWidget(page)
			item_width = self.list_widget.visualItemRect(item).width()

			if item_width > list_width:
				list_width = item_width
		
		list_width += 10
		self.list_widget.setFixedWidth(list_width)

	def save_settings(self):
		for i in range(self.stack_widget.count()):
			widget = self.stack_widget.widget(i)
			widget.write_settings()

		self.accept()

	def reset_settings(self):
		for i in range(self.stack_widget.count()):
			widget = self.stack_widget.widget(i)
			widget.reset_settings()
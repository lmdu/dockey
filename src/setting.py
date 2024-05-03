from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from param import *
from widget import *

__all__ = ['DockeyGlobalSettingDialog']

class AutodockSettingTabWidget(QTabWidget):
	params = AutodockParameter()

	def __init__(self, parent):
		super().__init__(parent)

		self.algorithms = {
			'global': "General",
			'LGA': "Lamarckian GA",
			'GA': "Genetic Algorithm",
			'SA': "Simulated Annealing",
			'LS': "Local Search"
		}

		self.create_algorithm_widgets()

	def update_parameter(self, cmd, value, idx=-1):
		self.params.set_value(cmd, value, idx)

	def create_algorithm_widgets(self, index=0):
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
		for p, m in self.params.items():
			if not m.user:
				continue

			if scope not in m.scope:
				continue

			if m.type is int:
				editor = QSpinBox(self)
				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)
				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)
				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)
				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif p in ('geometric_schedule', 'set_sw1'):
				pass

			elif p == 'linear_schedule':
				btn_layout = QHBoxLayout()
				btn_group = QButtonGroup()
				line_btn = QRadioButton("Linear", self)
				line_btn.setChecked(m.value)
				line_btn.toggled.connect(lambda x: self.update_parameter('linear_schedule', x))
				geome_btn = QRadioButton("Geometric", self)
				geome_btn.toggled.connect(lambda x: self.update_parameter('geometric_schedule', x))
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
				sw_btn.toggled.connect(lambda x: self.update_parameter('set_sw1', x))
				btn_group.addButton(sw_btn)
				btn_layout.addWidget(sw_btn)

				psw_btn = QRadioButton("pseudo", self)
				psw_btn.setChecked(m.value)
				psw_btn.toggled.connect(lambda x: self.update_parameter('set_psw1', x))
				btn_group.addButton(psw_btn)
				btn_layout.addWidget(psw_btn)
				btn_group.setExclusive(True)
				layout.addRow(m.comment, btn_layout)

			elif p == 'parameter_file':
				editor = BrowseLineEdit(self)
				editor.text_changed.connect(lambda x: self.update_parameter('parameter_file', x))
				layout.addRow(m.comment, editor)

			elif m.type is bool:
				editor = QCheckBox(self)
				editor.setChecked(m.value)
				editor.stateChanged.connect(lambda x,p=p: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif m.type is str:
				editor = QLineEdit(self)
				editor.setText(m.value)
				editor.textChanged.connect(lambda x,p=p: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif 'choices' in m:
				editor = QComboBox(self)
				editor.addItems(m.choices)
				idx = editor.findText(m.value)
				editor.setCurrentIndex(idx)
				editor.currentTextChanged.connect(lambda x,p=p: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif isinstance(m.type, list):
				for i, t in enumerate(m.type):
					if t is int:
						editor = QSpinBox(self)
						editor.valueChanged.connect(lambda x,p=p,i=i: self.update_parameter(p, x, i))
						editor.setValue(m.value[i])
					elif t is float:
						editor = QDoubleSpinBox(self)
						editor.valueChanged.connect(lambda x,p=p,i=i: self.update_parameter(p, x, i))
						editor.setValue(m.value[i])
					else:
						editor = QLineEdit(self)
						editor.textChanged.connect(lambda x,p=p,i=i: self.update_parameter(p, x, i))
						editor.setText(m.value[i])

					if 'range' in m:
						_min, _max = m.range
						editor.setRange(_min, _max)

					if i == 0:
						layout.addRow(m.comment, editor)
					else:
						layout.addRow('', editor)

			else:
				layout.addRow(QLabel(m.comment))

class AutodockVinaSettingWidget(QWidget):
	params = AutodockVinaParameter()

	def __init__(self, parent=None):
		super().__init__(parent)
		self.create_algorithm_widgets()

	def update_parameter(self, key, val):
		if key == 'scoring' and val == 'autodock4':
			val = 'ad4'

		self.params[key]['value'] = val

	def create_algorithm_widgets(self):
		param_layout = QFormLayout()
		self.setLayout(param_layout)
		self.create_parameter_widgets(param_layout)

	def create_parameter_widgets(self, param_layout):
		for p, m in self.params.items():
			if not m.user:
				continue

			if m.type is int:
				editor = QSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				param_layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				param_layout.addRow(m.comment, editor)

			elif p == 'scoring':
				editor = QComboBox(self)
				editor.addItems(m.choices)
				editor.setCurrentIndex(0)
				editor.currentTextChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				param_layout.addRow(m.comment, editor)

class QuickVinaSettingWidget(QWidget):
	params = QuickVinaParameter()

	def __init__(self, parent=None):
		super().__init__(parent)
		self.create_algorithm_widgets()

	def update_parameter(self, key, val):
		self.params[key]['value'] = val

	def create_algorithm_widgets(self):
		param_layout = QFormLayout()
		self.setLayout(param_layout)
		self.create_parameter_widgets(param_layout)

	def create_parameter_widgets(self, param_layout):
		for p, m in self.params.items():
			if not m.user:
				continue

			if m.type is int:
				editor = QSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				param_layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				param_layout.addRow(m.comment, editor)

class AutodockToolSettingTabWidget(QTabWidget):
	def __init__(self, parent=None):
		super().__init__(parent)

		receptor_widget = ReceptorPreparationConfigPage(self)
		ligand_widget = LigandPreparationConfigPage(self)
		self.addTab(receptor_widget , "Prepare Receptors")
		self.addTab(ligand_widget, "Prepare Ligands")

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
			('Meeko', QWidget),
			('OpenBabel', QWidget),
			('PDBFixer', QWidget),
			('PDB2PQR', QWidget)
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
			widget.restore_settings()
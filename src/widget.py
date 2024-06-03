import os
import json
import psutil
import multiprocessing

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from plip.basic.remote import VisualizerData
from plip.structure.preparation import PDBComplex

from utils import *
from table import *
from backend import *

__all__ = ['BrowseInput', 'CreateProjectDialog', 'AutodockConfigDialog',
	'AutodockVinaConfigDialog', 'ExportImageDialog', 'FeedbackBrowser',
	'PymolSettingDialog', 'InteractionTabWidget', 'PoseTabWidget',
	'TaskManagerDialog', 'PDBDownloadDialog', 'ZINCDownloadDialog',
	'QuickVinaConfigDialog', 'LigandFilterDialog', 'FlexResiduesDialog',
	'PubchemDownloadDialog', 'ChemblDownloadDialog', 'DockeyConfigDialog',
	'AcknowledgementDialog', 'CPUAndMemoryViewDialog', 'DockingToolConfigPage',
	'ReceptorPreparationConfigPage', 'LigandPreparationConfigPage',
	'DockeyConfigPage', 'DockeyRunAutodockDialog'
]

class QHLine(QFrame):
	def __init__(self):
		super().__init__()
		self.setFrameShape(QFrame.HLine)
		self.setFrameShadow(QFrame.Sunken)

class BrowseInput(QWidget):
	text_changed = Signal(str)

	def __init__(self, parent=None, is_file=True, is_save=False, _filter=""):
		super().__init__(parent)
		self.filter = _filter
		self.input = QLineEdit(self)
		self.input.setReadOnly(True)
		self.browse = QPushButton(self)
		self.browse.setFlat(True)
		self.browse.setIcon(QIcon(":/icons/folder.svg"))

		self.input.textChanged.connect(self.on_text_chanaged)

		if is_file:
			if is_save:
				self.browse.clicked.connect(self.select_save)
			else:
				self.browse.clicked.connect(self.select_file)
		else:
			self.browse.clicked.connect(self.select_folder)

		layout = QHBoxLayout()
		layout.setSpacing(0)
		layout.setContentsMargins(0,0,0,0)
		layout.addWidget(self.input)
		layout.addWidget(self.browse)

		self.setLayout(layout)

	@Slot()
	def on_text_chanaged(self, text):
		self.text_changed.emit(text)

	@Slot()
	def select_save(self):
		exe, _ = QFileDialog.getSaveFileName(self, filter=self.filter)

		if exe:
			self.input.setText(exe)

	@Slot()
	def select_file(self):
		exe, _ = QFileDialog.getOpenFileName(self, filter=self.filter)

		if exe:
			self.input.setText(exe)

	@Slot()
	def select_folder(self):
		folder = QFileDialog.getExistingDirectory(self)

		if folder:
			self.input.setText(folder)

	def get_text(self):
		return self.input.text()

	def set_text(self, text):
		self.input.setText(text)

class CreateProjectDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Create a New Project")

		self.name_input = QLineEdit(self)
		self.location_input = BrowseInput(self, False)
		tips = (
			'<font color="gray">'
			'This will create a new folder with input name in the '
			'location you selected'
			'</font>'
		)

		self.tip_label = QLabel(tips, self)

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.check_input_name)
		action_box.rejected.connect(self.reject)

		layout = QFormLayout()
		layout.addRow("Name", self.name_input)
		layout.addRow("Location", self.location_input)
		layout.addRow(self.tip_label)
		layout.addRow(action_box)

		self.setLayout(layout)

	@Slot()
	def check_input_name(self):
		name = self.name_input.text()
		folder = self.location_input.get_text()
		if not name:
			return QMessageBox.warning(self, "Input Warning",
				"Please input a project name"
			)

		if not folder:
			return QMessageBox.warning(self, "Input Warning",
				"Please select a directory where to create project"
			)

		if ' ' in name:
			return QMessageBox.critical(self, "Input Error",
				"White-spaces are not allowed in project name"
			)
		
		if any(c in name for c in '\\/:*?"<>|'):
			return QMessageBox.critical(self, "Input Error",
				"Your input name has illegal character"
			)

		if os.path.exists(os.path.join(folder, name)):
			return QMessageBox.critical(self, "Input Error",
				"Your input project name <b>{}</b> exists in folder <b>{}</b>. Please input a new project name".format(
					name, folder
				)
			)

		self.accept()

	@classmethod
	def get_project_folder(cls, parent):
		dlg = cls(parent)
		if dlg.exec():
			name = dlg.name_input.text()
			location = dlg.location_input.get_text()
			return os.path.join(location, name)

class DockingToolSettingDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.resize(QSize(550, 100))
		self.settings = QSettings()
		self.ad4 = self.settings.value('Tools/autodock_4', '')
		self.ag4 = self.settings.value('Tools/autogrid_4', '')
		self.vina = self.settings.value('Tools/autodock_vina', '')
		self.qvina = self.settings.value('Tools/quick_vina_w', '')
		self.autodock_input = BrowseInput(self)
		self.autodock_input.set_text(self.ad4)
		self.autogrid_input = BrowseInput(self)
		self.autogrid_input.set_text(self.ag4)
		self.vina_input = BrowseInput(self)
		self.vina_input.set_text(self.vina)
		self.qvina_input = BrowseInput(self)
		self.qvina_input.set_text(self.qvina)

		self.layout = QVBoxLayout()
		self.layout.addWidget(QLabel("Autodock4 executable", self))
		self.layout.addWidget(self.autodock_input)
		self.layout.addWidget(QLabel("Autogrid4 executable", self))
		self.layout.addWidget(self.autogrid_input)
		self.layout.addWidget(QLabel("Autodock vina executable", self))
		self.layout.addWidget(self.vina_input)
		self.layout.addWidget(QLabel("QuickVina-W executable", self))
		self.layout.addWidget(self.qvina_input)

		btn_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
		btn_box.accepted.connect(self.save_settings)
		btn_box.rejected.connect(self.reject)
		self.layout.addWidget(btn_box)
		self.setLayout(self.layout)

	@Slot()
	def save_settings(self):
		self.settings = QSettings()
		ad4 = self.autodock_input.get_text()
		ag4 = self.autogrid_input.get_text()
		vina = self.vina_input.get_text()
		qvina = self.qvina_input.get_text()

		if ad4 != self.ad4:
			self.settings.setValue('Tools/autodock_4', ad4)

		if ag4 != self.ag4:
			self.settings.setValue('Tools/autogrid_4', ag4)

		if vina != self.vina:
			self.settings.setValue('Tools/autodock_vina', vina)

		if qvina != self.qvina:
			self.settings.setValue('Tools/quick_vina_w', qvina)

		self.accept()

class ProgramConfigDialog(QDialog):
	title = ""
	progs = []

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle(self.title)
		self.resize(QSize(500, 10))
		self.settings = QSettings()

		layout = QVBoxLayout()
		self.setLayout(layout)

		self.inputs = []

		for i, prog in enumerate(self.progs):
			k, l = prog
			v = self.settings.value(k, '')

			if v:
				self.inputs.append(None)
			else:
				self.inputs.append(BrowseInput(self))
				layout.addWidget(QLabel(l, self))
				layout.addWidget(self.inputs[i])

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.save_config)
		action_box.rejected.connect(self.reject)

		layout.addWidget(action_box)

	@Slot()
	def save_config(self):
		for put in self.inputs:
			if put is None:
				continue

			v = put.get_text()

			if not v:
				return

		for i, put in enumerate(self.inputs):
			if put is None:
				continue

			v = put.get_text()
			k = self.progs[i][0]

			self.settings.setValue(k, v)

		self.accept()

	@classmethod
	def check_executable(cls, parent):
		settings = QSettings()
		configs = [settings.value(k) for k, l in cls.progs]

		if all(configs):
			return True

		dlg = cls(parent)
		return dlg.exec()

class AutodockConfigDialog(ProgramConfigDialog):
	title = "Autodock Config"
	progs = [
		('Tools/autodock_4', 'Autodock executable'),
		('Tools/autogrid_4', 'Autogrid executable'),
		#('Tools/mgltools', 'MGLTools directory')
	]

class QuickVinaConfigDialog(ProgramConfigDialog):
	title = "QuickVina Config"
	progs = [
		('Tools/quick_vina_w', 'QuickVina-W executable')
	]

class AutodockVinaConfigDialog(ProgramConfigDialog):
	title = "Autodock Vina Config"
	progs = [
		('Tools/autodock_vina', 'Autodock Vina executable'),
		#('Tools/mgltools', 'MGLTools directory')
	]

class ExportImageDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.options = None

		self.setWindowTitle("Export as PNG Image")
		self.resize(QSize(450, 10))
		layout = QFormLayout()
		layout.setVerticalSpacing(10)
		self.setLayout(layout)

		self.file_input = BrowseInput(self, is_save=True, _filter="Image (*.png)")
		self.width_input = QSpinBox(self)
		self.width_input.setRange(0, 100000)
		self.width_input.setValue(1000)
		self.height_input = QSpinBox(self)
		self.height_input.setRange(0, 100000)
		self.height_input.setValue(1000)
		self.dpi_input = QComboBox(self)
		self.dpi_input.addItems(['72', '96', '150', '200', '300', '600'])
		self.dpi_input.setCurrentIndex(1)

		layout.addRow("Image file: ", self.file_input)

		sub_layout = QHBoxLayout()
		sub_layout.addWidget(QLabel("Width: ", self))
		sub_layout.addWidget(self.width_input, 1)
		sub_layout.addSpacing(20)
		sub_layout.addWidget(QLabel("Height: ", self))
		sub_layout.addWidget(self.height_input, 1)
		sub_layout.addSpacing(20)
		sub_layout.addWidget(QLabel("DPI: ", self))
		sub_layout.addWidget(self.dpi_input, 1)

		layout.addRow(sub_layout)

		action_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.save_image)
		action_box.rejected.connect(self.reject)

		layout.addRow(action_box)

	@Slot()
	def save_image(self):
		image = self.file_input.get_text()

		if not image:
			return QMessageBox.warning(self, "Input Warning",
				"Please select a save file"
			)

		self.options = AttrDict({
			'image': self.file_input.get_text(),
			'width': self.width_input.value(),
			'height': self.height_input.value(),
			'dpi': int(self.dpi_input.currentText()),
		})

		self.accept()

	@classmethod
	def save_to_png(cls, parent):
		dlg = cls(parent)
		if dlg.exec():
			return dlg.options

class FeedbackBrowser(QPlainTextEdit):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setReadOnly(True)

	def sizeHint(self):
		return QSize(300, 50)

class GradientColorBar(QWidget):
	def __init__(self, parent):
		super().__init__(parent)

	def sizeHint(self):
		return QSize(200, 20)

	def paintEvent(self, event):
		painter = QPainter(self)
		gradient = QLinearGradient(event.rect().topLeft(), event.rect().topRight())
		gradient.setColorAt(0, Qt.white)
		gradient.setColorAt(1, Qt.black)
		painter.fillRect(event.rect(), gradient)

class PymolSettingDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.settings = QSettings()
		self.setWindowTitle("Pymol Settings")
		self.resize(QSize(400, 100))
		self.bar = GradientColorBar(self)
		self.slider = QSlider(Qt.Horizontal, self)
		self.slider.setRange(0, 100)
		self.slider.valueChanged.connect(self.on_color_changed)
		self.layout = QVBoxLayout()
		self.layout.addWidget(QLabel("Background", self))
		self.layout.addWidget(self.bar)
		self.layout.addWidget(self.slider)
		self.setLayout(self.layout)
		val = self.settings.value('Pymol/background', 100, int)
		self.slider.setValue(val)

	@Slot()
	def on_color_changed(self):
		val = self.slider.value()

		if val == 0:
			color = 'white'
		elif val == 100:
			color = 'black'
		else:
			color = 'grey{:0>2d}'.format(100-val)

		self.parent.cmd.bg_color(color)
		self.settings.setValue('Pymol/background', val)

#class JobConcurrentSettingDialog(QDialog):
#	def __init__(self, parent=None):
#		super().__init__(parent)
#		self.parent = parent
#		self.setWindowTitle("Concurrent Job Setting")
#		self.settings = QSettings()
#		self.job_num = self.settings.value('Job/concurrent', 1, int)
#		self.layout = QVBoxLayout()
#		self.label = QLabel("Number of concurrent running jobs", self)
#		self.number = QSpinBox(self)
#		self.number.setValue(self.job_num)
#		self.number.setRange(1, max(1, psutil.cpu_count()-2))
#		self.number.valueChanged.connect(self.on_job_number_changed)
#		self.layout.addWidget(self.label)
#		self.layout.addWidget(self.number)
#		self.setLayout(self.layout)
#
#	@Slot()
#	def on_job_number_changed(self, num):
#		self.settings.setValue('Job/concurrent', num)

class PoseTabWidget(QTabWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent

		self.create_best_table()
		self.create_pose_table()
		self.create_error_text()

		self.setUsesScrollButtons(False)
		self.setTabBarAutoHide(True)

	def create_best_table(self):
		self.best_table = BestTableView(self.parent)
		self.best_model = BestTableModel()
		self.best_table.setModel(self.best_model)
		self.best_table.hideColumn(1)
		self.best_table.hideColumn(2)
		self.best_table.verticalHeader().hide()
		self.best_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
		self.best_table.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.best_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.addTab(self.best_table, "Best poses")

	def create_pose_table(self):
		self.pose_table = PoseTableView(self.parent)
		self.pose_model = PoseTableModel()
		self.pose_table.setModel(self.pose_model)
		self.pose_table.hideColumn(0)
		self.pose_table.hideColumn(1)
		self.pose_table.verticalHeader().hide()
		self.pose_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
		self.pose_table.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.pose_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.addTab(self.pose_table, "Task poses")
		self.setTabVisible(1, False)

	def create_error_text(self):
		self.pose_error = DockeyTextBrowser(self)
		self.addTab(self.pose_error, "Error")
		self.setTabVisible(2, False)

	def show_job_pose(self, title, job_id):
		self.pose_model.set_job(job_id)
		self.setTabText(1, title)
		self.setTabVisible(2, False)
		self.setTabVisible(1, True)
		self.setCurrentIndex(1)

	def show_job_error(self, title, message):
		self.pose_error.setPlainText(message)
		self.setTabText(2, title)
		self.setTabVisible(1, False)
		self.setTabVisible(2, True)
		self.setCurrentIndex(2)

	def select(self):
		self.best_model.select()

	def reset(self):
		self.pose_model.reset()
		self.best_model.reset()
		self.pose_error.clear()

	def clear(self):
		self.pose_model.clear()
		self.best_model.clear()
		self.pose_error.clear()

class InteractionTabWidget(QTabWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.pose_id = 0
		self.binding_site = None
		self.complex_pdb = None
		self.complex_vis = None

		self.setUsesScrollButtons(False)

		self.create_hydrogen_table()
		self.create_halogen_table()
		self.create_hydrophobic_table()
		self.create_salt_table()
		self.create_water_table()
		self.create_stacking_table()
		self.create_cation_table()
		self.create_metal_table()
		self.create_site_select()

		for i in range(8):
			self.setTabVisible(i, False)
			self.widget(i).hideColumn(0)
			self.widget(i).hideColumn(1)
			self.widget(i).clicked.connect(self.on_interaction_changed)

	def sizeHint(self):
		return QSize(300, 100)

	def create_hydrogen_table(self):
		self.hydrogen_table = InteractionTableView(self.parent)
		self.hydrogen_model = HydrogenBondsModel()
		self.hydrogen_table.setModel(self.hydrogen_model)
		self.addTab(self.hydrogen_table, "Hydrogen bonds")
		self.hydrogen_model.row_count.connect(lambda x: self.setTabVisible(0, x))

	def create_halogen_table(self):
		self.halogen_table = InteractionTableView(self.parent)
		self.halogen_model = HalogenBondsModel()
		self.halogen_table.setModel(self.halogen_model)
		self.addTab(self.halogen_table, "Halogen bonds")
		self.halogen_model.row_count.connect(lambda x: self.setTabVisible(1, x))

	def create_hydrophobic_table(self):
		self.hydrophobic_table = InteractionTableView(self.parent)
		self.hydrophobic_model = HydrophobicInteractionModel()
		self.hydrophobic_table.setModel(self.hydrophobic_model)
		self.addTab(self.hydrophobic_table, "hydrophobic interactions")
		self.hydrophobic_model.row_count.connect(lambda x: self.setTabVisible(2, x))

	def create_salt_table(self):
		self.salt_table = InteractionTableView(self.parent)
		self.salt_model = SaltBridgesModel()
		self.salt_table.setModel(self.salt_model)
		self.addTab(self.salt_table, "Salt bridges")
		self.salt_model.row_count.connect(lambda x: self.setTabVisible(3, x))
		self.salt_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
		self.salt_table.horizontalHeader().setStretchLastSection(True)

	def create_water_table(self):
		self.water_table = InteractionTableView(self.parent)
		self.water_model = WaterBridgesModel()
		self.water_table.setModel(self.water_model)
		self.addTab(self.water_table, "Water bridges")
		self.water_model.row_count.connect(lambda x: self.setTabVisible(4, x))

	def create_stacking_table(self):
		self.stacking_table = InteractionTableView(self.parent)
		self.stacking_model = PiStackingModel()
		self.stacking_table.setModel(self.stacking_model)
		self.addTab(self.stacking_table, "π-stacking")
		self.stacking_model.row_count.connect(lambda x: self.setTabVisible(5, x))
		self.stacking_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
		self.stacking_table.horizontalHeader().setStretchLastSection(True)

	def create_cation_table(self):
		self.cation_table = InteractionTableView(self.parent)
		self.cation_model = PiCationModel()
		self.cation_table.setModel(self.cation_model)
		self.addTab(self.cation_table, "π-cation")
		self.cation_model.row_count.connect(lambda x: self.setTabVisible(6, x))
		self.cation_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
		self.cation_table.horizontalHeader().setStretchLastSection(True)

	def create_metal_table(self):
		self.metal_table = InteractionTableView(self.parent)
		self.metal_model = MetalComplexModel()
		self.metal_table.setModel(self.metal_model)
		self.addTab(self.metal_table, "Metal complexes")
		self.metal_model.row_count.connect(lambda x: self.setTabVisible(7, x))

	def create_site_select(self):
		layout = QHBoxLayout()
		layout.setAlignment(Qt.AlignTop)
		layout.setContentsMargins(0,0,0,0)
		spacer = QSpacerItem(10, 1, QSizePolicy.Expanding, QSizePolicy.Minimum)
		layout.addItem(spacer)
		self.site_select = QComboBox(self)
		self.site_model = BindingSiteModel()
		self.site_select.setModel(self.site_model)
		self.site_select.setModelColumn(2)
		self.site_select.currentTextChanged.connect(self.on_binding_site_changed)
		layout.addWidget(QLabel("Binding site: ", self))
		layout.addWidget(self.site_select)
		self.setLayout(layout)

	def change_pose(self, pose):
		if pose != self.pose_id:
			self.pose_id = pose
			self.binding_site = None
			sql = "SELECT complex FROM pose WHERE id=? LIMIT 1"
			self.complex_pdb = DB.get_one(sql, (self.pose_id,))

		self.site_model.change_pose(pose)

		if self.complex_pdb:
			self.parent.cmd.delete('all')
			self.parent.cmd.initialize()
			self.parent.cmd.read_pdbstr(self.complex_pdb, 'complex')
			self.parent.cmd.color('cyan', 'organic')
			self.parent.cmd.color('red', 'elem O and organic')
			self.parent.cmd.color('blue', 'elem N and organic')
			self.parent.cmd.color('yellow', 'elem S and organic')

	@Slot(str)
	def on_binding_site_changed(self, site):
		sql = "SELECT id FROM binding_site WHERE pid=? AND site=? LIMIT 1"
		site_id = DB.get_one(sql, (self.pose_id, site))

		for i in range(8):
			self.widget(i).model().change_binding_site(site_id)

		if site != self.binding_site:
			self.binding_site = site
			temp_dir = QTemporaryDir()
			temp_dir.setAutoRemove(False)

			if not temp_dir.isValid():
				self.parent.show_error_message(
					"Could not create temporary work directory, {}".format(
					temp_dir.errorString()
				))
				return

			try:
				mol = PDBComplex()
				mol.output_path = temp_dir.path()
				mol.load_pdb(self.complex_pdb, as_string=True)

				for ligand in mol.ligands:
					mol.characterize_complex(ligand)

				self.complex_vis = VisualizerData(mol, self.binding_site)
			finally:
				temp_dir.remove()

	@Slot()
	def on_interaction_changed(self):
		if self.complex_vis:
			interaction_visualize(self.complex_vis)

	def clear(self):
		self.site_model.clear()

		for i in range(8):
			self.widget(i).model().clear()

	def reset(self):
		self.site_model.reset()

		for i in range(8):
			self.widget(i).model().reset()

		self.pose_id = 0
		self.binding_site = None
		self.complex_pdb = None
		self.complex_vis = None

class DockeyConfigPage(QWidget):
	def __init__(self, parent):
		super().__init__(parent)
		self.input_widgets = []
		self.settings = QSettings()

		self.main_layout = QVBoxLayout()
		self.setLayout(self.main_layout)

		self.create_config_inputs()
		self.read_settings()

	def create_config_inputs(self):
		pass

	def read_settings(self):
		pass

	def register_widget(self, widget, wgtype, option, default, convert):
		self.input_widgets.append(AttrDict(
			widget = widget,
			wgtype = wgtype,
			option = option,
			default = default,
			convert = convert
		))


class DockingToolConfigPage(DockeyConfigPage):
	def create_config_inputs(self):
		self.autodock_input = BrowseInput(self)
		self.autogrid_input = BrowseInput(self)
		self.vina_input = BrowseInput(self)
		self.qvina_input = BrowseInput(self)

		self.register_widget(self.autodock_input, 'browse', 'Tools/autodock_4', '', str)
		self.register_widget(self.autogrid_input, 'browse', 'Tools/autogrid_4', '', str)
		self.register_widget(self.vina_input, 'browse', 'Tools/autodock_vina',  '', str)
		self.register_widget(self.qvina_input, 'browse', 'Tools/quick_vina_w', '', str)

		self.main_layout.addWidget(QLabel("<b>1.</b>", self))
		self.main_layout.addWidget(QLabel("Autogrid4 executable", self))
		self.main_layout.addWidget(self.autogrid_input)
		self.main_layout.addWidget(QLabel("Autodock4 executable", self))
		self.main_layout.addWidget(self.autodock_input)

		self.main_layout.addSpacing(20)
		self.main_layout.addWidget(QHLine())
		self.main_layout.addSpacing(20)

		self.main_layout.addWidget(QLabel("<b>2.</b>", self))
		self.main_layout.addWidget(QLabel("Autodock Vina executable", self))
		self.main_layout.addWidget(self.vina_input)

		self.main_layout.addSpacing(20)
		self.main_layout.addWidget(QHLine())
		self.main_layout.addSpacing(20)

		self.main_layout.addWidget(QLabel("<b>3.</b>", self))
		self.main_layout.addWidget(QLabel("QuickVina-W executable", self))
		self.main_layout.addWidget(self.qvina_input)

		self.main_layout.addStretch(1)

		self.main_layout.addWidget(QLabel("<font color='red'>Tips: you are allowed to specify location for any one of these tools or for all tools</font>"))

	def read_settings(self):
		for i in self.input_widgets:
			val = self.settings.value(i.option, i.default, i.convert)
			i.widget.set_text(val)

	def write_settings(self):
		for i in self.input_widgets:
			val = i.widget.get_text()
			self.settings.setValue(i.option, val)

class TaskManagerDialog(QDialog):
	def __init__(self, parent):
		super().__init__(parent)
		self.setWindowTitle("Concurrent task manager")
		self.settings = QSettings()

		self.main_layout = QVBoxLayout()
		self.main_layout.setContentsMargins(20, 20, 20, 20)
		self.main_layout.setSpacing(10)
		self.setLayout(self.main_layout)

		job_layout = QGridLayout()
		job_layout.setColumnStretch(1, 1)
		self.main_layout.addLayout(job_layout)
		job_label = QLabel("Number of concurrent running tasks:", self)
		job_layout.addWidget(job_label, 0, 0)

		sys_cpus = psutil.cpu_count()

		self.job_num = QSpinBox(self)
		self.job_num.setRange(1, max(1, sys_cpus))
		job_layout.addWidget(self.job_num, 1, 1)

		job_layout.addWidget(QLabel("<font color='red'>Total CPU numbers: {}</font>".format(sys_cpus), self), 1, 2)

		self.main_layout.addWidget(QLabel("<b>Tips:</b>", self))
		tips = (
			"<p>1. You can set the number of tasks that can run concurrently to improve the docking speed.</p>"
			"<p>2. Limiting the number of concurrent tasks helps you reduce CPU consumption.</p>"
			"<p>3. You can set the number of concurrent tasks to any number from 1 to maximum CPU threads.</p>"
			"<p>4. For AutoDock Vina or QuickVina-W, the total CPU consumed is task number x vina CPU number.</p>"
		)
		tips_label = QLabel(tips, self)

		self.main_layout.addWidget(tips_label)

		self.btn_widget = QDialogButtonBox(
			QDialogButtonBox.Ok |
			QDialogButtonBox.Cancel
		)

		self.btn_widget.accepted.connect(self.write_settings)
		self.btn_widget.rejected.connect(self.reject)
		self.main_layout.addWidget(QHLine())
		self.main_layout.addWidget(self.btn_widget)

		self.read_settings()

	def read_settings(self):
		val = self.settings.value('Job/concurrent', 1, int)
		self.job_num.setValue(val)

	def write_settings(self):
		val = self.job_num.value()
		self.settings.setValue('Job/concurrent', val)
		self.accept()

class ReceptorPreprocessConfigPage(DockeyConfigPage):
	def create_config_inputs(self):
		self.create_pdbfix_widgets()
		self.main_layout.addSpacing(20)
		self.create_pdbpqr_widgets()
		self.main_layout.addStretch(1)

	def create_pdbfix_widgets(self):
		pdbfix_grid = QGridLayout()
		pdbfix_grid.setColumnStretch(1, 1)
		pdbfix_grid.setColumnStretch(2, 1)
		self.main_layout.addLayout(pdbfix_grid)
		pdbfix_label = QLabel("<b>Use PDBFixer to fix receptor PDB file</b>", self)
		use_pdbfix = QCheckBox(self)
		self.register_widget(use_pdbfix, 'check', 'PDBFixer/use_pdbfix', False, bool)

		pdbfix_grid.addWidget(use_pdbfix, 0, 0)
		pdbfix_grid.addWidget(pdbfix_label, 0, 1)
		
		replace_nonres = QCheckBox("Replace non-standard residues", self)
		replace_nonres.setEnabled(False)
		self.register_widget(replace_nonres, 'check', 'PDBFixer/replace_nonres', True, bool)
		pdbfix_grid.addWidget(replace_nonres, 1, 1)
		use_pdbfix.stateChanged.connect(replace_nonres.setEnabled)

		remove_heterogen = QCheckBox("Remove all heterogens including water", self)
		remove_heterogen.setEnabled(False)
		self.register_widget(remove_heterogen, 'check', 'PDBFixer/remove_heterogen', True, bool)
		pdbfix_grid.addWidget(remove_heterogen, 2, 1)
		use_pdbfix.stateChanged.connect(remove_heterogen.setEnabled)

		add_misheavy = QCheckBox("Add missing heavy atoms", self)
		add_misheavy.setEnabled(False)
		self.register_widget(add_misheavy, 'check', 'PDBFixer/add_misheavy', True, bool)
		pdbfix_grid.addWidget(add_misheavy, 3, 1)
		use_pdbfix.stateChanged.connect(add_misheavy.setEnabled)

		add_mishydrogen = QCheckBox("Add missing hydrogen atoms, the pH to use for adding hydrogens:", self)
		add_mishydrogen.setEnabled(False)
		self.register_widget(add_mishydrogen, 'check', 'PDBFixer/add_mishydrogen', True, bool)
		pdbfix_grid.addWidget(add_mishydrogen, 4, 1)
		mishydrogen_ph = QDoubleSpinBox(self)
		mishydrogen_ph.setEnabled(False)
		mishydrogen_ph.setRange(0, 14)
		mishydrogen_ph.setDecimals(1)
		self.register_widget(mishydrogen_ph, 'spin', 'PDBFixer/mishydrogen_ph', 7.0, float)
		pdbfix_grid.addWidget(mishydrogen_ph, 4, 2)
		use_pdbfix.stateChanged.connect(add_mishydrogen.setEnabled)
		use_pdbfix.stateChanged.connect(mishydrogen_ph.setEnabled)

	def create_pdbpqr_widgets(self):
		pdbpqr_grid = QGridLayout()
		pdbpqr_grid.setColumnStretch(1, 1)
		pdbpqr_grid.setColumnStretch(2, 1)
		self.main_layout.addLayout(pdbpqr_grid)
		pdbpqr_label = QLabel("<b>Use PDB2PQR to convert receptor PDB file to PQR file</b>", self)
		use_pdbpqr = QCheckBox(self)
		self.register_widget(use_pdbpqr, 'check', 'PDB2PQR/use_pdbpqr', False, bool)

		pdbpqr_grid.addWidget(use_pdbpqr, 0, 0)
		pdbpqr_grid.addWidget(pdbpqr_label, 0, 1)

		pdbpqr_grid.addWidget(QLabel("choose a forcefield to use:", self), 1, 1)
		force_field = QComboBox(self)
		force_field.addItems(['AMBER', 'CHARMM', 'PEOEPB', 'PARSE', 'SWANSON', 'TYL06'])
		force_field.setCurrentIndex(3)
		force_field.setEnabled(False)
		pdbpqr_grid.addWidget(force_field, 1, 2)
		self.register_widget(force_field, 'select', 'PDB2PQR/force_field', 'PARSE', str)
		use_pdbpqr.stateChanged.connect(force_field.setEnabled)

		use_propka = QCheckBox("Use PROPKA to assign protonation states at provided pH:", self)
		use_propka.setEnabled(False)
		use_pdbpqr.stateChanged.connect(use_propka.setEnabled)
		pdbpqr_grid.addWidget(use_propka, 2, 1)
		self.register_widget(use_propka, 'check', 'PDB2PQR/use_propka', True, bool)

		propka_ph = QDoubleSpinBox(self)
		propka_ph.setRange(0, 14)
		propka_ph.setDecimals(1)
		propka_ph.setEnabled(False)
		use_pdbpqr.stateChanged.connect(propka_ph.setEnabled)
		self.register_widget(propka_ph, 'spin', 'PDB2PQR/propka_ph', 7.0, float)
		pdbpqr_grid.addWidget(propka_ph, 2, 2)

		node_bump = QCheckBox("Ensure that new atoms are not rebuilt too close to existing atoms", self)
		node_bump.setEnabled(False)
		use_pdbpqr.stateChanged.connect(node_bump.setEnabled)
		pdbpqr_grid.addWidget(node_bump, 3, 1, 1, 2)
		self.register_widget(node_bump, 'check', 'PDB2PQR/node_bump', False, bool)

		no_hopt = QCheckBox("Optimize the hydrogen bonding network", self)
		no_hopt.setEnabled(False)
		use_pdbpqr.stateChanged.connect(no_hopt.setEnabled)
		pdbpqr_grid.addWidget(no_hopt, 4, 1)
		self.register_widget(no_hopt, 'check', 'PDB2PQR/no_hopt', False, bool)

		remove_water = QCheckBox("Remove the waters from the output file", self)
		remove_water.setEnabled(False)
		use_pdbpqr.stateChanged.connect(remove_water.setEnabled)
		pdbpqr_grid.addWidget(remove_water, 5, 1)
		self.register_widget(remove_water, 'check', 'PDB2PQR/remove_water', True, bool)

		neutraln = QCheckBox("Make the N-terminus of a protein neutral (only for PARSE)", self)
		neutraln.setEnabled(False)
		use_pdbpqr.stateChanged.connect(neutraln.setEnabled)
		pdbpqr_grid.addWidget(neutraln, 6, 1)
		self.register_widget(neutraln, 'check', 'PDB2PQR/neutraln', False, bool)

		neutralc = QCheckBox("Make the C-terminus of a protein neutral (only for PARSE)", self)
		neutralc.setEnabled(False)
		use_pdbpqr.stateChanged.connect(neutralc.setEnabled)
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
		self.settings.remove('PDBFixer')
		self.settings.remove('PDB2PQR')
		self.read_settings()

class ReceptorPreparationConfigPage(DockeyConfigPage):
	def create_config_inputs(self):
		b_h_radio = QRadioButton('Build bonds and add hydrogens', self)
		b_radio = QRadioButton('Build a single bond from each atom with no bonds to its closest neighbor', self)
		h_radio = QRadioButton('Add hydrogens', self)
		c_h_radio = QRadioButton('Add hydrogens only if there are none already', self)
		non_radio = QRadioButton('Do not make any repairs', self)

		self.register_widget(b_h_radio, 'radio', 'Receptor/repairs', 'bonds_hydrogens', str)
		self.register_widget(b_radio, 'radio', 'Receptor/repairs', 'bonds', str)
		self.register_widget(h_radio, 'radio', 'Receptor/repairs', 'hydrogens', str)
		self.register_widget(c_h_radio, 'radio', 'Receptor/repairs', 'checkhydrogens', str)
		self.register_widget(non_radio, 'radio', 'Receptor/repairs', '', str)

		repair_group = QButtonGroup(self)
		repair_group.addButton(b_h_radio)
		repair_group.addButton(b_radio)
		repair_group.addButton(h_radio)
		repair_group.addButton(c_h_radio)
		repair_group.addButton(non_radio)

		self.main_layout.addWidget(QLabel("<b>Types of repairs to make:</b>", self))
		repair_layout = QGridLayout()
		repair_layout.addWidget(non_radio, 0, 0)
		repair_layout.addWidget(b_h_radio, 0, 1)
		repair_layout.addWidget(h_radio, 1, 0)
		repair_layout.addWidget(c_h_radio, 1, 1)
		repair_layout.addWidget(b_radio, 2, 0, 1, 2)
		self.main_layout.addLayout(repair_layout)

		n_c_radio = QRadioButton("Preserve all input charges, do not add new charges", self)
		a_g_radio = QRadioButton("Add gasteiger charges", self)

		self.register_widget(n_c_radio, 'radio', 'Receptor/charges_to_add', None, str)
		self.register_widget(a_g_radio, 'radio', 'Receptor/charges_to_add', 'gasteiger', str)

		charge_group = QButtonGroup(self)
		charge_group.addButton(a_g_radio)
		charge_group.addButton(n_c_radio)

		charge_layout = QHBoxLayout()
		charge_layout.addWidget(a_g_radio)
		charge_layout.addWidget(n_c_radio)
		self.main_layout.addWidget(QLabel("<b>Charges:</b>", self))
		self.main_layout.addLayout(charge_layout)

		nphs_check = QCheckBox("Merge charges and remove non-polar hydrogens")
		lps_check = QCheckBox("Merge charges and remove lone pairs")
		water_check = QCheckBox("Remove water residues")
		nonstd_check = QCheckBox("Remove chains composed entirely of residues of types other than the standard 20 amino acids")
		delalt_check = QCheckBox("Remove XX@B atoms and rename XX@A atoms->XX")

		self.register_widget(nphs_check, 'check', 'Receptor/cleanup', 'nphs', str)
		self.register_widget(lps_check, 'check', 'Receptor/cleanup', 'lps', str)
		self.register_widget(water_check, 'check', 'Receptor/cleanup', 'waters', str)
		self.register_widget(nonstd_check, 'check', 'Receptor/cleanup', 'nonstdres', str)
		self.register_widget(delalt_check, 'check', 'Receptor/cleanup', 'deleteAltB', str)

		self.main_layout.addWidget(QLabel("<b>Clean types:</b>", self))
		self.main_layout.addWidget(nphs_check)
		self.main_layout.addWidget(lps_check)
		self.main_layout.addWidget(water_check)
		self.main_layout.addWidget(nonstd_check)
		self.main_layout.addWidget(delalt_check)

		self.main_layout.addWidget(QLabel("<b>Residues:</b>", self))
		e_check = QCheckBox("Delete every non-standard residue from any chain", self)
		self.main_layout.addWidget(e_check)
		self.main_layout.addWidget(QLabel("<font color='gray'>Any residue whose name is not in below list will be deleted from any chain</font>", self))
		self.main_layout.addWidget(QLabel("<font color='gray'><small>[CYS,ILE,SER,VAL,GLN,LYS,ASN,PRO,THR,PHE,ALA,HIS,GLY,ASP,LEU,ARG,TRP,GLU,TYR,MET,HID,HSP,HIE,HIP,CYX,CSS]<small></font>", self))
		self.register_widget(e_check, 'check', 'Receptor/delete_single_nonstd_residues', False, bool)

		self.main_layout.addWidget(QLabel("<b>Atoms:</b>", self))
		w_check = QCheckBox("assign each receptor atom a unique name: newname is original name plus its index(1-based)", self)
		self.main_layout.addWidget(w_check)
		self.register_widget(w_check, 'check', 'Receptor/unique_atom_names', False, bool)

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
		self.settings.remove('Receptor')
		self.read_settings()


class LigandPreparationConfigPage(DockeyConfigPage):
	def create_config_inputs(self):
		tool_select = QComboBox(self)
		tool_select.addItems(["prepare_ligand4", "meeko"])
		stack_widget = QStackedWidget(self)
		tool_select.currentIndexChanged.connect(stack_widget.setCurrentIndex)

		self.register_widget(tool_select, 'select', 'Ligand/prepare_tool', 'meeko', str)

		tool_layout = QHBoxLayout()
		tool_layout.addWidget(QLabel("<b>Select ligand preparation tool:</b>"))
		tool_layout.addWidget(tool_select)

		self.main_layout.addLayout(tool_layout)
		self.main_layout.addWidget(QLabel("<font color='gray'>prepare_ligand4 recommended for AutoDock4 and meeko recommended for AutoDock Vina</font>", self))
		self.main_layout.addWidget(stack_widget)

		#prepare_ligand4.py settings
		prelig_page = QWidget(self)
		prelig_layout = QVBoxLayout()
		prelig_page.setLayout(prelig_layout)
		stack_widget.addWidget(prelig_page)

		b_h_radio = QRadioButton('Build bonds and add hydrogens', self)
		b_radio = QRadioButton('Build a single bond from each atom with no bonds to its closest neighbor', self)
		h_radio = QRadioButton('Add hydrogens', self)
		c_h_radio = QRadioButton('Add hydrogens only if there are none already', self)
		non_radio = QRadioButton('Do not make any repairs', self)

		self.register_widget(b_h_radio, 'radio', 'Ligand/repairs', 'bonds_hydrogens', str)
		self.register_widget(b_radio, 'radio', 'Ligand/repairs', 'bonds', str)
		self.register_widget(h_radio, 'radio', 'Ligand/repairs', 'hydrogens', str)
		self.register_widget(c_h_radio, 'radio', 'Ligand/repairs', 'checkhydrogens', str)
		self.register_widget(non_radio, 'radio', 'Ligand/repairs', '', str)

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

		self.register_widget(n_c_radio, 'radio', 'Ligand/charges_to_add', None, str)
		self.register_widget(a_g_radio, 'radio', 'Ligand/charges_to_add', 'gasteiger', str)

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

		self.register_widget(nphs_check, 'check', 'Ligand/cleanup', 'nphs', str)
		self.register_widget(lps_check, 'check', 'Ligand/cleanup', 'lps', str)
	
		prelig_layout.addWidget(QLabel("<b>Clean types:</b>", self))
		prelig_layout.addWidget(nphs_check)
		prelig_layout.addWidget(lps_check)

		bb_check = QCheckBox("backbone", self)
		am_check = QCheckBox("amide", self)
		gd_check = QCheckBox("guaninium", self)

		self.register_widget(bb_check, 'check', 'Ligand/allowed_bonds', 'backbone', str)
		self.register_widget(am_check, 'check', 'Ligand/allowed_bonds', 'amide', str)
		self.register_widget(gd_check, 'check', 'Ligand/allowed_bonds', 'guaninium', str)

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

		self.register_widget(f_check, 'check', 'Ligand/check_for_fragments', False, bool)
		self.register_widget(z_check, 'check', 'Ligand/inactivate_all_torsions', False, bool)
		self.register_widget(g_check, 'check', 'Ligand/attach_nonbonded_fragments', False, bool)
		self.register_widget(s_check, 'check', 'Ligand/attach_singletons', False, bool)

		bound_layout = QGridLayout()
		bound_layout.addWidget(f_check, 0, 0)
		bound_layout.addWidget(z_check, 0, 1)
		bound_layout.addWidget(g_check, 1, 0)
		bound_layout.addWidget(s_check, 1, 1)
		prelig_layout.addLayout(bound_layout)

		#meeko settings
		meeko_page = QWidget(self)
		meeko_layout = QVBoxLayout()
		meeko_page.setLayout(meeko_layout)
		stack_widget.addWidget(meeko_page)

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
		self.register_widget(rmc_check, 'check', 'Meeko/rigid_macrocycles', False, bool)
		self.register_widget(kcr_check, 'check', 'Meeko/keep_chorded_rings', False, bool)
		self.register_widget(ker_check, 'check', 'Meeko/keep_equivalent_rings', False, bool)
		self.register_widget(hyd_check, 'check', 'Meeko/hydrate', False, bool)
		#self.register_widget(knh_check, 'check', 'Meeko/keep_nonpolar_hydrogens', False, bool)
		self.register_widget(far_check, 'check', 'Meeko/flexible_amides', False, bool)
		self.register_widget(aim_check, 'check', 'Meeko/add_index_map', False, bool)
		self.register_widget(rms_check, 'check', 'Meeko/remove_smiles', False, bool)

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

		self.register_widget(mat_input, 'edit', 'Meeko/merge_these_atom_types', 'H', str)
		self.register_widget(rbs_input, 'edit', 'Meeko/rigidify_bonds_smarts', '', str)
		self.register_widget(rbi_input, 'edit', 'Meeko/rigidify_bonds_indices', '', str)
		self.register_widget(ats_input, 'edit', 'Meeko/atom_type_smarts', '', str)
		self.register_widget(dbp_spin, 'spin', 'Meeko/double_bond_penalty', 50, int)

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
		self.settings.remove('Ligand')
		self.settings.remove('Meeko')
		self.read_settings()

class DockeyConfigDialog(QDialog):
	def __init__(self, parent):
		super().__init__(parent)
		self.setWindowTitle("Settings")

		self.tab_widget = QTabWidget(self)
		self.docking_tool_tab = DockingToolConfigPage(self)
		self.receptor_preprocess_tab = ReceptorPreprocessConfigPage(self)
		self.receptor_prepare_tab = ReceptorPreparationConfigPage(self)
		self.ligand_prepare_tab = LigandPreparationConfigPage(self)
		self.tab_widget.addTab(self.docking_tool_tab, "Docking tools")
		self.tab_widget.addTab(self.receptor_preprocess_tab, "Receptor preprocessing")
		self.tab_widget.addTab(self.receptor_prepare_tab, "Receptor preparation")
		self.tab_widget.addTab(self.ligand_prepare_tab, "Ligand preparation")

		self.btn_box = QDialogButtonBox(QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		self.btn_box.accepted.connect(self.save_settings)
		self.btn_box.rejected.connect(self.reject)
		self.btn_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.reset_settings)

		main_layout = QVBoxLayout()
		main_layout.addWidget(self.tab_widget)
		main_layout.addWidget(self.btn_box)
		self.setLayout(main_layout)

	def write_settings(self):
		self.docking_tool_tab.write_settings()
		self.receptor_prepare_tab.write_settings()
		self.receptor_preprocess_tab.write_settings()
		self.ligand_prepare_tab.write_settings()

	def save_settings(self):
		self.write_settings()
		self.accept()

	def reset_settings(self):
		#self.docking_tool_tab.reset_settings()
		self.receptor_prepare_tab.reset_settings()
		self.receptor_preprocess_tab.reset_settings()
		self.ligand_prepare_tab.reset_settings()

class DownloaderDialog(QDialog):
	title = ''
	mol_id = ''
	mol_type = 1
	mol_fmt = 'pdb'
	id_label = ''

	def __init__(self, parent=None):
		super().__init__(parent)

		self.setWindowTitle(self.title)
		self.layout = QVBoxLayout()
		#self.layout.setVerticalSpacing(10)
		self.setLayout(self.layout)

		self.add_logo()
		self.create_widgets()

		btn_box = QDialogButtonBox(
			QDialogButtonBox.Ok |
			QDialogButtonBox.Cancel
		)
		btn_box.accepted.connect(self.accept)
		btn_box.rejected.connect(self.reject)
		self.layout.addWidget(btn_box)

	def sizeHint(self):
		return QSize(500, 10)

	def create_widgets(self):
		self.layout.addWidget(QLabel(self.id_label, self))
		self.text_widget = QPlainTextEdit(self)
		self.layout.addWidget(self.text_widget)
		tip_label = QLabel("<font color='gray'>Use comma to separate multiple IDs</font>", self) 
		self.layout.addWidget(tip_label)

	def add_logo(self):
		logo_label = QLabel(self)
		logo_image = QPixmap(self.logo)
		logo_label.setPixmap(logo_image)
		self.layout.addWidget(logo_label)

	@classmethod
	def get_urls(cls, parent):
		dlg = cls(parent)

		mol_urls = []

		if dlg.exec() == QDialog.Accepted:
			text = dlg.text_widget.toPlainText().strip()

			if not text:
				QMessageBox.warning(parent, 'Warning', 
					"There are no input ids to download"
				)
				return

			entries = text.split(',')

			for entry in entries:
				entry = entry.strip()
				mol_urls.append((entry, cls.base_url.format(entry), cls.mol_fmt))

		return mol_urls

class PDBDownloadDialog(DownloaderDialog):
	title = "Import Receptors from PDB Database"
	logo = ':/icons/pdb.png'
	id_label = "PDB ID:"
	base_url = "https://files.rcsb.org/download/{}.pdb"

class ZINCDownloadDialog(DownloaderDialog):
	title = "Import Ligands from ZINC Database"
	logo = ":/icons/zinc.png"
	id_label = "ZINC IDs:"
	mol_type = 2
	mol_fmt = 'sdf'
	base_url = "https://zinc.docking.org/substances/{}.sdf"

class PubchemDownloadDialog(DownloaderDialog):
	title = "Import Ligand from PubChem Database"
	logo = ":/icons/pubchem.png"
	id_label = "PubChem CIDs:"
	mol_type = 2
	mol_fmt = 'sdf'
	base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/sdf"

class ChemblDownloadDialog(DownloaderDialog):
	title = "Import Ligand from ChEMBL Database"
	logo = ":/icons/chembl.png"
	id_label = "ChEMBL IDs:"
	mol_type = 2
	mol_fmt = 'sdf'
	base_url = "https://www.ebi.ac.uk/chembl/api/data/molecule/{}.sdf"

class FlexResiduesDialog(QDialog):
	def __init__(self, parent, mol_id):
		super().__init__(parent)
		self.mol_id = mol_id
		self.residue_checked = True
		self.bond_checked = True
		self.flex_residues = {}

		self.setWindowTitle("Select flex residues")

		self.residue_tree = QTreeWidget(self)
		self.residue_tree.setColumnCount(5)
		self.residue_tree.setRootIsDecorated(False)
		self.residue_tree.setHeaderLabels(["Index", "Chain", "AA", "Residue", "Atoms"])
		self.residue_tree.itemChanged.connect(self.on_residue_checked)
		self.residue_tree.itemClicked.connect(self.on_residue_clicked)

		self.bond_tree = QTreeWidget(self)
		self.bond_tree.setColumnCount(5)
		self.bond_tree.setRootIsDecorated(False)
		self.bond_tree.setVisible(False)
		self.bond_tree.setHeaderLabels(["Atoms", "Length", "Aromatic", "Rotatable", "Amide"])
		self.bond_tree.itemChanged.connect(self.on_bond_checked)

		self.bond_check = QCheckBox("Select bonds to disallowed", self)
		self.bond_check.stateChanged.connect(self.bond_tree.setVisible)

		self.btns = QDialogButtonBox(
			QDialogButtonBox.Cancel |
			QDialogButtonBox.Ok
		)
		self.btns.accepted.connect(self.on_accept_clicked)
		self.btns.rejected.connect(self.reject)

		layout = QVBoxLayout()
		layout.addWidget(QLabel("Set residues as flexible:", self))
		layout.addWidget(self.residue_tree, 2)
		layout.addWidget(self.bond_check)
		layout.addWidget(self.bond_tree, 1)
		layout.addWidget(self.btns)
		self.setLayout(layout)

		self.get_molecule_and_flexres()
		self.display_molecule_residues()

	def sizeHint(self):
		return QSize(550, 350)

	def get_molecule_and_flexres(self):
		sql = "SELECT content,format FROM molecular WHERE id=? LIMIT 1"
		self.mol = DB.get_dict(sql, (self.mol_id,))

		sql = "SELECT * FROM flex WHERE rid=?"
		for row in DB.query(sql, (self.mol_id,)):
			self.flex_residues[row[2]] = AttrDict({
				'checked': True,
				'data': list(row),
				'bonds': set(row[-1].split(','))
			})

	def display_molecule_residues(self):
		self.residue_checked = False
		for res in get_molecule_residues(self.mol.content, self.mol.format):
			item = QTreeWidgetItem(self.residue_tree, res)
			item.setFlags(item.flags() | Qt.ItemIsUserCheckable)

			idx = int(res[0])
			
			if idx in self.flex_residues:
				item.setCheckState(0, Qt.Checked)
			else:
				item.setCheckState(0, Qt.Unchecked)

			self.residue_tree.addTopLevelItem(item)
		self.residue_checked = True

	@Slot(QTreeWidgetItem, int)
	def on_residue_checked(self, item, column):
		if not self.residue_checked:
			return

		idx = int(item.text(0))
		chain = item.text(1)
		aa = item.text(2)
		res = item.text(3)
		if item.checkState(0) == Qt.Checked:
			if idx in self.flex_residues:
				self.flex_residues[idx].checked = True
				self.flex_residues[idx].bonds = set()
			else:
				self.flex_residues[idx] = AttrDict({
					'checked': True,
					'data': [None, self.mol_id, idx, chain, aa, res, ''],
					'bonds': set()
				})
		else:
			if idx in self.flex_residues:
				self.flex_residues[idx].checked = False

	@Slot(QTreeWidgetItem, int)
	def on_residue_clicked(self, item, column):
		self.bond_checked = False
		self.bond_tree.clear()
		idx = int(item.text(0))
		for bond in get_residue_bonds(self.mol.content, self.mol.format, idx):
			item = QTreeWidgetItem(self.bond_tree, bond)
			item.setFlags(item.flags() | Qt.ItemIsUserCheckable)

			if idx in self.flex_residues and bond[0] in self.flex_residues[idx].bonds:
				item.setCheckState(0, Qt.Checked)
			else:
				item.setCheckState(0, Qt.Unchecked)

			for i in range(2, 5):
				if bond[i]:
					item.setIcon(i, QIcon(':/icons/yes.svg'))
				else:
					item.setIcon(i, QIcon(':/icons/no.svg'))

			self.bond_tree.addTopLevelItem(item)
		self.bond_checked = True

	@Slot(QTreeWidgetItem, int)
	def on_bond_checked(self, item, column):
		if not self.bond_checked:
			return

		res_idx = int(self.residue_tree.currentItem().text(0))

		if res_idx not in self.flex_residues:
			return

		if not self.flex_residues[res_idx].checked:
			return

		if item.checkState(0) == Qt.Checked:
			self.flex_residues[res_idx].bonds.add(item.text(0))
		else:
			if item.text(0) in self.flex_residues[res_idx].bonds:
				self.flex_residues[res_idx].bonds.remove(item.text(0))

	@Slot()
	def on_accept_clicked(self):
		for idx, flexres in self.flex_residues.items():

			#if bond checkbox is not check, clear bonds
			if not self.bond_check.isChecked():
				flexres.bonds = set()

			if flexres.checked:
				flexres.data[-1] = ','.join(flexres.bonds)

				if flexres.data[0] is None:
					sql = "INSERT INTO flex VALUES (?,?,?,?,?,?,?)"
					DB.query(sql, flexres.data)
				else:
					sql = "UPDATE flex SET bonds=? WHERE id=?"
					DB.query(sql, (flexres.data[-1], flexres.data[0]))

			else:
				if flexres.data[0] is not None:
					sql = "DELETE FROM flex WHERE id=?"
					DB.query(sql, (flexres.data[0],))

		self.accept()

class LigandFilterDialog(QDialog):
	def __init__(self, parent):
		super().__init__(parent)
		self.setWindowTitle("Filter and Remove Ligands")
		self.weight = QDoubleSpinBox(self)
		self.weight.setRange(0, 100000)
		self.weight.setDecimals(3)
		self.mwsign = QComboBox(self)
		self.mwsign.addItems(['>', '>=', '=', '<=', '<'])
		self.mcheck = QCheckBox("Molecular weight", self)

		self.rotator = QSpinBox(self)
		self.rotator.setRange(0, 100000)
		self.rbsign = QComboBox(self)
		self.rbsign.addItems(['>', '>=', '=', '<=', '<'])
		self.rcheck = QCheckBox("Rotatable bonds", self)
		self.rbjoin = QComboBox(self)
		self.rbjoin.addItems(['AND', 'OR'])

		self.pcheck = QCheckBox("Calculated logP", self)
		self.clogp = QDoubleSpinBox(self)
		self.clogp.setRange(0, 100000)
		self.clogp.setDecimals(3)
		self.cpsign = QComboBox(self)
		self.cpsign.addItems(['>', '>=', '=', '<=', '<'])
		self.cpjoin = QComboBox(self)
		self.cpjoin.addItems(['AND', 'OR'])

		self.btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		self.btns.accepted.connect(self.accept)
		self.btns.rejected.connect(self.reject)

		grid_layout = QGridLayout()
		grid_layout.setColumnStretch(2, 1)
		grid_layout.setColumnStretch(3, 2)
		grid_layout.addWidget(self.mcheck, 0, 1)
		grid_layout.addWidget(self.mwsign, 0, 2)
		grid_layout.addWidget(self.weight, 0, 3)
		grid_layout.addWidget(self.rbjoin, 1, 0)
		grid_layout.addWidget(self.rcheck, 1, 1)
		grid_layout.addWidget(self.rbsign, 1, 2)
		grid_layout.addWidget(self.rotator, 1, 3)
		grid_layout.addWidget(self.cpjoin, 2, 0)
		grid_layout.addWidget(self.pcheck, 2, 1)
		grid_layout.addWidget(self.cpsign, 2, 2)
		grid_layout.addWidget(self.clogp, 2, 3)

		self.layout = QVBoxLayout()
		self.setLayout(self.layout)
		self.layout.addWidget(QLabel("<b>Remove ligands that match these filters:</b>", self))
		self.layout.addLayout(grid_layout)
		self.layout.addSpacing(20)
		self.layout.addWidget(self.btns)

	def sizeHint(self):
		return QSize(450, 100)

	@classmethod
	def filter(cls, parent):
		dlg = cls(parent)

		if dlg.exec() == QDialog.Accepted:
			mwt = dlg.weight.value()
			mws = dlg.mwsign.currentText()

			rbn = dlg.rotator.value()
			rbs = dlg.rbsign.currentText()
			rbj = dlg.rbjoin.currentText()

			cpv = dlg.clogp.value()
			cps = dlg.cpsign.currentText()
			cpj = dlg.cpjoin.currentText()

			conditions = []

			if dlg.mcheck.isChecked():
				conditions.append(" weight {} {} ".format(mws, mwt))

			if dlg.rcheck.isChecked():
				if conditions:
					conditions.append(rbj)

				conditions.append(" rotors {} {} ".format(rbs, rbn))

			if dlg.pcheck.isChecked():
				if conditions:
					conditions.append(cpj)

				conditions.append(" logp {} {} ".format(cps, cpv))

			return ''.join(conditions)


class AcknowledgementDialog(QDialog):
	def __init__(self, parent, text):
		super().__init__(parent)
		self.setWindowTitle("Acknowledgements")

		layout = QVBoxLayout()
		self.setLayout(layout)

		info = QLabel(
			"Dockey has integrated several external useful tools "
			"and denpends on preinstallation of docking engines"
		, self)
		layout.addWidget(info)

		thanks = QTextBrowser(self)
		thanks.setOpenLinks(True)
		thanks.setOpenExternalLinks(True)
		thanks.setHtml(text)
		layout.addWidget(thanks)

		btnbox = QDialogButtonBox(QDialogButtonBox.Ok)
		btnbox.accepted.connect(self.accept)
		layout.addWidget(btnbox)

	def sizeHint(self):
		return QSize(630, 350)

	@classmethod
	def thank(cls, parent, text):
		dlg = cls(parent, text)
		dlg.exec()

class CPUAndMemoryViewDialog(QDialog):
	def __init__(self, parent):
		super().__init__(parent)
		self.setWindowTitle("Computing resource usage")
		self.proc = psutil.Process(os.getpid())

		main_layout = QFormLayout()
		main_layout.setVerticalSpacing(20)
		main_layout.addRow(QLabel("CPU and memory used by Dockey", self))

		self.cpu_usage = QProgressBar(self)
		self.memory_usage = QProgressBar(self)
		self.memory_text = QLabel(self)
		main_layout.addRow('CPU:', self.cpu_usage)
		main_layout.addRow('Memory:', self.memory_usage)
		main_layout.addRow('', self.memory_text)

		btnbox = QDialogButtonBox(QDialogButtonBox.Ok)
		btnbox.accepted.connect(self.accept)
		main_layout.addRow(btnbox)

		self.setLayout(main_layout)

		self.timer = QTimer(self)
		self.timer.timeout.connect(self.update_resource_usage)
		self.timer.start(1000)
		self.accepted.connect(self.timer.stop)
		self.rejected.connect(self.timer.stop)

	def sizeHint(self):
		return QSize(400, 10)

	@Slot()
	def update_resource_usage(self):
		process_count = 1
		cpu_percent = 0
		memory_size = 0
		memory_percent = 0

		try:
			with self.proc.oneshot():
				cpu_percent = self.proc.cpu_percent(interval=0.1)
				memory_size = self.proc.memory_info().rss
				memory_percent = self.proc.memory_percent()
		except psutil.NoSuchProcess:
			pass

		try:
			children = self.proc.children(recursive=True)
		except psutil.NoSuchProcess:
			children = []

		for child in children:
			try:
				with child.oneshot():
					process_count += 1
					cpu_percent += child.cpu_percent(interval=0.1)
					memory_size += child.memory_info().rss
					memory_percent += child.memory_percent()
			except psutil.NoSuchProcess:
				pass

		self.cpu_usage.setValue(int(cpu_percent/psutil.cpu_count()))
		self.memory_usage.setValue(int(memory_percent))
		self.memory_text.setText("Memory: {} \t Running processes: {}".format(
			memory_format(memory_size), process_count
		))

class DockeyRunTaskDialog(QDialog):
	title = "Run Tasks"

	def __init__(self, parent=None):
		super().__init__(parent)

		self.setWindowTitle(self.title)

		self.main_layout = QVBoxLayout()
		self.setLayout(self.main_layout)

		self.widget_layout = QGridLayout()
		self.widget_layout.setColumnStretch(1, 1)
		self.widget_rownum = -1
		self.main_layout.addLayout(self.widget_layout)

		self.create_customs()
		self.create_widgets()

		button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		button_box.accepted.connect(self.accept)
		button_box.rejected.connect(self.reject)
		self.main_layout.addWidget(button_box)

	def create_customs(self):
		pass

	def create_widgets(self):
		self.receptor_tool = QComboBox(self)
		self.receptor_tool.addItems(['AutoDockTools', 'OpenBabel'])

		self.ligand_tool = QComboBox(self)
		self.ligand_tool.addItems(['Meeko', 'AutoDockTools', 'OpenBabel'])

		self.widget_rownum += 1
		self.widget_layout.addWidget(QLabel("Select receptor preparation tool:"), self.widget_rownum, 0)
		self.widget_layout.addWidget(self.receptor_tool, self.widget_rownum, 1)

		self.widget_rownum += 1
		self.widget_layout.addWidget(QLabel("Select ligand preparation tool:"), self.widget_rownum, 0)
		self.widget_layout.addWidget(self.ligand_tool, self.widget_rownum, 1)

		self.use_pdbfix = QCheckBox("Use PDBFixer to fix PDB receptor file")
		self.use_pdbpqr = QCheckBox("Use PDB2PQR to convert PDB to PQR")
		advance_layout = QVBoxLayout()
		advance_layout.addWidget(self.use_pdbfix)
		advance_layout.addWidget(self.use_pdbpqr)

		self.advance_group = QGroupBox("Advanced")
		self.advance_group.setFlat(True)
		self.advance_group.toggled.connect(self.use_pdbfix.setVisible)

		self.advance_group.toggled.connect(self.use_pdbpqr.setVisible)
		self.advance_group.setCheckable(True)
		self.advance_group.setChecked(False)
		self.advance_group.setLayout(advance_layout)
		self.widget_rownum += 1
		self.main_layout.addWidget(self.advance_group)

	def get_parameters(self):
		pass

	@classmethod
	def start(cls, parent=None):
		dlg = cls(parent)

		if dlg.exec() == QDialog.Accepted:
			return self.get_parameters()

class DockeyRunAutodockDialog(DockeyRunTaskDialog):
	title = "Run AutoDock4"

	def create_customs(self):
		self.algorithm_select = QComboBox(self)
		self.algorithm_select.addItems(["Lamarckian GA", "Genetic Algorithm",
			"Simulated Annealing", "Local Search"])
		self.widget_rownum += 1
		self.widget_layout.addWidget(QLabel("Select algorithm for AutoDock4:"), self.widget_rownum, 0)
		self.widget_layout.addWidget(self.algorithm_select, self.widget_rownum, 1)

	def get_parameters(self):
		pass

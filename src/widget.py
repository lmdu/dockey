import os
import psutil

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtNetwork import *
from PySide6.QtWidgets import *
from plip.basic.remote import VisualizerData
from plip.structure.preparation import PDBComplex

from utils import *
from table import *
from backend import *

__all__ = ['BrowseInput', 'CreateProjectDialog', 'AutodockConfigDialog',
			'AutodockVinaConfigDialog', 'ExportImageDialog', 'FeedbackBrowser',
			'PymolSettingDialog', 'DockingToolSettingDialog', 'InteractionTabWidget',
			'JobConcurrentSettingDialog', 'DockeyConfigDialog', 'PDBDownloader',
			'ZINCDownloader'
			]

class BrowseInput(QWidget):
	def __init__(self, parent=None, is_file=True, is_save=False, _filter=""):
		super(BrowseInput, self).__init__(parent)
		self.filter = _filter
		self.input = QLineEdit(self)
		self.input.setReadOnly(True)
		self.browse = QPushButton(self)
		self.browse.setFlat(True)
		self.browse.setIcon(QIcon(":/icons/folder.svg"))

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
		super(CreateProjectDialog, self).__init__(parent)
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
		super(DockingToolSettingDialog, self).__init__(parent)
		self.resize(QSize(550, 100))
		self.settings = QSettings()
		self.ad4 = self.settings.value('Tools/autodock_4', '')
		self.ag4 = self.settings.value('Tools/autogrid_4', '')
		self.vina = self.settings.value('Tools/autodock_vina', '')
		#self.mglt = self.settings.value('Tools/mgltools', '')
		self.autodock_input = BrowseInput(self)
		self.autodock_input.set_text(self.ad4)
		self.autogrid_input = BrowseInput(self)
		self.autogrid_input.set_text(self.ag4)
		self.vina_input = BrowseInput(self)
		self.vina_input.set_text(self.vina)
		#self.mgltools_input = BrowseInput(self, False)
		#self.mgltools_input.set_text(self.mglt)

		self.layout = QVBoxLayout()
		self.layout.addWidget(QLabel("Autodock4 executable", self))
		self.layout.addWidget(self.autodock_input)
		self.layout.addWidget(QLabel("Autogrid4 executable", self))
		self.layout.addWidget(self.autogrid_input)
		self.layout.addWidget(QLabel("Autodock vina executable", self))
		self.layout.addWidget(self.vina_input)
		#self.layout.addWidget(QLabel("MGLTools directory", self))
		#self.layout.addWidget(self.mgltools_input)

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
		#mglt = self.mgltools_input.get_text()

		if ad4 != self.ad4:
			self.settings.setValue('Tools/autodock_4', ad4)

		if ag4 != self.ag4:
			self.settings.setValue('Tools/autogrid_4', ag4)

		if vina != self.vina:
			self.settings.setValue('Tools/autodock_vina', vina)

		#if mglt != self.mglt:
		#	self.settings.setValue('Tools/mgltools', mglt)

		self.accept()

class ProgramConfigDialog(QDialog):
	title = ""
	progs = []

	def __init__(self, parent=None):
		super(ProgramConfigDialog, self).__init__(parent)
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

class AutodockVinaConfigDialog(ProgramConfigDialog):
	title = "Autodock Vina Config"
	progs = [
		('Tools/autodock_vina', 'Autodock Vina executable'),
		#('Tools/mgltools', 'MGLTools directory')
	]

class ExportImageDialog(QDialog):
	def __init__(self, parent=None):
		super(ExportImageDialog, self).__init__(parent)
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
		super(FeedbackBrowser, self).__init__(parent)
		self.setReadOnly(True)

	def sizeHint(self):
		return QSize(300, 50)

class GradientColorBar(QWidget):
	def __init__(self, parent):
		super(GradientColorBar, self).__init__(parent)

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
		super(PymolSettingDialog, self).__init__(parent)
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

class JobConcurrentSettingDialog(QDialog):
	def __init__(self, parent=None):
		super(JobConcurrentSettingDialog, self).__init__(parent)
		self.parent = parent
		self.settings = QSettings()
		self.setWindowTitle("Concurrent Job Setting")
		self.layout = QVBoxLayout()
		self.label = QLabel("Number of concurrent running jobs", self)
		self.number = QSpinBox(self)
		self.number.setValue(self.settings.value('Job/concurrent', 1, int))
		self.number.setRange(1, max(1, psutil.cpu_count()-2))
		self.number.valueChanged.connect(self.on_job_number_changed)
		self.layout.addWidget(self.label)
		self.layout.addWidget(self.number)
		self.setLayout(self.layout)

	@Slot()
	def on_job_number_changed(self, num):
		self.parent.pool.setMaxThreadCount(num)
		self.settings.setValue('Job/concurrent', num)

class InteractionTabWidget(QTabWidget):
	def __init__(self, parent=None):
		super(InteractionTabWidget, self).__init__(parent)
		self.parent = parent
		self.pose_id = 0
		self.binding_site = None
		self.complex_pdb = None
		self.complex_vis = None

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

	@Slot()
	def on_binding_site_changed(self, site):
		sql = "SELECT id FROM binding_site WHERE pid=? AND site=? LIMIT 1"
		site_id = DB.get_one(sql, (self.pose_id, site))

		for i in range(8):
			self.widget(i).model().change_binding_site(site_id)

		if site != self.binding_site:
			self.binding_site = site
			mol = PDBComplex()
			mol.load_pdb(self.complex_pdb, as_string=True)

			for ligand in mol.ligands:
				mol.characterize_complex(ligand)

			self.complex_vis = VisualizerData(mol, self.binding_site)

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


class DockeyConfigDialog(QDialog):
	def __init__(self, parent):
		super(DockeyConfigDialog, self).__init__(parent)
		self.parent = parent
		self.setWindowTitle("Settings")
		self.resize(QSize(550, 100))
		self.input_widgets = []
		self.settings = QSettings()
		self.tab_widget = QTabWidget(self)

		self.btn_widget = QDialogButtonBox(
			QDialogButtonBox.Ok | QDialogButtonBox.Cancel
		)

		self.btn_widget.accepted.connect(self.write_settings)
		self.btn_widget.rejected.connect(self.reject)

		main_layout = QVBoxLayout()
		main_layout.addWidget(self.tab_widget)
		main_layout.addWidget(self.btn_widget)
		self.setLayout(main_layout)

		self.create_tool_tab()
		self.create_job_tab()
		self.read_settings()

	def register_widget(self, option, widget, wgtype, default='', convert=str):
		self.input_widgets.append(AttrDict(
			widget = widget,
			option = option,
			wgtype = wgtype,
			default = default,
			convert = convert
		))

	def create_tool_tab(self):
		page = QWidget(self)
		layout = QVBoxLayout()
		page.setLayout(layout)

		self.autodock_input = BrowseInput(self)
		self.autogrid_input = BrowseInput(self)
		self.vina_input = BrowseInput(self)

		self.register_widget('Tools/autodock_4', self.autodock_input, 'browse')
		self.register_widget('Tools/autogrid_4', self.autogrid_input, 'browse')
		self.register_widget('Tools/autodock_vina', self.vina_input, 'browse')
                          
		layout.addWidget(QLabel("Autogrid executable", self))
		layout.addWidget(self.autogrid_input)
		layout.addWidget(QLabel("Autodock executable", self))
		layout.addWidget(self.autodock_input)
		layout.addWidget(QLabel("Autodock Vina executable", self))
		layout.addWidget(self.vina_input)

		self.tab_widget.addTab(page, 'Docking Tools')

	def create_job_tab(self):
		page = QWidget(self)
		layout = QGridLayout()
		layout.setColumnStretch(0, 5)
		layout.setColumnStretch(1, 1)
		page.setLayout(layout)
		label = QLabel("Number of concurrent running jobs", self)
		layout.addWidget(label, 0, 0)
		self.thread_input = QSpinBox(self)
		self.thread_input.setRange(1, max(1, psutil.cpu_count()))
		self.register_widget('Job/concurrent', self.thread_input, 'spinbox', 1, int)

		tips = (
			"You can set the number of jobs that can run concurrently to improve the docking speed."
			"Limiting the number of concurrent jobs helps you reduce CPU consumption."
			"You can set the number of concurrent jobs to any number from 1 to maximum CPU threads."
		)
		tips_label = QLabel(tips, self)
		tips_label.setWordWrap(True)

		layout.addWidget(self.thread_input, 0, 1)
		layout.addWidget(tips_label, 1, 0)

		self.tab_widget.addTab(page, 'Job Manager')

	def read_settings(self):
		for i in self.input_widgets:
			val = self.settings.value(i.option, i.default, i.convert)

			if i.wgtype == 'browse':
				i.widget.set_text(val)
			elif i.wgtype == 'linedit':
				i.widget.setText(val)
			elif i.wgtype == 'spinbox':
				i.widget.setValue(val)

	def write_settings(self):
		self.change_concurrent_jobs()

		for i in self.input_widgets:
			if i.wgtype == 'browse':
				val = i.widget.get_text()
			elif i.wgtype == 'linedit':
				val = i.widget.text()
			elif i.wgtype == 'spinbox':
				val = i.widget.value()
			else:
				val = None

			self.settings.setValue(i.option, val)

		self.accept()

	def change_concurrent_jobs(self):
		val = self.settings.value('Job/concurrent', 1, int)
		num = self.thread_input.value()

		if num == val:
			return

		self.parent.pool.setMaxThreadCount(num)

class DownloaderDialog(QDialog):
	title = ''
	mol_id = ''
	mol_type = 1
	mol_fmt = 'pdb'

	def __init__(self, parent=None):
		super(DownloaderDialog, self).__init__(parent)
		self.parent = parent
		self.reply = None
		self.resize(QSize(350, 10))
		self.setWindowTitle(self.title)
		self.layout = QFormLayout()
		self.setLayout(self.layout)

		self.manager = QNetworkAccessManager(self.parent)

		self.create_widgets()

		self.progress_bar = QProgressBar(self)
		self.progress_bar.setAlignment(Qt.AlignCenter)
		self.layout.addRow("Progress", self.progress_bar)

		btn_box = QDialogButtonBox()
		self.start_btn = btn_box.addButton("Import", QDialogButtonBox.AcceptRole)
		self.abort_btn = btn_box.addButton("Cancel", QDialogButtonBox.RejectRole)
		self.layout.addRow(btn_box)

		self.start_btn.clicked.connect(self.on_start)
		self.abort_btn.clicked.connect(self.on_abort)

	def create_widgets(self):
		pass

	def get_url(self):
		pass

	@Slot()
	def on_start(self):
		url = self.get_url()

		if url:
			self.start_btn.setDisabled(True)
			self.reply = self.manager.get(QNetworkRequest(url))
			self.reply.downloadProgress.connect(self.on_progress)
			self.reply.finished.connect(self.on_finished)
			self.reply.errorOccurred.connect(self.on_error)

	@Slot()
	def on_abort(self):
		if self.reply:
			self.reply.abort()
			self.progress_bar.setValue(0)
		else:
			self.reject()

		self.start_btn.setDisabled(False)

	@Slot()
	def on_finished(self):
		if self.reply:
			if self.reply.error() == QNetworkReply.NoError:
				content = str(self.reply.readAll(), encoding='utf-8')

				self.save_molecular(content)

			self.reply.deleteLater()

		self.start_btn.setDisabled(False)
		self.accept()

	@Slot(int, int)
	def on_progress(self, received_bytes, total_bytes):
		self.progress_bar.setRange(0, total_bytes)
		self.progress_bar.setValue(received_bytes)

	@Slot(QNetworkReply.NetworkError)
	def on_error(self, code):
		self.reject()

		if self.reply:
			QMessageBox.critical(self.parent, "Error", self.reply.errorString())

	def save_molecular(self, mol):
		sql = "INSERT INTO molecular VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
		mi = get_molecule_information(mol, True, self.mol_id, self.mol_fmt)
		row = [None, mi.name, self.mol_type, mi.pdb, mi.atoms, mi.bonds,
			mi.hvyatoms, mi.residues, mi.rotors, mi.formula, mi.energy,
			mi.weight, mi.logp
		]
		DB.query(sql, row)
		mol_type = ['', 'receptor', 'ligand'][self.mol_type]
		self.parent.show_message("Import {} {}".format(mol_type, row[1]))
		self.parent.mol_model.select()

class PDBDownloader(DownloaderDialog):
	title = "Import Receptor from PDB"

	def create_widgets(self):
		self.text_widget = QLineEdit(self)
		self.layout.addRow("PDB ID:", self.text_widget)

	def get_url(self):
		base_url = "https://files.rcsb.org/download/{}.pdb"

		pdb_id = self.text_widget.text().strip().lower()

		if not pdb_id:
			return None

		self.mol_id = pdb_id

		return QUrl(base_url.format(pdb_id))

class ZINCDownloader(DownloaderDialog):
	title = "Import Ligand from ZINC"
	mol_type = 2
	mol_fmt = 'sdf'

	def create_widgets(self):
		self.text_widget = QLineEdit(self)
		self.version_widget = QComboBox(self)
		self.version_widget.addItems(['zinc15', 'zinc20'])
		self.layout.addRow("ZINC ID:", self.text_widget)
		self.layout.addRow("ZINC Version:", self.version_widget)

	def get_url(self):
		zinc_id = self.text_widget.text().strip().upper()

		if not zinc_id:
			return None

		if len(zinc_id) != 16:
			zinc_id = zinc_id.strip('ZINC')
			zinc_id = "ZINC{:0>12}".format(zinc_id)

		self.mol_id = zinc_id

		ver = self.version_widget.currentText()

		if ver == 'zinc20':
			base_url = "https://zinc20.docking.org/substances/{}.sdf"
		else:
			base_url = "https://zinc.docking.org/substances/{}.sdf"

		return QUrl(base_url.format(zinc_id))

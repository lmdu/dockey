import os

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from param import *
from config import *
from pymolview import *
from gridbox import *

__all__ = ['DockeyMainWindow']

class DockeyMainWindow(QMainWindow):
	def __init__(self):
		super(DockeyMainWindow, self).__init__()

		self.setWindowTitle("Dockey v{}".format(DOCKEY_VERSION))

		self.pymol_viewer = PymolGLWidget(self)
		self.setCentralWidget(self.pymol_viewer)
		self.cmd = self.pymol_viewer.cmd

		#grid box dock widget
		
		self.box_dock = QDockWidget("Gridbox setting panel", self)
		self.box_panel = GridBoxSettingPanel(self)
		GBD.updated.connect(self.box_panel.set_value)
		self.box_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.box_dock.setWidget(self.box_panel)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.box_dock)
		self.box_dock.setVisible(False)

		self.create_actions()
		self.create_menus()
		self.create_toolbar()
		self.create_statusbar()

		self.project_folder = None

	def create_actions(self):
		self.new_project_act = QAction("&New project", self,
			shortcut = QKeySequence.New,
			triggered = self.new_project
		)

		self.open_project_act = QAction("&Open project", self,
			shortcut = QKeySequence.Open,
			triggered = self.open_project
		)

		self.open_project_dir_act = QAction("&Show project in explorer", self,
			triggered = self.open_project_dir
		)

		self.close_project_act = QAction("&Close project", self,
			shortcut = QKeySequence.Close,
			triggered = self.close_project
		)

		self.import_act = QAction("&Import ligand", self,
			triggered = self.import_ligand
		)

		#view actions
		#pymol sidebar
		self.pymol_sidebar_act = QAction("Pymol control sidebar", self,
			checkable = True,
			checked = False
		)
		self.pymol_sidebar_act.toggled.connect(self.pymol_sidebar_toggle)

		#gridbox setting sidebar
		self.box_sidebar_act = self.box_dock.toggleViewAction()
		self.box_sidebar_act.setChecked(False)

		self.bounding_box_act = QAction("Bounding box", self)
		self.bounding_box_act.triggered.connect(self.pymol_viewer.draw_bounding_box)

		self.custom_box_act = QAction("Custom box", self)
		self.custom_box_act.triggered.connect(self.pymol_viewer.draw_custom_box)

		self.run_autodock_act = QAction("Autodock", self)
		self.run_autodock_act.triggered.connect(self.run_autodock)

	def create_menus(self):
		self.file_menu = self.menuBar().addMenu("&File")
		self.file_menu.addAction(self.new_project_act)
		self.file_menu.addAction(self.open_project_act)
		self.file_menu.addAction(self.close_project_act)
		self.file_menu.addSeparator()
		self.file_menu.addAction(self.import_act)

		self.view_menu = self.menuBar().addMenu("&View")
		self.view_menu.addAction(self.open_project_dir_act)
		self.view_menu.addSeparator()
		self.view_menu.addAction(self.box_sidebar_act)
		self.view_menu.addAction(self.pymol_sidebar_act)

		self.receptor_menu = self.menuBar().addMenu("&Receptor")
		self.ligand_menu = self.menuBar().addMenu("&Ligand")

		self.grid_menu = self.menuBar().addMenu("&Grid")
		self.grid_menu.addAction(self.bounding_box_act)
		self.grid_menu.addAction(self.custom_box_act)

		self.run_menu = self.menuBar().addMenu("&Run")
		self.run_menu.addAction(self.run_autodock_act)

	def create_toolbar(self):
		self.toolbar = self.addToolBar('')
		self.toolbar.setMovable(False)

		self.toolbar.addAction(self.import_act)

	def create_statusbar(self):
		self.statusbar = self.statusBar()
		self.statusbar.showMessage("Welcome to Dockey")

	def new_project(self):
		self.project_folder = CreateProjectDialog.get_project_folder(self)
		os.mkdir(self.project_folder)

		#make data subfolder
		#os.mkdir(os.path.join(self.project_folder, 'data'))

		#make param subfolder
		#os.mkdir(os.path.join(self.project_folder, 'param'))

		#make temp subfolder

	def open_project(self):
		folder = QFileDialog.getExistingDirectory(self)
		if folder:
			self.project_folder = folder

	def open_project_dir(self):
		QDesktopServices.openUrl(QUrl.fromLocalFile(self.project_folder))

	def close_project(self):
		pass

	def import_ligand(self):
		pdb_file, _ = QFileDialog.getOpenFileName(self,
			caption = "Select pdb file",
			filter = "PDB file (*.pdb)"
		)

		if not pdb_file:
			return

		file_info = QFileInfo(pdb_file)
		name = file_info.baseName()

		with open(pdb_file) as fh:
			content = fh.read()

		self.pymol_viewer.load_receptor(name, content)

	def pymol_sidebar_toggle(self, checked):
		self.pymol_viewer.sidebar_controler(checked)

	def run_autodock(self):
		dlg = AutodockParameterDialog(self)
		dlg.exec()

class CreateProjectDialog(QDialog):
	def __init__(self, parent=None):
		super(CreateProjectDialog, self).__init__(parent)
		self.parent = parent
		self.setWindowTitle("Create a New Project")

		self.name_input = QLineEdit(self)
		self.location_input = QLineEdit(self)
		self.location_input.setReadOnly(True)
		self.browse_btn = QPushButton(self)
		self.browse_btn.setFlat(True)
		self.browse_btn.clicked.connect(self.select_location)
		self.browse_btn.setIcon(QIcon("icons/folder.svg"))
		self.tip_label = QLabel('<font color="gray">This will create a new folder with input name in the location you selected</font>', self)

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.check_input_name)
		action_box.rejected.connect(self.reject)

		layout = QGridLayout()
		layout.addWidget(QLabel("Name"), 0, 0)
		layout.addWidget(self.name_input, 0, 1)
		layout.addWidget(QLabel("Loacation"), 1, 0)
		layout.addWidget(self.location_input, 1, 1)
		layout.addWidget(self.browse_btn, 1, 2)
		layout.addWidget(self.tip_label, 2, 0, 1, 3)
		layout.addWidget(action_box, 3, 1, 1, 2)

		self.setLayout(layout)

	@Slot()
	def check_input_name(self):
		name = self.name_input.text()
		folder = self.location_input.text()
		if not name:
			return QMessageBox.warning(self, "Input Warning",
				"Please input a project name"
			)

		if not folder:
			return QMessageBox.warning(self, "Input Warning",
				"Please select a directory to hold project"
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

	@Slot()
	def select_location(self):
		folder = QFileDialog.getExistingDirectory(self)

		if folder:
			self.location_input.setText(folder)

	@classmethod
	def get_project_folder(cls, parent):
		dlg = cls(parent)
		if dlg.exec():
			name = dlg.name_input.text()
			location = dlg.location_input.text()
			return os.path.join(location, name)

import os

from PySide6.QtGui import *
from PySide6.QtSql import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from db import *
from param import *
from config import *
from pymolview import *
from gridbox import *

__all__ = ['DockeyMainWindow']

class DockeyMainWindow(QMainWindow):
	def __init__(self):
		super(DockeyMainWindow, self).__init__()

		self.setWindowTitle("Dockey v{}".format(DOCKEY_VERSION))
		self.setWindowIcon(QIcon('icons/logo.svg'))
		self.resize(QSize(900, 600))

		self.db = DockeyDatabase()

		self.create_pymol_viewer()
		self.create_molecular_viewer()
		self.create_gridbox_adjuster()
		self.create_task_table()

		self.create_actions()
		self.create_menus()
		self.create_toolbar()
		self.create_statusbar()

		self.project_folder = None

	def create_pymol_viewer(self):
		self.pymol_viewer = PymolGLWidget(self)
		self.cmd = self.pymol_viewer.cmd
		self.setCentralWidget(self.pymol_viewer)

	def create_molecular_viewer(self):
		#self.mol_viewer = QTreeWidget(self)
		#self.mol_viewer.setColumnCount(1)
		#self.mol_viewer.header().hide()
		self.mol_viewer = QListView(self)
		self.mol_viewer.clicked.connect(self.on_molecular_changed)
		self.mol_dock = QDockWidget("Moleculars", self)
		self.mol_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.mol_dock.setWidget(self.mol_viewer)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.mol_dock)
		self.mol_dock.setVisible(True)

	@Slot()
	def on_molecular_changed(self, index):
		self.cmd.delete('all')
		self.cmd.reinitialize()

		name = index.data(Qt.DisplayRole)
		content = self.db.get_content(name)
		
		if content:
			self.cmd.read_pdbstr(content, name)

	def create_gridbox_adjuster(self):
		self.box_adjuster = GridBoxSettingPanel(self)
		GBD.updated.connect(self.box_adjuster.set_value)

		self.box_dock = QDockWidget("Gridbox", self)
		self.box_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.box_dock.setWidget(self.box_adjuster)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.box_dock)
		self.box_dock.setVisible(False)

	def create_task_table(self):
		self.task_table = QTableView(self)
		self.task_dock = QDockWidget("Docks", self)
		self.task_dock.setWidget(self.task_table)
		self.task_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.task_dock)
		self.task_dock.setVisible(True)

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

		self.import_receptor_act = QAction("&Import Receptors", self,
			triggered = self.import_receptors
		)

		self.import_ligand_act = QAction("&Import Ligands", self,
			triggered = self.import_ligands
		)

		#view actions
		#pymol sidebar
		self.pymol_sidebar_act = QAction("Show pymol sidebar", self,
			checkable = True,
			checked = False
		)
		self.pymol_sidebar_act.toggled.connect(self.pymol_sidebar_toggle)

		#gridbox setting sidebar
		self.box_sidebar_act = self.box_dock.toggleViewAction()
		self.box_sidebar_act.setText("Show gridbox")
		self.box_sidebar_act.setChecked(False)

		self.mol_list_act = self.mol_dock.toggleViewAction()
		self.mol_list_act.setText("Show molecular list")

		self.task_table_act = self.task_dock.toggleViewAction()
		self.task_table_act.setText("Show task table")

		self.bounding_box_act = QAction("Bounding box", self)
		self.bounding_box_act.triggered.connect(self.pymol_viewer.draw_bounding_box)

		self.custom_box_act = QAction("Custom box", self)
		self.custom_box_act.triggered.connect(self.pymol_viewer.draw_custom_box)

		self.run_autodock_act = QAction("Autodock", self)
		self.run_autodock_act.triggered.connect(self.run_autodock)

		#help actions
		self.about_act = QAction("&About", self,
			triggered = self.open_about
		)

		self.doc_act = QAction("Documentation", self,
			triggered = self.open_documentation
		)

		self.issue_act = QAction("Report Issue", self,
			triggered = self.report_issue
		)

		self.update_act = QAction("Check for updates", self,
			triggered = self.check_update
		)

	def create_menus(self):
		self.file_menu = self.menuBar().addMenu("&File")
		self.file_menu.addAction(self.new_project_act)
		self.file_menu.addAction(self.open_project_act)
		self.file_menu.addAction(self.close_project_act)
		self.file_menu.addSeparator()
		self.file_menu.addAction(self.import_receptor_act)
		self.file_menu.addAction(self.import_ligand_act)

		self.edit_menu = self.menuBar().addMenu("&Edit")

		self.view_menu = self.menuBar().addMenu("&View")
		#self.view_menu.addAction(self.open_project_dir_act)
		#self.view_menu.addSeparator()
		self.view_menu.addAction(self.mol_list_act)
		self.view_menu.addAction(self.task_table_act)
		self.view_menu.addSeparator()
		self.view_menu.addAction(self.box_sidebar_act)
		self.view_menu.addAction(self.pymol_sidebar_act)

		self.grid_menu = self.menuBar().addMenu("&Grid")
		self.grid_menu.addAction(self.bounding_box_act)
		self.grid_menu.addAction(self.custom_box_act)

		self.run_menu = self.menuBar().addMenu("&Run")
		self.run_menu.addAction(self.run_autodock_act)

		self.help_menu = self.menuBar().addMenu("&Help")
		self.help_menu.addAction(self.about_act)
		self.help_menu.addAction(self.doc_act)
		self.help_menu.addAction(self.issue_act)
		self.help_menu.addAction(self.update_act)

	def create_toolbar(self):
		self.toolbar = self.addToolBar('')
		self.toolbar.setMovable(False)

		#self.toolbar.addAction(self.import_act)

	def create_statusbar(self):
		self.statusbar = self.statusBar()
		self.statusbar.showMessage("Welcome to Dockey")

	def show_message(self, msg):
		self.statusbar.showMessage(msg)

	def create_molecular_model(self):
		model = MolecularListModel()
		self.mol_viewer.setModel(model)

	def update_molecular_model(self):
		model = self.mol_viewer.model()
		model.beginResetModel()
		model.endResetModel()

	def create_tasks_model(self):
		model = QSqlTableModel()
		model.setTable('tasks')
		model.select()
		self.task_table.setModel(model)

	def create_db_connect(self, db_file):
		if not self.db.connect(db_file):
			QMessageBox.critical(self, "Error",
				"Could not connect to dockey.db file")

		self.create_molecular_model()

	def new_project(self):
		folder = CreateProjectDialog.get_project_folder(self)

		if not folder:
			return

		try:
			os.mkdir(folder)
		except:
			QMessageBox.critical(self, "Error",
				"Could not create project directory")

		db_file = os.path.join(folder, 'dockey.db')
		self.create_db_connect(db_file)
		self.project_folder = folder

	def open_project(self):
		folder = QFileDialog.getExistingDirectory(self)

		if not folder:
			return

		db_file = os.path.join(folder, 'dockey.db')

		if not os.path.exists(db_file):
			QMessageBox.critical(self, "Error",
				"Could not find dockey.db file in {}".format(folder))

		self.create_db_connect(db_file)
		self.project_folder = folder

		self.update_molecular_model()

	def open_project_dir(self):
		QDesktopServices.openUrl(QUrl.fromLocalFile(self.project_folder))

	def close_project(self):
		pass

	def import_receptors(self):
		receptors, _ = QFileDialog.getOpenFileNames(self,
			caption = "Select receptor files",
			filter = (
				"Receptors (*.pdb *.pdbqt *.mol2 *.mol *.smi *.sdf);;"
				"All files (*.*)"
			)
		)

		if not receptors:
			return

		self.db.prepare('molecular')

		for receptor in receptors:
			name = QFileInfo(receptor).baseName()

			with open(receptor) as fh:
				content = fh.read()

			self.db.insert(name, 1, content)

			self.show_message("Import recepotr {}".format(name))

		self.update_molecular_model()

	def import_ligands(self):
		ligands, _ = QFileDialog.getOpenFileNames(self,
			caption = "Select ligand files",
			filter = (
				"Receptors (*.pdb *.pdbqt *.mol2 *.mol *.smi *.sdf);;"
				"All files (*.*)"
			)
		)

		if not ligands:
			return

		self.db.prepare('molecular')

		for ligand in ligands:
			name = QFileInfo(ligand).baseName()

			with open(ligand) as fh:
				content = fh.read()

			self.db.insert(name, 2, content)

			self.show_message("Import ligand {}".format(name))

		self.update_molecular_model()

	def pymol_sidebar_toggle(self, checked):
		self.pymol_viewer.sidebar_controler(checked)

	def run_autodock(self):
		#dlg = AutodockParameterDialog(self)
		dlg = AutodockParamterWizard(self)
		if not dlg.exec():
			return

	def open_about(self):
		QMessageBox.about(self, "About dockey", DOCKEY_ABOUT)

	def open_documentation(self):
		QDesktopServices.openUrl(QUrl("https://dockey.readthedocs.io"))

	def report_issue(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/dockey/issues"))

	def check_update(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/dockey/releases"))

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

class MolecularListModel(QSqlQueryModel):
	def __init__(self, parent=None):
		super(MolecularListModel, self).__init__(parent)
		self.setQuery("SELECT name,type,id FROM molecular ORDER BY type,id")

	def data(self, index, role):
		if role == Qt.DecorationRole:
			r = self.record(index.row())
			v = r.value('type')

			if v == 1:
				return QIcon("icons/receptor.svg")
			elif v == 2:
				return QIcon("icons/ligand.svg")

		return super(MolecularListModel, self).data(index, role)

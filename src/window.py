import os
import psutil

from PySide6.QtGui import *
from PySide6.QtSql import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from view import *
from table import *
from utils import *
from param import *
from widget import *
from worker import *
from config import *
from backend import *
from gridbox import *


__all__ = ['DockeyMainWindow']

class DockeyMainWindow(QMainWindow):
	def __init__(self):
		super(DockeyMainWindow, self).__init__()

		self.setWindowTitle("Dockey v{}".format(DOCKEY_VERSION))
		self.setWindowIcon(QIcon('icons/logo.svg'))
		self.setDockOptions(QMainWindow.AnimatedDocks | QMainWindow.AllowTabbedDocks)

		self.create_pymol_viewer()
		self.create_molecular_viewer()
		self.create_gridbox_adjuster()
		self.create_job_table()
		self.create_pose_table()

		self.create_molecular_model()
		self.create_job_model()
		self.create_pose_model()

		self.create_actions()
		self.create_menus()
		self.create_toolbar()
		self.create_statusbar()

		self.project_folder = None
		self.current_molecular = None
		self.dock_engine = None
		self.dock_params = None
		self.job_num = 0
		self.job_id = 0

		self.pool = QThreadPool(self)
		self.pool.setMaxThreadCount(3)

		self.read_settings()

	def closeEvent(self, event):
		self.write_settings()

		#remove all jobs
		self.pool.clear()

		#kill all subprocess
		parent = psutil.Process(os.getpid())
		for child in parent.children():
			child.kill()

	def read_settings(self):
		settings = QSettings()
		settings.beginGroup('Window')
		self.resize(settings.value('size', QSize(900, 600)))
		self.move(settings.value('pos', QPoint(200, 200)))
		settings.endGroup()

	def write_settings(self):
		settings = QSettings()

		if not self.isMaximized():
			settings.beginGroup('Window')
			settings.setValue('size', self.size())
			settings.setValue('pos', self.pos())
			settings.endGroup()

	def create_pymol_viewer(self):
		self.pymol_viewer = PymolGLWidget(self)
		self.cmd = self.pymol_viewer.cmd
		self.setCentralWidget(self.pymol_viewer)

	def create_molecular_viewer(self):
		#self.mol_viewer = QTreeWidget(self)
		#self.mol_viewer.setColumnCount(1)
		#self.mol_viewer.header().hide()
		self.mol_viewer = DockeyListView(self)
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
		sql = "SELECT * FROM molecular WHERE id=?"
		mol = DB.get_dict(sql, (index.row()+1,))

		self.cmd.load(mol.path)
		self.current_molecular = mol

		sql = "SELECT * FROM grid WHERE rid=?"
		data = DB.get_dict(sql, (mol.id,))

		if data:
			self.box_adjuster.params.update_grid(data)
			draw_gridbox(self.cmd, self.box_adjuster.params)

	def create_gridbox_adjuster(self):
		self.box_adjuster = GridBoxSettingPanel(self)
		self.box_dock = QDockWidget("Gridbox", self)
		self.box_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.box_dock.setWidget(self.box_adjuster)
		self.addDockWidget(Qt.RightDockWidgetArea, self.box_dock)
		self.box_dock.setVisible(False)

	def create_job_table(self):
		self.job_table = DockeyTableView(self)
		self.job_table.clicked.connect(self.on_job_changed)
		self.job_table.setItemDelegateForColumn(4, JobsTableDelegate(self))
		self.job_table.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.job_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.job_table.horizontalHeader().setStretchLastSection(True)
		#self.job_table.verticalHeader().hide()
		#self.job_table.setAlternatingRowColors(True)
		self.job_dock = QDockWidget("Jobs", self)
		self.job_dock.setWidget(self.job_table)
		self.job_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea | Qt.LeftDockWidgetArea)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.job_dock, Qt.Vertical)
		self.job_dock.setVisible(True)

	@Slot()
	def on_job_changed(self, index):
		job_id = index.row() + 1
		sql = "SELECT status,message FROM jobs WHERE id=? LIMIT 1"
		status, message = DB.get_row(sql, (job_id,))

		if status == 1:
			self.pose_dock.setWidget(self.pose_table)
			self.pose_error.hide()
			self.pose_model.set_job(job_id)
			self.pose_table.show()

		elif status == 3:
			self.pose_error.setPlainText(message)
			self.pose_dock.setWidget(self.pose_error)
			self.pose_table.hide()
			self.pose_error.show()	

	def create_pose_table(self):
		self.pose_table = DockeyTableView(self)
		self.pose_table.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.pose_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.pose_table.clicked.connect(self.on_pose_changed)
		self.pose_error = DockeyTextBrowser(self)

		self.pose_dock = QDockWidget("Poses", self)
		self.pose_table.hide()
		self.pose_dock.setWidget(self.pose_error)
		self.pose_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.pose_dock, Qt.Vertical)
		self.pose_dock.setVisible(True)

	@Slot()
	def on_pose_changed(self, index):
		new_index = self.pose_model.index(index.row(), 0)
		pid = new_index.data(Qt.DisplayRole)

		sql = "SELECT jid, ligand FROM pose WHERE id=? LIMIT 1"
		jid, pose = DB.get_row(sql, (pid,))

		sql = (
			"SELECT m1.name,m2.name,m1.path FROM jobs AS j "
			"LEFT JOIN molecular AS m1 ON m1.id=j.rid "
			"LEFT JOIN molecular AS m2 ON m2.id=j.lid "
			"WHERE j.id=? LIMIT 1"
		)

		receptor, ligand, rpath = DB.get_row(sql, (jid,))

		self.cmd.delete('all')
		self.cmd.reinitialize()
		self.cmd.load(rpath, receptor)
		self.cmd.read_pdbstr(pose, ligand)
		#self.cmd.orient()
		self.cmd.zoom()

	def create_molecular_model(self):
		self.mol_model = MolecularTableModel()
		self.mol_viewer.setModel(self.mol_model)
		self.mol_viewer.setModelColumn(1)

	def create_job_model(self):
		self.job_model = JobsTableModel()
		self.job_table.setModel(self.job_model)
		#self.job_table.hideColumn(0)
		#self.job_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

	
	def create_pose_model(self):
		self.pose_model = PoseTableModel()
		self.pose_table.setModel(self.pose_model)

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
		self.box_sidebar_act.setText("Show grid box")
		self.box_sidebar_act.setChecked(False)

		self.mol_list_act = self.mol_dock.toggleViewAction()
		self.mol_list_act.setText("Show molecular list")

		self.job_table_act = self.job_dock.toggleViewAction()
		self.job_table_act.setText("Show task table")

		self.pose_table_act = self.pose_dock.toggleViewAction()
		self.pose_table_act.setText("Show pose table")

		self.bounding_box_act = QAction("Bounding box", self,
			triggered = self.draw_bounding_box
		)

		self.custom_box_act = QAction("Custom box", self,
			triggered = self.draw_custom_box
		)

		self.delete_box_act = QAction("Delete box", self,
			triggered = self.delete_grid_box
		)

		self.run_autodock_act = QAction("Autodock", self,
			triggered = self.run_autodock
		)

		self.run_vina_act = QAction("Autodock Vina", self,
			triggered = self.run_autodock_vina
		)

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
		self.view_menu.addAction(self.job_table_act)
		self.view_menu.addAction(self.pose_table_act)
		self.view_menu.addSeparator()
		self.view_menu.addAction(self.box_sidebar_act)
		self.view_menu.addAction(self.pymol_sidebar_act)

		self.grid_menu = self.menuBar().addMenu("&Grid")
		self.grid_menu.addAction(self.bounding_box_act)
		self.grid_menu.addAction(self.custom_box_act)
		self.grid_menu.addSeparator()
		self.grid_menu.addAction(self.delete_box_act)

		#self.tool_menu = self.menuBar().addMenu("&Tool")

		self.run_menu = self.menuBar().addMenu("&Run")
		self.run_menu.addAction(self.run_autodock_act)
		self.run_menu.addAction(self.run_vina_act)

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

	@Slot()
	def show_error_message(self, msg):
		QMessageBox.critical(self, "Error", msg)

	def create_db_connect(self, db_file):
		#if not DB.connect(db_file):
		#	QMessageBox.critical(self, "Error",
		#		"Could not connect to dockey.db file")
		DB.connect(db_file)

	def new_project(self):
		folder = CreateProjectDialog.get_project_folder(self)

		if not folder:
			return

		try:
			#make project folder
			os.mkdir(folder)

			#make data folder
			data_dir = os.path.join(folder, 'data')
			os.mkdir(data_dir)

			#make jobs folder
			jobs_dir = os.path.join(folder, 'jobs')
			os.mkdir(jobs_dir)

		except:
			return QMessageBox.critical(self, "Error",
				"Could not create project directory")

		db_file = os.path.join(folder, 'dockey.db')
		self.create_db_connect(db_file)
		self.project_folder = folder

		DB.insert_rows("INSERT INTO option VALUES (?,?,?)", [
			(None, 'root_dir', folder),
			(None, 'data_dir', data_dir),
			(None, 'jobs_dir', jobs_dir)
		])

	def open_project(self):
		folder = QFileDialog.getExistingDirectory(self)

		if not folder:
			return

		db_file = os.path.join(folder, 'dockey.db')

		if not os.path.exists(db_file):
			return QMessageBox.critical(self, "Error",
				"Could not find dockey.db file in {}".format(folder))

		self.create_db_connect(db_file)
		self.project_folder = folder

		self.mol_model.select()
		self.job_model.select()

	def open_project_dir(self):
		QDesktopServices.openUrl(QUrl.fromLocalFile(self.project_folder))

	def close_project(self):
		pass

	def import_moleculars(self, mols, _type=1):
		sql = "INSERT INTO molecular VALUES (?,?,?,?,?)"
		rows = []

		for mol in mols:
			qfi = QFileInfo(mol)
			base_name = qfi.baseName()
			file_name = qfi.fileName()
			suffix = qfi.suffix()

			new_path = os.path.join(self.project_folder, 'data', file_name)
			QFile(mol).copy(new_path)

			rows.append([None, base_name, _type, suffix, new_path])

		DB.insert_rows(sql, rows)

		molname = ['','receptor', 'ligand'][_type]

		if len(rows) == 1:
			self.show_message("Import {} {}".format(molname, rows[0][1]))
		else:
			self.show_message("Import {} {}s".format(len(rows), molname))

		self.mol_model.select()

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

		self.import_moleculars(receptors, 1)

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

		self.import_moleculars(ligands, 2)

	def pymol_sidebar_toggle(self, checked):
		self.pymol_viewer.sidebar_controler(checked)

	def draw_bounding_box(self):
		if not self.current_molecular:
			return

		if self.current_molecular.type == 2:
			return

		points = self.cmd.get_extent(self.current_molecular.name)
		self.box_adjuster.params.update_dimension(points)

		draw_gridbox(self.cmd, self.box_adjuster.params)

		self.box_dock.setVisible(True)

	def draw_custom_box(self):
		if not self.current_molecular:
			return

		if self.current_molecular.type == 2:
			return

		self.box_adjuster.params.custom()

		draw_gridbox(self.cmd, self.box_adjuster.params)

		self.box_dock.setVisible(True)

	def delete_grid_box(self):
		self.cmd.delete('gridbox')
		self.box_adjuster.params.reset()
		sql = "DELETE FROM grid WHERE id=?"
		DB.query(sql, (self.current_molecular.id,))

	def get_job(self):
		if not self.job_num:
			self.job_num = DB.get_count('jobs')

		self.job_id += 1

		if self.job_id > self.job_num:
			return None

		sql = (
			"SELECT j.id,j.rid,j.lid,m1.name,m1.path,m1.format,m2.name,m2.path,m2.format "
			"FROM jobs AS j LEFT JOIN molecular AS m1 ON m1.id=j.rid "
			"LEFT JOIN molecular AS m2 ON m2.id=j.lid WHERE j.id=?"
		)
		job = DB.get_row(sql, (self.job_id,))

		sql = "SELECT * FROM grid WHERE rid=?"
		grid = DB.get_row(sql, (job[1],))

		return AttrDict({
			'id': job[0],
			'ri': job[1],
			'rn': job[3],
			'rp': job[4],
			'rf': job[5],
			'li': job[2],
			'ln': job[6],
			'lp': job[7],
			'lf': job[8],
			'x': grid[2],
			'y': grid[3],
			'z': grid[4],
			'cx': grid[5],
			'cy': grid[6],
			'cz': grid[7],
			'spacing': grid[8],
			'root': self.project_folder
		})

	def generate_job_list(self):
		receptors = DB.get_column("SELECT id FROM molecular WHERE type=1")
		ligands = DB.get_column("SELECT id FROM molecular WHERE type=2")

		# check grid box of receptors
		grids = DB.get_column("SELECT rid FROM grid")

		# receptors with no grid box
		rng = set(receptors) - set(grids)

		if rng:
			sql = "SELECT name FROM molecular WHERE id IN ({})".format(
				','.join(map(str,rng))
			)
			rs = DB.get_column(sql)
			QMessageBox.critical(self, 'Error',
				"No grid box was set for receptor {}".format(','.join(rs))
			)
			return False



		rows = []
		for r in receptors:
			for l in ligands:
				rows.append((None, r, l, 0, 0, None, None, None))

		sql = "INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, rows)

		self.job_model.select()

		return len(rows)

	#@Slot()
	#def new_job(self):
	#	self.submit_job()

	#@Slot()
	def submit_job(self):
		job = self.get_job()

		if job is None:
			return

		if self.dock_engine == 'autodock':
			worker = AutodockWorker(job, self.dock_params)
		elif self.dock_engine == 'vina':
			worker = AutodockVinaWorker(job)

		self.pool.start(worker)

		#worker.signals.finished.connect(self.new_job)
		#worker.signals.error.connect(self.show_error_message)
		worker.signals.refresh.connect(self.job_model.update_row)
		

	def check_tool_before_run(self, tool):
		settings = QSettings()
		folder = settings.value('Tools/{}'.format(tool))

		if not folder:
			dlg = ConfigDialog(self)
			dlg.exec()

	def check_jobs(self):
		jobs = DB.get_one("SELECT COUNT(1) FROM jobs LIMIT 1")

		if jobs:
			ret = QMessageBox.warning(self, "Warning",
				"Are you sure you want to delete previous jobs and submit new jobs?",
				QMessageBox.Yes | QMessageBox.No
			)

			if ret == QMessageBox.No:
				return False

		DB.query("DELETE FROM jobs")
		DB.query("DELETE FROM pose")
		self.job_model.select()
		self.pose_model.select()

		return True

	def run_autodock(self):
		if not AutodockConfigDialog.check_executable(self):
			return

		if not self.check_jobs():
			return

		self.check_tool_before_run('autodock_4')

		dlg = AutodockParamterWizard(self)
		if not dlg.exec():
			return

		self.dock_engine = 'autodock'
		self.dock_params = dlg.params

		num = self.generate_job_list()
		for i in range(num):
			self.submit_job()

	def run_autodock_vina(self):
		pass

	def open_about(self):
		QMessageBox.about(self, "About dockey", DOCKEY_ABOUT)

	def open_documentation(self):
		QDesktopServices.openUrl(QUrl("https://dockey.readthedocs.io"))

	def report_issue(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/dockey/issues"))

	def check_update(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/dockey/releases"))


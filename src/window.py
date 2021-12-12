import os

from PySide6.QtGui import *
from PySide6.QtSql import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from param import *
from worker import *
from config import *
from backend import *
from gridbox import *
from pymolview import *

__all__ = ['DockeyMainWindow']

class DockeyMainWindow(QMainWindow):
	def __init__(self):
		super(DockeyMainWindow, self).__init__()

		self.setWindowTitle("Dockey v{}".format(DOCKEY_VERSION))
		self.setWindowIcon(QIcon('icons/logo.svg'))

		self.create_pymol_viewer()
		self.create_molecular_viewer()
		self.create_gridbox_adjuster()
		self.create_job_table()

		self.create_molecular_model()
		self.create_job_model()

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
		self.job_table = QTableView(self)
		self.job_table.setItemDelegateForColumn(4, JobsTableDelegate(self))
		self.job_table.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.job_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.job_table.horizontalHeader().setStretchLastSection(True)
		self.job_table.verticalHeader().hide()
		self.job_table.setAlternatingRowColors(True)
		self.job_dock = QDockWidget("Jobs", self)
		self.job_dock.setWidget(self.job_table)
		self.job_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.job_dock)
		self.job_dock.setVisible(True)

	def create_molecular_model(self):
		self.mol_model = MolecularTableModel()
		self.mol_viewer.setModel(self.mol_model)
		self.mol_viewer.setModelColumn(1)

	def create_job_model(self):
		self.job_model = JobsTableModel()
		self.job_table.setModel(self.job_model)
		#self.job_table.hideColumn(0)
		#self.job_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
		

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

		self.bounding_box_act = QAction("Bounding box", self,
			triggered = self.draw_bounding_box
		)

		self.custom_box_act = QAction("Custom box", self,
			triggered = self.draw_custom_box
		)

		self.delete_box_act = QAction("Delete box", self,
			triggered = self.delete_grid_box
		)

		self.run_autodock_act = QAction("Autodock", self)
		self.run_autodock_act.triggered.connect(self.run_autodock)

		#tool actions
		self.dock_config_act = QAction("Docking tools config", self,
			triggered = self.dock_tools_config
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
		self.view_menu.addSeparator()
		self.view_menu.addAction(self.box_sidebar_act)
		self.view_menu.addAction(self.pymol_sidebar_act)

		self.grid_menu = self.menuBar().addMenu("&Grid")
		self.grid_menu.addAction(self.bounding_box_act)
		self.grid_menu.addAction(self.custom_box_act)
		self.grid_menu.addSeparator()
		self.grid_menu.addAction(self.delete_box_act)

		self.tool_menu = self.menuBar().addMenu("&Tool")
		self.tool_menu.addAction(self.dock_config_act)

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

	@Slot()
	def show_error_message(self, msg):
		QMessageBox.critical(self, "Error", msg)

	def dock_tools_config(self):
		dlg = ConfigDialog(self)
		dlg.exec()

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
		worker.signals.refresh.connect(self.job_model.update)
		

	def check_tool_before_run(self, tool):
		settings = QSettings()
		folder = settings.value('Tools/{}'.format(tool))

		if not folder:
			dlg = ConfigDialog(self)
			dlg.exec()

	def check_jobs(self):
		jobs = DB.get_one("SELECT COUNT(1) FROM jobs LIMIT 1")

		if jobs:
			QMessageBox.warning(self, "Warning", "ddd")

		DB.query("DELETE FROM jobs")
		self.job_model.select()


	def run_autodock(self):
		self.check_jobs()
		self.check_tool_before_run('autodock_4')

		dlg = AutodockParamterWizard(self)
		if not dlg.exec():
			return

		self.dock_engine = 'autodock'
		self.dock_params = dlg.params

		num = self.generate_job_list()
		for i in range(num):
			self.submit_job()

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

class DockeyTableModel(QAbstractTableModel):
	table = None
	custom_headers = []
	row_count = Signal(int)

	def __init__(self, parent=None):
		super(DockeyTableModel, self).__init__(parent)

		#store ids of displayed row
		self.displayed = []

		#cache the current row
		self.cache_row = [-1, None]

		#total row counts
		self.total_count = 0

		#readed row counts
		self.read_count = 0

		#number of readed rows once time
		self._reads = 100

	def rowCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return self.read_count

	def columnCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return len(self.custom_headers)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			return self.get_value(row, col)

		elif role == Qt.BackgroundRole:
			pass

	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if orientation == Qt.Horizontal and role == Qt.DisplayRole:
			return self.custom_headers[section]

		elif orientation == Qt.Vertical and role == Qt.DisplayRole:
			return section+1

	def canFetchMore(self, parent):
		if parent.isValid():
			return False

		if self.read_count < self.total_count:
			return True

		return False

	def fetchMore(self, parent):
		if parent.isValid():
			return

		ids = DB.get_column(self.read_sql)
		fetch_count = len(ids)
		self.beginInsertRows(QModelIndex(), self.read_count, self.read_count+fetch_count-1)
		self.displayed.extend(ids)
		self.read_count += fetch_count
		self.endInsertRows()

	@property
	def count_sql(self):
		return "SELECT COUNT(1) FROM {} LIMIT 1".format(self.table)

	@property
	def read_sql(self):
		remainder = self.total_count - self.read_count
		fetch_count = min(self._reads, remainder)

		return "SELECT id FROM {} LIMIT {},{}".format(
			self.table,
			self.read_count,
			fetch_count
		)

	@property
	def all_sql(self):
		return "SELECT id FROM {}".format(self.table)

	@property
	def get_sql(self):
		return "SELECT * FROM {} WHERE id=? LIMIT 1".format(self.table)

	def get_value(self, row, col):
		if row != self.cache_row[0]:
			self.update_cache(row)

		return self.cache_row[1][col]

	def update_cache(self, row):
		_id = self.displayed[row]
		self.cache_row[0] = row
		self.cache_row[1] = DB.get_row(self.get_sql, (_id,))

	def select(self):
		self.read_count = 0
		self.cache_row = [-1, None]

		self.beginResetModel()
		self.total_count = DB.get_one(self.count_sql)
		self.displayed = DB.get_column(self.read_sql)
		self.read_count = len(self.displayed)
		self.endResetModel()

		self.row_count.emit(self.total_count)

class MolecularTableModel(DockeyTableModel):
	table = 'molecular'
	custom_headers = ['ID','Name', 'Type']

	def data(self, index, role):
		if role == Qt.DecorationRole:
			row = index.row()
			v = self.get_value(row, 2)

			if v == 1:
				return QIcon("icons/receptor.svg")
			elif v == 2:
				return QIcon("icons/ligand.svg")

		return super(MolecularTableModel, self).data(index, role)

class JobsTableModel(DockeyTableModel):
	table = 'jobs'
	custom_headers = ['Job ID', 'Receptor', 'Ligand', 'Status', 'Progress',
		'Start Time', 'Finish Time', 'Message']

	status = {
		0: 'Pending',
		1: 'Success',
		2: 'Running',
		3: 'Error'
	}
	status_colors = {
		0: QColor(255, 193, 7),
		1: QColor(25, 135, 84),
		2: QColor(13, 110, 253),
		3: QColor(220, 53, 69)
	}
	status_icons = {
		0: QIcon('icons/pend.svg'),
		1: QIcon('icons/success.svg'),
		2: QIcon('icons/run.svg'),
		3: QIcon('icons/error.svg')
	}

	@property
	def get_sql(self):
		return (
			"SELECT j.id,m1.name,m2.name,j.status,j.progress,"
			"j.started,j.finished,j.message FROM jobs AS j "
			"LEFT JOIN molecular AS m1 ON m1.id=j.rid "
			"LEFT JOIN molecular AS m2 ON m2.id=j.lid "
			"WHERE j.id=? LIMIT 1"
		)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()
		val = self.get_value(row, col)

		if role == Qt.DisplayRole:
			if col == 3:
				return self.status[val]

			elif col == 5 or col == 6:
				return time_format(val)

			return val

		elif role == Qt.ForegroundRole:
			if col == 3:
				return self.status_colors[val]

		elif role == Qt.DecorationRole:
			if col == 3:
				return self.status_icons[val]

	@Slot()
	def update(self, rowid):
		rowid -= 1
		colid = len(self.custom_headers)-1
		self.update_cache(rowid)
		index1 = self.index(rowid, 0)
		index2 = self.index(rowid, colid)
		self.dataChanged.emit(index1, index2)

class JobsTableDelegate(QStyledItemDelegate):
	def paint(self, painter, option, index):
		percent = index.model().data(index, Qt.DisplayRole)
		percent = percent if percent else 0
		bar = QStyleOptionProgressBar(2)
		bar.rect = option.rect.adjusted(0, 5, 0, -5)
		#bar.rect = QRect(option.rect.x(), option.rect.y()+5, option.rect.width(), option.rect.height()/1.5)
		#bar.rect = option.rect
		bar.minimum = 0
		bar.maximum = 100
		bar.progress = percent
		bar.text = "{}%".format(percent)
		bar.textVisible = True
		bar.textAlignment = Qt.AlignCenter
		bar.state |= QStyle.StateFlag.State_Horizontal

		QApplication.style().drawControl(QStyle.CE_ProgressBar, bar, painter)

class BrowseInput(QWidget):
	def __init__(self, parent=None, is_file=True):
		super(BrowseInput, self).__init__(parent)
		self.input = QLineEdit(self)
		self.input.setReadOnly(True)
		self.browse = QPushButton(self)
		self.browse.setFlat(True)
		self.browse.setIcon(QIcon("icons/folder.svg"))

		if is_file:
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
	def select_file(self):
		exe, _ = QFileDialog.getOpenFileName(self)

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

class ConfigDialog(QDialog):
	def __init__(self, parent=None):
		super(ConfigDialog, self).__init__(parent)

		self.setWindowTitle("Configuration")
		self.resize(QSize(500, -1))

		self.autodock_4_input = BrowseInput(self)
		self.autogrid_4_input = BrowseInput(self)
		self.autodock_vina_input = BrowseInput(self)
		self.autodock_gpu_input = BrowseInput(self)
		self.mgltools_input = BrowseInput(self, False)

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.save_config)
		action_box.rejected.connect(self.reject)

		tips = (
			"Specify the location of tools ("
			"<font color='green'>#</font> any one of them, "
			"<font color='red'>*</font> required)"
		)

		layout = QVBoxLayout()
		layout.addWidget(QLabel(tips, self))
		layout.addWidget(QLabel("Autodock 4 (<font color='green'>#</font>)", self))
		layout.addWidget(self.autodock_4_input)
		layout.addWidget(QLabel("Autogrid 4", self))
		layout.addWidget(self.autogrid_4_input)
		layout.addWidget(QLabel("Autodock Vina (<font color='green'>#</font>)", self))
		layout.addWidget(self.autodock_vina_input)
		layout.addWidget(QLabel("Autodock GPU (<font color='green'>#</font>)", self))
		layout.addWidget(self.autodock_gpu_input)
		layout.addWidget(QLabel("MGLTools Directory (<font color='red'>*</font>)", self))
		layout.addWidget(self.mgltools_input)
		layout.addWidget(action_box)
		self.setLayout(layout)

		self.read_config()

	def save_config(self):
		autodock = self.autodock_4_input.get_text()
		autogrid = self.autogrid_4_input.get_text()
		vina = self.autodock_vina_input.get_text()
		adgpu = self.autodock_gpu_input.get_text()
		mgltools = self.mgltools_input.get_text()

		if not any((autodock, vina, adgpu)):
			return QMessageBox.critical(self, "Config Error",
				"Please specify location for one of Autodock 4, Autodock Vina and Autodock GPU"
			)

		if not mgltools:
			return QMessageBox.critical(self, "Config Error",
				"You must specify install location for MGLTools"
			)

		settings = QSettings()
		settings.beginGroup('Tools')
		settings.setValue('autodock_4', autodock)
		settings.setValue('autogrid_4', autogrid)
		settings.setValue('autodock_vina', vina)
		settings.setValue('autodock_gpu', adgpu)
		settings.setValue('mgltools', mgltools)
		settings.endGroup()

		self.accept()

	def read_config(self):
		settings = QSettings()
		settings.beginGroup('Tools')
		autodock = settings.value('autodock_4', '')
		autogrid = settings.value('autogrid_4', '')
		vina = settings.value('autodock_vina', '')
		adgpu = settings.value('autodock_gpu', '')
		mgltools = settings.value('mgltools', '')
		settings.endGroup()

		self.autodock_4_input.set_text(autodock)
		self.autogrid_4_input.set_text(autogrid)
		self.autodock_vina_input.set_text(vina)
		self.autodock_gpu_input.set_text(adgpu)
		self.mgltools_input.set_text(mgltools)

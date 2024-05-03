import os
import time
import psutil
import multiprocessing

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from pymol._gui import PyMOLDesktopGUI
from pymol.importing import read_mol2str

from view import *
from table import *
from utils import *
from param import *
from widget import *
from worker import *
from config import *
from backend import *
from gridbox import *
from setting import *

__all__ = ['DockeyMainWindow']

class DockeyMainWindow(QMainWindow, PyMOLDesktopGUI):
	project_ready = Signal(bool)

	def __init__(self):
		super().__init__()
		#self.job_query = None
		self.job_params = None
		self.job_engine = None
		#self.job_mutex = QMutex()
		self.job_worker = None

		self.pymol_actions = {}
		#self.pool = QThreadPool(self)
		#self.pool.setMaxThreadCount(3)

		self.setWindowTitle("Dockey v{}".format(DOCKEY_VERSION))
		self.setWindowIcon(QIcon(':/icons/logo.svg'))
		self.setAcceptDrops(True)
		self.setDockOptions(QMainWindow.AnimatedDocks | QMainWindow.AllowTabbedDocks)

		self.create_pymol_viewer()
		self.create_molecular_viewer()
		self.create_gridbox_adjuster()
		self.create_job_table()
		self.create_pose_tab()
		self.create_interaction_tables()
		self.create_pymol_feedback()

		self.create_molecular_model()
		self.create_job_model()
		#self.create_pose_model()

		self.create_actions()
		self.create_menus()
		self.create_toolbar()
		self.create_statusbar()
		self.create_pymol_menus()

		self.read_settings()

		self.show()

	def closeEvent(self, event):
		self.write_settings()

		if DB.active():
			sql = "UPDATE jobs SET status=? WHERE status>2"
			DB.query(sql, (2,))

		if self.job_worker:
			self.job_worker.stop_jobs()

		event.accept()

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

	def create_pymol_menus(self):
		#tb = self.addToolBar('PyMOL')

		#settingaction and _addmenu function was extracted from
		#https://github.com/schrodinger/pymol-open-source/blob/master/modules/pmg_qt/pymol_qt_gui.py
		actiongroups = {}
		cmd = self.cmd

		def SettingAction(parent, cmd, name, label='', true_value=1, false_value=0, command=None):
			'''
			Menu toggle action for a PyMOL setting
			parent: parent QObject
			cmd: PyMOL instance
			name: setting name
			label: menu item text
			'''
			if not label:
				label = name

			index = cmd.setting._get_index(name)
			type_, values = cmd.get_setting_tuple(index)
			action = QAction(label, parent)

			if not command:
				command = lambda: cmd.set(
					index,
					true_value if action.isChecked() else false_value,
					log=1,
					quiet=0)

			parent.pymol_actions.setdefault(index, []).append(lambda v: action.setChecked(v != false_value))

			if type_ in (1, 2, 3, 5, 6):
				action.setCheckable(True)
				if values[0] == true_value:
					action.setChecked(True)
			else:
				print('TODO', type_, name)

			action.triggered.connect(command)
			return action

		def _addmenu(data, menu):
			#menu.setTearOffEnabled(True)
			#menu.setWindowTitle(menu.title())  # needed for Windows
			for item in data:
				if item[0] == 'separator':
					menu.addSeparator()
				elif item[0] == 'menu':
					_addmenu(item[2], menu.addMenu(item[1].replace('&', '&&')))
				elif item[0] == 'command':
					command = item[2]
					if command is None:
						print('warning: skipping', item)
					else:
						if isinstance(command, str):
							command = lambda c=command: cmd.do(c)
						menu.addAction(item[1], command)
				elif item[0] == 'check':
					if len(item) > 4:
						menu.addAction(SettingAction(self, cmd, item[2], item[1], item[3], item[4]))
					else:
						menu.addAction(SettingAction(self, cmd, item[2], item[1]))
				elif item[0] == 'radio':
					label, name, value = item[1:4]
					try:
						group, type_, values = actiongroups[item[2]]
					except KeyError:
						group = QActionGroup(self)
						type_, values = cmd.get_setting_tuple(name)
						actiongroups[item[2]] = group, type_, values
					action = QAction(label, self)
					action.triggered.connect(lambda _=0, args=(name, value): cmd.set(*args, log=1, quiet=0))

					index = cmd.setting._get_index(name)
					self.pymol_actions.setdefault(index, []).append(
						lambda v, V=value, a=action: a.setChecked(v == V)
					)

					group.addAction(action)
					menu.addAction(action)
					action.setCheckable(True)

					if values[0] == value:
						action.setChecked(True)
				elif item[0] == 'open_recent_menu':
					self.open_recent_menu = menu.addMenu('Open Recent...')
				else:
					print('error:', item)

		#fix scene panel dialog
		self.scene_panel_menu_dialog = None
		for _, label, data in self.get_menudata(self.cmd):
			assert _ == 'menu'

			if label in ['File', 'Edit', 'Plugin', 'Help']:
				continue

			menu = self.pymol_menu.addMenu(label)
			_addmenu(data, menu)

	def create_pymol_feedback(self):
		self.feedback_edit = QLineEdit(self)
		self.feedback_edit.returnPressed.connect(self.execute_pymol_cmd)
		self.feedback_edit.setPlaceholderText("PyMOL>")
		self.feedback_browser = FeedbackBrowser(self)
		self.feedback_dock = QDockWidget("Feedbacks", self)
		self.feedback_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
		self.feedback_dock.setWidget(self.feedback_browser)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.feedback_dock, Qt.Vertical)
		self.feedback_dock.setVisible(False)

		#timer
		self.feedback_timer = QTimer(self)
		self.feedback_timer.setSingleShot(True)
		self.feedback_timer.timeout.connect(self.update_feedback)
		self.feedback_timer.start(100)

	def create_interaction_tables(self):
		self.interaction_tab = InteractionTabWidget(self)
		#self.interaction_title = InteractionTitleWidget(self)
		self.interaction_dock = QDockWidget("Interactions", self)
		#self.interaction_dock.setTitleBarWidget(self.interaction_title)
		self.interaction_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
		self.interaction_dock.setWidget(self.interaction_tab)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.interaction_dock, Qt.Vertical)

	@Slot()
	def execute_pymol_cmd(self):
		cmd = self.feedback_edit.text().strip()

		if cmd:
			self.cmd.do(cmd)
			self.feedback_edit.clear()
			self.feedback_timer.start(0)

	@Slot()
	def update_feedback(self):
		feedbacks = self.cmd._get_feedback()
		if feedbacks:
			feedbacks = ["<p>{}</p>".format(feedback) for feedback in feedbacks]
			self.feedback_browser.appendHtml('\n'.join(feedbacks))

		settings = self.cmd.get_setting_updates() or ()
		for setting in settings:
			if setting in self.pymol_actions:
				val = self.cmd.get_setting_tuple(setting)[1][0]
				#self.pymol_actions[setting].setChecked(val)
				for callback in self.pymol_actions[setting]:
					callback(val)

		self.feedback_timer.start(500)

	def create_molecular_viewer(self):
		#self.mol_viewer = QTreeWidget(self)
		#self.mol_viewer.setColumnCount(1)
		#self.mol_viewer.header().hide()
		self.mol_viewer = DockeyListView(self)
		self.mol_viewer.clicked.connect(self.on_molecular_changed)
		self.mol_dock = QDockWidget("Molecules", self)
		self.mol_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.mol_dock.setWidget(self.mol_viewer)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.mol_dock)
		self.mol_dock.setVisible(True)

	def draw_grid_box(self, rid):
		sql = "SELECT * FROM grid WHERE rid=? LIMIT 1"
		data = DB.get_dict(sql, (rid,))

		if data:
			self.box_adjuster.params.update_grid(data)
			draw_gridbox(self.cmd, self.box_adjuster.params)

	def load_mol_to_view(self, mol_name, mol_content, mol_format):
		if mol_format == 'pdb':
			self.cmd.read_pdbstr(mol_content, mol_name)

		elif mol_format == 'mol':
			self.cmd.read_molstr(mol_content, mol_name)

		elif mol_format == 'mol2':
			#content = convert_string_to_pdb(mol_content, mol_format)
			#self.cmd.read_pdbstr(content, mol_name)
			read_mol2str(mol_content, mol_name)

		elif mol_format == 'sdf':
			self.cmd.read_sdfstr(mol_content, mol_name)

	@Slot(QModelIndex)
	def on_molecular_changed(self, index):
		self.cmd.delete('all')
		self.cmd.initialize()

		#name = index.data(Qt.DisplayRole)
		mol_id = index.siblingAtColumn(0).data(Qt.DisplayRole)

		sql = "SELECT id,name,type,content,format FROM molecular WHERE id=? LIMIT 1"
		mol = DB.get_dict(sql, (mol_id,))

		self.load_mol_to_view(mol.name, mol.content, mol.format)

		if mol.type == 1 and self.box_sidebar_act.isChecked():
			self.draw_grid_box(mol.id)

	@Slot(bool)
	def display_grid_box(self, flag):
		objects = self.cmd.get_names('objects')

		if flag:
			sql = "SELECT id FROM molecular WHERE name=? AND type=1 LIMIT 1"
			for name in objects:
				rid = DB.get_one(sql, (name,))

				if rid:
					self.draw_grid_box(rid)
		else:
			if 'gridbox' in objects:
				self.cmd.delete('gridbox')

	def create_gridbox_adjuster(self):
		self.box_adjuster = GridBoxSettingPanel(self)
		self.box_dock = QDockWidget("Grid", self)
		self.box_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.box_dock.setWidget(self.box_adjuster)
		self.addDockWidget(Qt.RightDockWidgetArea, self.box_dock)
		self.box_dock.setVisible(False)

	def create_job_table(self):
		self.job_table = JobTableView(self)
		self.job_table.clicked.connect(self.on_job_changed)
		self.job_table.setItemDelegateForColumn(4, JobsTableDelegate(self))
		self.job_table.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.job_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
		#self.job_table.horizontalHeader().setStretchLastSection(True)
		#self.job_table.verticalHeader().hide()
		#self.job_table.setAlternatingRowColors(True)
		self.job_dock = QDockWidget("Tasks", self)
		self.job_dock.setWidget(self.job_table)
		self.job_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea
			| Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.job_dock, Qt.Vertical)
		self.job_dock.setVisible(True)

	@Slot(QModelIndex)
	def on_job_changed(self, index):
		job_id = index.row() + 1
		sql = (
			"SELECT m1.name,m2.name,j.status,j.message FROM jobs AS j "
			"LEFT JOIN molecular AS m1 ON m1.id=j.lid "
			"LEFT JOIN molecular AS m2 ON m2.id=j.rid "
			"WHERE j.id=? LIMIT 1"
		)
		ligand, receptor, status, message = DB.get_row(sql, (job_id,))
		title = "{} vs {}".format(ligand, receptor)

		if status == 1:
			self.pose_tab.show_job_pose(title, job_id)

		elif status == 0:
			self.pose_tab.show_job_error(title, message)

	def create_pose_tab(self):
		self.pose_tab = PoseTabWidget(self)
		self.pose_tab.pose_table.clicked.connect(self.on_task_pose_changed)
		self.pose_tab.best_table.clicked.connect(self.on_best_pose_changed)

		self.pose_dock = QDockWidget("Poses", self)
		self.pose_dock.setWidget(self.pose_tab)
		self.pose_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.BottomDockWidgetArea)
		self.addDockWidget(Qt.BottomDockWidgetArea, self.pose_dock, Qt.Vertical)
		self.pose_dock.setVisible(True)

	@Slot(QModelIndex)
	def on_task_pose_changed(self, index):
		pid = index.siblingAtColumn(0).data(Qt.DisplayRole)
		self.interaction_tab.change_pose(pid)

	@Slot(QModelIndex)
	def on_best_pose_changed(self, index):
		pid = index.siblingAtColumn(1).data(Qt.DisplayRole)
		self.interaction_tab.change_pose(pid)

	def create_molecular_model(self):
		self.mol_model = MolecularTableModel()
		self.mol_viewer.setModel(self.mol_model)
		self.mol_viewer.setModelColumn(1)

	def create_job_model(self):
		self.job_model = JobsTableModel()
		self.job_table.setModel(self.job_model)
		#fm = self.fontMetrics()
		header = self.job_table.horizontalHeader()
		header.setStretchLastSection(False)
		#header.resizeSection(0, fm.maxWidth()*3)
		#header.resizeSection(3, fm.maxWidth()*5)
		#self.job_table.hideColumn(0)
		#self.job_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
		header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
		header.setSectionResizeMode(1, QHeaderView.Stretch)
		header.setSectionResizeMode(2, QHeaderView.Stretch)
		header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
		header.setSectionResizeMode(4, QHeaderView.ResizeToContents)
	
	#def create_pose_model(self):
	#	self.pose_model = PoseTableModel()
	#	self.pose_table.setModel(self.pose_model)
	#	self.pose_table.hideColumn(0)
	#	self.pose_table.hideColumn(1)

	def create_actions(self):
		self.new_project_act = QAction(QIcon(':/icons/new.svg'), "&New Project", self)
		self.new_project_act.setShortcut(QKeySequence.New)
		self.new_project_act.triggered.connect(self.new_project)
		self.new_project_act.setIconVisibleInMenu(False)

		self.open_project_act = QAction(QIcon(':/icons/open.svg'), "&Open Project", self)
		self.open_project_act.setShortcut(QKeySequence.Open)
		self.open_project_act.triggered.connect(self.open_project)
		self.open_project_act.setIconVisibleInMenu(False)

		self.save_project_act = QAction("&Save Project", self)
		self.save_project_act.setShortcut(QKeySequence.Save)
		self.save_project_act.setDisabled(True)
		self.save_project_act.triggered.connect(self.save_project)
		self.project_ready.connect(self.save_project_act.setEnabled)

		self.save_project_as_act = QAction("&Save Project As", self)
		self.save_project_as_act.setDisabled(True)
		self.save_project_as_act.setShortcut(QKeySequence.SaveAs)
		self.save_project_as_act.triggered.connect(self.save_project_as)
		self.project_ready.connect(self.save_project_as_act.setEnabled)

		self.close_project_act = QAction("&Close Project", self)
		self.close_project_act.setDisabled(True)
		self.close_project_act.setShortcut(QKeySequence.Close)
		self.close_project_act.triggered.connect(self.close_project)
		self.project_ready.connect(self.close_project_act.setEnabled)

		self.import_receptor_act = QAction("&Import Receptors", self)
		self.import_receptor_act.setDisabled(True)
		self.import_receptor_act.triggered.connect(self.import_receptors)
		self.project_ready.connect(self.import_receptor_act.setEnabled)

		self.import_pdb_act = QAction("&PDB", self)
		self.import_pdb_act.setDisabled(True)
		self.import_pdb_act.triggered.connect(self.import_receptor_from_pdb)
		self.project_ready.connect(self.import_pdb_act.setEnabled)

		self.import_ligand_act = QAction("&Import Ligands", self)
		self.import_ligand_act.setDisabled(True)
		self.import_ligand_act.triggered.connect(self.import_ligands)
		self.project_ready.connect(self.import_ligand_act.setEnabled)

		self.import_ligand_sdf_act = QAction("&SDF", self)
		self.import_ligand_sdf_act.setDisabled(True)
		self.import_ligand_sdf_act.triggered.connect(self.import_ligands_from_sdf)
		self.project_ready.connect(self.import_ligand_sdf_act.setEnabled)

		self.import_ligand_dir_act = QAction("&Import Ligands from Directory", self)
		self.import_ligand_dir_act.setDisabled(True)
		self.import_ligand_dir_act.triggered.connect(self.import_ligands_from_dir)
		self.project_ready.connect(self.import_ligand_dir_act.setEnabled)

		self.import_zinc_act = QAction("&ZINC", self)
		self.import_zinc_act.setDisabled(True)
		self.import_zinc_act.triggered.connect(self.import_ligand_from_zinc)
		self.project_ready.connect(self.import_zinc_act.setEnabled)

		self.import_pubchem_act = QAction("&PubChem", self)
		self.import_pubchem_act.setDisabled(True)
		self.import_pubchem_act.triggered.connect(self.import_ligand_from_pubchem)
		self.project_ready.connect(self.import_pubchem_act.setEnabled)

		self.import_coconut_act = QAction("&COCONUT", self)
		self.import_coconut_act.setDisabled(True)
		self.import_coconut_act.triggered.connect(self.import_ligand_from_coconut)
		self.project_ready.connect(self.import_coconut_act.setEnabled)

		self.import_chembl_act = QAction("&ChEMBL", self)
		self.import_chembl_act.setDisabled(True)
		self.import_chembl_act.triggered.connect(self.import_ligand_from_chembl)
		self.project_ready.connect(self.import_chembl_act.setEnabled)

		self.export_image_act = QAction(QIcon(':/icons/image.svg'), "&Export As Image", self)
		self.export_image_act.triggered.connect(self.export_as_image)
		self.export_image_act.setIconVisibleInMenu(False)

		self.export_file_act = QAction("&Export As File", self)
		self.export_file_act.triggered.connect(self.export_as_file)

		self.exit_act = QAction("&Exit", self)
		self.exit_act.setShortcut(QKeySequence.Quit)
		self.exit_act.triggered.connect(self.close)
		#triggered = lambda : print(GridBoxSettingPanel.params.spacing),

		#edit actions
		self.pymol_undo_act = QAction("Undo", self)
		self.pymol_undo_act.setShortcut(QKeySequence.Undo)
		self.pymol_undo_act.triggered.connect(self.cmd.undo)
	
		self.pymol_redo_act = QAction("Redo", self)
		self.pymol_redo_act.setShortcut(QKeySequence.Redo)
		self.pymol_redo_act.triggered.connect(self.cmd.redo)

		self.all_hydro_act = QAction("Add All Hydrogens", self)
		self.all_hydro_act.triggered.connect(self.pymol_add_all_hydrogens)

		self.polar_hydro_act = QAction("Add Polar Hydrogens", self)
		self.polar_hydro_act.triggered.connect(self.pymol_add_polar_hydrogens)

		self.del_water_act = QAction("Remove Water", self)
		self.del_water_act.triggered.connect(self.pymol_remove_water)

		self.del_solvent_act = QAction("Remove Solvent", self)
		self.del_solvent_act.triggered.connect(self.pymol_remove_solvent)

		self.del_organic_act = QAction("Remove Organic", self)
		self.del_organic_act.triggered.connect(self.pymol_remove_organic)

		self.del_chain_act = QAction("Remove Chain", self)
		self.del_chain_act.triggered.connect(self.pymol_remove_chain)

		#preferences actions
		self.setting_act = QAction(QIcon(':/icons/setting.svg'), "Settings", self)
		self.setting_act.setIconVisibleInMenu(False)
		self.setting_act.triggered.connect(self.open_setting_dialog)
		self.setting_act.setShortcut(QKeySequence.Preferences)

		#self.dock_tool_act = QAction("Docking Tools", self)
		#self.dock_tool_act.triggered.connect(self.docking_tool_settings)

		#self.job_num_act = QAction("Concurrent execution", self)
		#self.job_num_act.triggered.connect(self.concurrent_job_number)

		self.task_act = QAction("Concurrent Task Manager", self)
		self.task_act.triggered.connect(self.open_job_manager_dialog)

		#view actions
		#pymol sidebar
		self.pymol_sidebar_act = QAction("Show Pymol Sidebar", self)
		self.pymol_sidebar_act.setCheckable(True)
		self.pymol_sidebar_act.setChecked(False)
		self.pymol_sidebar_act.toggled.connect(self.pymol_sidebar_toggle)

		index = self.cmd.setting._get_index('internal_gui')
		self.pymol_actions.setdefault(index, []).append(lambda v: self.pymol_sidebar_act.setChecked(v))

		self.pymol_feedback_act = self.feedback_dock.toggleViewAction()
		self.pymol_feedback_act.setText("Show Pymol Feedback")

		self.pymol_cmd_act = QAction(QIcon(':/icons/cmd.svg'), "Show Pymol Command", self)
		self.pymol_cmd_act.setCheckable(True)
		self.pymol_cmd_act.setChecked(False)
		self.pymol_cmd_act.setIconVisibleInMenu(False)

		self.check_cpu_act = QAction(QIcon(':/icons/cpu.svg'), "CPU and Memory Usage", self)
		self.check_cpu_act.setIconVisibleInMenu(False)
		self.check_cpu_act.triggered.connect(self.view_cpu_consumption)

		#gridbox setting sidebar
		self.box_sidebar_act = self.box_dock.toggleViewAction()
		self.box_sidebar_act.setIcon(QIcon(':/icons/box.svg'))
		self.box_sidebar_act.setIconVisibleInMenu(False)
		self.box_sidebar_act.setText("Show Grid Box")
		self.box_sidebar_act.setChecked(False)
		self.box_sidebar_act.toggled.connect(self.display_grid_box)

		self.mol_list_act = self.mol_dock.toggleViewAction()
		self.mol_list_act.setIcon(QIcon(':/icons/molecular.svg'))
		self.mol_list_act.setIconVisibleInMenu(False)
		self.mol_list_act.setText("Show Molecular List")

		self.job_table_act = self.job_dock.toggleViewAction()
		self.job_table_act.setIcon(QIcon(":/icons/job.svg"))
		self.job_table_act.setIconVisibleInMenu(False)
		self.job_table_act.setText("Show Job Table")

		self.pose_table_act = self.pose_dock.toggleViewAction()
		self.pose_table_act.setIcon(QIcon(':/icons/pose.svg'))
		self.pose_table_act.setIconVisibleInMenu(False)
		self.pose_table_act.setText("Show Pose Table")

		self.interaction_table_act = self.interaction_dock.toggleViewAction()
		self.interaction_table_act.setIcon(QIcon(':/icons/share.svg'))
		self.interaction_table_act.setIconVisibleInMenu(False)
		self.interaction_table_act.setText("Show Interaction Table")

		self.bounding_box_act = QAction(QIcon(':/icons/bounding.svg'), "Draw Bounding Box", self)
		self.bounding_box_act.triggered.connect(self.draw_bounding_box)
		self.bounding_box_act.setIconVisibleInMenu(False)

		self.centered_box_act = QAction("Draw Residue Centered Box", self)
		self.centered_box_act.triggered.connect(self.draw_centered_box)

		self.custom_box_act = QAction("Draw Custom Box", self)
		self.custom_box_act.triggered.connect(self.draw_custom_box)

		self.delete_box_act = QAction("Delete Grid Box", self)
		self.delete_box_act.triggered.connect(self.delete_grid_box)

		self.run_autodock_act = QAction("Autodock", self)
		self.run_autodock_act.setDisabled(True)
		self.run_autodock_act.triggered.connect(self.run_autodock)
		self.project_ready.connect(self.run_autodock_act.setEnabled)

		self.run_vina_act = QAction("Autodock Vina", self)
		self.run_vina_act.setDisabled(True)
		self.run_vina_act.triggered.connect(self.run_autodock_vina)
		self.project_ready.connect(self.run_vina_act.setEnabled)

		self.run_qvina_act = QAction("QuickVina-W", self)
		self.run_qvina_act.setDisabled(True)
		self.run_qvina_act.triggered.connect(self.run_quick_vina_w)
		self.project_ready.connect(self.run_qvina_act.setEnabled)

		#help actions
		self.about_act = QAction("&About", self)
		self.about_act.triggered.connect(self.open_about)

		self.thank_act = QAction("Acknowledgements", self)
		self.thank_act.triggered.connect(self.open_thanks)

		self.doc_act = QAction("Documentation", self)
		self.doc_act.triggered.connect(self.open_documentation)

		self.issue_act = QAction("Report Issue", self)
		self.issue_act.triggered.connect(self.report_issue)

		self.update_act = QAction("Check for updates", self)
		self.update_act.triggered.connect(self.check_update)

	def create_menus(self):
		self.file_menu = self.menuBar().addMenu("&File")
		self.file_menu.addAction(self.new_project_act)
		self.file_menu.addAction(self.open_project_act)
		#self.file_menu.addAction(self.save_project_act)
		self.file_menu.addAction(self.save_project_as_act)
		self.file_menu.addAction(self.close_project_act)
		self.file_menu.addSeparator()
		self.file_menu.addAction(self.import_receptor_act)
		receptor_menu = self.file_menu.addMenu("&Import Receptors from Database")
		receptor_menu.addAction(self.import_pdb_act)
		#self.file_menu.addAction(self.import_receptors_act)
		#self.file_menu.addAction(self.import_pdb_act)
		self.file_menu.addAction(self.import_ligand_act)
		self.file_menu.addAction(self.import_ligand_dir_act)
		ligand_menu = self.file_menu.addMenu("&Import Ligands from Database")
		ligand_menu.addAction(self.import_ligand_sdf_act)
		ligand_menu.addAction(self.import_zinc_act)
		ligand_menu.addAction(self.import_pubchem_act)
		ligand_menu.addAction(self.import_chembl_act)
		#ligand_menu.addAction(self.import_coconut_act)

		self.file_menu.addSeparator()
		self.file_menu.addAction(self.export_image_act)
		self.file_menu.addAction(self.export_file_act)
		self.file_menu.addSeparator()
		self.file_menu.addAction(self.exit_act)

		self.edit_menu = self.menuBar().addMenu("&Edit")
		#self.edit_menu.addAction(self.pymol_undo_act)
		#self.edit_menu.addAction(self.pymol_redo_act)
		#self.edit_menu.addSeparator()
		#bg_menu = self.edit_menu.addMenu("Background")
		#bg_menu.addAction(self.pymol_color_act)
		#bg_menu.addAction(self.pymol_opaque_act)
		self.edit_menu.addAction(self.all_hydro_act)
		self.edit_menu.addAction(self.polar_hydro_act)
		self.edit_menu.addSeparator()
		self.edit_menu.addAction(self.del_water_act)
		self.edit_menu.addAction(self.del_solvent_act)
		self.edit_menu.addAction(self.del_organic_act)
		self.edit_menu.addAction(self.del_chain_act)
		self.edit_menu.addSeparator()
		self.edit_menu.addAction(self.setting_act)

		self.view_menu = self.menuBar().addMenu("&View")
		#self.view_menu.addAction(self.open_project_dir_act)
		#self.view_menu.addSeparator()
		self.view_menu.addAction(self.mol_list_act)
		self.view_menu.addAction(self.job_table_act)
		self.view_menu.addAction(self.pose_table_act)
		self.view_menu.addAction(self.interaction_table_act)
		self.view_menu.addSeparator()
		self.view_menu.addAction(self.box_sidebar_act)
		self.view_menu.addSeparator()
		self.view_menu.addAction(self.pymol_sidebar_act)
		self.view_menu.addAction(self.pymol_cmd_act)
		self.view_menu.addAction(self.pymol_feedback_act)

		self.grid_menu = self.menuBar().addMenu("&Grid")
		self.grid_menu.addAction(self.bounding_box_act)
		self.grid_menu.addAction(self.centered_box_act)
		self.grid_menu.addAction(self.custom_box_act)
		self.grid_menu.addSeparator()
		self.grid_menu.addAction(self.delete_box_act)

		#self.tool_menu = self.menuBar().addMenu("&Tool")

		self.run_menu = self.menuBar().addMenu("&Run")
		self.run_menu.addAction(self.run_autodock_act)
		self.run_menu.addAction(self.run_vina_act)
		self.run_menu.addAction(self.run_qvina_act)

		self.task_menu = self.menuBar().addMenu("&Task")
		self.task_menu.addAction(self.task_act)
		self.task_menu.addAction(self.check_cpu_act)

		self.pymol_menu = self.menuBar().addMenu("&PyMOL")

		self.help_menu = self.menuBar().addMenu("&Help")
		self.help_menu.addAction(self.about_act)
		self.help_menu.addAction(self.thank_act)
		self.help_menu.addAction(self.doc_act)
		self.help_menu.addAction(self.issue_act)
		self.help_menu.addAction(self.update_act)

	def create_toolbar(self):
		self.toolbar = self.addToolBar('')
		#self.toolbar.setIconSize(QSize(28,28))
		self.toolbar.setMovable(False)
		#self.toolbar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
		self.toolbar.addAction(self.new_project_act)
		self.toolbar.addAction(self.open_project_act)
		self.toolbar.addSeparator()
		self.toolbar.addAction(self.mol_list_act)
		self.toolbar.addAction(self.job_table_act)
		self.toolbar.addAction(self.pose_table_act)
		self.toolbar.addAction(self.interaction_table_act)
		self.toolbar.addSeparator()
		self.toolbar.addAction(self.box_sidebar_act)
		self.toolbar.addAction(self.bounding_box_act)
		self.toolbar.addSeparator()
		self.toolbar.addAction(self.export_image_act)
		self.toolbar.addSeparator()
		self.toolbar.addAction(self.setting_act)
		self.toolbar.addAction(self.check_cpu_act)
		self.toolbar.addAction(self.pymol_cmd_act)
		self.feedback_act = self.toolbar.addWidget(self.feedback_edit)
		self.feedback_act.setVisible(False)
		self.pymol_cmd_act.toggled.connect(self.feedback_act.setVisible)

	def create_statusbar(self):
		self.statusbar = self.statusBar()
		self.statusbar.showMessage("Welcome to Dockey")

		self.progressbar = QProgressBar(self)
		self.progressbar.setStyleSheet("""
			QProgressBar {
				text-align: right;
				max-width: 100px;
				min-width: 100px;
				max-height: 15px;
				min-height: 15px;
			}
		""")
		self.statusbar.addPermanentWidget(self.progressbar)

	def show_message(self, msg):
		self.statusbar.showMessage(msg)

	@Slot(str)
	def show_error_message(self, msg):
		QMessageBox.critical(self, "Error", msg)

	@Slot(str)
	def show_popup_message(self, msg):
		QMessageBox.information(self, "Message", msg)

	def create_db_connect(self, db_file):
		try:
			DB.connect(db_file)
		except:
			return QMessageBox.critical(self, "Error",
				"Could not open the project file, it's not a dockey project file"
			)

		self.mol_model.select()
		self.job_model.select()
		self.pose_tab.select()
		self.project_ready.emit(True)

		self.setWindowTitle("{} - Dockey".format(db_file))

	def dragEnterEvent(self, event):
		event.acceptProposedAction()

	def dragMoveEvent(self, event):
		event.acceptProposedAction()

	def dropEvent(self, event):
		if event.mimeData().hasUrls():
			url = event.mimeData().urls()[0]
			project_file = url.toLocalFile()

			if project_file.endswith('.dock'):
				if DB.active():
					ret = QMessageBox.warning(self, "Warning",
						"Are you sure you want to close the current project and open the dragged project?",
						QMessageBox.Yes | QMessageBox.No
					)

					if ret == QMessageBox.No:
						return
					else:
						self.close_project()

				self.create_db_connect(project_file)
				
			else:
				QMessageBox.critical(self, "Error",
					"{} is not a dockey project file".format(project_file)
				)

	def open_setting_dialog(self):
		#dlg = DockeyConfigDialog(self)
		dlg = DockeyGlobalSettingDialog(self)
		dlg.exec()

	def open_job_manager_dialog(self):
		dlg = TaskManagerDialog(self)

		if dlg.exec() == QDialog.Accepted:
			settings = QSettings()
			num = settings.value('Job/concurrent', 1, int)

			if self.job_worker:
				self.job_worker.set_thread(num)

	def new_project(self):
		if DB.active():
			ret = QMessageBox.warning(self, "Warning",
				"Are you sure you want to close the current project and create a new project?",
				QMessageBox.Yes | QMessageBox.No
			)

			if ret == QMessageBox.No:
				return
			else:
				self.close_project()

		project_file, _ = QFileDialog.getSaveFileName(self,
			caption = "New Project File",
			filter = "Dockey Project File (*.dock)"
		)

		if not project_file:
			return

		if not project_file.endswith('.dock'):
			project_file += '.dock'

		if os.path.isfile(project_file):
			os.remove(project_file)

		self.create_db_connect(project_file)

	def open_project(self):
		if DB.active():
			ret = QMessageBox.warning(self, "Warning",
				"Are you sure you want to close the current project and open a new project?",
				QMessageBox.Yes | QMessageBox.No
			)

			if ret == QMessageBox.No:
				return
			else:
				self.close_project()

		project_file, _ = QFileDialog.getOpenFileName(self, filter="Dockey Project File (*.dock)")

		if not project_file:
			return

		self.create_db_connect(project_file)

	def save_project(self):
		if DB.changed():
			DB.commit()

	def save_project_as(self):
		project_file, _ = QFileDialog.getSaveFileName(self, filter="Dockey Project File (*.dock)")

		if not project_file:
			return

		if not project_file.endswith('.dock'):
			project_file += '.dock'

		#if os.path.isfile(project_file):
		#	os.remove(project_file)

		DB.save(project_file)

	def close_project(self):
		if DB.changed():
			ret = QMessageBox.question(self, "Confirmation",
				"Would you like to save project before closing project",
				QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
			)

			if ret == QMessageBox.Cancel:
				return

			elif ret == QMessageBox.Yes:
				self.save_project()

		self.mol_model.reset()
		self.job_model.reset()
		self.pose_tab.reset()
		self.interaction_tab.reset()

		#reset pymol
		self.cmd.delete('all')
		self.cmd.initialize()

		DB.close()
		self.project_ready.emit(False)

		self.setWindowTitle("Dockey v{}".format(DOCKEY_VERSION))

	def import_moleculars(self, mols, _type=1):
		sql = "INSERT INTO molecular VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
		rows = []

		for mol in mols:
			mi = get_molecule_information(mol)

			if 'error' in mi:
				return QMessageBox.critical(self, "Error", mi.error)

			rows.append([None, mi.name, _type, mi.pdb, mi.atoms, mi.bonds,
				mi.hvyatoms, mi.residues, mi.rotors, mi.formula, mi.energy,
				mi.weight, mi.logp
			])

		DB.insert_rows(sql, rows)

		mol_type = ['', 'receptor', 'ligand'][_type]

		if len(rows) == 1:
			self.show_message("Import {} {}".format(mol_type, rows[0][1]))
		else:
			self.show_message("Import {} {}s".format(len(rows), mol_type))

		self.mol_model.select()

	def start_import_thread(self, thread):
		pool = QThreadPool.globalInstance()

		if pool.activeThreadCount() > 0:
			return QMessageBox.warning(self, "Warning", "There is already a task running")

		thread.signals.message.connect(self.show_message)
		thread.signals.failure.connect(self.show_error_message)
		thread.signals.success.connect(self.show_popup_message)
		thread.signals.finished.connect(self.mol_model.select)
		thread.signals.progress.connect(self.progressbar.setValue)
		pool.start(thread)

	def import_receptors(self):
		receptors, _ = QFileDialog.getOpenFileNames(self,
			caption = "Select receptor files",
			filter = (
				"Receptors (*.pdb *.mol2 *.mol *.sdf);;"
				"All files (*.*)"
			)
		)

		if not receptors:
			return

		importor = ImportFileWorker(receptors, 1)
		self.start_import_thread(importor)

	def import_receptor_from_pdb(self):
		urls = PDBDownloadDialog.get_urls(self)
		
		if not urls:
			return

		importor = ImportURLWorker(urls, 1)
		self.start_import_thread(importor)

	def import_ligands(self):
		ligands, _ = QFileDialog.getOpenFileNames(self,
			caption = "Select ligand files",
			filter = (
				"Receptors (*.pdb *.mol2 *.mol *.sdf);;"
				"All files (*.*)"
			)
		)

		if not ligands:
			return

		importor = ImportFileWorker(ligands, 2)
		self.start_import_thread(importor)

	def import_ligands_from_sdf(self):
		sdf_file, _ = QFileDialog.getOpenFileName(self,
			caption = "Select a SDF file",
			filter = ("SDF (*.sdf);;All files (*.*)")
		)

		if not sdf_file:
			return

		props = get_sdf_props(sdf_file)
		prop, ret = QInputDialog.getItem(self, "Select property", 
			"Select a property as molecular name", props)

		if not ret:
			return

		importor = ImportSDFWorker([sdf_file, prop], 2)
		self.start_import_thread(importor)

	def import_ligands_from_dir(self):
		folder = QFileDialog.getExistingDirectory(self)

		if not folder:
			return

		importor = ImportFolderWorker(folder, 2)
		self.start_import_thread(importor)

	def import_ligand_from_zinc(self):
		urls = ZINCDownloadDialog.get_urls(self)
		
		if not urls:
			return

		importor = ImportURLWorker(urls, 2)
		self.start_import_thread(importor)

	def import_ligand_from_pubchem(self):
		urls = PubchemDownloadDialog.get_urls(self)

		if not urls:
			return

		importor = ImportURLWorker(urls, 2)
		self.start_import_thread(importor)

	def import_ligand_from_chembl(self):
		urls = ChemblDownloadDialog.get_urls(self)

		if not urls:
			return

		importor = ImportURLWorker(urls, 2)
		self.start_import_thread(importor)

	def import_ligand_from_coconut(self):
		pass

	def export_as_image(self):
		opt = ExportImageDialog.save_to_png(self)

		if not opt: return

		self.cmd.png(opt.image, opt.width, opt.height, opt.dpi, ray=1)

	def export_as_file(self):
		outfile, _ = QFileDialog.getSaveFileName(self, filter="PDB (*.pdb)")

		if not outfile: return

		if not outfile.endswith('.pdb'):
			outfile += '.pdb'

		self.cmd.save(outfile)

	#@Slot()
	#def pymol_undo(self):
	#	self.cmd.undo()

	#@Slot()
	#def pymol_redo(self):
	#	self.cmd.redo()

	@Slot()
	def pymol_settings(self):
		dlg = PymolSettingDialog(self)
		dlg.exec()

	@Slot()
	def docking_tool_settings(self):
		dlg = DockingToolSettingDialog(self)
		dlg.exec()

	@Slot()
	def concurrent_job_number(self):
		dlg = JobConcurrentSettingDialog(self)
		dlg.exec()

	def pymol_sidebar_toggle(self, checked):
		self.pymol_viewer.sidebar_controler(checked)

	def pymol_opaque_background(self, checked):
		self.cmd.set('ray_opaque_background', checked)

	def pymol_auto_save_mol(self):
		index = self.mol_viewer.currentIndex()
		mol_id = index.siblingAtColumn(0).data()
		mol_name = index.siblingAtColumn(1).data()
		sql = "SELECT format FROM molecular WHERE id=? LIMIT 1"
		mol_format = DB.get_one(sql, (mol_id,))

		if mol_name not in self.cmd.get_names('all'):
			return

		content = self.cmd.get_str(mol_format, mol_name)

		sql = "UPDATE molecular SET content=? WHERE id=?"
		DB.query(sql, (content, mol_id))

	def pymol_add_all_hydrogens(self):
		if not self.cmd.get_object_list():
			return

		self.cmd.h_add("(all)")
		#self.cmd.sort("(all) extend 1")
		self.pymol_auto_save_mol()

	def pymol_add_polar_hydrogens(self):
		if not self.cmd.get_object_list():
			return

		self.cmd.h_add("(all) & (don.|acc.)")
		#self.cmd.sort("(all) extend 1")
		self.pymol_auto_save_mol()

	def pymol_remove_water(self):
		self.cmd.remove('resn hoh')
		self.pymol_auto_save_mol()

	def pymol_remove_solvent(self):
		self.cmd.remove('solvent')
		self.pymol_auto_save_mol()

	def pymol_remove_organic(self):
		self.cmd.remove('organic')
		self.pymol_auto_save_mol()

	def pymol_remove_chain(self):
		chains = self.cmd.get_chains("(all)")

		if not chains:
			return

		#chain, ok = QInputDialog.getItem(self, "Remove Chain",
		#	"Select a chain to remove",
		#	chains, 0, False
		#)

		dlg = QInputDialog(self)
		dlg.setWindowTitle("Remove chain")
		dlg.setLabelText("Select a chain:")
		dlg.setComboBoxItems(chains)
		dlg.textValueChanged.connect(lambda c: self.cmd.select('chain {}'.format(c)))

		self.cmd.select('chain {}'.format(chains[0]))

		if dlg.exec() == QDialog.Accepted:
			chain = dlg.textValue()
			self.cmd.remove('chain {}'.format(chain))
			self.pymol_auto_save_mol()

	def get_current_receptor(self):
		objs = self.cmd.get_names('objects')
		sql = "SELECT * FROM molecular WHERE name=? LIMIT 1"

		for obj in objs:
			if obj == 'gridbox':
				continue

			mol = DB.get_dict(sql, (objs[0],))
			if mol.type == 1:
				return mol

		QMessageBox.critical(self, "Error", "There is no receptor in pymol viewer")

	def draw_bounding_box(self):
		r = self.get_current_receptor()

		if not r: return

		points = self.cmd.get_extent(r.name)
		self.box_adjuster.params.update_dimension(points)
		draw_gridbox(self.cmd, self.box_adjuster.params)
		self.box_dock.setVisible(True)

	def draw_centered_box(self):
		r = self.get_current_receptor()

		if not r: return

		items = []
		for res in get_molecule_residues(r.content, r.format):
			items.append("Chain {}, {}-{}".format(res[1], res[2].capitalize(), res[3]))

		res, ok = QInputDialog.getItem(self, "Select a residue", "Draw grid box centered on a selected residue:", items)

		if not ok: return

		sel = res.replace(', ', ' and resn ').replace('-', ' and resi ')
		

		points = self.cmd.get_extent(sel)
		self.box_adjuster.params.update_dimension(points)
		draw_gridbox(self.cmd, self.box_adjuster.params)
		self.box_dock.setVisible(True)

	def draw_custom_box(self):
		r = self.get_current_receptor()

		if not r: return

		self.box_adjuster.params.custom()
		draw_gridbox(self.cmd, self.box_adjuster.params)
		self.box_dock.setVisible(True)

	def delete_grid_box(self):
		r = self.get_current_receptor()

		if not r:
			return

		self.cmd.delete('gridbox')
		self.box_adjuster.params.reset()
		DB.query("DELETE FROM grid WHERE id=?", (r.id,))

	def generate_job_list(self):
		pool = QThreadPool.globalInstance() 

		if pool.activeThreadCount() > 0:
			return QMessageBox.warning(self, "Warning", "There is already a task running")

		DB.set_option('tool', self.job_engine)

		thread = JobListGenerator()
		thread.signals.message.connect(self.show_message)
		thread.signals.finished.connect(self.job_model.select)
		thread.signals.finished.connect(self.start_jobs)
		thread.signals.progress.connect(self.progressbar.setValue)
		pool.start(thread)

	def stop_job(self, job):
		ret = QMessageBox.warning(self, "Warning",
			"Are you sure you want to stop this job?",
			QMessageBox.Yes | QMessageBox.No
		)

		if ret == QMessageBox.No:
			return False

		if self.job_worker is not None:
			self.job_worker.stop_job(job)

		#sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		#DB.query(sql, (4, 'Stopped', job))
		#self.job_model.update_row(job)

	def start_jobs(self):
		ret = DB.get_one("SELECT 1 FROM jobs LIMIT 1")

		if not ret: return

		params = self.job_params.deep_copy()

		if self.job_engine == 'autodock':
			#self.job_worker = AutodockWorker(params)
			self.job_worker = WorkerManager(self, AutodockWorker, params)

		elif self.job_engine == 'vina':
			#self.job_worker = AutodockVinaWorker(params)
			self.job_worker = WorkerManager(self, AutodockVinaWorker, params)

		elif self.job_engine == 'qvina':
			#self.job_worker = QuickVinaWorker(params)
			self.job_worker = WorkerManager(self, QuickVinaWorker, params)

		#self.job_worker.signals.refresh.connect(self.job_model.update_row)
		#self.job_worker.signals.failure.connect(self.show_error_message)
		#self.job_worker.signals.finished.connect(self.pose_tab.select)
		#self.job_worker.signals.progress.connect(self.progressbar.setValue)

		#self.job_worker.signals.updated.connect(self.job_model.update_row)
		self.job_worker.signals.finished.connect(self.pose_tab.select)
		self.job_worker.signals.progress.connect(self.progressbar.setValue)
		self.job_worker.signals.message.connect(self.show_message)
		#self.job_thread = QThread(self)
		#self.job_thread.started.connect(self.job_worker.run)
		#self.job_worker.finished.connect(self.job_thread.quit)
		#self.job_worker.finished.connect(self.job_worker.deleteLater)
		#self.job_worker.moveToThread(self.job_thread)
		#self.job_thread.start()
		QThreadPool.globalInstance().start(self.job_worker)

	def check_mols(self):
		receptors = int(DB.get_option('receptor_count'))
		ligands = int(DB.get_option('ligand_count'))

		if not (receptors and ligands):
			QMessageBox.critical(self, "Error", "There are no receptors and ligands")
			return False

		if not receptors:
			QMessageBox.critical(self, "Error", "There are no receptors")
			return False

		if not ligands:
			QMessageBox.critical(self, "Error", "There are no ligands")
			return False

		return True

	def check_jobs(self):
		jobs = int(DB.get_option("job_count") or 0)

		if jobs:
			ret = QMessageBox.warning(self, "Warning",
				"Are you sure you want to delete the previous jobs and analysis results?",
				QMessageBox.Yes | QMessageBox.No
			)

			if ret == QMessageBox.No:
				return False

		self.job_model.clear()
		self.pose_tab.clear()
		self.interaction_tab.clear()

		#clear logs
		DB.query("DELETE FROM logs")
		#delete tool option
		DB.set_option('tool', '')
		DB.set_option('job_count', 0)
		DB.set_option('best_count', 0)

		return True

	def run_autodock(self):
		if not self.check_mols():
			return

		if not AutodockConfigDialog.check_executable(self):
			return

		if not self.check_jobs():
			return

		dlg = AutodockParameterWizard(self)
		if not dlg.exec():
			return

		self.job_params = dlg.params
		self.job_engine = 'autodock'
		self.generate_job_list()

		#self.submit_jobs('autodock', params)

	def run_autodock_vina(self):
		if not self.check_mols():
			return

		if not AutodockVinaConfigDialog.check_executable(self):
			return

		if not self.check_jobs():
			return

		dlg = AutodockVinaParameterWizard(self)
		if not dlg.exec():
			return

		self.job_params = dlg.params
		self.job_engine = 'vina'
		self.generate_job_list()
		#self.submit_jobs('vina', params)

	def run_quick_vina_w(self):
		if not self.check_mols():
			return

		if not QuickVinaConfigDialog.check_executable(self):
			return

		if not self.check_jobs():
			return

		dlg = QuickVinaParameterWizard(self)
		if not dlg.exec():
			return

		self.job_params = dlg.params
		self.job_engine = 'qvina'
		self.generate_job_list()
		#self.submit_jobs('qvina', params)

	def view_cpu_consumption(self):
		dlg = CPUAndMemoryViewDialog(self)
		dlg.exec()

	def open_about(self):
		QMessageBox.about(self, "About dockey", DOCKEY_ABOUT)

	def open_thanks(self):
		AcknowledgementDialog.thank(self, DOCKEY_THANKS)

	def open_documentation(self):
		QDesktopServices.openUrl(QUrl("https://dockey.readthedocs.io"))

	def report_issue(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/dockey/issues"))

	def check_update(self):
		QDesktopServices.openUrl(QUrl("https://github.com/lmdu/dockey/releases"))

	def set_receptor_flexres(self, mol_id):
		dlg = FlexResiduesDialog(self, mol_id)
		dlg.exec()

	def filter_ligands(self):
		condition = LigandFilterDialog.filter(self)

		if condition:
			sql = "SELECT COUNT(1) FROM molecular WHERE type=2 AND {}".format(condition)
			count = DB.get_one(sql)

			if count == 0:
				QMessageBox.warning(self, "Warning", "No molecules match these filters")

			else:
				ret = QMessageBox.question(self, "Comfirmation",
					"Are you sure you want to remove {} ligand(s) that match these filters?".format(count))

				if ret == QMessageBox.Yes:
					DB.query("DELETE FROM molecular WHERE type=2 AND {}".format(condition))

					ligand_count = int(DB.get_option('ligand_count'))
					ligand_count -= count
					DB.set_option('ligand_count', ligand_count)

					total_count = int(DB.get_option('molecule_count'))
					total_count -= count
					DB.set_option('molecule_count', total_count)

					self.mol_model.select()


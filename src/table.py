import os
import sys
import csv
import psutil

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from worker import *
from backend import *

__all__ = ['DockeyListView', 'JobTableView', 'PoseTableView',
	'DockeyTextBrowser', 'DockeyTableModel', 'MolecularTableModel',
	'JobsTableModel', 'JobsTableDelegate', 'PoseTableModel',
	'HydrogenBondsModel', 'HydrophobicInteractionModel', 'HalogenBondsModel',
	'SaltBridgesModel', 'WaterBridgesModel', 'PiStackingModel', 'PiCationModel',
	'MetalComplexModel', 'BindingSiteModel', 'InteractionTableView',
	'BestTableView', 'BestTableModel'
]

class MoleculeDetailDialog(QDialog):
	def __init__(self, parent, info):
		super(MoleculeDetailDialog, self).__init__(parent)
		self.setWindowTitle("Molecule counts")
		layout = QVBoxLayout()
		self.setLayout(layout)

		label = QLabel(info, self)
		label.setWordWrap(True)
		layout.addWidget(label)

		btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
		btn_box.accepted.connect(self.accept)
		btn_box.rejected.connect(self.reject)
		layout.addWidget(btn_box)

class JobStatusesDialog(QDialog):
	def __init__(self, parent, info):
		super(JobStatusesDialog, self).__init__(parent)
		self.setWindowTitle("Statuses Overview")
		layout = QVBoxLayout()
		self.setLayout(layout)

		layout.addWidget(QLabel("Number of Jobs:", self))

		label = QLabel(info, self)
		label.setWordWrap(True)
		layout.addWidget(label)

		btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
		btn_box.accepted.connect(self.accept)
		btn_box.rejected.connect(self.reject)
		layout.addWidget(btn_box)

class JobDetailDialog(QDialog):
	def __init__(self, parent, info, pid):
		super(JobDetailDialog, self).__init__(parent)
		self.parent = parent
		self.procs = []

		self.setWindowTitle("Job Details")
		layout = QVBoxLayout()
		self.setLayout(layout)

		info_label = QLabel(info, self)
		info_label.setWordWrap(True)
		layout.addWidget(info_label)

		self.usage_label = QLabel("<table cellspacing='10'><tr><td>CPU: 0 | Memory: 0</td></tr></table>", self)
		self.usage_label.setAlignment(Qt.AlignCenter)
		layout.addWidget(self.usage_label)

		btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
		btn_box.accepted.connect(self.accept)
		btn_box.rejected.connect(self.reject)
		layout.addWidget(btn_box)

		if pid and psutil.pid_exists(pid):
			proc = psutil.Process(pid)

			for child in proc.children():
				self.procs.append(child)

			self.procs.append(proc)

		else:
			self.procs = None

		self.timer = QTimer(self)
		self.timer.timeout.connect(self.update_usage)
		self.timer.start(1000)
		self.accepted.connect(self.timer.stop)

	@Slot()
	def update_usage(self):
		if self.procs is None:
			self.usage_label.setText("<table cellspacing='10'><tr><td>CPU: 0 | Memory: 0</td></tr></table>")
			return

		cpu = 0
		mem = 0
		memp = 0

		try:
			for proc in self.procs:
				cpu += proc.cpu_percent()
				mem += proc.memory_info().rss
				memp += proc.memory_percent()
		except:
			self.procs = None
			self.usage_label.setText("<table cellspacing='10'><tr><td>CPU: 0 | Memory: 0</td></tr></table>")
			return
			
		self.usage_label.setText(
			"<table cellspacing='10'><tr><td>CPU: {:.2f}% | Memory: {}, {:.2f}%</td></tr></table>".format(
				cpu/psutil.cpu_count(False),
				memory_format(mem),
				memp
			)
		)

class JobLogDialog(QDialog):
	def __init__(self, parent, job):
		super(JobLogDialog, self).__init__(parent)
		self.job = job
		self.resize(QSize(650, 400))
		self.tab_bar = QTabBar(self)
		self.log_viewer = QTextBrowser(self)
		self.log_names = self.get_log_names()

		for name in self.log_names:
			self.tab_bar.addTab(name)

		self.tab_bar.currentChanged.connect(self.on_tab_changed)
		self.on_tab_changed(0)

		btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
		btn_box.accepted.connect(self.accept)

		self.layout = QVBoxLayout()
		self.layout.setSpacing(0)
		self.layout.addWidget(self.tab_bar)
		self.layout.addWidget(self.log_viewer)
		self.layout.addItem(QSpacerItem(5,10))
		self.layout.addWidget(btn_box)
		self.setLayout(self.layout)

	def get_log_names(self):
		sql = "SELECT name FROM logs WHERE jid=?"
		return DB.get_column(sql, (self.job,))

	@Slot(int)
	def on_tab_changed(self, index):
		name = self.log_names[index]
		sql = "SELECT content FROM logs WHERE jid=? AND name=? LIMIT 1"
		content = DB.get_one(sql, (self.job, name))
		content = content.replace('\0', '')

		self.log_viewer.setText(content)

	@classmethod
	def view_log(cls, parent, job):
		dlg = cls(parent, job)
		dlg.exec()

class DockeyListView(QListView):
	def __init__(self, parent=None):
		super(DockeyListView, self).__init__(parent)
		self.parent = parent
		self.setContextMenuPolicy(Qt.CustomContextMenu)
		self.customContextMenuRequested.connect(self.on_custom_menu)

	def sizeHint(self):
		return QSize(180, 100)

	@Slot(QPoint)
	def on_custom_menu(self, pos):
		self.current_index = self.indexAt(pos)

		#add_r_act = QAction("Add Receptors", self,
		#	triggered = self.add_receptors
		#)

		#add_l_act = QAction("Add Ligands", self,
		#	triggered = self.add_ligands
		#)

		rep_f_act = QAction("Specify Flexible Residues", self)
		rep_f_act.triggered.connect(self.set_flex_residules)
		rep_f_act.setEnabled(self.current_index.isValid() and self.current_index.siblingAtColumn(2).data() == 1)

		rep_a_act = QAction("Preset Active Binding Sites", self)
		rep_a_act.triggered.connect(self.set_active_sites)
		rep_a_act.setEnabled(self.current_index.isValid() and self.current_index.siblingAtColumn(2).data() == 1)

		lig_f_act = QAction("Filter and Remove Ligands", self)
		lig_f_act.triggered.connect(self.ligand_filter)
		lig_f_act.setEnabled(DB.active())

		del_m_act = QAction("Delete Current Molecule", self)
		del_m_act.triggered.connect(self.delete_molecular)
		del_m_act.setEnabled(self.current_index.isValid())

		clr_r_act = QAction("Delete All Receptors", self)
		clr_r_act.setEnabled(DB.active())
		clr_r_act.triggered.connect(self.delete_receptors)

		clr_l_act = QAction("Delete All Ligands", self)
		clr_l_act.triggered.connect(self.delete_ligands)
		clr_l_act.setEnabled(DB.active())

		clr_m_act = QAction("Delete All Molecules", self)
		clr_m_act.triggered.connect(self.delete_all)
		clr_m_act.setEnabled(DB.active())

		view_act = QAction("View Molecule Details", self)
		view_act.triggered.connect(self.view_details)
		view_act.setEnabled(self.current_index.isValid())

		stat_act = QAction("View Molecule Counts", self)
		stat_act.triggered.connect(self.view_stats)
		stat_act.setEnabled(DB.active())

		menu = QMenu(self)
		#menu.addAction(add_r_act)
		#menu.addAction(add_l_act)
		#menu.addAction(self.parent.import_receptor_act)
		#menu.addAction(self.parent.import_ligand_act)
		menu.addAction(rep_f_act)
		menu.addAction(rep_a_act)
		menu.addSeparator()
		menu.addAction(lig_f_act)
		menu.addSeparator()
		menu.addAction(del_m_act)
		menu.addAction(clr_r_act)
		menu.addAction(clr_l_act)
		menu.addAction(clr_m_act)
		menu.addSeparator()
		menu.addAction(view_act)
		menu.addSeparator()
		menu.addAction(stat_act)
		menu.popup(self.mapToGlobal(pos))

	@Slot()
	def add_receptors(self):
		self.parent.import_receptors()

	@Slot()
	def add_ligands(self):
		self.parent.import_ligands()

	@Slot()
	def set_flex_residules(self):
		mol_id = self.current_index.siblingAtColumn(0).data()
		self.parent.set_receptor_flexres(mol_id)

	@Slot()
	def set_active_sites(self):
		mol_id = self.current_index.siblingAtColumn(0).data()
		self.parent.set_active_sites(mol_id)

	@Slot()
	def ligand_filter(self):
		self.parent.filter_ligands()

	@Slot()
	def delete_molecular(self):
		ret = QMessageBox.question(self.parent, "Comfirmation",
			"Are you sure you want to delete select moleculue?")

		if ret == QMessageBox.No:
			return

		if not self.current_index.isValid():
			return

		mol_id = self.current_index.siblingAtColumn(0).data()
		mol_type = self.current_index.siblingAtColumn(2).data()

		if mol_type == 1:
			option = 'receptor_count'
		else:
			option = 'ligand_count'

		count = int(DB.get_option(option))
		count -= 1
		DB.set_option(option, count)

		total = int(DB.get_option('molecule_count'))
		total -= 1
		DB.set_option('molecule_count', total)

		self.parent.mol_model.remove(self.current_index.row())
		DB.query("DELETE FROM molecular WHERE id=?", (mol_id,))

	@Slot()
	def delete_ligands(self):
		ret = QMessageBox.question(self.parent, "Comfirmation",
			"Are you sure you want to delete all ligands?")

		if ret == QMessageBox.No:
			return

		DB.query("DELETE FROM molecular WHERE type=2")

		count = DB.get_option('ligand_count')

		if count:
			count = int(count)
			total = int(DB.get_option('molecule_count'))
			total -= count
			DB.set_option('ligand_count', 0)
			DB.set_option('molecule_count', total)
			self.parent.mol_model.select()

	@Slot()
	def delete_receptors(self):
		ret = QMessageBox.question(self.parent, "Comfirmation",
			"Are you sure you want to delete all receptors?")

		if ret == QMessageBox.No:
			return

		DB.query("DELETE FROM molecular WHERE type=1")
		count = DB.get_option('receptor_count')

		if count:
			count = int(count)
			total = int(DB.get_option('molecule_count'))
			total -= count
			DB.set_option('receptor_count', 0)
			DB.set_option('molecule_count', total)
			self.parent.mol_model.select()

	@Slot()
	def delete_all(self):
		ret = QMessageBox.question(self.parent, "Comfirmation",
			"Are you sure you want to delete all moleculues?")

		if ret == QMessageBox.No:
			return

		DB.query("DELETE FROM molecular")
		DB.set_option('receptor_count', 0)
		DB.set_option('ligand_count', 0)
		DB.set_option('molecule_count', 0)

		self.parent.mol_model.select()

	@Slot()
	def view_stats(self):
		rep_count = DB.get_option('receptor_count') or 0
		lig_count = DB.get_option('ligand_count') or 0
		total_count = DB.get_option('molecule_count') or 0

		info = (
			"<table cellspacing='10'>"
			"<tr><td>Number of Receptors: </td><td>{}</td></tr>"
			"<tr><td>Number of Ligands: </td><td>{}</td></tr>"
			"<tr><td>Total Number of Molecules: </td><td>{}</td></tr>"
			"</table>"
		)

		dlg = MoleculeDetailDialog(self.parent, info.format(
			rep_count, lig_count, total_count))
		dlg.exec()

	@Slot()
	def view_details(self):
		if not self.current_index.isValid():
			return

		mid = self.current_index.siblingAtColumn(0).data()
		sql = "SELECT * FROM molecular WHERE id=? LIMIT 1"
		mol = DB.get_dict(sql, (mid,))

		if mol:
			info = (
				"<table cellspacing='10'>"
				"<tr><td>{} name: </td><td>{}</td></tr>"
				"<tr><td>{} format: </td><td>{}</td></tr>"
				"<tr><td>Number of atoms: </td><td>{}</td></tr>"
				"<tr><td>Number of bonds: </td><td>{}</td></tr>"
				"<tr><td>Number of heavy atoms: </td><td>{}</td></tr>"
				"<tr><td>Number of residues: </td><td>{}</td></tr>"
				"<tr><td>Number of rotors: </td><td>{}</td></tr>"
				"<tr><td>Chemical formula: </td><td>{}</td></tr>"
				#"<tr><td>Heat of formation: </td><td>{} kcal/mol</td></tr>"
				"<tr><td>Molecular weight: </td><td>{}</td></tr>"
				"<tr><td>Calculated logP: </td><td>{}</td></tr>"
				"</table>"
			)

			dlg = MoleculeDetailDialog(self.parent, info.format(
				'Receptor' if mol.type == 1 else 'Ligand',
				mol.name,
				'Receptor' if mol.type == 1 else 'Ligand',
				mol.format,
				mol.atoms,
				mol.bonds,
				mol.hvyatoms,
				mol.residues,
				mol.rotors,
				self.format_formula(mol.formula),
				#mol.energy,
				round(mol.weight, 3),
				round(mol.logp, 3)
			))
			dlg.exec()

	def format_formula(self, formula):
		nums = []
		chars = []

		formulas = []

		for char in formula:
			if char.isdigit():
				nums.append(char)
			else:
				if chars:
					formulas.append("{}<sub>{}</sub>".format(
						''.join(chars),
						''.join(nums)
					))

					nums = []
					chars = []
				chars.append(char)

		if chars:
			formulas.append("{}<sub>{}</sub>".format(
				''.join(chars),
				''.join(nums)
			))

		return ''.join(formulas)

class DockeyTableModel(QAbstractTableModel):
	table = None
	custom_headers = []
	row_count = Signal(int)

	def __init__(self, parent=None):
		super(DockeyTableModel, self).__init__(parent)
		self.parent = parent

		#store ids of displayed row
		self.displayed = []

		#cache the current row
		self.cache_row = [-1, None]

		#total row counts
		self.total_count = 0

		#readed row counts
		self.read_count = 0

		#number of readed rows once time
		self._reads = 200

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

	def get_total(self):
		return DB.get_one(self.count_sql)

	def select(self):
		self.beginResetModel()
		self.read_count = 0
		self.total_count = self.get_total()
		self.displayed = DB.get_column(self.read_sql)
		self.read_count = len(self.displayed)
		self.cache_row = [-1, None]
		self.endResetModel()
		self.row_count.emit(self.total_count)

	def reset(self):
		self.beginResetModel()
		self.cache_row = [-1, None]
		self.read_count = 0
		self.displayed = []
		self.total_count = 0
		self.endResetModel()

	def clear(self):
		DB.query("DELETE FROM {}".format(self.table))
		self.reset()

	def remove(self, row):
		self.beginRemoveRows(QModelIndex(), row, row)
		self.displayed.pop(row)
		self.total_count -= 1
		self.read_count -= 1
		self.endRemoveRows()

class MolecularTableModel(DockeyTableModel):
	table = 'molecular'
	custom_headers = ['ID', 'Name', 'Type']

	def get_total(self):
		total = DB.get_option('molecule_count')

		if total:
			return int(total)
		else:
			return DB.get_one(self.count_sql)

	@property
	def read_sql(self):
		remainder = self.total_count - self.read_count
		fetch_count = min(self._reads, remainder)

		return "SELECT id FROM {} ORDER BY type LIMIT {},{}".format(
			self.table,
			self.read_count,
			fetch_count
		)

	def data(self, index, role):
		if role == Qt.DecorationRole:
			row = index.row()
			v = self.get_value(row, 2)

			if v == 1:
				return QIcon(":/icons/receptor.svg")

			elif v == 2:
				return QIcon(":/icons/ligand.svg")

		return super().data(index, role)

class JobsTableModel(DockeyTableModel):
	table = 'jobs'
	custom_headers = [' ID ', 'Receptor', 'Ligand', ' Status ', ' Progress ']

	statuses = {
		0: 'Failure',
		1: 'Success',
		2: 'Stopped',
		3: 'Running',
		4: 'Waiting'
		
	}
	status_colors = {
		0: QColor(220, 53, 69),
		1: QColor(25, 135, 84),
		2: QColor(188, 188, 188),
		3: QColor(13, 110, 253),
		4: QColor(255, 193, 7)
	}
	status_icons = {
		0: QIcon(':/icons/error.svg'),
		1: QIcon(':/icons/success.svg'),
		2: QIcon(':/icons/queue.svg'),
		3: QIcon(':/icons/run.svg'),
		4: QIcon(':/icons/pend.svg')
	}

	def get_total(self):
		total = DB.get_option('job_count')

		if total:
			return int(total)
		else:
			return DB.get_one(self.count_sql)

	@property
	def get_sql(self):
		return (
			"SELECT j.id,r.name,l.name,j.status,j.progress,"
			"j.started,j.finished,j.message FROM jobs AS j "
			"LEFT JOIN molecular AS r ON r.id=j.rid "
			"LEFT JOIN molecular AS l ON l.id=j.lid "
			"WHERE j.id=? LIMIT 1"
		)

	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if orientation == Qt.Horizontal and role == Qt.DisplayRole:
			return self.custom_headers[section]

		elif orientation == Qt.Vertical and role == Qt.DecorationRole:
			status = self.get_value(section, 3)
			return self.status_icons[status]

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()
		val = self.get_value(row, col)

		if role == Qt.DisplayRole:
			if col == 3:
				return self.statuses[val]

			elif col == 5 or col == 6:
				return time_format(val)

			return val

		elif role == Qt.ForegroundRole:
			if col == 3:
				return self.status_colors[val]

		#elif role == Qt.DecorationRole:
		#	if col == 3:
		#		return self.status_icons[val]

		#elif role == Qt.TextAlignmentRole:
		#	if col != 7:
		#		return Qt.AlignCenter

	@Slot(int)
	def update_row(self, rowid):
		if rowid not in self.displayed:
			return

		rowidx = rowid - 1
		colidx = len(self.custom_headers)-1
		self.update_cache(rowidx)
		index1 = self.index(rowidx, 3)
		index2 = self.index(rowidx, colidx)
		self.dataChanged.emit(index1, index2)
		self.headerDataChanged.emit(Qt.Vertical, rowidx, rowidx)

class JobsTableDelegate(QStyledItemDelegate):
	def paint(self, painter, option, index):
		percent = index.data() or 0
		#bar = QStyleOptionProgressBar(2)
		bar = QStyleOptionProgressBar()
		bar.rect = option.rect.adjusted(0, 5, 0, -5)
		bar.minimum = 0
		bar.maximum = 100
		bar.progress = int(percent)
		bar.text = "{}%".format(percent)
		bar.textVisible = True
		bar.textAlignment = Qt.AlignCenter
		bar.state |= QStyle.State_Horizontal
		painter.save()

		#fix incorrect progressbar position on MacOS
		if sys.platform == 'darwin':
			painter.translate(bar.rect.x(), bar.rect.y())
		
		QApplication.style().drawControl(QStyle.CE_ProgressBar, bar, painter)
		painter.restore()

class PoseTableModel(DockeyTableModel):
	table = 'pose'
	job = 0
	custom_headers = ['ID', 'Job']

	@property
	def count_sql(self):
		return "SELECT COUNT(1) FROM {} WHERE jid={} LIMIT 1".format(self.table, self.job)

	@property
	def get_sql(self):
		return "SELECT * FROM {} WHERE jid={} AND id=? LIMIT 1".format(self.table, self.job)

	@property
	def read_sql(self):
		remainder = self.total_count - self.read_count
		fetch_count = min(self._reads, remainder)

		return "SELECT id FROM {} WHERE jid={} LIMIT {},{}".format(
			self.table,
			self.job,
			self.read_count,
			fetch_count
		)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			return self.get_value(row, col)

		elif role == Qt.BackgroundRole:
			val = self.get_value(row, -1)

			if val:
				return QColor(169, 223, 191)

	#def data(self, index, role=Qt.DisplayRole):
	#	if not index.isValid():
	#		return None

	#	if role == Qt.TextAlignmentRole:
	#		return Qt.AlignCenter

	#	return super(PoseTableModel, self).data(index, role)

	def switch_table(self, tool=None):
		if tool == 'autodock':
			self.custom_headers = ['ID', 'Job', 'Run', 'Energy', 'cRMSD', 'rRMSD',
				'logKi', 'LE', 'FQ', 'SILE', 'LLE', 'LELP', 'Ki'
			]
		elif tool in ['vina', 'qvina']:
			self.custom_headers = ['ID', 'Job', 'Mode', 'Affinity', 'lRMSD', 'uRMSD',
				'logKi', ' LE', 'FQ', 'SILE', 'LLE', 'LELP', 'Ki'
			]
		else:
			self.custom_headers = ['ID', 'Job']

	def select(self):
		self.read_count = 0
		self.cache_row = [-1, None]

		self.beginResetModel()
		self.total_count = DB.get_one(self.count_sql)
		self.displayed = DB.get_column(self.read_sql)
		self.read_count = len(self.displayed)
		self.switch_table(DB.get_option('tool'))
		self.column_count = len(self.custom_headers)
		self.endResetModel()

		self.row_count.emit(self.total_count)

	def set_job(self, job_id):
		self.job = job_id
		self.select()

	def clear(self):
		DB.query("DELETE FROM {}".format(self.table))
		self.switch_table()
		self.reset()

class BestTableModel(PoseTableModel):
	table = 'best'
	sorter = ''
	custom_headers = ['', 'PID', 'Job']

	def get_total(self):
		total = DB.get_option('best_count')

		if total:
			return int(total)
		else:
			return DB.get_one(self.count_sql)

	@property
	def count_sql(self):
		return "SELECT COUNT(1) FROM {} LIMIT 1".format(self.table)

	@property
	def get_sql(self):
		return (
			"SELECT b.id,b.pid,p.jid,r.name,l.name,p.energy,p.rmsd1,"
			"p.rmsd2,p.logki,p.le,p.sile,p.fq,p.lle,p.lelp,p.ki,"
			"p.mode,p.complex,p.actives FROM best AS b "
			"LEFT JOIN pose AS p ON p.id=b.pid "
			"LEFT JOIN jobs AS j ON j.id=p.jid "
			"LEFT JOIN molecular AS r ON r.id=j.rid "
			"LEFT JOIN molecular AS l ON l.id=j.lid "
			"WHERE b.id=?"
		)

	@property
	def read_sql(self):
		remainder = self.total_count - self.read_count
		fetch_count = min(self._reads, remainder)

		if self.sorter:
			sql = "SELECT b.id FROM {} AS b LEFT JOIN pose \
			AS p ON p.id=b.pid ORDER BY {}".format(
				self.table,
				self.sorter
			)
		else:
			sql = "SELECT id FROM {} LIMIT {},{}".format(
				self.table,
				self.read_count,
				fetch_count
			)
		return sql

	def switch_table(self, tool=None):
		if tool == 'autodock':
			self.custom_headers = ['ID', 'PID', 'Job', 'Receptor', 'Ligand', 'Energy',
				'cRMSD', 'rRMSD','logKi', 'LE', 'SILE', 'FQ', 'LLE', 'LELP', 'Ki'
			]
		elif tool in ['vina', 'qvina']:
			self.custom_headers = ['ID', 'PID', 'Job', 'Receptor', 'Ligand', 'Affinity',
				'lRMSD', 'uRMSD', 'logKi', ' LE', 'SILE', 'FQ', 'LLE', 'LELP', 'Ki'
			]
		else:
			self.custom_headers = ['', 'PID', 'Job']

	def sort(self, column, order):
		fields = ['b.id', '', '', '', '', 'p.energy', 'p.rmsd1', 'p.rmsd2', 'p.logki', 'p.le',
			'p.sile', 'p.fq', 'p.lle', 'p.lelp'
		]

		if column in [1, 2, 3, 4, 14]:
			return

		if order == Qt.DescendingOrder:
			self.sorter = "{} DESC".format(fields[column])

		elif order == Qt.AscendingOrder:
			self.sorter = fields[column]

		else:
			self.sorter = ''

		self.select()


class BindingSiteModel(DockeyTableModel):
	table = 'binding_site'
	pose_id = 0
	custom_headers = ['ID', 'PID', 'Site']

	@property
	def count_sql(self):
		return "SELECT COUNT(1) FROM {} WHERE pid={} LIMIT 1".format(self.table, self.pose_id)

	@property
	def read_sql(self):
		remainder = self.total_count - self.read_count
		fetch_count = min(self._reads, remainder)

		return "SELECT id FROM {} WHERE pid={} LIMIT {},{}".format(
			self.table,
			self.pose_id,
			self.read_count,
			fetch_count
		)

	@property
	def all_sql(self):
		return "SELECT id FROM {} WHERE pid={}".format(self.table, self.pose_id)

	@property
	def get_sql(self):
		return "SELECT * FROM {} WHERE pid={} AND id=? LIMIT 1".format(self.table, self.pose_id)

	def change_pose(self, pose_id):
		self.pose_id = pose_id
		self.select()

class InteractionModel(DockeyTableModel):
	binding_site = 0
	custom_headers = ['ID', 'BID', 'Residue']

	@property
	def count_sql(self):
		return "SELECT COUNT(1) FROM {} WHERE bid={} LIMIT 1".format(self.table, self.binding_site)

	@property
	def read_sql(self):
		remainder = self.total_count - self.read_count
		fetch_count = min(self._reads, remainder)

		return "SELECT id FROM {} WHERE bid={} LIMIT {},{}".format(
			self.table,
			self.binding_site,
			self.read_count,
			fetch_count
		)

	@property
	def all_sql(self):
		return "SELECT id FROM {} WHERE bid={}".format(self.table, self.binding_site)

	@property
	def get_sql(self):
		return "SELECT * FROM {} WHERE bid={} AND id=? LIMIT 1".format(self.table, self.binding_site)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			return self.get_value(row, col)

		elif role == Qt.BackgroundRole:
			val = self.get_value(row, -1)

			if val:
				return QColor(169, 223, 191)

	def change_binding_site(self, site_id):
		self.binding_site = site_id
		self.select()

	def clear(self):
		DB.query("DELETE FROM {}".format(self.table))
		self.change_binding_site(0)

class HydrogenBondsModel(InteractionModel):
	table = 'hydrogen_bond'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance HA',
		'Distance DA', 'Donor Angle', 'Protein Donor', 'Side Chain', 'Donor Atom',
		'Acceptor Atom'
	]

class HalogenBondsModel(InteractionModel):
	table = 'halogen_bond'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance',
		'Donor Angle', 'Acceptor Angle', 'Donor Atom', 'Acceptor Atom'
	]

class HydrophobicInteractionModel(InteractionModel):
	table = 'hydrophobic_interaction'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance',
		'Ligand Atom', 'Protein Atom'
	]

class SaltBridgesModel(InteractionModel):
	table = 'salt_bridge'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance',
		'Protein Positive', 'Ligand Group', 'Ligand Atoms'
	]

class WaterBridgesModel(InteractionModel):
	table = 'water_bridge'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance AW',
		'Distance DW', 'Donor Angle', 'Water Angle', 'Protein Donor', 'Donor Atom',
		'Acceptor Atom', 'Water Atom'
	]

class PiStackingModel(InteractionModel):
	table = 'pi_stacking'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance',
		'Angle', 'Offset', 'Stacking Type', 'Ligand Atoms'
	]

class PiCationModel(InteractionModel):
	table = 'pi_cation'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Distance',
		'Offset', 'Protein Charged', 'Ligand Group', 'Ligand Atoms'
	]

class MetalComplexModel(InteractionModel):
	table = 'metal_complex'
	custom_headers = ['ID', 'Site', 'Chain', 'Residue', 'Amino Acid', 'Metal',
		'Target', 'Distance', 'Location'
	]

class JobTableView(QTableView):
	def __init__(self, parent=None):
		super(JobTableView, self).__init__(parent)
		self.parent = parent
		self.setContextMenuPolicy(Qt.CustomContextMenu)
		self.customContextMenuRequested.connect(self.on_custom_menu)

	def sizeHint(self):
		return QSize(300, 150)

	def get_pid(self, job):
		if self.parent.job_worker:
			return self.parent.job_worker.get_job_pid(job)

	@Slot(QPoint)
	def on_custom_menu(self, pos):
		self.current_index = self.indexAt(pos)

		view_detail_act = QAction("View Current Task", self)
		view_detail_act.triggered.connect(self.view_job_details)
		view_detail_act.setEnabled(self.current_index.isValid())

		view_status_act = QAction("View Task Counts", self)
		view_status_act.triggered.connect(self.view_job_statuses)
		view_status_act.setDisabled(not DB.active())

		#view_logs_act = QAction("View Log Files", self,
		#	triggered = self.view_logs
		#)
		#view_logs_act.setEnabled(self.current_index.isValid())
		stop_job_act = QAction("Stop Current Task", self)
		stop_job_act.triggered.connect(self.stop_select_job)
		stop_job_act.setEnabled(False)

		if self.current_index.isValid():
			status = self.current_index.siblingAtColumn(3).data(Qt.DisplayRole)

			if status == 'Running':
				stop_job_act.setEnabled(True)

		export_act = QAction("Export Task Table", self)
		export_act.triggered.connect(self.export_table)

		menu = QMenu(self)
		menu.addAction(view_detail_act)
		menu.addAction(stop_job_act)
		menu.addSeparator()
		menu.addAction(view_status_act)

		#menu.addAction(view_logs_act)
		menu.addSeparator()
		menu.addAction(export_act)
		menu.popup(self.viewport().mapToGlobal(pos))

	@Slot()
	def view_job_details(self):
		if not self.current_index.isValid():
			return

		jid = self.current_index.row() + 1
		sql = "SELECT * FROM jobs WHERE id=? LIMIT 1"
		job = DB.get_dict(sql, (jid,))
		sql = "SELECT name FROM molecular WHERE id=? LIMIT 1"
		receptor = DB.get_one(sql, (job.rid,))
		sql = "SELECT name FROM molecular WHERE id=? LIMIT 1"
		ligand = DB.get_one(sql, (job.lid,))

		info = (
			"<table cellspacing='10'>"
			"<tr><td>Receptor name: </td><td>{}</td></tr>"
			"<tr><td>Ligand name: </td><td>{}</td></tr>"
			"<tr><td>Job status: </td><td>{}</td></tr>"
			"<tr><td>Start time: </td><td>{}</td></tr>"
			"<tr><td>Finish time: </td><td>{}</td></tr>"
			"<tr><td>Elapsed time: </d><td>{}</td></tr>"
			"</table>"
		)

		statuses = [
			'<font color="red">Failure</font>',
			'<font color="green">Success</font>',
			'<font color="gray">Stopped</font>',
			'<font color="blue">Running</font>',
			'<font color="orange">Waiting</font>'
		]

		pid = self.get_pid(jid)

		dlg = JobDetailDialog(self.parent, info.format(
			receptor,
			ligand,
			statuses[job.status],
			time_format(job.started),
			time_format(job.finished),
			time_elapse(job.started, job.finished)
		), pid)
		dlg.exec()

	@Slot()
	def view_job_statuses(self):
		statuses = ['Failure', 'Success', 'Stopped', 'Running', 'Waiting']
		sql = "SELECT status, COUNT(1) from jobs group by status"

		total = 0
		rows = ["<table cellspacing='10'>"]
		for s, c in DB.query(sql):
			rows.append("<tr><td>{}:</td><td>{}</td></tr>".format(
				statuses[s], c
			))
			total += c
		rows.append("<tr><td>Total:</td><td>{}</td></tr>".format(total))
		info = ''.join(rows)

		dlg = JobStatusesDialog(self.parent, info)
		dlg.exec()

	@Slot()
	def stop_select_job(self):
		if not self.current_index.isValid():
			return

		jid = self.current_index.row() + 1
		self.parent.stop_job(jid)

	@Slot()
	def view_logs(self):
		if not self.current_index.isValid():
			return

		jid = self.current_index.row() + 1
		JobLogDialog.view_log(self.parent, jid)

	@Slot()
	def export_table(self):
		out, _ = QFileDialog.getSaveFileName(self.parent, filter="CSV file (*.csv)")

		if not out:
			return

		sql = self.model().get_sql.split('WHERE')[0]
		statuses = self.model().statuses

		with open(out, 'w') as fw:
			fw.write("ID,Receptor,Ligand,Status,Started,Finished\n")

			for row in DB.query(sql):
				fw.write("{},{},{},{},{},{}\n".format(
					row[0], row[1], row[2],
					statuses[row[3]],
					time_format(row[5]),
					time_format(row[6])
				)
			)

		self.parent.show_message("Export table to {}".format(out))


class PoseTableView(QTableView):
	def __init__(self, parent=None):
		super(PoseTableView, self).__init__(parent)
		self.parent = parent
		self.setContextMenuPolicy(Qt.CustomContextMenu)
		self.customContextMenuRequested.connect(self.on_custom_menu)

	def sizeHint(self):
		return QSize(300, 150)

	@Slot(QPoint)
	def on_custom_menu(self, pos):
		self.current_index = self.indexAt(pos)
		
		save_pose_act = QAction("Save Selected Pose", self)
		save_pose_act.triggered.connect(self.save_pose)
		save_pose_act.setDisabled(not self.current_index.isValid())
		save_all_act = QAction("Save All Poses", self)
		save_all_act.triggered.connect(self.save_all)

		save_complex_act = QAction("Save Receptor Ligand Complex", self)
		save_complex_act.triggered.connect(self.save_complex)
		save_complex_act.setDisabled(not self.current_index.isValid())

		export_act = QAction("Export Table", self)
		export_act.triggered.connect(self.export_table)

		view_act = QAction("View Details", self)
		view_act.triggered.connect(self.view_details)
		view_act.setDisabled(not self.current_index.isValid())

		menu = QMenu(self)
		
		menu.addAction(save_pose_act)
		menu.addAction(save_all_act)
		menu.addAction(save_complex_act)
		menu.addSeparator()
		menu.addAction(export_act)
		menu.addSeparator()
		menu.addAction(view_act)
		menu.popup(self.viewport().mapToGlobal(pos))

	def get_pose_ids(self):
		if not self.current_index.isValid():
			return None, None

		row = self.current_index.row()
		index = self.current_index.model().createIndex(row, 0)
		pid = index.data(role=Qt.DisplayRole)
		index = self.current_index.model().createIndex(row, 1)
		jid = index.data(role=Qt.DisplayRole)
		return pid, jid

	def export_table(self):
		out, _ = QFileDialog.getSaveFileName(self.parent, filter="CSV file (*.csv)")
		if not out:
			return

		_, jid = self.get_pose_ids()

		if not jid:
			return

		headers = self.model().custom_headers
		sql = "SELECT * FROM pose WHERE jid=?"

		with open(out, 'w') as fw:
			fw.write("{}\n".format(','.join(headers)))

			for row in DB.query(sql, (jid,)):
				fw.write("{}\n".format(','.join(map(str,row[0:len(headers)]))))

		self.parent.show_message("Export table to {}".format(out))

	def save_all(self):
		pdb, _ = QFileDialog.getSaveFileName(self.parent, filter="PDB file (*.pdb)")

		if not pdb:
			return

		_, jid = self.get_pose_ids()

		if not jid:
			return

		sql = "SELECT mode FROM pose WHERE jid=?"

		with open(pdb, 'w') as fw:
			for row in DB.query(sql, (jid,)):
				if row[0]:
					fw.write(row[0])
					fw.write('\n')

		self.parent.show_message("Save all poses to {}".format(pdb))

	def save_pose(self):
		pdb, _ = QFileDialog.getSaveFileName(self.parent, filter="PDB file (*.pdb)")

		if not pdb:
			return

		pid, _ = self.get_pose_ids()
		sql = "SELECT mode FROM pose WHERE id=? LIMIT 1"
		pose = DB.get_one(sql, (pid,))

		with open(pdb, 'w') as fw:
			fw.write(pose)

		self.parent.show_message("Save pose to {}".format(pdb))

	def save_complex(self):
		pdb_file, _ = QFileDialog.getSaveFileName(self.parent, filter="PDB file (*.pdb)")

		if not pdb_file:
			return

		pid, jid = self.get_pose_ids()
		sql = "SELECT complex FROM pose WHERE id=? LIMIT 1"
		content = DB.get_one(sql, (pid,))

		with open(pdb_file, 'w') as fw:
			fw.write(content)

	def view_details(self):
		pid, _ = self.get_pose_ids()
		sql = "SELECT * FROM pose WHERE id=? LIMIT 1"
		pose = DB.get_dict(sql, (pid,))

		tool = DB.get_option('tool')
		if tool == 'autodock':
			titles = ['Run', 'Free energy of binding', 'Cluster RMSD', 'Reference RMSD']

		elif tool in ['vina', 'qvina']:
			titles = ['Mode', 'Binding affinity', 'RMSD l.b.', 'RMSD u.b.']

		info = (
			"<table cellspacing='10'>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>Ki: </td><td>{}</td></tr>"
			"<tr><td>logKi: </td><td>{}</td></tr>"
			"<tr><td>Ligand efficiency (LE): </d><td>{}</td></tr>"
			"<tr><td>Size-independent ligand efficiency (SILE): </td><td>{}</td></tr>"
			"<tr><td>Fit Quality (FQ): </td><td>{}</td></tr>"
			"<tr><td>Lipophilic ligand efficiency (LLE): </td><td>{}</td></tr>"
			"<tr><td>Ligand efficiency lipophilic price (LELP): </td><td>{}</td></tr>"
			"</table>"
		)

		dock_info = info.format(
			titles[0],
			pose.run,
			titles[1],
			pose.energy,
			titles[2],
			pose.rmsd1,
			titles[3],
			pose.rmsd2,
			pose.ki,
			pose.logki,
			pose.le,
			pose.sile,
			pose.fq,
			pose.lle,
			pose.lelp
		)

		if pose.actives:
			rows = []
			for active in pose.actives.split(','):
				itype, chain, res, num = active.split(':')
				rows.append("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
					itype, chain, res, num
				))

			interact_info = """
			Interaction with preset active binding sites:
			<table>
			<tr><td>Interaction</td><td>Chain</td><td>Residue</td><td>Number</td></tr>
			{}
			</table>
			""".format('\n'.join(rows))
		else:
			interact_info = ""

		dlg = MoleculeDetailDialog(self.parent, "{}{}".format(dock_info, interact_info))
		dlg.exec()

class BestTableView(PoseTableView):
	def __init__(self, parent):
		super(BestTableView, self).__init__(parent)
		self.setSortingEnabled(True)

	@Slot(QPoint)
	def on_custom_menu(self, pos):
		self.current_index = self.indexAt(pos)
		
		save_pose_act = QAction("Save Selected Pose", self)
		save_pose_act.triggered.connect(self.save_pose)
		save_pose_act.setDisabled(not self.current_index.isValid())
		save_all_act = QAction("Save All Poses", self)
		save_all_act.triggered.connect(self.save_all)

		save_complex_act = QAction("Save Receptor Ligand Complex", self)
		save_complex_act.triggered.connect(self.save_complex)
		save_complex_act.setDisabled(not self.current_index.isValid())

		export_act = QAction("Export Table", self)
		export_act.triggered.connect(self.export_table)

		view_act = QAction("View Details", self)
		view_act.triggered.connect(self.view_details)
		view_act.setDisabled(not self.current_index.isValid())

		menu = QMenu(self)
		
		menu.addAction(save_pose_act)
		menu.addAction(save_all_act)
		menu.addAction(save_complex_act)
		menu.addSeparator()
		menu.addAction(export_act)
		menu.addSeparator()
		menu.addAction(view_act)
		menu.popup(self.viewport().mapToGlobal(pos))

	def get_pose_ids(self):
		if not self.current_index.isValid():
			return None, None

		row = self.current_index.row()
		index = self.current_index.model().createIndex(row, 1)
		pid = index.data(role=Qt.DisplayRole)
		index = self.current_index.model().createIndex(row, 2)
		jid = index.data(role=Qt.DisplayRole)
		return pid, jid

	def export_table(self):
		out, _ = QFileDialog.getSaveFileName(self.parent, filter="CSV file (*.csv)")

		if not out:
			return

		headers = self.model().custom_headers
		sql = self.model().get_sql.split('WHERE')[0]

		with open(out, 'w') as fw:
			fw.write("{}\n".format(','.join(headers)))

			for row in DB.query(sql):
				fw.write("{}\n".format(','.join(map(str,row[0:len(headers)]))))

		self.parent.show_message("Export table to {}".format(out))

	def save_all(self):
		pdb, _ = QFileDialog.getSaveFileName(self.parent, filter="PDB file (*.pdb)")

		if not pdb:
			return

		sql = "SELECT mode FROM pose,best WHERE pose.id=best.pid"

		with open(pdb, 'w') as fw:
			for row in DB.query(sql):
				if row[0]:
					fw.write(row[0])
					fw.write('\n')

		self.parent.show_message("Save all poses to {}".format(pdb))

	def save_pose(self):
		pdb, _ = QFileDialog.getSaveFileName(self.parent, filter="PDB file (*.pdb)")

		if not pdb:
			return

		pid, _ = self.get_pose_ids()
		sql = "SELECT mode FROM pose WHERE id=? LIMIT 1"
		pose = DB.get_one(sql, (pid,))

		with open(pdb, 'w') as fw:
			fw.write(pose)

		self.parent.show_message("Save pose to {}".format(pdb))

	def save_complex(self):
		pdb_file, _ = QFileDialog.getSaveFileName(self.parent, filter="PDB file (*.pdb)")

		if not pdb_file:
			return

		pid, jid = self.get_pose_ids()
		sql = "SELECT complex FROM pose WHERE id=? LIMIT 1"
		content = DB.get_one(sql, (pid,))

		with open(pdb_file, 'w') as fw:
			fw.write(content)

	def view_details(self):
		pid, _ = self.get_pose_ids()
		sql = "SELECT * FROM pose WHERE id=? LIMIT 1"
		pose = DB.get_dict(sql, (pid,))

		tool = DB.get_option('tool')
		if tool == 'autodock':
			titles = ['Run', 'Free energy of binding', 'Cluster RMSD', 'Reference RMSD']

		elif tool in ['vina', 'qvina']:
			titles = ['Mode', 'Binding affinity', 'RMSD l.b.', 'RMSD u.b.']

		info = (
			"<table cellspacing='10'>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>{}: </td><td>{}</td></tr>"
			"<tr><td>Ki: </td><td>{}</td></tr>"
			"<tr><td>logKi: </td><td>{}</td></tr>"
			"<tr><td>Ligand efficiency (LE): </d><td>{}</td></tr>"
			"<tr><td>Size-independent ligand efficiency (SILE): </td><td>{}</td></tr>"
			"<tr><td>Fit Quality (FQ): </td><td>{}</td></tr>"
			"<tr><td>Lipophilic ligand efficiency (LLE): </td><td>{}</td></tr>"
			"<tr><td>Ligand efficiency lipophilic price (LELP): </td><td>{}</td></tr>"
			"</table>"
		)

		dlg = MoleculeDetailDialog(self.parent, info.format(
			titles[0],
			pose.run,
			titles[1],
			pose.energy,
			titles[2],
			pose.rmsd1,
			titles[3],
			pose.rmsd2,
			pose.ki,
			pose.logki,
			pose.le,
			pose.sile,
			pose.fq,
			pose.lle,
			pose.lelp
		))
		dlg.exec()

class DockeyTextBrowser(QTextBrowser):
	def sizeHint(self):
		return QSize(300, 150)

class InteractionTableView(QTableView):
	def __init__(self, parent):
		super(InteractionTableView, self).__init__(parent)
		#self.verticalHeader().hide()
		self.parent = parent
		self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
		self.setSelectionBehavior(QAbstractItemView.SelectRows)
		self.setEditTriggers(QAbstractItemView.NoEditTriggers)
		self.doubleClicked.connect(self.on_double_clicked)
		self.setContextMenuPolicy(Qt.CustomContextMenu)
		self.customContextMenuRequested.connect(self.on_custom_menu)

	def sizeHint(self):
		return QSize(300, 150)

	def on_double_clicked(self, index):
		chain = index.siblingAtColumn(2).data()
		residue = index.siblingAtColumn(3).data()
		self.parent.cmd.zoom('{}/{}/'.format(chain, residue), animate=0.5)

	def on_custom_menu(self, pos):
		self.table_models = {
			'hydrogen_bond': HydrogenBondsModel, 
			'halogen_bond': HalogenBondsModel,
			'hydrophobic_interaction': HydrophobicInteractionModel,
			'salt_bridge': SaltBridgesModel,
			'water_bridge': WaterBridgesModel,
			'pi_stacking': PiStackingModel,
			'pi_cation': PiCationModel,
			'metal_complex': MetalComplexModel
		}

		exp_cur_act = QAction("Export Current Table", self)
		exp_cur_act.triggered.connect(self.export_current_table)
		exp_all_act = QAction("Export All Tables", self)
		exp_all_act.triggered.connect(self.export_all_tables)
		exp_best_act = QAction("Export Interactions for Best Poses", self)
		exp_best_act.triggered.connect(self.export_best_poses_interactions)
		exp_pose_act = QAction("Export Interactions for All Poses", self)
		exp_pose_act.triggered.connect(self.export_all_poses_interactions)

		menu = QMenu(self)
		menu.addAction(exp_cur_act)
		menu.addAction(exp_all_act)
		menu.addSeparator()
		menu.addAction(exp_best_act)
		menu.addAction(exp_pose_act)
		menu.popup(self.mapToGlobal(pos))

	def export_current_table(self):
		out, _ = QFileDialog.getSaveFileName(self.parent, filter="CSV file (*.csv)")

		if not out:
			return

		table = self.model().table
		bid = self.model().binding_site
		sql = "SELECT site FROM binding_site WHERE id=? LIMIT 1"
		site = DB.get_one(sql, (bid,))
		sql = "SELECT * FROM {} WHERE bid=?".format(table)

		with open(out, 'w', newline='') as fw:
			writer = csv.writer(fw)
			writer.writerow(self.model().custom_headers)
			for row in DB.query(sql, (bid,)):
				row = list(row)
				row[1] = site
				writer.writerow(row)

		self.parent.show_popup_message("Export table to {}".format(out))

	def export_all_tables(self):
		out_dir = QFileDialog.getExistingDirectory(self.parent)

		if not out_dir:
			return

		bid = self.model().binding_site
		sql = "SELECT site FROM binding_site WHERE id=? LIMIT 1"
		site = DB.get_one(sql, (bid,))

		for table in self.table_models:
			rows = []
			sql = "SELECT * FROM {} WHERE bid={}".format(table, bid)
			headers = self.table_models[table].custom_headers

			for row in DB.query(sql):
				row = list(row)
				row[1] = site
				rows.append(row)

			if rows:
				out_file = os.path.join(out_dir, '{}.csv'.format(table))
				with open(out_file, 'w', newline='') as fw:
					writer = csv.writer(fw)
					writer.writerow(headers)
					writer.writerows(rows)

		self.parent.show_popup_message("Successfully export all tables")

	def export_best_poses_interactions(self):
		out_dir = QFileDialog.getExistingDirectory(self.parent)

		if not out_dir:
			return

		threader = BestInteractionExportWorker(out_dir, self.table_models)
		threader.signals.message.connect(self.parent.show_message)
		message = "Successfully export interactions for best poses"
		threader.signals.finished.connect(lambda: self.parent.show_popup_message(message))
		QThreadPool.globalInstance().start(threader)

	def export_all_poses_interactions(self):
		out_dir = QFileDialog.getExistingDirectory(self.parent)

		if not out_dir:
			return

		threader = PoseInteractionExportWorker(out_dir, self.table_models)
		threader.signals.message.connect(self.parent.show_message)
		message = "Successfully export interactions for all poses"
		threader.signals.finished.connect(lambda: self.parent.show_popup_message(message))
		QThreadPool.globalInstance().start(threader)


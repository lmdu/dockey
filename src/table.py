from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from backend import *

__all__ = ['DockeyListView', 'DockeyTableView', 'DockeyTextBrowser',
	'DockeyTableModel', 'MolecularTableModel', 'JobsTableModel',
	'JobsTableDelegate', 'PoseTableModel'
]

class MoleculeDetailDialog(QDialog):
	def __init__(self, parent, info):
		super(MoleculeDetailDialog, self).__init__(parent)
		self.setWindowTitle("Molecular Details")
		layout = QVBoxLayout()
		self.setLayout(layout)

		label = QLabel(info, self)
		layout.addWidget(label)

		btn_box = QDialogButtonBox(QDialogButtonBox.Ok)
		btn_box.accepted.connect(self.accept)
		btn_box.rejected.connect(self.reject)
		layout.addWidget(btn_box)

class DockeyListView(QListView):
	def __init__(self, parent=None):
		super(DockeyListView, self).__init__(parent)
		self.parent = parent
		self.setContextMenuPolicy(Qt.CustomContextMenu)
		self.customContextMenuRequested.connect(self.on_custom_menu)

	def sizeHint(self):
		return QSize(180, 100)

	@Slot()
	def on_custom_menu(self, pos):
		self.current_index = self.indexAt(pos)

		menu = QMenu(self)

		add_r_act = QAction("Add Receptors", self,
			triggered = self.add_receptors
		)
		menu.addAction(add_r_act)
		add_l_act = QAction("Add Ligands", self,
			triggered = self.add_ligands
		)
		menu.addAction(add_l_act)

		menu.addSeparator()

		del_m_act = QAction("Delete", self,
			triggered = self.delete_molecular
		)
		
		if not self.current_index.isValid():
			del_m_act.setDisabled(True)
		
		menu.addAction(del_m_act)
		
		clr_m_act = QAction("Delete All", self,
			triggered = self.delete_all
		)
		menu.addAction(clr_m_act)

		view_act = QAction("View Details", self,
			triggered = self.view_details
		)
		menu.addAction(view_act)

		if not self.current_index.isValid():
			view_act.setDisabled(True)

		action = menu.exec(self.mapToGlobal(pos))

	@Slot()
	def add_receptors(self):
		self.parent.import_receptors()

	@Slot()
	def add_ligands(self):
		self.parent.import_ligands()

	@Slot()
	def delete_molecular(self):
		if not self.current_index.isValid():
			return

		name = self.current_index.data(role=Qt.DisplayRole)
		DB.query("DELETE FROM molecular WHERE name=?",(name,))
		self.parent.mol_model.select()

	@Slot()
	def delete_all(self):
		DB.query("DELETE FROM molecular")
		self.parent.mol_model.select()

	@Slot()
	def view_details(self):
		if not self.current_index.isValid():
			return

		name = self.current_index.data(role=Qt.DisplayRole)
		sql = "SELECT * FROM molecular WHERE name=? LIMIT 1"
		mol = DB.get_dict(sql, (name,))

		if mol:
			info = (
				"<table cellspacing='10'>"
				"<tr><td>{} name: </td><td>{}</td></tr>"
				"<tr><td>Number of atoms: </td><td>{}</td></tr>"
				"<tr><td>Number of bonds: </td><td>{}</td></tr>"
				"<tr><td>Number of heavy atoms: </td><td>{}</td></tr>"
				"<tr><td>Number of residues: </td><td>{}</td></tr>"
				"<tr><td>Number of rotors: </td><td>{}</td></tr>"
				"<tr><td>Stochoimetric formula: </td><td>{}</td></tr>"
				"<tr><td>Heat of formation: </td><td>{} kcal/mol</td></tr>"
				"<tr><td>Molecular Weight: </td><td>{}</td></tr>"
				"<tr><td>Calculated logP: </td><td>{}</td></tr>"
				"</table>"
			)

			dlg = MoleculeDetailDialog(self.parent, info.format(
				'Receptor' if mol.type == 1 else 'Ligand',
				mol.name,
				mol.atoms,
				mol.bonds,
				mol.hvyatoms,
				mol.residues,
				mol.rotors,
				mol.formula,
				mol.energy,
				round(mol.weight, 3),
				round(mol.logp, 3)
			))
			dlg.exec()


class DockeyTableView(QTableView):
	def sizeHint(self):
		return QSize(300, 150)

class DockeyTextBrowser(QTextBrowser):
	def sizeHint(self):
		return QSize(300, 150)

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
	custom_headers = [' ID ', 'Receptor', 'Ligand', ' Status ', ' Progress ']

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
				return self.status[val]

			elif col == 5 or col == 6:
				return time_format(val)

			return val

		elif role == Qt.ForegroundRole:
			if col == 3:
				return self.status_colors[val]

		#elif role == Qt.DecorationRole:
		#	if col == 3:
		#		return self.status_icons[val]

		elif role == Qt.TextAlignmentRole:
			if col != 7:
				return Qt.AlignCenter

	@Slot()
	def update_row(self, rowid):
		rowid -= 1
		colid = len(self.custom_headers)-1
		self.update_cache(rowid)
		index1 = self.index(rowid, 0)
		index2 = self.index(rowid, colid)
		self.dataChanged.emit(index1, index2)
		self.headerDataChanged.emit(Qt.Vertical, rowid, rowid)

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

class PoseTableModel(DockeyTableModel):
	table = 'pose'
	job = 0
	tool = None
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

	def switch_table(self, tool):
		if tool == 'autodock4':
			self.custom_headers = ['ID', 'Job', 'Run', 'Energy', 'Cluster RMSD', 'Reference RMSD',
				'logKi', 'LE', 'LLE', 'FQ', 'LELP'
			]
		elif tool == 'vina':
			self.custom_headers = ['ID', 'Job', 'Mode', 'Affinity', 'RMSD l.b.', 'RMSD u.b.',
				'logKi', ' Ligand efficiency', 'Lipophilic\nligand efficiency', 'Fit quality', 'Ligand efficiency\nlipophilic price'
			]

	def set_job(self, job_id):
		if self.tool is None:
			self.tool = DB.get_option('tool')

		self.switch_table(self.tool)
		self.job = job_id
		self.select()

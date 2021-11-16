import os

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

__all__ = ['AutogridParameter', 'AutodockParameterDialog']

class AutogridParameter:
	comments = {
		'npts': "# num.grid points in xyz",
		'gridfld': "# grid_data_file",
		'spacing': "# spacing(A)",
		'receptor_types': "# receptor atom types",
		'ligand_types': "# ligand atom types",
		'receptor': "# macromolecule",
		'gridcenter': "# xyz-coordinates or auto",
		'smooth': "# store minimum energy w/in rad(A)",
		'map': "# atom-specific affinity map",
		'elecmap': "# electrostatic potential map",
		'dsolvmap': "# desolvation potential map",
		'dielectric': "# <0, AD4 distance-dep.diel;>0, constant"
	}

	def __init__(self, receptor, ligand, size, center, spacing=0.375,
				 smooth=0.5, dielectric=-0.1465):
		self.receptor = receptor
		self.ligand = ligand
		self.size = size
		self.center = center
		self.spacing = spacing
		self.smooth = smooth
		self.dielectric = dielectric

		self.receptor_name = os.path.splitext(self.receptor)[0]
		self.receptor_types = self.get_atom_types(self.receptor)
		self.ligand_types = self.get_atom_types(self.ligand)
		self.gpf_file = self.receptor.replace('.pdbqt', '.gpf')

	def get_atom_types(self, pdbqt_file):
		atom_types = set()
		with open(pdbqt_file) as fh:
			for line in fh:
				if not line.startswith(('ATOM', 'HETATM')):
					continue

				atom_types.add(line[78:80].strip())

		return sorted(atom_types)

	def make_gpf_file(self):
		self.rows = [
			"npts {0[0]} {0[1]} {0[2]}".format(self.size),
			"gridfld {}.maps.fld".format(self.receptor_name),
			"spacing {}".format(self.spacing),
			"receptor_types {}".format(' '.join(self.receptor_types)),
			"ligand_types {}".format(' '.join(self.ligand_types)),
			"receptor {}".format(self.receptor),
			"gridcenter {0[0]} {0[1]} {0[2]}".format(self.center),
			"smooth {}".format(self.smooth),
		]

		for _type in self.ligand_types:
			self.rows.append("map {}.{}.map".format(self.receptor_name, _type))

		self.rows.append("elecmap {}.e.map".format(self.receptor_name))
		self.rows.append("dsolvmap {}.d.map".format(self.receptor_name))
		self.rows.append("dielectric {}".format(self.dielectric))

		max_len = max([len(row) for row in self.rows]) + 6

		with open(self.gpf_file, 'w') as fw:
			for row in self.rows:
				field = row.split()[0]
				fw.write("{:<{}}{}".format(row, max_len, self.comments[field]))

class AutodockParameter:
	def __init__(self):
		self.num = 0
		self.params = {
			'autodock_parameter_version': {
				'type': str,
				'default': '4.2',
				'value': '4.2',
				'comment': 'used by autodock to validate parameter set',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'outlev': {
				'type': int,
				'default': 1,
				'value': 1,
				'comment': 'diagnostic output level',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'parameter_file': {
				'type': str,
				'default': 'AD4.1_bound.dat',
				'value': '',
				'comment': 'parameter library filename',
				'scope': 'global',
				'required': False,
				'user': True,
				'order': self.order
			},
			'intelec': {
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'calculate internal electrostatics',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'intnbp_r_eps': {
				'type': [float,float,int,int,str,str],
				'default': [],
				'value': [],
				'comment': 'internal energy potential for a given class of interactions',
				'required': False,
				'user': False,
				'order': self.order
			},
			'torsdof': {
				'type': int,
				'default': 0,
				'value': 0,
				'comment': 'torsional degrees of freedom',
				'required': True,
				'user': False,
				'order': self.order
			},
			'seed': {
				'type': list,
				'default': ['pid', 'time'],
				'value': ['pid', 'time'],
				'comment': 'seeds for random generator',
				'scope': 'global',
				'required': True,
				'user': True,
				'order': self.order
			},
			'ligand_types': {
				'type': list,
				'default': [],
				'value': [],
				'comment': 'atoms types in ligand',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'fld': {
				'type': str,
				'default': '',
				'value': '',
				'comment': 'grid_data_file',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'map': {
				'type': iter,
				'default': [],
				'value': [],
				'comment': 'atom-specific affinity map',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'elecmap': {
				'type': str,
				'default': '',
				'value': '',
				'comment': 'electrostatics map',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'desolvmap': {
				'type': str,
				'default': '',
				'value': '',
				'comment': 'desolvation map',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'move': {
				'type': str,
				'default': '',
				'value': '',
				'comment': 'small molecule',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'flexres': {
				'type': str,
				'default': '',
				'value': '',
				'comment': 'file containing flexible residues',
				'scope': 'global',
				'required': False,
				'user': True,
				'order': self.order,
			},
			'about': {
				'type': list,
				'default': [],
				'value': [],
				'comment': 'small molecule center',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'tran0': {
				'type': list,
				'default': ['random'],
				'value': ['random'],
				'comment': 'initial coordinates/A or random',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'quaternion0': {
				'type': list,
				'default': ['random'],
				'value': ['random'],
				'comment': 'initial orientation',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'dihe0': {
				'type': list,
				'default': ['random'],
				'value': ['random'],
				'comment': 'initial dihedrals (relative) or random',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'tstep': {
				'type': float,
				'default': 0.2,
				'value': 0.2,
				'comment': 'translation step/A',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'qstep': {
				'type': float,
				'default': 5.0,
				'value': 5.0,
				'comment': 'quaternion step/deg',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'dstep': {
				'type': float,
				'default': 5.0,
				'value': 5.0,
				'comment': 'torsion step/deg',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'rmstol': {
				'type': float,
				'default': 2.0,
				'value': 2.0,
				'comment': 'cluster_tolerance/A',
				'scope': 'global',
				'required': True,
				'user': True,
				'order': self.order
			},
			'epdb': {
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'evaluate ligand specified with move command',
				'scope': 'global',
				'required': True,
				'user': True,
				'order': self.order
			},
			'e0max': {
				'type': [float, int],
				'default': [0.0, 10000],
				'value': [0.0, 10000],
				'comment': 'max initial energy; max number of retries',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'rt0': {
				'type': float,
				'default': 616.0,
				'value': 616.0,
				'comment': 'initial annealing temperature (times gas constant)',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'line_schedule': {
				'type': bool,
				'default': False,
				'value': True,
				'comment': 'use linear, arithmetic temperature reduction',
				'scope': ['SA'],
				'required': False,
				'user': True,
				'order': self.order
			},
			'geometric_schedule': {
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'use geometric, arithmetic temperature reduction',
				'scope': ['SA'],
				'required': False,
				'user': True,
				'order': self.order
			},
			'rtrf': {
				'type': float,
				'default': 0.95,
				'value': 0.95,
				'comment': 'annealing temperature reduction factor',
				'scope': ['SA'],
				'required': False,
				'user': True,
				'order': self.order
			},
			'runs': {
				'type': int,
				'default': 10,
				'value': 10,
				'comment': 'number of automated docking runs',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'cycles': {
				'type': int,
				'default': 50,
				'value': 50,
				'comment': 'number of temperature reduction cycles',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'accs': {
				'type': int,
				'default': 25000,
				'value': 25000,
				'comment': 'maximum number of accepted steps per cycle',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'rejs': {
				'type': int,
				'default': 25000,
				'value': 25000,
				'comment': 'maximum number of rejected steps per cycle',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'rejs': {
				'type': int,
				'default': 25000,
				'value': 25000,
				'comment': 'maximum number of rejected steps per cycle',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'select': {
				'type': str,
				'default': 'm',
				'value': 'm',
				'comment': 'state selection flag: (m)inimum or (l)ast state',
				'choices': ['m', 'l'],
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'trnrf': {
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'comment': 'per cycle reduction factor for translation',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'quarf': {
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'comment': 'per cycle reduction factor for quaternions',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'dihrf': {
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'comment': 'per cycle reduction factor for dihedrals',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_pop_size': {
				'type': int,
				'default': 150,
				'value': 150,
				'comment': 'number of individuals in population',
				'scope': ['GA', 'LGA', 'LS'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_num_evals': {
				'type': int,
				'default': 2500000,
				'value': 2500000,
				'comment': 'maximum number of energy evaluations',
				'scope': ['GA', 'LGA', 'LS'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_num_generations': {
				'type': int,
				'default': 27000,
				'value': 27000,
				'comment': 'maximum number of energy evaluations',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_elitism': {
				'type': int,
				'default': 1,
				'value': 1,
				'comment': 'number of top individuals to survive to next generation',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_mutation_rate': {
				'type': float,
				'default': 0.02,
				'value': 0.02,
				'comment': 'rate of gene mutation',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_crossover_rate': {
				'type': float,
				'default': 0.8,
				'value': 0.8,
				'comment': 'rate of crossover',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_window_size': {
				'type': int,
				'default': 10,
				'value': 10,
				'comment': 'number of preceding generations',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ga_cauchy_alpha': {
				'type': float,
				'default': 0.0,
				'value': 0.0,
				'comment': 'Alpha parameter of Cauchy distribution',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'ga_cauchy_beta': {
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'comment': 'Beta parameter Cauchy distribution',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'set_ga': {
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'set the above parameters for GA or LGA',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'set_ga': {
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'set the above parameters for GA or LGA',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'sw_max_its': {
				'type': int,
				'default': 300,
				'value': 300,
				'comment': 'iterations of Solis & Wets local search',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'sw_max_succ': {
				'type': int,
				'default': 4,
				'value': 4,
				'comment': 'consecutive successes before changing rho',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'sw_max_fail': {
				'type': int,
				'default': 4,
				'value': 4,
				'comment': 'consecutive failures before changing rho',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'sw_rho': {
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'comment': 'size of local search space to sample',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'sw_lb_rho': {
				'type': float,
				'default': 0.01,
				'value': 0.01,
				'comment': 'lower bound on rho',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'ls_search_freq': {
				'type': float,
				'default': 0.06,
				'value': 0.06,
				'comment': 'lower bound on rho',
				'scope': ['LGA'],
				'required': True,
				'user': True,
				'order': self.order
			},
			'set_sw1': {
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'set the above classical Solis & Wets parameters',
				'scope': ['LGA', 'LS'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'set_psw1': {
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'set the above pseudo-Solis & Wets parameters',
				'scope': ['LGA', 'LS'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'unbound_model': {
				'type': str,
				'default': 'bound',
				'value': 'bound',
				'comment': 'state of unbound ligand',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			},
			'simanneal': {
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'do as many SA runs as set by runs keyword above',
				'scope': ['SA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'do_local_only': {
				'type': int,
				'default': 50,
				'value': 50,
				'comment': 'do this many LS runs',
				'scope': ['LS'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'do_global_only': {
				'type': int,
				'default': 50,
				'value': 50,
				'comment': 'do this many GA runs',
				'scope': ['GA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'ga_run': {
				'type': int,
				'default': 50,
				'value': 50,
				'comment': 'do this many hybrid GA-LS runs',
				'scope': ['LGA'],
				'required': True,
				'user': False,
				'order': self.order
			},
			'analysis': {
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'perform a ranked cluster analysis',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}
		}

	@property
	def order(self):
		self.num += 1
		return self.num

	def set_value(self, k, v):
		self.params[k]['value'] = v

	def get_value(self, k):
		return self.params[k]['value']

	def get_param(self, k):
		return self.params[k]

	def get_count(self, algorithm='LGA'):
		count = 0
		for param in self.params:
			if param['user'] and (param['scope'] == 'global' or algorithm in param['scope']):
				count += 1

		return count

	def get_scope_items(self, algorithm='LGA'):
		params = {}
		for k in self.params:
			if self.params[k]['user'] and algorithm in self.params[k]['scope']:
				params[k] = self.params[k]

		return sorted(params.items(), key=lambda x: x[1]['order'])

	def get_ordered_items(self, algorithm='LGA'):
		params = {}
		for k in self.params:
			if self.params[k]['user'] and (self.params[k]['scope'] == 'global' or algorithm in self.params[k]['scope']):
				params[k] = self.params[k]

		return sorted(params.items(), key=lambda x: x[1]['order'])

'''
class ParameterTreeModel(QAbstractItemModel):
	def __init__(self, parent, params):
		super(ParameterTreeModel, self).__init__(parent)
		self.params = params
		self.items = self.params.get_ordered_items()

	def index(self, row, column, parent):
		return self.createIndex(row, column)

	def parent(self, child):
		return QModelIndex()

	def rowCount(self, parent):
		return 1

	def columnCount(self, parent):
		return 3

	def data(self, index, role):
		if not index.isValid():
			return None

		if role == Qt.DisplayRole:
			return 1
'''

class ParameterTreeModel(QStandardItemModel):
	def __init__(self, parent, params):
		super(ParameterTreeModel, self).__init__(parent)
		self.params = params

	def data(self, index, role):
		param = self.params.get_param(self.get_key(index))
		column = index.column()

		if role == Qt.DisplayRole:
			if column < 2:
				return super(ParameterTreeModel, self).data(index, role)
			else:
				if param['type'] is not bool:
					return param['value']

		elif role == Qt.CheckStateRole:
			if column == 2 and param['type'] is bool:
				if param['value']:
					return Qt.Checked
				else:
					return Qt.Unchecked


	def get_key(self, index):
		item = self.item(index.row())
		return item.text()

class ParameterTreeDelegate(QStyledItemDelegate):
	def __init__(self, parent, params):
		super(ParameterTreeDelegate, self).__init__(parent)
		self.parent = parent
		self.params = params

	def paint(self, painter, option, index):
		QStyledItemDelegate.paint(self, painter, option, index)

	def createEditor(self, parent, option, index):
		col = index.column()
		row = index.row()

		if col == 2:
			k = self.parent.model.get_key(index)
			i = self.params.params[k]
			t = i['type']

			if t is float:
				e = QDoubleSpinBox(parent)
				e.setMaximum(100000000)
				return e
			elif t is int:
				e = QSpinBox(parent)
				e.setMaximum(100000000)
				return e
			elif 'choices' in i:
				e = QComboBox(parent)
				e.addItems(i['choices'])
				return e
			else:
				return QLineEdit(parent)
		
		return QStyledItemDelegate.createEditor(self, parent, option, index)

	def setEditorData(self, editor, index):
		col = index.column()
		row = index.row()

		if col == 2:
			k = self.parent.model.get_key(index)
			i = self.params.params[k]
			t = i['type']
			v = i['value']

			if t is float or t is int:
				editor.setValue(v)

			elif 'choices' in i:
				idx = editor.findText(v)
				editor.setCurrentIndex(idx)
			else:
				editor.setText(v)

		#return QStyledItemDelegate.setEditorData(self, editor, index)

	def setModelData(self, editor, model, index):
		QStyledItemDelegate.setModelData(self, editor, model, index)

class ParameterTreeView(QTreeView):
	def __init__(self, parent, params):
		super(ParameterTreeView, self).__init__(parent)
		self.params = params
		self.delegate = ParameterTreeDelegate(self, self.params)
		self.model = ParameterTreeModel(self, self.params)
		self.model.setHorizontalHeaderLabels(["Command", "Comment", "Value"])
		self.root = self.model.invisibleRootItem()
		self.setModel(self.model)
		self.setItemDelegate(self.delegate)

	def add_item(self, *args):
		item = [QStandardItem(str(it)) for it in args]
		self.root.appendRow(item)
		return item[0]

	def clear(self):
		for i in range(self.root.rowCount()):
			self.root.removeRow(i)

	def sizeHint(self):
		return QSize(500, 200)

class AutodockParameterDialog(QDialog):
	def __init__(self, parent):
		super(AutodockParameterDialog, self).__init__(parent)
		self.params = AutodockParameter()

		self.dock_tree = ParameterTreeView(self, self.params)
		self.dock_tree.header().setStretchLastSection(False)
		self.dock_tree.header().setSectionResizeMode(0, QHeaderView.ResizeToContents)
		self.dock_tree.header().setSectionResizeMode(1, QHeaderView.Stretch)
		self.dock_tree.header().setSectionResizeMode(2, QHeaderView.ResizeToContents)
		self.search_tree = ParameterTreeView(self, self.params)
		self.search_tree.header().setStretchLastSection(False)
		self.search_tree.header().setSectionResizeMode(0, QHeaderView.ResizeToContents)
		self.search_tree.header().setSectionResizeMode(1, QHeaderView.Stretch)
		self.search_tree.header().setSectionResizeMode(2, QHeaderView.ResizeToContents)

		self.algorithm_select = QComboBox(self)
		self.algorithm_select.addItems([
			"Lamarckian GA",
			"Genetic Algorithm",
			"Simulated Annealing",
			"Local Search"
		])
		self.algorithm_select.currentIndexChanged.connect(self.change_search_params)

		layout = QFormLayout()
		self.setLayout(layout)

		layout.addRow("Search algorithm", self.algorithm_select)
		layout.addRow(QLabel("Search parameters"))
		layout.addRow(self.search_tree)
		layout.addRow(QLabel("Docking parameters"))
		layout.addRow(self.dock_tree)

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.accept)
		action_box.rejected.connect(self.reject)

		layout.addRow(action_box)

		items = self.params.get_ordered_items()

		for param, meta in items:
			if meta['scope'] == 'global':
				self.dock_tree.add_item(param, meta['comment'], meta['value'])

			else:
				self.search_tree.add_item(param, meta['comment'], meta['value'])

		#self.dock_tree.resizeColumnToContents(0)
		#self.search_tree.resizeColumnToContents(0)
		#self.dock_tree.resizeColumnToContents(2)
		#self.search_tree.resizeColumnToContents(2)

	@Slot()
	def change_search_params(self, index):
		algorithm = ['LGA', 'GA', 'SA', 'LS'][index]
		items = self.params.get_scope_items(algorithm)

		self.search_tree.clear()

		for param, meta in items:
			parent = self.search_tree.add_item(param, meta['comment'], meta['value'])
			if type(meta['value']) == list and len(meta['value']) > 1:
				for v in meta['value']:
					parent.appendRow([QStandardItem(''), QStandardItem(''), QStandardItem(str(v))])

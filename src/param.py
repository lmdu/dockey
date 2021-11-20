import os
import sys

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *

__all__ = ['AutogridParameter', 'AutodockParameterDialog', 'AutodockParamterWizard']

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
		self.algorithm = 0
		self.params = {
			'autodock_parameter_version': AttrDict({
				'type': str,
				'default': '4.2',
				'value': '4.2',
				'comment': 'used by autodock to validate parameter set',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'outlev': AttrDict({
				'type': int,
				'default': 1,
				'value': 1,
				'comment': 'diagnostic output level',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'parameter_file': AttrDict({
				'type': str,
				'default': 'AD4.1_bound.dat',
				'value': '',
				'comment': 'parameter library filename',
				'scope': 'global',
				'required': False,
				'user': True,
				'order': self.order
			}),
			'intelec': AttrDict({
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'calculate internal electrostatics',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'intnbp_r_eps': AttrDict({
				'type': [float,float,int,int,str,str],
				'default': [],
				'value': [],
				'comment': 'internal energy potential for a given class of interactions',
				'required': False,
				'user': False,
				'order': self.order
			}),
			'torsdof': AttrDict({
				'type': int,
				'default': 0,
				'value': 0,
				'comment': 'torsional degrees of freedom',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'seed': AttrDict({
				'type': [str, str],
				'default': ['pid', 'time'],
				'value': ['pid', 'time'],
				'comment': 'seeds for random generator',
				'scope': 'global',
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ligand_types': AttrDict({
				'type': list,
				'default': [],
				'value': [],
				'comment': 'atoms types in ligand',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'fld': AttrDict({
				'type': str,
				'default': '',
				'value': '',
				'comment': 'grid_data_file',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'map': AttrDict({
				'type': iter,
				'default': [],
				'value': [],
				'comment': 'atom-specific affinity map',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'elecmap': AttrDict({
				'type': str,
				'default': '',
				'value': '',
				'comment': 'electrostatics map',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'desolvmap': AttrDict({
				'type': str,
				'default': '',
				'value': '',
				'comment': 'desolvation map',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'move': AttrDict({
				'type': str,
				'default': '',
				'value': '',
				'comment': 'small molecule',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'flexres': AttrDict({
				'type': str,
				'default': '',
				'value': '',
				'comment': 'file containing flexible residues',
				'scope': 'global',
				'required': False,
				'user': True,
				'order': self.order,
			}),
			'about': AttrDict({
				'type': list,
				'default': [],
				'value': [],
				'comment': 'small molecule center',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'tran0': AttrDict({
				'type': list,
				'default': ['random'],
				'value': ['random'],
				'comment': 'initial coordinates/A or random',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'quaternion0': AttrDict({
				'type': list,
				'default': ['random'],
				'value': ['random'],
				'comment': 'initial orientation',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'dihe0': AttrDict({
				'type': list,
				'default': ['random'],
				'value': ['random'],
				'comment': 'initial dihedrals (relative) or random',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'tstep': AttrDict({
				'type': float,
				'default': 0.2,
				'value': 0.2,
				'range': [0, 1000000000],
				'comment': 'translation step/A',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'qstep': AttrDict({
				'type': float,
				'default': 5.0,
				'value': 5.0,
				'range': [0, 1000000000],
				'comment': 'quaternion step/deg',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'dstep': AttrDict({
				'type': float,
				'default': 5.0,
				'value': 5.0,
				'range': [0, 1000000000],
				'comment': 'torsion step/deg',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'rmstol': AttrDict({
				'type': float,
				'default': 2.0,
				'value': 2.0,
				'range': [0, 1000000000],
				'comment': 'cluster_tolerance/A',
				'scope': 'global',
				'required': True,
				'user': True,
				'order': self.order
			}),
			'epdb': AttrDict({
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'evaluate ligand specified with move command',
				'scope': 'global',
				'required': True,
				'user': True,
				'order': self.order
			}),
			'e0max': AttrDict({
				'type': [float, int],
				'default': [0.0, 10000],
				'value': [0.0, 10000],
				'range': [0, 1000000000],
				'comment': 'max initial energy; max number of retries',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'rt0': AttrDict({
				'type': float,
				'default': 616.0,
				'value': 616.0,
				'range': [0, 1000000000],
				'comment': 'initial annealing temperature (times gas constant)',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'linear_schedule': AttrDict({
				'type': bool,
				'default': False,
				'value': True,
				'comment': 'use linear, arithmetic temperature reduction',
				'scope': ['SA'],
				'required': False,
				'user': True,
				'order': self.order
			}),
			'geometric_schedule': AttrDict({
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'use geometric, arithmetic temperature reduction',
				'scope': ['SA'],
				'required': False,
				'user': True,
				'order': self.order
			}),
			'rtrf': AttrDict({
				'type': float,
				'default': 0.95,
				'value': 0.95,
				'range': [0, 1],
				'comment': 'annealing temperature reduction factor',
				'scope': ['SA'],
				'required': False,
				'user': True,
				'order': self.order
			}),
			'runs': AttrDict({
				'type': int,
				'default': 10,
				'value': 10,
				'range': [0, 1000000000],
				'comment': 'number of automated docking runs',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'cycles': AttrDict({
				'type': int,
				'default': 50,
				'value': 50,
				'range': [0, 1000000000],
				'comment': 'number of temperature reduction cycles',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'accs': AttrDict({
				'type': int,
				'default': 25000,
				'value': 25000,
				'range': [0, 1000000000],
				'comment': 'maximum number of accepted steps per cycle',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'rejs': AttrDict({
				'type': int,
				'default': 25000,
				'value': 25000,
				'range': [0, 1000000000],
				'comment': 'maximum number of rejected steps per cycle',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'rejs': AttrDict({
				'type': int,
				'default': 25000,
				'value': 25000,
				'range': [0, 1000000000],
				'comment': 'maximum number of rejected steps per cycle',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'select': AttrDict({
				'type': str,
				'default': 'm',
				'value': 'm',
				'comment': 'state selection flag: (m)inimum or (l)ast state',
				'choices': ['m', 'l'],
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'trnrf': AttrDict({
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'range': [0, 1000000000],
				'comment': 'per cycle reduction factor for translation',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'quarf': AttrDict({
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'range': [0, 1000000000],
				'comment': 'per cycle reduction factor for quaternions',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'dihrf': AttrDict({
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'range': [0, 1000000000],
				'comment': 'per cycle reduction factor for dihedrals',
				'scope': ['SA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_pop_size': AttrDict({
				'type': int,
				'default': 150,
				'value': 150,
				'range': [50, 200],
				'comment': 'No. of individuals in population',
				'scope': ['GA', 'LGA', 'LS'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_num_evals': AttrDict({
				'type': int,
				'default': 2500000,
				'value': 2500000,
				'range': [0, 1000000000],
				'comment': 'Max number of energy evaluations',
				'scope': ['GA', 'LGA', 'LS'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_num_generations': AttrDict({
				'type': int,
				'default': 27000,
				'value': 27000,
				'range': [0, 1000000000],
				'comment': 'maximum number of energy evaluations',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_elitism': AttrDict({
				'type': int,
				'default': 1,
				'value': 1,
				'range': [0, 200],
				'comment': 'number of top individuals to survive to next generation',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_mutation_rate': AttrDict({
				'type': float,
				'default': 0.02,
				'value': 0.02,
				'range': [0, 1],
				'comment': 'Rate of gene mutation',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_crossover_rate': AttrDict({
				'type': float,
				'default': 0.8,
				'value': 0.8,
				'range': [0, 1],
				'comment': 'Rate of crossover',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_window_size': AttrDict({
				'type': int,
				'default': 10,
				'value': 10,
				'range': [0, 1000000000],
				'comment': 'Number of preceding generations',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ga_cauchy_alpha': AttrDict({
				'type': float,
				'default': 0.0,
				'value': 0.0,
				'comment': 'Alpha parameter of Cauchy distribution',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'ga_cauchy_beta': AttrDict({
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'comment': 'Beta parameter Cauchy distribution',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'set_ga': AttrDict({
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'set the above parameters for GA or LGA',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'set_ga': AttrDict({
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'set the above parameters for GA or LGA',
				'scope': ['GA', 'LGA'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'sw_max_its': AttrDict({
				'type': int,
				'default': 300,
				'value': 300,
				'range': [0, 1000000000],
				'comment': 'iterations of Solis & Wets local search',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'sw_max_succ': AttrDict({
				'type': int,
				'default': 4,
				'value': 4,
				'range': [0, 10],
				'comment': 'consecutive successes before changing rho',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'sw_max_fail': AttrDict({
				'type': int,
				'default': 4,
				'value': 4,
				'range': [0, 10],
				'comment': 'consecutive failures before changing rho',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'sw_rho': AttrDict({
				'type': float,
				'default': 1.0,
				'value': 1.0,
				'range': [0, 1000000000],
				'comment': 'size of local search space to sample',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'sw_lb_rho': AttrDict({
				'type': float,
				'default': 0.01,
				'value': 0.01,
				'range': [0, 1],
				'comment': 'lower bound on rho',
				'scope': ['LS', 'LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'ls_search_freq': AttrDict({
				'type': float,
				'default': 0.06,
				'value': 0.06,
				'range': [0, 1],
				'comment': 'lower bound on rho',
				'scope': ['LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'set_sw1': AttrDict({
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'set the above classical Solis & Wets parameters',
				'scope': ['LGA', 'LS'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'set_psw1': AttrDict({
				'type': bool,
				'default': False,
				'value': False,
				'comment': 'set the above pseudo-Solis & Wets parameters',
				'scope': ['LGA', 'LS'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'unbound_model': AttrDict({
				'type': str,
				'default': 'bound',
				'value': 'bound',
				'comment': 'state of unbound ligand',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			}),
			'simanneal': AttrDict({
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'do as many SA runs as set by runs keyword above',
				'scope': ['SA'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'do_local_only': AttrDict({
				'type': int,
				'default': 50,
				'value': 50,
				'comment': 'do this many LS runs',
				'scope': ['LS'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'do_global_only': AttrDict({
				'type': int,
				'default': 50,
				'value': 50,
				'comment': 'do this many GA runs',
				'scope': ['GA'],
				'required': True,
				'user': False,
				'order': self.order
			}),
			'ga_run': AttrDict({
				'type': int,
				'default': 10,
				'value': 10,
				'range': [0, 1000000000],
				'comment': 'do this many hybrid GA-LS runs',
				'scope': ['LGA'],
				'required': True,
				'user': True,
				'order': self.order
			}),
			'analysis': AttrDict({
				'type': bool,
				'default': True,
				'value': True,
				'comment': 'perform a ranked cluster analysis',
				'scope': 'global',
				'required': True,
				'user': False,
				'order': self.order
			})
		}

	@property
	def order(self):
		self.num += 1
		return self.num

	def get_items(self):
		return self.params.items()

	def set_value(self, k, v, idx=-1):
		if idx >= 0:
			self.params[k]['value'][idx] = v
		else:
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
'''
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
'''

class AutodockParameterDialog(QDialog):
	def __init__(self, parent):
		super(AutodockParameterDialog, self).__init__(parent)
		self.setWindowTitle("Parameter setting for AutoDock")
		self.params = AutodockParameter()
		self.titles = ["Lamarckian GA", "Genetic Algorithm", "Simulated Annealing", "Local Search"]
		self.algorithms = ['LGA', 'GA', 'SA', 'LS']
		self.algorithm_select = QComboBox(self)
		self.algorithm_select.addItems(self.titles)
		self.algorithm_select.currentIndexChanged.connect(self.change_tab_page)

		layout = QVBoxLayout()
		self.setLayout(layout)

		select_layout = QHBoxLayout()
		select_layout.addWidget(QLabel("Search algorithm", self))
		select_layout.addWidget(self.algorithm_select)
		layout.addLayout(select_layout)
		layout.addWidget(QLabel("Search parameters", self))

		self.tab = QTabWidget(self)
		self.tab.setDocumentMode(True)
		self.tab.setTabBarAutoHide(True)
		layout.addWidget(self.tab)

		btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		btn_box.accepted.connect(self.accept)
		btn_box.rejected.connect(self.reject)
		layout.addWidget(btn_box)

		self.create_tab_page()

	def update_parameter(self, cmd, value, idx=-1):
		self.params.set_value(cmd, value, idx)

	def change_tab_page(self, index):
		self.tab.clear()
		self.create_tab_page(index)

	def create_tab_page(self, index=0):
		self.page = QWidget(self)
		self.tab.addTab(self.page, self.titles[index])
		page_layout = QFormLayout()
		self.page.setLayout(page_layout)

		for cmd, meta in self.params.get_items():
			if not meta.user:
				continue

			if self.algorithms[index] not in meta.scope:
				continue

			if meta.type is int:
				editor = QSpinBox(self)
				if 'range' in meta:
					_min, _max = meta.range
					editor.setRange(_min, _max)
				editor.setValue(meta.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(cmd, x))
				page_layout.addRow(meta.comment, editor)

			elif meta.type is float:
				editor = QDoubleSpinBox(self)
				if 'range' in meta:
					_min, _max = meta.range
					editor.setRange(_min, _max)
				editor.setValue(meta.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(cmd, x))
				page_layout.addRow(meta.comment, editor)

			elif cmd == 'geometric_schedule':
				pass

			elif cmd == 'linear_schedule':
				slay = QHBoxLayout()
				bgrp = QButtonGroup()
				lbtn = QRadioButton("Linear", self)
				lbtn.toggled.connect(lambda x: self.update_parameter('linear_schedule', x))
				gbtn = QRadioButton("Geometric", self)
				gbtn.toggled.connect(lambda x: self.update_parameter('geometric_schedule', x))
				slay.addWidget(lbtn)
				bgrp.addButton(lbtn)
				slay.addWidget(gbtn)
				bgrp.addButton(gbtn)
				bgrp.setExclusive(True)
				page_layout.addRow(meta.comment, slay)

			elif meta.type is bool:
				editor = QCheckBox(self)
				editor.setChecked(meta.value)
				editor.stateChanged.connect(lambda x: self.update_parameter(cmd, x))
				page_layout.addRow(meta.comment, editor)

			elif 'choices' in meta:
				editor = QComboBox(self)
				editor.addItems(meta.choices)
				idx = editor.findText(meta.value)
				editor.setCurrentIndex(idx)
				editor.currentTextChanged.connect(lambda x: self.update_parameter(cmd, x))
				page_layout.addRow(meta.comment, editor)

			elif isinstance(meta.type, list):
				for i, t in enumerate(meta.type):
					if t is int:
						editor = QSpinBox(self)
					elif t is float:
						editor = QDoubleSpinBox(self)
					else:
						editor = QLineEdit(self)

					if 'range' in meta:
						_min, _max = meta.range
						editor.setRange(_min, _max)

					editor.setValue(meta.value[i])

					if i == 0:
						page_layout.addRow(meta.comment, editor)
					else:
						page_layout.addRow('', editor)

			else:
				page_layout.addRow(QLabel(meta.comment))

class AutodockParamterWizard(QWizard):
	def __init__(self, parent):
		super(AutodockParamterWizard, self).__init__(parent)
		self.setWindowTitle("Parameter settings for AutoDock")
		self.setWizardStyle(QWizard.ClassicStyle)
		self.params = AutodockParameter()

		self.algorithms = ["Lamarckian GA", "Genetic Algorithm",
							"Simulated Annealing", "Local Search"]

		self.create_select_page()
		self.create_algorithm_page()
		self.create_docking_page()
		self.create_finish_page()

	def update_parameter(self, cmd, value, idx=-1):
		self.params.set_value(cmd, value, idx)

	@Slot()
	def on_algorithm_changed(self, index):
		#remove widgets
		for i in range(self.algorithm_layout.rowCount()):
			self.algorithm_layout.removeRow(0)

		self.params.algorithm = index

		self.create_algorithm_widgets(index)

	def create_select_page(self):
		tips = QLabel((
			"<p><b>Lamarckian Genetic Algorithm</b> provides the most efficient "
			"search for general applications, and in most cases will be the technique "
			"used. It is typically effective for systems with about 10 rotatable bonds "
			"in the ligand.</p><p><b>Genetic Algorithm</b> may also be run without the "
			"local search, but this is typically less efficient than the Lamarckian GA-LS "
			"combination.</p><p><b>Simulated Annealing</b> is also less efficient that the "
			"Lamarckian Genetic Algorithm, but it can be useful in applications where search "
			"starting from a given point is desired.</p><p><b>Local Search</b> may be used "
			"to optimize a molecule in its local environment.</p>"
		), self)
		tips.setWordWrap(True)

		self.select_page = QWizardPage(self)
		self.select_page.setTitle("Select a property search algorithm")
		self.addPage(self.select_page)
		self.select_layout = QVBoxLayout()
		self.select_page.setLayout(self.select_layout)
		self.select_layout.addWidget(tips)
		self.select_layout.addWidget(QLabel("<b>Select search algorithm</b>", self))
		self.select_input = QComboBox(self)
		self.select_input.addItems(self.algorithms)
		self.select_layout.addWidget(self.select_input)
		self.select_input.currentIndexChanged.connect(self.on_algorithm_changed)

	def create_algorithm_page(self):
		self.algorithm_page = QWizardPage(self)
		self.addPage(self.algorithm_page)
		self.algorithm_layout = QFormLayout()
		self.algorithm_page.setLayout(self.algorithm_layout)
		self.create_algorithm_widgets()

	def create_algorithm_widgets(self, index=0):
		self.algorithm_page.setTitle("Set parameters for {} algorithm".format(self.algorithms[index]))
		scope = ['LGA', 'GA', 'SA', 'LS'][index]
		self.create_parameter_widgets(self.algorithm_layout, scope)

	def create_parameter_widgets(self, layout, scope='global'):
		for cmd, meta in self.params.get_items():
			if not meta.user:
				continue

			if scope not in meta.scope:
				continue

			if meta.type is int:
				editor = QSpinBox(self)
				if 'range' in meta:
					_min, _max = meta.range
					editor.setRange(_min, _max)
				editor.setValue(meta.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(cmd, x))
				layout.addRow(meta.comment, editor)

			elif meta.type is float:
				editor = QDoubleSpinBox(self)
				if 'range' in meta:
					_min, _max = meta.range
					editor.setRange(_min, _max)
				editor.setValue(meta.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(cmd, x))
				layout.addRow(meta.comment, editor)

			elif cmd == 'geometric_schedule':
				pass

			elif cmd == 'linear_schedule':
				slay = QHBoxLayout()
				bgrp = QButtonGroup()
				lbtn = QRadioButton("Linear", self)
				lbtn.setChecked(meta.value)
				lbtn.toggled.connect(lambda x: self.update_parameter('linear_schedule', x))
				gbtn = QRadioButton("Geometric", self)
				gbtn.toggled.connect(lambda x: self.update_parameter('geometric_schedule', x))
				slay.addWidget(lbtn)
				bgrp.addButton(lbtn)
				slay.addWidget(gbtn)
				bgrp.addButton(gbtn)
				bgrp.setExclusive(True)
				layout.addRow(meta.comment, slay)

			elif meta.type is bool:
				editor = QCheckBox(self)
				editor.setChecked(meta.value)
				editor.stateChanged.connect(lambda x: self.update_parameter(cmd, x))
				layout.addRow(meta.comment, editor)

			elif meta.type is str:
				editor = QLineEdit(self)
				editor.setText(meta.value)
				editor.textChanged.connect(lambda x: self.update_parameter(cmd, x))
				layout.addRow(meta.comment, editor)

			elif 'choices' in meta:
				editor = QComboBox(self)
				editor.addItems(meta.choices)
				idx = editor.findText(meta.value)
				editor.setCurrentIndex(idx)
				editor.currentTextChanged.connect(lambda x: self.update_parameter(cmd, x))
				layout.addRow(meta.comment, editor)

			elif isinstance(meta.type, list):
				for i, t in enumerate(meta.type):
					if t is int:
						editor = QSpinBox(self)
						editor.valueChanged.connect(lambda x: self.update_parameter(cmd, x, i))
						editor.setValue(meta.value[i])
					elif t is float:
						editor = QDoubleSpinBox(self)
						editor.valueChanged.connect(lambda x: self.update_parameter(cmd, x, i))
						editor.setValue(meta.value[i])
					else:
						editor = QLineEdit(self)
						editor.textChanged.connect(lambda x: self.update_parameter(cmd, x, i))
						editor.setText(meta.value[i])

					if 'range' in meta:
						_min, _max = meta.range
						editor.setRange(_min, _max)

					if i == 0:
						layout.addRow(meta.comment, editor)
					else:
						layout.addRow('', editor)

			else:
				layout.addRow(QLabel(meta.comment))

	def create_docking_page(self):
		self.docking_page = QWizardPage(self)
		self.docking_page.setTitle("Set other docking parameters")
		self.addPage(self.docking_page)
		self.docking_layout = QFormLayout()
		self.docking_page.setLayout(self.docking_layout)
		self.create_parameter_widgets(self.docking_layout)

	def create_finish_page(self):
		self.finish_page = QWizardPage(self)
		self.finish_page.setTitle("Start AutoDock")
		self.addPage(self.finish_page)
		self.finish_layout = QVBoxLayout()
		self.finish_page.setLayout(self.finish_layout)
		self.finish_layout.addWidget(
			QLabel((
				"<p>Everything is ready, please confirm the docking tasks<p>"
				"<p>Receptors: <b>{}</b></p>"
				"<p>Ligands: <b>{}</b></p>"
				"<p>Search algorithm: <b>{}</b></p>"
				"<p>Click <b>Finish</b> button to start docking tasks</p>".format(
					1, 20, self.algorithms[self.params.algorithm])
			), self)
		)

import os
import sys
import psutil

from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from backend import *

__all__ = ['AutogridParameter', 'AutodockParameterWizard', 'AutodockVinaParameterWizard']

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

	def __init__(self, receptor_file, ligand_file, size, center, spacing=0.375,
				 smooth=0.5, dielectric=-0.1465):
		self.receptor_file = receptor_file
		self.ligand_file = ligand_file
		self.size = size
		self.center = center
		self.spacing = spacing
		self.smooth = smooth
		self.dielectric = dielectric

		self.receptor_name = os.path.splitext(os.path.basename(self.receptor_file))[0]
		self.receptor_types = get_atom_types_from_pdbqt(self.receptor_file)
		self.ligand_types = get_atom_types_from_pdbqt(self.ligand_file)
		self.gpf_file = self.receptor_file.replace('.pdbqt', '.gpf')

	def make_gpf_file(self):
		self.rows = [
			"npts {0[0]} {0[1]} {0[2]}".format(self.size),
			"gridfld {}.maps.fld".format(self.receptor_name),
			"spacing {}".format(self.spacing),
			"receptor_types {}".format(' '.join(self.receptor_types)),
			"ligand_types {}".format(' '.join(self.ligand_types)),
			"receptor {}".format(os.path.basename(self.receptor_file)),
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
				fw.write("{:<{}}{}\n".format(row, max_len, self.comments[field]))

		return self.gpf_file

class Parameter(dict):
	num = 0

	def __setattr__(self, attr, val):
		self[attr] = val

	def __getattr__(self, attr):
		try:
			return self[attr]
		except KeyError:
			return getattr(self, attr)

	def make_param(self, **kwargs):
		p = AttrDict(kwargs)
		p['order'] = self.order
		return p

	@property
	def order(self):
		return self.get_order()

	@classmethod
	def get_order(cls):
		cls.num += 1
		return cls.num


class AutodockParameter(Parameter):
	algorithm = 0

	def __init__(self):
		self.autodock_parameter_version = self.make_param(
			type = str,
			default = '4.2',
			value = '4.2',
			comment = 'used by autodock to validate parameter set',
			scope = 'global',
			required = True,
			user = False
		)
		self.outlev = self.make_param(
			type = int,
			default = 1,
			value = 1,
			comment = 'diagnostic output level',
			scope = 'global',
			required = True,
			user = False
		)
		self.parameter_file = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'Parameter library filename',
			scope = 'global',
			required = False,
			user = True
		)
		self.intelec = self.make_param(
			type = bool,
			default = True,
			value = True,
			comment = 'calculate internal electrostatics',
			scope = 'global',
			required = True,
			user = False
		)
		self.intnbp_r_eps = self.make_param(
			type = [float,float,int,int,str,str],
			default = [],
			value = [],
			comment = 'internal energy potential for a given class of interactions',
			scope = 'global',
			required = False,
			user = False
		)
		self.torsdof = self.make_param(
			type = int,
			default = 0,
			value = 0,
			comment = 'torsional degrees of freedom',
			scope = 'global',
			required = False,
			user = False
		)
		self.seed = self.make_param(
			type = [str, str],
			default = ['pid', 'time'],
			value = ['pid', 'time'],
			comment = 'Seeds for random generator',
			scope = 'global',
			required = True,
			user = True
		)
		self.ligand_types = self.make_param(
			type = list,
			default = [],
			value = [],
			comment = 'atoms types in ligand',
			scope = 'global',
			required = True,
			user = False
		)
		self.fld = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'grid_data_file',
			scope = 'global',
			required = True,
			user = False
		)
		self.map = self.make_param(
			type = iter,
			default = [],
			value = [],
			comment = 'atom-specific affinity map',
			scope = 'global',
			required = True,
			user = False
		)
		self.elecmap = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'electrostatics map',
			scope = 'global',
			required = True,
			user = False
		)
		self.desolvmap = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'desolvation map',
			scope = 'global',
			required = True,
			user = False
		)
		self.move = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'small molecule',
			scope = 'global',
			required = True,
			user = False
		)
		self.flexres = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'File containing flexible residues',
			scope = 'global',
			required = False,
			user = True,
		)
		self.about = self.make_param(
			type = list,
			default = [],
			value = [],
			comment = 'small molecule center',
			scope = 'global',
			required = False,
			user = False
		)
		self.tran0 = self.make_param(
			type = list,
			default = ['random'],
			value = ['random'],
			comment = 'initial coordinates/A or random',
			scope = 'global',
			required = True,
			user = False
		)
		self.quaternion0 = self.make_param(
			type = list,
			default = ['random'],
			value = ['random'],
			comment = 'initial orientation',
			scope = 'global',
			required = True,
			user = False
		)
		self.dihe0 = self.make_param(
			type = list,
			default = ['random'],
			value = ['random'],
			comment = 'initial dihedrals (relative) or random',
			scope = 'global',
			required = True,
			user = False
		)
		self.tstep = self.make_param(
			type = float,
			default = 0.2,
			value = 0.2,
			range = [0, 1000000000],
			comment = 'Translation step/A',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.qstep = self.make_param(
			type = float,
			default = 5.0,
			value = 5.0,
			range = [0, 1000000000],
			comment = 'Quaternion step/deg',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.dstep = self.make_param(
			type = float,
			default = 5.0,
			value = 5.0,
			range = [0, 1000000000],
			comment = 'Torsion step/deg',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.rmstol = self.make_param(
			type = float,
			default = 2.0,
			value = 2.0,
			range = [0, 1000000000],
			comment = 'Cluster tolerance (Å)',
			scope = 'global',
			required = True,
			user = True
		)
		self.epdb = self.make_param(
			type = bool,
			default = False,
			value = False,
			comment = 'Evaluate ligand specified with move command',
			scope = 'global',
			required = True,
			user = True
		)
		self.e0max = self.make_param(
			type = [float, int],
			default = [0.0, 10000],
			value = [0.0, 10000],
			range = [0, 1000000000],
			comment = 'Max initial energy and max number of retries',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.rt0 = self.make_param(
			type = float,
			default = 616.0,
			value = 616.0,
			range = [0, 1000000000],
			comment = 'Initial annealing temperature (times gas constant)',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.linear_schedule = self.make_param(
			type = bool,
			default = False,
			value = True,
			comment = 'use linear, arithmetic temperature reduction',
			scope = ['SA'],
			required = False,
			user = True
		)
		self.geometric_schedule = self.make_param(
			type = bool,
			default = False,
			value = False,
			comment = 'use geometric, arithmetic temperature reduction',
			scope = ['SA'],
			required = False,
			user = True
		)
		self.rtrf = self.make_param(
			type = float,
			default = 0.95,
			value = 0.95,
			range = [0, 1],
			comment = 'Annealing temperature reduction factor',
			scope = ['SA'],
			required = False,
			user = True
		)
		self.runs = self.make_param(
			type = int,
			default = 10,
			value = 10,
			range = [0, 1000000000],
			comment = 'Number of automated docking runs',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.cycles = self.make_param(
			type = int,
			default = 50,
			value = 50,
			range = [0, 1000000000],
			comment = 'Number of temperature reduction cycles',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.accs = self.make_param(
			type = int,
			default = 25000,
			value = 25000,
			range = [0, 1000000000],
			comment = 'Maximum number of accepted steps per cycle',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.rejs = self.make_param(
			type = int,
			default = 25000,
			value = 25000,
			range = [0, 1000000000],
			comment = 'Maximum number of rejected steps per cycle',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.select = self.make_param(
			type = str,
			default = 'm',
			value = 'm',
			comment = 'State selection: flag = (m)inimum or (l)ast state',
			choics = ['m', 'l'],
			scope = ['SA'],
			requred = True,
			user = True
		)
		self.trnrf = self.make_param(
			type = float,
			default = 1.0,
			value = 1.0,
			range = [0, 1000000000],
			comment = 'Per cycle reduction factor for translation',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.quarf = self.make_param(
			type = float,
			default = 1.0,
			value = 1.0,
			range = [0, 1000000000],
			comment = 'Per cycle reduction factor for quaternions',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.dihrf = self.make_param(
			type = float,
			default = 1.0,
			value = 1.0,
			range = [0, 1000000000],
			comment = 'Per cycle reduction factor for dihedrals',
			scope = ['SA'],
			required = True,
			user = True
		)
		self.ga_pop_size = self.make_param(
			type = int,
			default = 150,
			value = 150,
			range = [50, 200],
			comment = 'Number of individuals in population',
			scope = ['GA', 'LGA', 'LS'],
			required = True,
			user = True
		)
		self.ga_num_evals = self.make_param(
			type = int,
			default = 2500000,
			value = 2500000,
			range = [0, 1000000000],
			comment = 'Maximum number of energy evaluations',
			scope = ['GA', 'LGA', 'LS'],
			required = True,
			user = True
		)
		self.ga_num_generations = self.make_param(
			type = int,
			default = 27000,
			value = 27000,
			range = [0, 1000000000],
			comment = 'Maximum number of generations',
			scope = ['GA', 'LGA'],
			required = True,
			user = True
		)
		self.ga_elitism = self.make_param(
			type = int,
			default = 1,
			value = 1,
			range = [0, 200],
			comment = 'Number of top individuals to survive to next generation',
			scope = ['GA', 'LGA'],
			required = True,
			user = True
		)
		self.ga_mutation_rate = self.make_param(
			type = float,
			default = 0.02,
			value = 0.02,
			range = [0, 1],
			comment = 'Rate of gene mutation',
			scope = ['GA', 'LGA'],
			required = True,
			user = True
		)
		self.ga_crossover_rate = self.make_param(
			type = float,
			default = 0.8,
			value = 0.8,
			range = [0, 1],
			comment = 'Rate of crossover',
			scope = ['GA', 'LGA'],
			required = True,
			user = True
		)
		self.ga_window_size = self.make_param(
			type = int,
			default = 10,
			value = 10,
			range = [0, 1000000000],
			comment = 'Number of preceding generations',
			scope = ['GA', 'LGA'],
			required = True,
			user = True
		)
		self.ga_cauchy_alpha = self.make_param(
			type = float,
			default = 0.0,
			value = 0.0,
			comment = 'Alpha parameter of Cauchy distribution',
			scope = ['GA', 'LGA'],
			required = True,
			user = False
		)
		self.ga_cauchy_beta = self.make_param(
			type = float,
			default = 1.0,
			value = 1.0,
			comment = 'Beta parameter Cauchy distribution',
			scope = ['GA', 'LGA'],
			required = True,
			user = False
		)	
		self.set_ga = self.make_param(
			type = bool,
			default = True,
			value = True,
			comment = 'set the above parameters for GA or LGA',
			scope = ['GA', 'LGA'],
			required = True,
			user = False
		)
		self.set_ga = self.make_param(
			type = bool,
			default = True,
			value = True,
			comment = 'set the above parameters for GA or LGA',
			scope = ['GA', 'LGA'],
			required = True,
			user = False
		)
		self.sw_max_its = self.make_param(
			type = int,
			default = 300,
			value = 300,
			range = [0, 1000000000],
			comment = 'Iterations of Solis and Wets local search',
			scope = ['LS', 'LGA'],
			required = True,
			user = True
		)
		self.sw_max_succ = self.make_param(
			type = int,
			default = 4,
			value = 4,
			range = [0, 10],
			comment = 'Consecutive successes before changing rho',
			scope = ['LS', 'LGA'],
			required = True,
			user = True
		)
		self.sw_max_fail = self.make_param(
			type = int,
			default = 4,
			value = 4,
			range = [0, 10],
			comment = 'Consecutive failures before changing rho',
			scope = ['LS', 'LGA'],
			required = True,
			user = True
		)
		self.sw_rho = self.make_param(
			type = float,
			default = 1.0,
			value = 1.0,
			range = [0, 1000000000],
			comment = 'Size of local search space to sample',
			scope = ['LS', 'LGA'],
			required = True,
			user = True
		)
		self.sw_lb_rho = self.make_param(
			type = float,
			default = 0.01,
			value = 0.01,
			range = [0, 1],
			comment = 'Lower bound on rho',
			scope = ['LS', 'LGA'],
			required = True,
			user = True
		)
		self.ls_search_freq = self.make_param(
			type = float,
			default = 0.06,
			value = 0.06,
			range = [0, 1],
			comment = 'Probability of particular phenotype for local research',
			scope = ['LGA'],
			required = True,
			user = True
		)
		self.set_sw1 = self.make_param(
			type = bool,
			default = False,
			value = False,
			comment = 'Classic Solis and Wets with uniform variances',
			scope = ['LGA', 'LS'],
			required = True,
			user = True
		)
		self.set_psw1 = self.make_param(
			type = bool,
			default = True,
			value = True,
			comment = 'Solis and Wets local searcher',
			scope = ['LGA', 'LS'],
			required = True,
			user = True
		)
		self.unbound_model = self.make_param(
			type = str,
			default = 'bound',
			value = 'bound',
			comment = 'state of unbound ligand',
			scope = 'global',
			required = True,
			user = False
		)
		self.simanneal = self.make_param(
			type = bool,
			default = True,
			value = True,
			comment = 'do as many SA runs as set by runs keyword above',
			scope = ['SA'],
			required = True,
			user = False
		)
		self.do_local_only = self.make_param(
			type = int,
			default = 50,
			value = 50,
			comment = 'do this many LS runs',
			scope = ['LS'],
			required = True,
			user = False
		)
		self.do_global_only = self.make_param(
			type = int,
			default = 50,
			value = 50,
			comment = 'do this many GA runs',
			scope = ['GA'],
			required = True,
			user = False
		)
		self.ga_run = self.make_param(
			type = int,
			default = 10,
			value = 10,
			range = [0, 1000000000],
			comment = 'Perform the requested number of dockings',
			scope = ['LGA'],
			required = True,
			user = True
		)
		self.analysis = self.make_param(
			type = bool,
			default = True,
			value = True,
			comment = 'Perform a ranked cluster analysis',
			scope = 'global',
			required = True,
			user = False
		)

	@classmethod
	def set_algorithm(cls, index):
		cls.algorithm = index

	def set_value(self, k, v, idx=-1):
		if idx >= 0:
			self[k]['value'][idx] = v
		else:
			self[k]['value'] = v

	def get_value(self, k):
		return self[k]['value']

	def get_param(self, k):
		return self[k]

	def get_count(self, algorithm='LGA'):
		count = 0
		for p in self:
			if self[p].user and (self[p].scope == 'global' or algorithm in self[p].scope):
				count += 1

		return count

	def get_scope_items(self, algorithm='LGA'):
		params = {}
		for p in self:
			if self[p].user and algorithm in self[p].scope:
				params[p] = self[p]

		return sorted(params.items(), key=lambda x: x[1].order)

	def get_ordered_items(self, algorithm='LGA'):
		params = {}
		for p in self:
			if self[p].user and (self[p].scope == 'global' or algorithm in self[p].scope):
				params[p] = self[p]

		return sorted(params.items(), key=lambda x: x[1].order)

	def make_dpf_file(self, receptor_file, ligand_file):
		ligand_types = get_atom_types_from_pdbqt(ligand_file)
		self.ligand_types.value = ligand_types

		receptor_name = QFileInfo(receptor_file).baseName()
		ligand_name = QFileInfo(ligand_file).fileName()

		self.fld.value = "{}.maps.fld".format(receptor_name)
		self.map.value = ["{}.{}.map".format(receptor_name, ligand_type) for ligand_type in ligand_types]
		self.elecmap.value = "{}.e.map".format(receptor_name)
		self.desolvmap.value = "{}.d.map".format(receptor_name)
		self.move.value = ligand_name
		self.about.value = get_molecule_center_from_pdbqt(ligand_file)

		algorithm = ['LGA', 'GA', 'SA', 'LS'][self.algorithm]
		
		params = {}
		for p in self:
			if self[p].scope == 'global' or algorithm in self[p].scope:
				params[p] = self[p]
		params = sorted(params.items(), key=lambda x: x[1].order)

		rows = []
		for k, v in params:
			if not v.required and v.value == v.default:
				continue

			if k == 'map':
				for m in v.value:
					rows.append("map {}".format(m))

			elif v.type in (int, str):
				rows.append("{} {}".format(k, v.value))

			elif v.type is list:
				rows.append("{} {}".format(k, ' '.join(map(str, v.value))))

			elif v.type is bool:
				if v.value:
					rows.append(k)

		max_len = max([len(row) for row in rows]) + 6

		gdf_file = receptor_file.replace('.pdbqt', '.dpf')

		with open(gdf_file, 'w') as fw:
			for row in rows:
				field = row.split()[0]
				fw.write("{:<{}}#{}\n".format(row, max_len, self[field].comment))

		return gdf_file

class AutodockParameterWizard(QWizard):
	params = AutodockParameter()

	def __init__(self, parent):
		super(AutodockParameterWizard, self).__init__(parent)
		self.setWindowTitle("Parameter settings for AutoDock")
		self.setWizardStyle(QWizard.ClassicStyle)

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

		#self.params.algorithm = index
		self.params.set_algorithm(index)

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
			"to optimize a molecule in its local environment.</p><p></p>"
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
		for p, m in self.params.items():
			if not m.user:
				continue

			if scope not in m.scope:
				continue

			if m.type is int:
				editor = QSpinBox(self)
				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)
				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)
				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)
				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif p in ('geometric_schedule', 'set_sw1'):
				pass

			elif p == 'linear_schedule':
				btn_layout = QHBoxLayout()
				btn_group = QButtonGroup()
				line_btn = QRadioButton("Linear", self)
				line_btn.setChecked(m.value)
				line_btn.toggled.connect(lambda x: self.update_parameter('linear_schedule', x))
				geome_btn = QRadioButton("Geometric", self)
				geome_btn.toggled.connect(lambda x: self.update_parameter('geometric_schedule', x))
				btn_layout.addWidget(line_btn)
				btn_group.addButton(line_btn)
				btn_layout.addWidget(geome_btn)
				btn_group.addButton(geome_btn)
				btn_group.setExclusive(True)
				layout.addRow(m.comment, btn_layout)

			elif p == 'set_psw1':
				btn_layout = QHBoxLayout()
				btn_group = QButtonGroup()
				sw_btn = QRadioButton("classic", self)
				sw_btn.toggled.connect(lambda x: self.update_parameter('set_sw1', x))
				btn_group.addButton(sw_btn)
				btn_layout.addWidget(sw_btn)

				psw_btn = QRadioButton("pseudo", self)
				psw_btn.setChecked(m.value)
				psw_btn.toggled.connect(lambda x: self.update_parameter('set_psw1', x))
				btn_group.addButton(psw_btn)
				btn_layout.addWidget(psw_btn)
				btn_group.setExclusive(True)
				layout.addRow(m.comment, btn_layout)

			elif m.type is bool:
				editor = QCheckBox(self)
				editor.setChecked(m.value)
				editor.stateChanged.connect(lambda x,p=p: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif m.type is str:
				editor = QLineEdit(self)
				editor.setText(m.value)
				editor.textChanged.connect(lambda x,p=p: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif 'choices' in m:
				editor = QComboBox(self)
				editor.addItems(m.choices)
				idx = editor.findText(m.value)
				editor.setCurrentIndex(idx)
				editor.currentTextChanged.connect(lambda x,p=p: self.update_parameter(p, x))
				layout.addRow(m.comment, editor)

			elif isinstance(m.type, list):
				for i, t in enumerate(m.type):
					if t is int:
						editor = QSpinBox(self)
						editor.valueChanged.connect(lambda x,p=p,i=i: self.update_parameter(p, x, i))
						editor.setValue(m.value[i])
					elif t is float:
						editor = QDoubleSpinBox(self)
						editor.valueChanged.connect(lambda x,p=p,i=i: self.update_parameter(p, x, i))
						editor.setValue(m.value[i])
					else:
						editor = QLineEdit(self)
						editor.textChanged.connect(lambda x,p=p,i=i: self.update_parameter(p, x, i))
						editor.setText(m.value[i])

					if 'range' in m:
						_min, _max = m.range
						editor.setRange(_min, _max)

					if i == 0:
						layout.addRow(m.comment, editor)
					else:
						layout.addRow('', editor)

			else:
				layout.addRow(QLabel(m.comment))

	def create_docking_page(self):
		self.docking_page = QWizardPage(self)
		self.docking_page.setTitle("Set other docking parameters")
		self.addPage(self.docking_page)
		self.docking_layout = QFormLayout()
		self.docking_page.setLayout(self.docking_layout)
		self.create_parameter_widgets(self.docking_layout)

	def create_finish_page(self):
		self.finish_page = QWizardPage(self)
		self.finish_page.setTitle("Start Autodock")
		self.addPage(self.finish_page)
		self.finish_layout = QVBoxLayout()
		self.finish_page.setLayout(self.finish_layout)
		rnum = DB.get_one("SELECT COUNT(1) FROM molecular WHERE type=1 LIMIT 1")
		lnum = DB.get_one("SELECT COUNT(1) FROM molecular WHERE type=2 LIMIT 1")
		info_wdg = QLabel((
			"<p>Everything is ready, please confirm the docking jobs<p>"
			"<p>The number of receptors: <b>{}</b></p>"
			"<p>The number of ligands: <b>{}</b></p>"
			"<p>The number of jobs will be generated: <b>{}</b></p>"
			"<p>Selected search algorithm: <b>{}</b></p>"
			"<p>Click <b>Finish</b> button to submit and start docking jobs</p>".format(
				rnum, lnum, rnum*lnum, self.algorithms[self.params.algorithm]
			)
		), self)
		self.finish_layout.addWidget(info_wdg)

class AutodockVinaParameter(Parameter):
	def __init__(self):
		self.receptor = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'rigid part of the receptor',
			required = False,
			user = False
		)
		self.flex = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'flexible side chains, if any',
			required = False,
			user = False
		)
		self.ligand = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'ligand',
			required = True,
			user = False
		)
		self.scoring = self.make_param(
			type = str,
			default = 'vina',
			value = 'vina',
			choices = ['vina', 'autodock4', 'vinardo'],
			comment = 'Scoring function',
			required = True,
			user = True
		)
		self.maps = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'affinity maps for the autodock4.2 (ad4) or vina',
			required = False,
			user = False
		)
		self.center_x = self.make_param(
			type = float,
			default = 0,
			value = 0,
			comment = 'X coordinate of the center',
			required = True,
			user = False
		)
		self.center_y = self.make_param(
			type = float,
			default = 0,
			value = 0,
			comment = 'Y coordinate of the center',
			required = True,
			user = False
		)
		self.center_z = self.make_param(
			type = float,
			default = 0,
			value = 0,
			comment = 'Z coordinate of the center',
			required = True,
			user = False
		)
		self.size_x = self.make_param(
			type = float,
			default = 0,
			value = 0,
			comment = 'size in the X dimension',
			required = True,
			user = False
		)
		self.size_y = self.make_param(
			type = float,
			default = 0,
			value = 0,
			comment = 'size in the Y dimension',
			required = True,
			user = False
		)
		self.size_z = self.make_param(
			type = float,
			default = 0,
			value = 0,
			comment = 'size in the Z dimension',
			required = True,
			user = False
		)
		self.out = self.make_param(
			type = str,
			default = '',
			value = '',
			comment = 'output models',
			required = True,
			user = False
		)
		self.cpu = self.make_param(
			type = int,
			default = 1,
			value = 1,
			range = [1, psutil.cpu_count()],
			comment = 'The number of CPUs to use',
			required = True,
			user = True
		)
		self.seed = self.make_param(
			type = int,
			default = 0,
			value = 0,
			range = [0, 1000000000],
			comment = 'Explicit random seed',
			required = False,
			user = True
		)
		self.exhaustiveness = self.make_param(
			type = int,
			default = 8,
			value = 8,
			range = [1, 100],
			comment = 'Exhaustiveness of the global search',
			required = False,
			user = True
		)
		self.max_evals = self.make_param(
			type = int,
			default = 0,
			value = 0,
			range = [0, 1000000000],
			comment = 'Number of evaluations in each MC run',
			required = False,
			user = True
		)
		self.num_modes = self.make_param(
			type = int,
			default = 9,
			value = 9,
			range = [0, 100],
			comment = 'Maximum number of binding modes to generate',
			required = False,
			user = True
		)
		self.min_rmsd = self.make_param(
			type = float,
			default = 1.0,
			value = 1.0,
			range = [0, 100],
			comment = 'Minimum RMSD between output poses',
			required = False,
			user = True
		)
		self.energy_range = self.make_param(
			type = float,
			default = 3.0,
			value = 3.0,
			range = [0, 100],
			comment = 'Maximum energy difference between the best and worst mode',
			required = False,
			user = True
		)
		self.spacing = self.make_param(
			type = float,
			default = 0.375,
			value = 0.375,
			comment = 'grid spacing',
			required = False,
			user = False
		)
		self.verbosity = self.make_param(
			type = int,
			default = 1,
			value = 1,
			choices = [0, 1, 2],
			comment = 'verbosity (0=no output, 1=normal, 2=verbose)',
			required = False,
			user = False
		)

	def get_sorted_items(self):
		return sorted(self.items(), key=lambda x: x[1].order)

	def make_config_file(self, config_file):
		items =  self.get_sorted_items()

		with open(config_file, 'w') as fw:
			for p, m in items:
				if m.required or m.default != m.value:
					fw.write("{} = {}\n".format(p, m.value))

class AutodockVinaParameterWizard(QWizard):
	params = AutodockVinaParameter()

	def __init__(self, parent=None):
		super(AutodockVinaParameterWizard, self).__init__(parent)
		self.setWindowTitle("Run Autodock Vina")
		self.setWizardStyle(QWizard.ClassicStyle)

		self.create_parameter_page()
		self.create_finish_page()

	def update_parameter(self, key, val):
		if key == 'scoring' and val == 'autodock4':
			val = 'ad4'

		self.params[key]['value'] = val

	def create_parameter_page(self):
		self.parameter_page = QWizardPage(self)
		self.parameter_page.setTitle("Parameter settings for Autodock Vina")
		self.addPage(self.parameter_page)
		self.parameter_layout = QFormLayout()
		self.parameter_page.setLayout(self.parameter_layout)
		self.create_parameter_widgets()

	def create_parameter_widgets(self):
		for p, m in self.params.get_sorted_items():
			if not m.user:
				continue

			if m.type is int:
				editor = QSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				self.parameter_layout.addRow(m.comment, editor)

			elif m.type is float:
				editor = QDoubleSpinBox(self)

				if 'range' in m:
					_min, _max = m.range
					editor.setRange(_min, _max)

				editor.setValue(m.value)
				editor.valueChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				self.parameter_layout.addRow(m.comment, editor)

			elif p == 'scoring':
				editor = QComboBox(self)
				editor.addItems(m.choices)
				editor.setCurrentIndex(0)
				editor.currentTextChanged.connect(lambda x, p=p: self.update_parameter(p, x))
				self.parameter_layout.addRow(m.comment, editor)

	def create_finish_page(self):
		self.finish_page = QWizardPage(self)
		self.finish_page.setTitle("Start Autodock Vina")
		self.addPage(self.finish_page)
		self.finish_layout = QVBoxLayout()
		self.finish_page.setLayout(self.finish_layout)
		rnum = DB.get_one("SELECT COUNT(1) FROM molecular WHERE type=1 LIMIT 1")
		lnum = DB.get_one("SELECT COUNT(1) FROM molecular WHERE type=2 LIMIT 1")
		info_wdg = QLabel((
			"<p>Everything is ready, please confirm the docking jobs<p>"
			"<p>The number of receptors: <b>{}</b></p>"
			"<p>The number of ligands: <b>{}</b></p>"
			"<p>The number of jobs will be generated: <b>{}</b></p>"
			"<p>Selected scoring function: <b>{}</b></p>"
			"<p>Click <b>Finish</b> button to submit and start docking jobs</p>".format(
				rnum, lnum, rnum*lnum, self.params.scoring.value
			)
		), self)
		self.finish_layout.addWidget(info_wdg)





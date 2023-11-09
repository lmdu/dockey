import os
import csv
import time
import psutil
import traceback
import subprocess
import multiprocessing

from PySide6.QtCore import *

from param import *
from utils import *
from process import *
from gridbox import *
from backend import *
from prepare import *

__all__ = ['AutodockWorker', 'AutodockVinaWorker', 'QuickVinaWorker',
	'ImportFileWorker', 'ImportFolderWorker', 'JobListGenerator',
	'PoseInteractionExportWorker', 'BestInteractionExportWorker',
	'ImportSDFWorker', 'ImportURLWorker', 'WorkerManager'
]

class ImportSignals(QObject):
	success = Signal(str)
	failure = Signal(str)
	message = Signal(str)
	finished = Signal()
	progress = Signal(int)

class ImportWorker(QRunnable):
	processer = None

	def __init__(self, mol_files, mol_type=0):
		super().__init__()
		self.setAutoDelete(True)
		self.signals = ImportSignals()
		self.mol_files = mol_files
		self.mol_type = mol_type
		self.mol_count = 0
		self.pipe = multiprocessing.SimpleQueue()
		self.mol_text = ['', 'receptors', 'ligands'][self.mol_type]

	def update_molecule_count(self):
		if self.mol_type == 1:
			option = 'receptor_count'
		else:
			option = 'ligand_count'

		count = DB.get_option(option)

		if count:
			count = int(count) + self.mol_count
		else:
			count = self.mol_count

		DB.set_option(option, count)

		total = DB.get_option('molecule_count')

		if total:
			total = int(total) + self.mol_count
		else:
			total = self.mol_count

		DB.set_option('molecule_count', total)

	def write_molecules(self, data):
		sql = "INSERT INTO molecular VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, data['rows'])
		self.signals.message.emit("Imported {} {}".format(data['total'], self.mol_text))

	def call_response(self, data):
		if data['action'] == 'insert':
			self.write_molecules(data)
			self.mol_count = data['total']

			if data['progress'] > 0:
				self.signals.progress.emit(data['progress'])

		elif data['action'] == 'failure':
			self.signals.failure.emit(data['message'])

		elif data['action'] == 'success':
			if self.mol_count > 1:
				self.signals.success.emit("Successfully imported {} {}".format(
					self.mol_count, self.mol_text))

		elif data['action'] == 'finish':
			while True:
				if self.pipe.empty():
					self.pipe.close()
					break

				QThread.msleep(10)

	@Slot()
	def run(self):
		self.signals.message.emit("Starting import molecules ...")
		self.signals.progress.emit(0)

		proc = self.processer(self.mol_files, self.mol_type, self.pipe)
		proc.start()

		while True:
			try:
				data = self.pipe.get()
				self.call_response(data)
			except:
				break

		if self.mol_count > 0:
			self.update_molecule_count()

		self.signals.progress.emit(100)
		self.signals.finished.emit()

class ImportFileWorker(ImportWorker):
	processer = ImportFileProcess

class ImportSDFWorker(ImportWorker):
	processer = ImportSDFProcess

class ImportFolderWorker(ImportWorker):
	processer = ImportFolderProcess

class ImportURLWorker(ImportWorker):
	processer = ImportURLProcess

class JobListSignals(QObject):
	failure = Signal(str)
	message = Signal(str)
	success = Signal()
	finished = Signal()
	progress = Signal(int)

class JobListGenerator(QRunnable):
	def __init__(self):
		super(JobListGenerator, self).__init__()
		self.setAutoDelete(True)
		self.signals = JobListSignals()
		self.pipe = multiprocessing.SimpleQueue()
		self.job_insert = 0

	def write_jobs(self, data):
		sql = "INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, data['rows'])
		self.signals.message.emit("Generated {} tasks".format(data['total']))
		self.job_insert += len(data['rows'])

		if data['progress']:
			self.signals.progress.emit(data['progress'])

	def call_response(self, data):
		if data['action'] == 'insert':
			self.write_jobs(data)

		elif data['action'] == 'failure':
			self.signals.failure.emit(data['message'])

		elif data['action'] == 'success':
			if self.job_insert > 0:
				DB.set_option('job_count', self.job_insert)
				#self.signals.success.emit()

		elif data['action'] == 'finish':
			while True:
				if self.pipe.empty():
					self.pipe.close()
					break

				QThread.msleep(10)

	def run(self):
		self.signals.progress.emit(0)

		rids = DB.get_column("SELECT id FROM molecular WHERE type=1")
		lids = DB.get_column("SELECT id FROM molecular WHERE type=2")

		proc = JobListProcess(rids, lids, self.pipe)
		proc.start()

		while True:
			try:
				data = self.pipe.get()
				self.call_response(data)
			except:
				break

		self.signals.progress.emit(100)
		self.signals.finished.emit()

class ManagerSignals(QObject):
	finished = Signal()
	feedback = Signal(int)
	updated = Signal(int)
	progress = Signal(int)
	message = Signal(str)

class WorkerManager(QRunnable):
	def __init__(self, parent, dock_worker, dock_params):
		super().__init__()
		self.parent = parent
		self.setAutoDelete(True)
		self.pool = QThreadPool()
		self.settings = QSettings()
		self.signals = ManagerSignals()
		self.dock_worker = dock_worker
		self.dock_params = dock_params
		self.job_query = None
		self.job_thread = self.settings.value('Job/concurrent', 1, int)
		self.job_list = {}
		self.job_count = 0
		self.job_params = self.get_all_params()

		self.signals.feedback.connect(self.remove_job)

	def before_run(self):
		self.job_total = int(DB.get_option('job_count'))
		self.job_query = DB.query("SELECT id FROM jobs")

	def get_job(self, jid):
		job = DB.get_dict("SELECT * FROM jobs WHERE id=? LIMIT 1", (jid,))
		rep = DB.get_dict("SELECT * FROM molecular WHERE id=? LIMIT 1", (job.rid,))
		lig = DB.get_dict("SELECT * FROM molecular WHERE id=? LIMIT 1", (job.lid,))
		grid = DB.get_dict("SELECT * FROM grid WHERE rid=? LIMIT 1", (job.rid,))

		if not grid:
			grid = [0, job.rid]
			sql = "SELECT content FROM molecular WHERE id=? LIMIT 1"
			pdb_str = DB.get_one(sql, (job.rid,))
			grid_space = self.settings.value('Grid/spacing', 0.375, float)
			dims = get_dimension_from_pdb(pdb_str, grid_space)
			grid = AttrDict({
				'x': dims[0],
				'y': dims[1],
				'z': dims[2],
				'cx': dims[3],
				'cy': dims[4],
				'cz': dims[5],
				'spacing': grid_space
			})

		# get receptor flex residues
		residues = {}
		for row in DB.query("SELECT * FROM flex WHERE rid=?", (job.rid,)):
			if row[3] not in residues:
				residues[row[3]] = {}

			aa = "{}{}".format(row[4], row[5])

			residues[row[3]][aa] = row[6].split(',')

		return AttrDict({
			'id': job.id,
			'rid': rep.id,
			'rn': rep.name,
			'rc': rep.content,
			'rf': rep.format,
			'lid': lig.id,
			'ln': lig.name,
			'lc': lig.content,
			'lf': lig.format,
			'x': grid.x,
			'y': grid.y,
			'z': grid.z,
			'cx': grid.cx,
			'cy': grid.cy,
			'cz': grid.cy,
			'spacing': grid.spacing,
			'flex': residues
		})

	def get_prepare_params(self):
		tool = self.settings.value('Ligand/prepare_tool', 'meeko')

		if tool == 'prepare_ligand4':
			self.settings.beginGroup('Ligand')
			lig_params = dict(
				tool = 'prepare_ligand4',
				repairs = self.settings.value('repairs', 'bonds_hydrogens'),
				charges_to_add = self.settings.value('charges_to_add', 'gasteiger'),
				#preserve_charge_types = self.settings.value('preserve_charge_types', ''),
				cleanup = self.settings.value('cleanup', 'nphs_lps'),
				allowed_bonds = self.settings.value('allowed_bonds', 'backbone'),
				check_for_fragments = self.settings.value('check_for_fragments', False, bool),
				inactivate_all_torsions = self.settings.value('inactivate_all_torsions', False, bool),
				attach_nonbonded_fragments = self.settings.value('attach_nonbonded_fragments', False, bool),
				attach_singletons = self.settings.value('attach_singletons', False, bool)
			)
			self.settings.endGroup()

		elif tool == 'meeko':
			self.settings.beginGroup('Meeko')
			lig_params = dict(
				tool = 'meeko',
				merge_these_atom_types = self.settings.value('merge_these_atom_types', 'H', str),
				rigid_macrocycles = self.settings.value('rigid_macrocycles', False, bool),
				keep_chorded_rings = self.settings.value('keep_chorded_rings', False, bool),
				keep_equivalent_rings = self.settings.value('keep_equivalent_rings', False, bool),
				hydrate = self.settings.value('hydrate', False, bool),
				flexible_amides = self.settings.value('flexible_amides', False, bool),
				add_index_map = self.settings.value('add_index_map', False, bool),
				remove_smiles = self.settings.value('remove_smiles', False, bool),
				rigidify_bonds_smarts = self.settings.value('rigidify_bonds_smarts', ''),
				rigidify_bonds_indices = self.settings.value('rigidify_bonds_indices', ''),
				atom_type_smarts = self.settings.value('atom_type_smarts', ''),
				double_bond_penalty = self.settings.value('double_bond_penalty', 50, int)
			)
			self.settings.endGroup()

		self.settings.beginGroup('Receptor')
		rep_params = dict(
			repairs = self.settings.value('repairs', 'bonds_hydrogens'),
			charges_to_add = self.settings.value('charges_to_add', 'gasteiger'),
			#preserve_charge_types = self.settings.value('preserve_charge_types', ''),
			cleanup = self.settings.value('cleanup', 'nphs_lps_waters_nonstdres'),
			delete_single_nonstd_residues = self.settings.value('delete_single_nonstd_residues', False, bool),
			unique_atom_names = self.settings.value('unique_atom_names', False, bool)
		)
		self.settings.endGroup()

		self.settings.beginGroup('PDBFixer')
		pdbfix_params = dict(
			use_pdbfix = self.settings.value('use_pdbfix', True, bool),
			replace_nonres = self.settings.value('replace_nonres', True, bool),
			remove_heterogen = self.settings.value('remove_heterogen', True, bool),
			add_misheavy = self.settings.value('add_misheavy', True, bool),
			add_mishydrogen = self.settings.value('add_mishydrogen', True, bool),
			mishydrogen_ph = self.settings.value('mishydrogen_ph', 7.0, float)
		)
		self.settings.endGroup()

		self.settings.beginGroup('PDB2PQR')
		pdbpqr_params = dict(
			use_pdbpqr = self.settings.value('use_pdbpqr', False, bool),
			force_field = self.settings.value('force_field', 'PARSE'),
			use_propka = self.settings.value('use_propka', True, bool),
			propka_ph = self.settings.value('propka_ph', 7.0, float),
			node_bump = self.settings.value('node_bump', False, bool),
			no_hopt = self.settings.value('no_hopt', False, bool),
			remove_water = self.settings.value('remove_water', True, bool),
			neutraln = self.settings.value('neutraln', False, bool),
			neutralc = self.settings.value('neutralc', False, bool)
		)
		self.settings.endGroup()

		pre_params = {}
		pre_params.update(pdbfix_params)
		pre_params.update(pdbpqr_params)

		return lig_params, rep_params, pre_params

	def get_all_params(self):
		lig_params, rep_params, pre_params = self.get_prepare_params()
		return [self.dock_params, lig_params, rep_params, pre_params]

	def start_job(self, jid):
		job = self.get_job(jid)
		self.job_list[jid] = self.dock_worker(self.job_params, job)
		self.job_list[jid].signals.finished.connect(self.parent.job_worker.signals.feedback)
		self.job_list[jid].signals.changed.connect(self.parent.job_model.update_row)
		#QThreadPool.globalInstance().start(worker)
		self.pool.start(self.job_list[jid])

	#@Slot(int)
	#def update_job(self, jid):
	#	self.signals.updated.emit(jid)

	@Slot(int)
	def remove_job(self, jid):
		if jid in self.job_list:
			self.job_list.pop(jid)

			self.job_count += 1
			p = int(self.job_count/self.job_total*100)
			self.signals.progress.emit(p)

	def stop_job(self, jid):
		if jid in self.job_list:
			self.job_list[jid].stop()

	def stop_jobs(self):
		self.job_thread = 0
		self.job_query = None

		for jid in list(self.job_list.keys()):
			if jid in self.job_list:
				self.job_list[jid].stop()

	def set_thread(self, num):
		self.job_thread = num

	def run(self):
		self.signals.progress.emit(0)
		self.signals.message.emit("Docking running ...")
		self.before_run()

		while True:
			if self.job_query is None:
				break

			if len(self.job_list) >= self.job_thread:
				QThread.msleep(10)
				continue

			row = self.job_query.fetchone()

			if row:
				self.start_job(row[0])

			if not self.job_list:
				break

		self.signals.progress.emit(100)
		self.signals.finished.emit()
		self.signals.message.emit("Docking completed")

class WorkerSignals(QObject):
	finished = Signal(int)
	changed = Signal(int)

class BaseWorker(QRunnable):
	processer = None

	def __init__(self, params, job):
		super().__init__()
		self.setAutoDelete(True)
		self.job = job
		self.params = params
		self.signals = WorkerSignals()
		self.settings = QSettings()
		self.pipe = multiprocessing.SimpleQueue()
		self.tempdir = None
		self.process = None

	def update_progress(self, data):
		sql = "UPDATE jobs SET progress=? WHERE id=?"
		DB.query(sql, (data['message'], self.job.id))
		self.signals.changed.emit(self.job.id)

	def update_started(self):
		sql = "UPDATE jobs SET status=?,progress=?,started=? WHERE id=?"
		DB.query(sql, (3, 0, int(time.time()), self.job.id))
		self.signals.changed.emit(self.job.id)

	def update_stopped(self):
		sql = "UPDATE jobs SET status=? WHERE id=?"
		DB.query(sql, (2, self.job.id))
		self.signals.changed.emit(self.job.id)

	def update_finished(self):
		sql = "UPDATE jobs SET progress=?,finished=? WHERE id=?"
		DB.query(sql, (100, int(time.time()), self.job.id))
		self.signals.changed.emit(self.job.id)
		self.signals.finished.emit(self.job.id)

	def update_success(self):
		sql = "UPDATE jobs SET status=? WHERE id=?"
		DB.query(sql, (1, self.job.id))
		self.signals.changed.emit(self.job.id)

	def update_error(self, data):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (0, data['message'], self.job.id))
		self.signals.changed.emit(self.job.id)

	def write_pose_interactions(self, data):
		jid = data['job']
		best = data['best']
		poses = data['poses']
		interactions = data['interactions']

		pose_sql = "INSERT INTO pose VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
		DB.insert_rows(pose_sql, poses)

		pid_sql = "SELECT id FROM pose WHERE jid=?"
		pose_ids = DB.get_column(pid_sql, (jid,))

		best_sql = "INSERT INTO best VALUES (?,?)"
		DB.query(best_sql, (None, pose_ids[best]))

		best_count = DB.get_option('best_count')
		if best_count:
			best_count = int(best_count) + 1
		else:
			best_count = 1
		DB.set_option('best_count', best_count)

		site_sql = "INSERT INTO binding_site VALUES (?,?,?)"
		sites = interactions['binding_site']

		#convert index to real id 
		for site in sites:
			site[1] = pose_ids[site[1]]

		DB.insert_rows(site_sql, sites)

		bs_sql = "SELECT * from binding_site WHERE pid IN ({})".format(
			','.join(map(str, pose_ids))
		)
		site_mapping = {"{}:{}".format(pose_ids.index(row[1]), row[2]) : row[0] for row in DB.query(bs_sql)}

		cols_mapping = {
			'hydrogen_bond': 12,
			'halogen_bond': 10,
			'hydrophobic_interaction': 8,
			'water_bridge': 13,
			'salt_bridge': 9,
			'pi_stacking': 10,
			'pi_cation': 10,
			'metal_complex': 9
		}

		for k in interactions:
			if k == 'binding_site':
				continue

			if not interactions[k]:
				continue

			for i in range(len(interactions[k])):
				site = interactions[k][i][1]
				interactions[k][i][1] = site_mapping[site]

			sql = "INSERT INTO {} VALUES ({})".format(k, ','.join(['?']*cols_mapping[k]))
			DB.insert_rows(sql, interactions[k])

	def make_temp_dir(self):
		self.tempdir = QTemporaryDir()
		#open for test
		self.tempdir.setAutoRemove(False)

		if not self.tempdir.isValid():
			raise Exception(
				"Could not create temporary work directory, {}".format(
					self.tempdir.errorString()
				)
			)

		return self.tempdir.path()

	def get_commands(self):
		pass

	def start_process(self):
		tmpdir = self.make_temp_dir()
		cmds = self.get_commands()
		self.process = self.processer(self.job, self.params, cmds, tmpdir, self.pipe)
		self.process.start()

	def stop(self):
		proc = psutil.Process(self.process.pid)

		for child in proc.children(recursive=True):
			child.kill()

		self.process.kill()
		self.update_stopped()
		self.pipe.close()

	def call_response(self, data):
		if data['action'] == 'progress':
			self.update_progress(data)

		elif data['action'] == 'result':
			self.write_pose_interactions(data)

		elif data['action'] == 'error':
			self.update_error(data)

		elif data['action'] == 'success':
			self.update_success()

		elif data['action'] == 'finish':
			while True:
				if self.pipe.empty():
					self.pipe.close()
					break

				QThread.msleep(10)

	@Slot()
	def run(self):
		self.start_process()
		self.update_started()

		while True:
			try:
				data = self.pipe.get()
				self.call_response(data)

			except:
				break

		self.update_finished()

class AutodockWorker(BaseWorker):
	processer = AutodockProcess

	def get_commands(self):
		autodock = self.settings.value('Tools/autodock_4')
		autogrid = self.settings.value('Tools/autogrid_4')

		return autodock, autogrid

class AutodockVinaWorker(BaseWorker):
	processer = AutodockVinaProcess

	def get_commands(self):
		vina = self.settings.value('Tools/autodock_vina')
		autogrid = self.settings.value('Tools/autogrid_4')

		return vina, autogrid

class QuickVinaWorker(BaseWorker):
	processer = QuickVinaProcess

	def get_commands(self):
		vina = self.settings.value('Tools/quick_vina_w')
		return vina

class InteractionExportSignal(QObject):
	message = Signal(str)
	finished = Signal()

class InteractionExportWorker(QRunnable):
	def __init__(self, out_dir, table_models):
		super(InteractionExportWorker, self).__init__()
		self.setAutoDelete(True)
		self.out_dir = out_dir
		self.table_models = table_models
		self.signals = InteractionExportSignal()

	@property
	def sql(self):
		pass

	def run(self):
		for site in DB.query(self.sql):
			for table in self.table_models:
				rows = []
				sql = "SELECT * FROM {} WHERE bid={}".format(table, site[0])

				for row in DB.query(sql):
					res = [site[2], site[3], site[4], site[1]]
					res.extend(row[2:])
					rows.append(res)

				if rows:
					out_file = os.path.join(self.out_dir, '{}.csv'.format(table))

					if not os.path.exists(out_file):
						headers = ['Receptor', 'Ligand', 'Mode', 'Binding Site']
						headers.extend(self.table_models[table].custom_headers[2:])

						with open(out_file, 'w', newline='') as fw:
							writer = csv.writer(fw)
							writer.writerow(headers)
							writer.writerows(rows)
					else:
						with open(out_file, 'a', newline='') as fw:
							writer = csv.writer(fw)
							writer.writerows(rows)

			self.signals.message.emit("Export tables for {} vs {}".format(site[2], site[3]))

		self.signals.finished.emit()

class BestInteractionExportWorker(InteractionExportWorker):
	@property
	def sql(self):
		return (
			"SELECT s.id,s.site,r.name,l.name,p.run FROM binding_site AS s "
			"INNER JOIN best AS b ON b.pid=s.pid "
			"INNER JOIN pose AS p ON p.id=b.pid "
			"INNER JOIN jobs AS j ON j.id=p.jid "
			"LEFT JOIN molecular AS r ON r.id=j.rid "
			"LEFT JOIN molecular AS l ON l.id=j.lid"
		)

class PoseInteractionExportWorker(InteractionExportWorker):
	@property
	def sql(self):
		return (
			"SELECT s.id,s.site,r.name,l.name,p.run FROM binding_site AS s "
			"INNER JOIN pose AS p ON p.id=s.pid "
			"INNER JOIN jobs AS j ON j.id=p.jid "
			"LEFT JOIN molecular AS r ON r.id=j.rid "
			"LEFT JOIN molecular AS l ON l.id=j.lid"
		)

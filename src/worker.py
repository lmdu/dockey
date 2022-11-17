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
	'ImportSDFWorker'
]

class ImportSignals(QObject):
	success = Signal()
	failure = Signal(str)
	message = Signal(str)

class ImportWorker(QRunnable):
	processer = None

	def __init__(self, mol_files, mol_type=0):
		super(ImportWorker, self).__init__()
		self.setAutoDelete(True)
		self.signals = ImportSignals()
		self.mol_files = mol_files
		self.mol_type = mol_type
		self.consumer, self.producer = multiprocessing.Pipe(duplex=False)

	def delete_imported_molecules(self):
		DB.query("DELETE FROM molecular")

	def write_molecules(self, data):
		mol_type = ['', 'receptors', 'ligands'][self.mol_type]
		sql = "INSERT INTO molecular VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, data['rows'])
		self.signals.message.emit("Imported {} {}".format(data['total'], mol_type))

	def run(self):
		proc = self.processer(self.mol_files, self.mol_type, self.producer)
		proc.start()

		self.producer.close()

		while 1:
			try:
				data = self.consumer.recv()

				if data['action'] == 'insert':
					self.write_molecules(data)

				elif data['action'] == 'failure':
					self.signals.failure.emit(data['message'])
					#self.delete_imported_molecules()

			except EOFError:
				break

			except:
				self.signals.failure.emit(traceback.format_exc())
				break

			else:
				self.signals.success.emit()

class ImportFileWorker(ImportWorker):
	processer = ImportFileProcess

class ImportSDFWorker(ImportWorker):
	processer = ImportSDFProcess

class ImportFolderWorker(ImportWorker):
	processer = ImportFolderProcess

class JobListSignals(QObject):
	failure = Signal(str)
	message = Signal(str)
	success = Signal()
	finished = Signal()

class JobListGenerator(QRunnable):
	def __init__(self):
		super(JobListGenerator, self).__init__()
		self.setAutoDelete(True)
		self.signals = JobListSignals()
		self.consumer, self.producer = multiprocessing.Pipe(duplex=False)

	def write_jobs(self, data):
		sql = "INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, data['rows'])
		self.signals.message.emit("Generate {} jobs".format(data['total']))

	def run(self):
		rids = DB.get_column("SELECT id FROM molecular WHERE type=1")
		lids = DB.get_column("SELECT id FROM molecular WHERE type=2")

		proc = JobListProcess(rids, lids, self.producer)
		proc.start()

		self.producer.close()

		while 1:
			try:
				data = self.consumer.recv()

				if data['action'] == 'insert':
					self.write_jobs(data)

				if data['action'] == 'failure':
					self.signals.failure.emit(data['message'])

			except EOFError:
				break

			except:
				error = traceback.format_exc()
				self.signals.failure.emit(error)
				print(error)
				break

			else:
				self.signals.success.emit()

		self.signals.finished.emit()

class WorkerSignals(QObject):
	finished = Signal()
	refresh = Signal(int)
	message = Signal(str)
	failure = Signal(str)
	threads = Signal(int)
	stopjob = Signal(int)

class BaseWorker(QRunnable):
	processer = None

	def __init__(self, params):
		super(BaseWorker, self).__init__()
		self.jobs = {}
		self.params = params
		self.signals = WorkerSignals()
		self.settings = QSettings()
		self.consumer, self.producer = multiprocessing.Pipe(duplex=False)
		self.job_query = None
		self.job_num = self.settings.value('Job/concurrent', 1, int)
		self.signals.threads.connect(self.change_job_numbers)
		self.signals.stopjob.connect(self.stop_job)
		self.setAutoDelete(True)

	def exit(self):
		self.job_query = None
		self.producer.close()

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
			if row[2] not in residues:
				residues[row[2]] = []

			residues[row[2]].append("{}{}".format(row[3], row[4]))


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

	@Slot()
	def change_job_numbers(self, num):
		self.job_num = num
		num = len(self.jobs)

		if self.job_num > num:
			for i in range(self.job_num - num):
				self.submit_job()

	def update_progress(self, data):
		sql = "UPDATE jobs SET progress=? WHERE id=?"
		DB.query(sql, (data['message'], data['job']))
		self.signals.refresh.emit(data['job'])

	def update_started(self, data):
		sql = "UPDATE jobs SET status=?,progress=?,started=? WHERE id=?"
		DB.query(sql, (2, 0, data['message'], data['job']))
		self.signals.refresh.emit(data['job'])

	def update_finished(self, data):
		sql = "UPDATE jobs SET progress=?,finished=? WHERE id=?"
		DB.query(sql, (100, data['message'], data['job']))
		self.signals.refresh.emit(data['job'])
		self.signals.finished.emit()

	def update_success(self, data):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (1, data['message'], data['job']))
		self.signals.refresh.emit(data['job'])

	def update_error(self, data):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (0, data['message'], data['job']))
		self.signals.refresh.emit(data['job'])

	def get_prepare_params(self):
		tool = self.settings.value('Ligand/prepare_tool', 'prepare_ligand4')

		if tool == 'prepare_ligand4':
			self.settings.beginGroup('Ligand')
			lig_params = dict(
				tool = 'prepare_ligand4',
				repairs = self.settings.value('repairs', ''),
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
				rigid_macrocycles = self.settings.value('rigid_macrocycles', False, bool),
				keep_chorded_rings = self.settings.value('keep_chorded_rings', False, bool),
				keep_equivalent_rings = self.settings.value('keep_equivalent_rings', False, bool),
				hydrate = self.settings.value('hydrate', False, bool),
				keep_nonpolar_hydrogens = self.settings.value('keep_nonpolar_hydrogens', False, bool),
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
			repairs = self.settings.value('repairs', ''),
			charges_to_add = self.settings.value('charges_to_add', 'gasteiger'),
			#preserve_charge_types = self.settings.value('preserve_charge_types', ''),
			cleanup = self.settings.value('cleanup', 'nphs_lps_waters_nonstdres'),
			delete_single_nonstd_residues = self.settings.value('delete_single_nonstd_residues', False, bool),
			unique_atom_names = self.settings.value('unique_atom_names', False, bool)
		)
		self. settings.endGroup()

		return lig_params, rep_params

	def get_all_params(self):
		lig_params, rep_params = self.get_prepare_params()
		return [self.params, lig_params, rep_params]

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
		temp_dir = QTemporaryDir()
		temp_dir.setAutoRemove(False)

		if not temp_dir.isValid():
			raise Exception(
				"Could not create temporary work directory, {}".format(
					temp_dir.errorString()
				)
			)

		return temp_dir

	def start_process(self, jid):
		temp_dir = self.make_temp_dir()
		params = self.get_all_params()
		cmds = self.get_commands()
		job = self.get_job(jid)
		proc = self.processer(job, params, cmds, temp_dir.path(), self.producer)
		proc.start()

		self.jobs[job.id] = AttrDict(
			tempdir = temp_dir,
			process = proc
		)

	def submit_job(self):
		if self.job_query is None:
			return

		if len(self.jobs) > self.job_num:
			return

		row = self.job_query.fetchone()

		if row:
			self.start_process(row[0])

	def delete_job(self, job):
		obj = self.jobs.pop(job)
		obj.tempdir.remove()

	def stop_job(self, job):
		if job not in self.jobs:
			return

		obj = self.jobs.pop(job)

		parent = psutil.Process(obj.process.pid)
		children = parent.children(recursive=True)

		for child in children:
			child.kill()

		obj.process.join()
		obj.tempdir.remove()

		if len(self.jobs) < self.job_num:
			self.submit_job()

	def get_job_pid(self, job):
		if job not in self.jobs:
			return

		return self.jobs[job].process.pid

	@Slot()
	def run(self):
		self.job_query = DB.query("SELECT id FROM jobs")

		try:
			for i in range(self.job_num):
				self.submit_job()

			while 1:
				try:
					data = self.consumer.recv()

					if data['action'] == 'start':
						self.update_started(data)

					elif data['action'] == 'progress':
						self.update_progress(data)

					elif data['action'] == 'error':
						self.update_error(data)

					elif data['action'] == 'success':
						self.update_success(data)

					elif data['action'] == 'result':
						self.write_pose_interactions(data)

					elif data['action'] == 'finish':
						self.update_finished(data)
						self.delete_job(data['job'])
						self.submit_job()

				except EOFError:
					break

				except:
					error = traceback.format_exc()
					raise Exception(error)
					break

		except:
			error = traceback.format_exc()
			self.signals.message.emit(error)
			print(error)

		finally:
			for job in self.jobs:
				self.jobs[job].tempdir.remove()

			self.signals.finished.emit()

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

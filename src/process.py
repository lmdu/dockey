import os
import time
import psutil
import traceback
import subprocess
import multiprocessing

from PySide6.QtCore import *

from utils import *
from param import *
from prepare import *
from backend import *

__all__ = ['ImportFileProcess', 'ImportFolderProcess', 'JobListProcess',
	'AutodockProcess', 'AutodockVinaProcess', 'QuickVinaProcess'
]

class ImportProcess(multiprocessing.Process):
	def __init__(self, mol_files, mol_type, producer):
		super(ImportProcess, self).__init__()
		self.daemon = True
		self.mol_files = mol_files
		self.mol_type = mol_type
		self.producer = producer
		self.mol_list = []
		self.mol_count = 0

	def write(self):
		mol_type = ['', 'receptors', 'ligands'][self.mol_type]
		sql = "INSERT INTO molecular VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

		if self.mol_list:
			self.mol_count += len(self.mol_list)
			self.producer.send({
				'action': 'insert',
				'rows': self.mol_list,
				'total': self.mol_count
			})		
			self.mol_list = []

	def process(self):
		pass

	def run(self):
		try:
			self.process()
		
		except:
			self.producer.send({
				'action': 'failure',
				'message': traceback.format_exc()
			})
			
		finally:
			self.producer.close()

class ImportFileProcess(ImportProcess):
	def process(self):
		for mol_file in self.mol_files:
			m = get_molecule_information(mol_file)

			if 'error' in m:
				self.producer.send({
					'action': 'failure',
					'message': "{}: {}".format(mol_file, m.error)
				})
				return

			self.mol_list.append([None, m.name, self.mol_type, m.content, m.format,
				m.atoms, m.bonds, m.hvyatoms, m.residues, m.rotors, m.formula,
				m.energy, m.weight, m.logp
			])

			if len(self.mol_list) == 200:
				self.write()

		self.write()

class ImportFolderProcess(ImportProcess):
	def process(self):
		mol_files = QDirIterator(self.mol_files, QDir.Files)

		while mol_files.hasNext():
			mol_file = mol_files.next()

			m = get_molecule_information(mol_file)

			if 'error' in m:
				self.producer.send({
					'action': 'failure',
					'message': "{}: {}".format(mol_file, m.error)
				})
				return

			self.mol_list.append([None, m.name, self.mol_type, m.content, m.format,
				m.atoms, m.bonds, m.hvyatoms, m.residues, m.rotors, m.formula,
				m.energy, m.weight, m.logp
			])

			if len(self.mol_list) == 200:
				self.write()

		self.write()

class JobListProcess(multiprocessing.Process):
	def __init__(self, rids, lids, producer):
		super(JobListProcess, self).__init__()
		self.daemon = True
		self.rids = rids
		self.lids = lids
		self.job_list = []
		self.job_count = 0
		self.producer = producer

	def write(self):
		if self.job_list:
			self.job_count += len(self.job_list)
			self.producer.send({
				'action': 'insert',
				'rows': self.job_list,
				'total': self.job_count
			})
			self.job_list = []

	def process(self):
		for r in self.rids:
			for l in self.lids:
				self.job_list.append((None, r, l, 3, 0, 0, 0, ''))

				if len(self.job_list) == 200:
					self.write()

		self.write()

	def run(self):
		try:
			self.process()

		except:
			error = traceback.format_exc()
			self.producer.send({
				'action': 'failure',
				'message': error
			})
			print(error)

		finally:
			self.producer.close()

class BaseProcess(multiprocessing.Process):
	def __init__(self, job, params, cmds, work_dir, producer):
		super(BaseProcess, self).__init__()
		self.daemon = True
		self.job = job
		self.producer = producer
		self.cmds = cmds
		self.work_dir = work_dir
		self.dock_params = params[0]
		self.ligp_params = params[1]
		self.repp_params = params[2]
		#self.grid_spacing = params[3]
		
		#if os.name == 'nt':
		#	self.creationflags = 0x08000000
		#else:
		#	self.creationflags = 0

	def send_result(self, poses):
		interactions = get_complex_interactions(poses)

		#get min energy
		e = 0
		best = 0
		for i, pose in enumerate(poses):
			if pose[3] < e:
				e = pose[3]
				best = i

		self.producer.send({
			'job': self.job.id,
			'action': 'result',
			'best': best,
			'poses': poses,
			'interactions': interactions
		})

	def send_message(self, mtype, message):
		self.producer.send({
			'job': self.job.id,
			'action': mtype,
			'message': message
		})

	def update_progress(self, progress):
		self.send_message('progress', progress)

	def update_started(self):
		started = int(time.time())
		self.send_message('start', started)

	def update_finished(self):
		finished = int(time.time())
		self.send_message('finish', finished)

	def update_error(self, error):
		self.send_message('error', error)

	def update_success(self):
		self.send_message('success', 'Done')

	def prepare_receptor(self):
		rpdbqt = os.path.join(self.work_dir, "{}.pdbqt".format(self.job.rn))

		if self.job.rf in ['pdb', 'mol2']:
			content = self.job.rc
			rfile = os.path.join(self.work_dir, "{}.{}".format(self.job.rn, self.job.rf))

		else:
			content = convert_string_to_pdb(self.job.rc, self.job.rf)
			rfile = os.path.join(self.work_dir, "{}.pdb".format(self.job.rn))

		with open(rfile, 'w') as fw:
			fw.write(content)

		prepare_autodock_receptor(rfile, rpdbqt, self.repp_params)

		return rpdbqt

	def prepare_ligand(self):
		lpdbqt = os.path.join(self.work_dir, "{}.pdbqt".format(self.job.ln))

		if self.ligp_params['tool'] == 'meeko' and self.job.lf not in ['pdb', 'mol', 'mol2', 'sdf']:
			content = convert_string_to_pdb(self.job.lc, self.job.lf)
			lfile = os.path.join(self.work_dir, "{}.pdb".format(self.job.ln))

		elif self.ligp_params['tool'] == 'prepare_ligand4' and self.job.lf not in ['pdb', 'mol2']:
			content = convert_string_to_pdb(self.job.lc, self.job.lf)
			lfile = os.path.join(self.work_dir, "{}.pdb".format(self.job.ln))

		else:
			content = self.job.lc
			lfile = os.path.join(self.work_dir, "{}.{}".format(self.job.ln, self.job.lf))

		with open(lfile, 'w') as fw:
			fw.write(content)

		prepare_ligand(lfile, lpdbqt, self.ligp_params)

		return lpdbqt

	def run(self):
		self.update_started()

		try:
			self.do()

		except:
			error = traceback.format_exc()
			self.update_error(error)
			print(error)

		else:
			self.update_success()

		finally:
			self.update_finished()
			self.producer.close()

class AutodockProcess(BaseProcess):
	def do(self):
		#convert receptor and ligand to pdbqt format
		rpdbqt = self.prepare_receptor()
		lpdbqt = self.prepare_ligand()

		#set autogrid parameter and create gpf parameter file
		ag_param = AutogridParameter(rpdbqt, lpdbqt, (self.job.x, self.job.y, self.job.z),
			(self.job.cx, self.job.cy, self.job.cz), self.job.spacing
		)
		gpf_file = ag_param.make_gpf_file()
		glg_file = gpf_file.replace('.gpf', '.glg')
		self.log_file = glg_file
		#self.save_log('gpf_file', gpf_file)

		autodock, autogrid = self.cmds

		#self.update_message("Running autogrid")
		proc = psutil.Popen([autogrid, '-p', gpf_file, '-l', glg_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8'
			#creationflags = self.creationflags
		)

		read_start = 0

		while proc.poll() is None:
			if QFileInfo(glg_file).size() > read_start:
				p = 0

				with open(glg_file) as log:
					log.seek(read_start)

					for line in log:
						if '%' in line:
							p = round(float(line.strip().split()[2].strip('%'))*0.05, 2)

					read_start = log.tell()

				if p > 0:
					self.update_progress(p)
			else:
				time.sleep(1)

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		#self.save_log('glg_file', glg_file)

		#set autodock4 parameters and make dpf parameter file
		dpf_file = self.dock_params.make_dpf_file(rpdbqt, lpdbqt)
		dlg_file = dpf_file.replace('.dpf', '.dlg')

		#run autodock4
		#self.update_message("Running autodock")
		proc = psutil.Popen([autodock, '-p', dpf_file, '-l', dlg_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8',
			#creationflags = self.creationflags
		)

		read_start = 0

		while proc.poll() is None:
			if QFileInfo(dlg_file).size() > read_start:
				p = 0

				with open(dlg_file) as log:
					log.seek(read_start)

					for line in log:
						if line.startswith('Run:'):
							cols = line.strip().split()
							p = round(5+int(cols[1])/int(cols[9])*90, 2)

					read_start = log.tell()

				if p > 0:
					self.update_progress(p)
			else:
				time.sleep(1)

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		#self.save_log('dlg_file', dlg_file)

		#parse autodock dlg log file
		#self.update_message("Analyzing docking results")
		runs = {}
		rows = []

		receptor = convert_pdbqt_to_pdb(rpdbqt, as_string=False)

		with open(dlg_file) as fh:
			for line in fh:
				if 'RMSD TABLE' in line:
					break

				if line.startswith("DOCKED: MODEL"):
					rid = int(line.strip().split()[2])
					runs[rid] = ['', []]

				elif line.startswith("DOCKED: USER    Estimated Inhibition Constant"):
					runs[rid][0] = line.split('=')[1].strip().split('(')[0].strip()

				if line.startswith("DOCKED:"):
					runs[rid][1].append(line[8:])

			for line in fh:
				if 'INFORMATION ENTROPY ANALYSIS' in line:
					break

				cols = line.strip().split()

				if len(cols) == 7:
					rid = int(cols[2])
					#rank = int(cols[0])
					energy = float(cols[3])
					ki = runs[rid][0]
					crmsd = float(cols[4])
					rrmsd = float(cols[5])
					pose = convert_pdbqt_to_pdb(''.join(runs[rid][1]))
					comp = generate_complex_pdb(receptor, pose)
					row = [None, self.job.id, rid, energy, crmsd, rrmsd]
					lea = ligand_efficiency_assessment(pose, energy, ki)
					row.extend(lea)
					row.append(pose)
					row.append(comp)
					rows.append(row)

		#pids = self.save_pose(rows)
		#self.save_interaction(pids, rows)
		self.send_result(rows)

class AutodockVinaProcess(BaseProcess):
	def do(self):
		#convert receptor and ligand to pdbqt format
		rpdbqt = self.prepare_receptor()
		lpdbqt = self.prepare_ligand()

		#get commands
		vina, autogrid = self.cmds

		#if use autodock scoring function
		#set autogrid parameter and create gpf parameter file
		if self.dock_params.scoring.value == 'ad4':
			ag_param = AutogridParameter(rpdbqt, lpdbqt, (self.job.x, self.job.y, self.job.z),
				(self.job.cx, self.job.cy, self.job.cz), self.job.spacing
			)
			gpf_file = ag_param.make_gpf_file()
			glg_file = gpf_file.replace('.gpf', '.glg')
			self.log_file = glg_file

			#self.save_log('gpf_file', gpf_file)

			#self.update_message("Running autogrid")
			proc = psutil.Popen([autogrid, '-p', gpf_file, '-l', glg_file],
				stdin = subprocess.PIPE,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				cwd = self.work_dir,
				encoding = 'utf8',
				#creationflags = self.creationflags
			)

			read_start = 0

			while proc.poll() is None:
				if QFileInfo(glg_file).size() > read_start:
					p = 0

					with open(glg_file) as log:
						log.seek(read_start)

						for line in log:
							if '%' in line:
								p = round(float(line.strip().split()[2].replace('%', ''))*0.05, 2)

						read_start = log.tell()

					if p > 0:
						self.update_progress(p)
				else:
					time.sleep(1)

			if proc.returncode != 0:
				raise Exception(proc.stderr.read())

			#self.save_log('glg_file', glg_file)

			self.dock_params.maps.value = self.job.rn
		else:
			#using autodock score function no need --receptor arg
			self.dock_params.receptor.value = os.path.basename(rpdbqt)

		#set autodock vina parameters and make config file
		config_file = os.path.join(self.work_dir, "config.txt".format(self.job.rn))
		out_file = os.path.join(self.work_dir, "out.pdbqt")
		log_file = os.path.join(self.work_dir, "out.log")

		self.dock_params.ligand.value = os.path.basename(lpdbqt)
		self.dock_params.center_x.value = self.job.cx
		self.dock_params.center_y.value = self.job.cy
		self.dock_params.center_z.value = self.job.cz
		self.dock_params.size_x.value = self.job.x
		self.dock_params.size_y.value = self.job.y
		self.dock_params.size_z.value = self.job.z
		self.dock_params.spacing.value = self.job.spacing
		self.dock_params.out.value = os.path.basename(out_file)
		self.dock_params.make_config_file(config_file)

		#self.save_log('config_file', config_file)

		#run autodock vina
		#self.update_message("Running autodock vina")
		proc = psutil.Popen([vina, '--config', config_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8',
			#creationflags = self.creationflags
		)

		star = 0
		log_fw = open(log_file, 'w')

		while proc.poll() is None:
			chars = proc.stdout.read(5)

			for char in chars:
				if char == '*':
					star += 1
					
			p = round(star/51*90, 2)
			self.update_progress(p)

			if chars:
				log_fw.write(chars)
			else:
				time.sleep(0.5)

		log_fw.write(proc.stdout.read())
		log_fw.close()

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		#self.save_log('log_file', log_file)
		#self.save_log('out_file', out_file)

		#self.update_message("Analyzing docking results")
		rows = []
		with open(log_file) as fh:
			for line in fh:
				if line.startswith('-----+'):
					break

			for line in fh:
				cols = line.strip().split()
				run = int(cols[0])
				energy = float(cols[1])
				rmsd1 = float(cols[2])
				rmsd2 = float(cols[3])

				rows.append([None, self.job.id, run, energy, rmsd1, rmsd2])

		modes = {}
		with open(out_file) as fh:
			for line in fh:
				if line.startswith('MODEL'):
					lines = [line]
					idx = int(line.strip().split()[1]) - 1

				elif line.startswith('ENDMDL'):
					lines.append(line)
					mode = convert_pdbqt_to_pdb(''.join(lines))
					modes[idx] = mode
				else:
					lines.append(line)

		receptor = convert_pdbqt_to_pdb(rpdbqt, as_string=False)

		for i, row in enumerate(rows):
			mode = modes.get(i, '')

			if mode:
				comp = generate_complex_pdb(receptor, mode)
			else:
				comp = ''

			lea = ligand_efficiency_assessment(mode, row[3])
			rows[i].extend(lea)
			rows[i].append(mode)
			rows[i].append(comp)

		#pids = self.save_pose(rows)
		#self.save_interaction(pids, rows)
		self.send_result(rows)

class QuickVinaProcess(BaseProcess):
	def do(self):
		#convert receptor and ligand to pdbqt format
		rpdbqt = self.prepare_receptor()
		lpdbqt = self.prepare_ligand()

		#get commands
		vina = self.cmds
		self.dock_params.receptor.value = os.path.basename(rpdbqt)

		#set autodock vina parameters and make config file
		config_file = os.path.join(self.work_dir, "config.txt".format(self.job.rn))
		out_file = os.path.join(self.work_dir, "out.pdbqt")
		log_file = os.path.join(self.work_dir, "out.log")

		self.dock_params.ligand.value = os.path.basename(lpdbqt)
		self.dock_params.center_x.value = self.job.cx
		self.dock_params.center_y.value = self.job.cy
		self.dock_params.center_z.value = self.job.cz
		self.dock_params.size_x.value = self.job.x
		self.dock_params.size_y.value = self.job.y
		self.dock_params.size_z.value = self.job.z
		self.dock_params.out.value = os.path.basename(out_file)
		self.dock_params.make_config_file(config_file)

		#self.save_log('config_file', config_file)

		#run quick vina
		#self.update_message("Running quick vina")
		proc = psutil.Popen([vina, '--config', config_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8',
			#creationflags = self.creationflags
		)

		star = 0
		log_fw = open(log_file, 'w')

		while proc.poll() is None:
			chars = proc.stdout.read(5)

			for char in chars:
				if char == '*':
					star += 1
					
			p = round(star/51*90, 2)
			self.update_progress(p)

			if chars:
				log_fw.write(chars)
			else:
				time.sleep(0.5)

		log_fw.write(proc.stdout.read())
		log_fw.close()

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		#self.save_log('log_file', log_file)
		#self.save_log('out_file', out_file)

		#self.update_message("Analyzing docking results")
		rows = []
		with open(log_file) as fh:
			for line in fh:
				if line.startswith('-----+'):
					break

			for line in fh:
				if line.startswith('Writing output'):
					break

				cols = line.strip().split()
				run = int(cols[0])
				energy = float(cols[1])
				rmsd1 = float(cols[2])
				rmsd2 = float(cols[3])

				rows.append([None, self.job.id, run, energy, rmsd1, rmsd2])

		modes = {}
		with open(out_file) as fh:
			for line in fh:
				if line.startswith('MODEL'):
					lines = [line]
					idx = int(line.strip().split()[1]) - 1

				elif line.startswith('ENDMDL'):
					lines.append(line)
					mode = convert_pdbqt_to_pdb(''.join(lines))
					modes[idx] = mode
				else:
					lines.append(line)

		receptor = convert_pdbqt_to_pdb(rpdbqt, as_string=False)

		for i, row in enumerate(rows):
			mode = modes.get(i, '')

			if mode:
				comp = generate_complex_pdb(receptor, mode)
			else:
				comp = ''

			lea = ligand_efficiency_assessment(mode, row[3])
			rows[i].extend(lea)
			rows[i].append(mode)
			rows[i].append(comp)

		#pids = self.save_pose(rows)
		#self.save_interaction(pids, rows)
		self.send_result(rows)

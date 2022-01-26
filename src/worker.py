import os
import time
import psutil
import traceback
import subprocess

from PySide6.QtCore import *

from param import *
from utils import *
from backend import *
from prepare import *

__all__ = ['AutodockWorker', 'AutodockVinaWorker']

class WorkerSignals(QObject):
	#finished = Signal()
	#error = Signal(str)
	#progress = Signal(int)
	refresh = Signal(int)

class BaseWorker(QRunnable):
	def __init__(self, job, params):
		super(BaseWorker, self).__init__()
		self.job = self.get_job(job)
		self.params = params
		self.signals = WorkerSignals()
		self.settings = QSettings()

	def get_job(self, _id):
		sql = (
			"SELECT j.id,j.rid,j.lid,m1.name,m1.pdb,m2.name,m2.pdb "
			"FROM jobs AS j LEFT JOIN molecular AS m1 ON m1.id=j.rid "
			"LEFT JOIN molecular AS m2 ON m2.id=j.lid WHERE j.id=?"
		)
		job = DB.get_row(sql, (_id,))
		grid = DB.get_row("SELECT * FROM grid WHERE rid=?", (job[1],))

		return AttrDict({
			'id': job[0],
			'rid': job[1],
			'lid': job[2],
			'rname': job[3],
			'rpdb': job[4],
			'lname': job[5],
			'lpdb': job[6],
			'x': grid[2],
			'y': grid[3],
			'z': grid[4],
			'cx': grid[5],
			'cy': grid[6],
			'cz': grid[7],
			'spacing': grid[8]
		})

	def save_pose(self, poses):
		sql = "INSERT INTO pose VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, poses)

	def save_log(self, name, log_file):
		sql = "INSERT INTO logs VALUES (?,?,?,?)"

		with open(log_file) as fh:
			content = fh.read()

		logs = [(None, self.job.id, name, content)]
		DB.insert_rows(sql, logs)

	def update_status(self, status):
		sql = "UPDATE jobs SET status=? WHERE id=?"
		DB.query(sql, (status, self.job.id))
		self.signals.refresh.emit(self.job.id)

	def update_progress(self, progress):
		sql = "UPDATE jobs SET progress=? WHERE id=?"
		DB.query(sql, (progress, self.job.id))
		self.signals.refresh.emit(self.job.id)

	def update_message(self, message):
		sql = "UPDATE jobs SET message=? WHERE id=?"
		DB.query(sql, (message, self.job.id))
		self.signals.refresh.emit(self.job.id)

	def update_started(self):
		started = int(time.time())
		sql = "UPDATE jobs SET status=?,progress=?,started=? WHERE id=?"
		DB.query(sql, (2, 0, started, self.job.id))
		self.signals.refresh.emit(self.job.id)

	def update_finished(self):
		finished = int(time.time())
		sql = "UPDATE jobs SET progress=?,finished=? WHERE id=?"
		DB.query(sql, (100, finished, self.job.id))
		self.signals.refresh.emit(self.job.id)

	def update_error(self, error):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (3, error, self.job.id))
		self.signals.refresh.emit(self.job.id)

	def update_success(self):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (1, "The job was successfully finished", self.job.id))
		self.signals.refresh.emit(self.job.id)

	def pipline(self):
		pass


	@Slot()
	def run(self):
		self.update_started()

		try:
			#make a temp work dir
			self.temp_dir = QTemporaryDir()
			if not self.temp_dir.isValid():
				raise Exception(
					"Could not create temporary work directory, {}".format(
						self.temp_dir.errorString()
					)
				)
			self.work_dir = self.temp_dir.path()

			#run worker
			self.pipline()
		except:
			self.update_error(traceback.format_exc())
		else:
			self.update_success()
		finally:
			self.temp_dir.remove()
			self.update_finished()

	def get_mgltools(self):
		mgltools_path = self.settings.value('Tools/MGLTools')
		py27 = os.path.join(mgltools_path, 'python.exe')
		script_path = os.path.join(mgltools_path, 'Lib', 'site-packages',
			'AutoDockTools', 'Utilities24')
		return py27, script_path

	def execute(self, cmd, args, wkdir):
		proc = QProcess()
		proc.setWorkingDirectory(wkdir)
		proc.start(cmd, args)
		proc.waitForFinished(-1)

	def prepare_receptor(self, infile, outfile):
		#py27, folder = self.get_mgltools()
		#script = os.path.join(folder, 'prepare_receptor4.py')
		#args = [script, '-r', infile, '-o', outfile, '-A', 'bonds_hydrogens']
		#self.execute(py27, args, os.path.dirname(infile))
		prepare_autodock_receptor(infile, outfile)

	def prepare_ligand(self, infile, outfile):
		#py27, folder = self.get_mgltools()
		#script = os.path.join(folder, 'prepare_ligand4.py')
		#args = [script, '-l', infile, '-o', outfile, '-A', 'bonds_hydrogens']
		#self.execute(py27, args, os.path.dirname(infile))
		prepare_autodock_ligand(infile, outfile)

class AutodockWorker(BaseWorker):
	def __init__(self, job, params):
		super(AutodockWorker, self).__init__(job, params)

	def get_commands(self):
		autodock = self.settings.value('Tools/autodock_4')
		autogrid = self.settings.value('Tools/autogrid_4')

		return autodock, autogrid

	def pipline(self):
		#convert receptor and ligand to pdbqt format
		rpdb = os.path.join(self.work_dir, "{}.pdb".format(self.job.rname))
		rpdbqt = os.path.join(self.work_dir, "{}.pdbqt".format(self.job.rname))

		with open(rpdb, 'w') as fw:
			fw.write(self.job.rpdb)

		self.prepare_receptor(rpdb, rpdbqt)

		lpdb = os.path.join(self.work_dir, "{}.pdb".format(self.job.lname))
		lpdbqt = os.path.join(self.work_dir, "{}.pdbqt".format(self.job.lname))

		with open(lpdb, 'w') as fw:
			fw.write(self.job.lpdb)
		
		self.prepare_ligand(lpdb, lpdbqt)

		#set autogrid parameter and create gpf parameter file
		ag_param = AutogridParameter(rpdbqt, lpdbqt, (self.job.x, self.job.y, self.job.z),
			(self.job.cx, self.job.cy, self.job.cz), self.job.spacing
		)
		gpf_file = ag_param.make_gpf_file()
		glg_file = gpf_file.replace('.gpf', '.glg')
		self.log_file = glg_file
		self.save_log('gpf_file', gpf_file)

		autodock, autogrid = self.get_commands()

		self.update_message("Running autogrid")
		proc = psutil.Popen([autogrid, '-p', gpf_file, '-l', glg_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8',
			creationflags = 0x08000000
		)

		read_start = 0

		while proc.poll() is None:
			if QFileInfo(glg_file).size() > read_start:
				p = 0

				with open(glg_file) as log:
					log.seek(read_start)

					for line in log:
						if '%' in line:
							p = round(float(line.strip().split()[2].replace('%', ''))*0.09, 2)

					read_start = log.tell()

				if p > 0:
					self.update_progress(p)
			else:
				QThread.sleep(1)

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		self.save_log('glg_file', glg_file)

		#set autodock4 parameters and make dpf parameter file
		dpf_file = self.params.make_dpf_file(rpdbqt, lpdbqt)
		dlg_file = dpf_file.replace('.dpf', '.dlg')

		#run autodock4
		self.update_message("Running autodock")
		proc = psutil.Popen([autodock, '-p', dpf_file, '-l', dlg_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8',
			creationflags = 0x08000000
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
							p = round(int(cols[1])/int(cols[9])*90, 2)

					read_start = log.tell()

				if p > 0:
					self.update_progress(p)
			else:
				QThread.sleep(1)

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		self.save_log('dlg_file', dlg_file)

		#parse autodock dlg log file
		self.update_message("Analyzing docking results")
		runs = {}
		rows = []

		with open(dlg_file) as fh:
			for line in fh:
				if 'RMSD TABLE' in line:
					break

				if line.startswith("DOCKED: MODEL"):
					rid = int(line.strip().split()[2])
					runs[rid] = [0, []]

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
					logki = convert_ki_to_log(ki)
					crmsd = float(cols[4])
					rrmsd = float(cols[5])
					pose = convert_pdbqt_to_pdb(''.join(runs[rid][1]))

					row = [None, self.job.id, rid, energy, crmsd, rrmsd]
					lea = ligand_efficiency_assessment(pose, energy, ki)
					row.extend(lea)
					row.append(pose)
					rows.append(row)

		self.save_pose(rows)

class AutodockVinaWorker(BaseWorker):
	def __init__(self, job, params):
		super(AutodockVinaWorker, self).__init__(job, params)

	def get_commands(self):
		vina = self.settings.value('Tools/autodock_vina')
		autogrid = self.settings.value('Tools/autogrid_4')

		return vina, autogrid

	def pipline(self):
		#convert receptor and ligand to pdbqt format
		rpdb = os.path.join(self.work_dir, "{}.pdb".format(self.job.rname))
		rpdbqt = os.path.join(self.work_dir, "{}.pdbqt".format(self.job.rname))

		with open(rpdb, 'w') as fw:
			fw.write(self.job.rpdb)

		self.prepare_receptor(rpdb, rpdbqt)

		lpdb = os.path.join(self.work_dir, "{}.pdb".format(self.job.lname))
		lpdbqt = os.path.join(self.work_dir, "{}.pdbqt".format(self.job.lname))

		with open(lpdb, 'w') as fw:
			fw.write(self.job.lpdb)

		self.prepare_ligand(lpdb, lpdbqt)

		#get commands
		vina, autogrid = self.get_commands()

		#if use autodock scoring function
		#set autogrid parameter and create gpf parameter file
		if self.params.scoring.value == 'ad4':
			ag_param = AutogridParameter(rpdbqt, lpdbqt, (self.job.x, self.job.y, self.job.z),
				(self.job.cx, self.job.cy, self.job.cz), self.job.spacing
			)
			gpf_file = ag_param.make_gpf_file()
			glg_file = gpf_file.replace('.gpf', '.glg')
			self.log_file = glg_file

			self.save_log('gpf_file', gpf_file)

			self.update_message("Running autogrid")
			proc = psutil.Popen([autogrid, '-p', gpf_file, '-l', glg_file],
				stdin = subprocess.PIPE,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				cwd = self.work_dir,
				encoding = 'utf8',
				creationflags = 0x08000000
			)

			read_start = 0

			while proc.poll() is None:
				if QFileInfo(glg_file).size() > read_start:
					p = 0

					with open(glg_file) as log:
						log.seek(read_start)

						for line in log:
							if '%' in line:
								p = round(float(line.strip().split()[2].replace('%', ''))*0.09, 2)

						read_start = log.tell()

					if p > 0:
						self.update_progress(p)
				else:
					QThread.sleep(1)

			if proc.returncode != 0:
				raise Exception(proc.stderr.read())

			self.save_log('glg_file', glg_file)

			self.params.maps.value = self.job.rname
		else:
			#using autodock score function no need --receptor arg
			self.params.receptor.value = os.path.basename(rpdbqt)

		#set autodock vina parameters and make config file
		config_file = os.path.join(self.work_dir, "config.txt".format(self.job.rname))
		out_file = os.path.join(self.work_dir, "out.pdbqt")
		log_file = os.path.join(self.work_dir, "out.log")

		self.params.ligand.value = os.path.basename(lpdbqt)
		self.params.center_x.value = self.job.cx
		self.params.center_y.value = self.job.cy
		self.params.center_z.value = self.job.cz
		self.params.size_x.value = self.job.x
		self.params.size_y.value = self.job.y
		self.params.size_z.value = self.job.z
		self.params.spacing.value = self.job.spacing
		self.params.out.value = os.path.basename(out_file)
		self.params.make_config_file(config_file)

		self.save_log('config_file', config_file)

		#run autodock vina
		self.update_message("Running autodock vina")
		proc = psutil.Popen([vina, '--config', config_file],
			stdin = subprocess.PIPE,
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = self.work_dir,
			encoding = 'utf8',
			creationflags = 0x08000000
		)

		star = 0
		log_fw = open(log_file, 'w')

		while proc.poll() is None:
			chars = proc.stdout.read(5)

			for char in chars:
				if char == '*':
					star += 1
					
			p = round(star/51*100, 2)
			self.update_progress(p)

			if chars:
				log_fw.write(chars)
			else:
				QThread.msleep(500)

		log_fw.write(proc.stdout.read())
		log_fw.close()

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		self.save_log('log_file', log_file)
		self.save_log('out_file', out_file)

		self.update_message("Analyzing docking results")
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

		for i, row in enumerate(rows):
			mode = modes.get(i, '')
			lea = ligand_efficiency_assessment(mode, row[3])
			rows[i].extend(lea)
			rows[i].append(mode)

		self.save_pose(rows)


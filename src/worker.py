import os
import time
import psutil
import traceback
import subprocess

from PySide6.QtSql import *
from PySide6.QtCore import *

from param import *
from utils import *
from backend import *

__all__ = ['AutodockWorker', 'AutodockVinaWorker']

class WorkerSignals(QObject):
	#finished = Signal()
	#error = Signal(str)
	#progress = Signal(int)
	refresh = Signal(int)

class BaseWorker(QRunnable):
	def __init__(self, *args, **kwargs):
		super(BaseWorker, self).__init__()
		self.args = args
		self.kwargs = kwargs
		self.signals = WorkerSignals()
		self.settings = QSettings()

		self.job_id = args[0].id

	def update_status(self, status):
		sql = "UPDATE jobs SET status=? WHERE id=?"
		DB.query(sql, (status, self.job_id))
		self.signals.refresh.emit(self.job_id)

	def update_progress(self, progress):
		sql = "UPDATE jobs SET progress=? WHERE id=?"
		DB.query(sql, (progress, self.job_id))
		self.signals.refresh.emit(self.job_id)

	def update_message(self, message):
		sql = "UPDATE jobs SET message=? WHERE id=?"
		DB.query(sql, (message, self.job_id))
		self.signals.refresh.emit(self.job_id)

	def update_started(self):
		started = int(time.time())
		sql = "UPDATE jobs SET status=?,progress=?,started=? WHERE id=?"
		DB.query(sql, (2, 0, started, self.job_id))
		self.signals.refresh.emit(self.job_id)

	def update_finished(self):
		finished = int(time.time())
		sql = "UPDATE jobs SET progress=?,finished=? WHERE id=?"
		DB.query(sql, (100, finished, self.job_id))
		self.signals.refresh.emit(self.job_id)

	def update_error(self, error):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (3, error, self.job_id))
		self.signals.refresh.emit(self.job_id)

	def update_success(self):
		sql = "UPDATE jobs SET status=?,message=? WHERE id=?"
		DB.query(sql, (1, "The job was successfully finished", self.job_id))
		self.signals.refresh.emit(self.job_id)


	@Slot()
	def run(self):
		self.update_started()

		try:
			self.pipline(*self.args, **self.kwargs)
		except:
			self.update_error(traceback.format_exc())
		else:
			self.update_success()
		finally:
			pass

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
		py27, folder = self.get_mgltools()
		script = os.path.join(folder, 'prepare_receptor4.py')
		args = [script, '-r', infile, '-o', outfile, '-A', 'bonds_hydrogens']
		self.execute(py27, args, os.path.dirname(infile))

	def prepare_ligand(self, infile, outfile):
		py27, folder = self.get_mgltools()
		script = os.path.join(folder, 'prepare_ligand4.py')
		args = [script, '-l', infile, '-o', outfile, '-A', 'bonds_hydrogens']
		self.execute(py27, args, os.path.dirname(infile))

class AutodockWorker(BaseWorker):
	def __init__(self, *args, **kwargs):
		super(AutodockWorker, self).__init__(*args, **kwargs)

	def get_commands(self):
		autodock = self.settings.value('Tools/autodock_4')
		autogrid = self.settings.value('Tools/autogrid_4')

		return autodock, autogrid

	def pipline(self, job, params):
		work_dir = os.path.join(job.root, 'jobs', "job.{}.{}.vs.{}".format(
			job.id, job.rn, job.ln
		))

		#create receptor and ligand dock dir
		if not os.path.exists(work_dir):
			os.mkdir(work_dir)

		#convert receptor and ligand to pdbqt format
		rpdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.rn))
		self.prepare_receptor(job.rp, rpdbqt)

		lpdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.ln))
		self.prepare_ligand(job.lp, lpdbqt)

		#set autogrid parameter and create gpf parameter file
		ag_param = AutogridParameter(rpdbqt, lpdbqt, (job.x, job.y, job.z),
			(job.cx, job.cy, job.cz), job.spacing
		)
		gpf_file = ag_param.make_gpf_file()
		glg_file = gpf_file.replace('.gpf', '.glg')
		self.log_file = glg_file

		autodock, autogrid = self.get_commands()

		#run autogrid4
		#delete the glg log file before new start
		if QFile.exists(glg_file):
			QFile.remove(glg_file)

		self.update_message("Running autogrid")
		proc = psutil.Popen([autogrid, '-p', gpf_file, '-l', glg_file],
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = work_dir
		)

		read_start = 0

		while proc.poll() is None:
			QThread.sleep(1)

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


		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		#set autodock4 parameters and make dpf parameter file
		dpf_file = params.make_dpf_file(rpdbqt, lpdbqt)
		dlg_file = dpf_file.replace('.dpf', '.dlg')

		#run autodock4
		self.update_message("Running autodock")
		proc = psutil.Popen([autodock, '-p', dpf_file, '-l', dlg_file],
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = work_dir
		)

		read_start = 0

		while proc.poll() is None:
			QThread.sleep(1)

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

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

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
					rank = int(cols[0])
					energy = float(cols[3])
					ki = runs[rid][0]
					crmsd = float(cols[4])
					rrmsd = float(cols[5])
					ligand = convert_pdbqt_to_pdb(''.join(runs[rid][1]))

					rows.append((None, job.id, rid, energy, ki, crmsd, rrmsd, rank, ligand))

		sql = "INSERT INTO pose VALUES (?,?,?,?,?,?,?,?,?)"
		DB.insert_rows(sql, rows)

class AutodockVinaWorker(BaseWorker):
	def __init__(self, *args, **kwargs):
		super(AutodockVinaWorker, self).__init__(*args, **kwargs)

	def get_commands(self):
		vina = self.settings.value('Tools/autodock_vina')
		autogrid = self.settings.value('Tools/autogrid_4')

		return vina, autogrid

	def pipline(self, job, params):
		work_dir = os.path.join(job.root, 'jobs', "job.{}.{}.vs.{}".format(
			job.id, job.rn, job.ln
		))

		#create receptor and ligand dock dir
		if not os.path.exists(work_dir):
			os.mkdir(work_dir)

		#convert receptor and ligand to pdbqt format
		rpdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.rn))
		self.prepare_receptor(job.rp, rpdbqt)

		lpdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.ln))
		self.prepare_ligand(job.lp, lpdbqt)

		#get commands
		vina, autogrid = self.get_commands()

		#if use autodock scoring function
		#set autogrid parameter and create gpf parameter file
		if params.scoring == 'ad4':
			ag_param = AutogridParameter(rpdbqt, lpdbqt, (job.x, job.y, job.z),
				(job.cx, job.cy, job.cz), job.spacing
			)
			gpf_file = ag_param.make_gpf_file()
			glg_file = gpf_file.replace('.gpf', '.glg')
			self.log_file = glg_file

			#run autogrid4
			#delete the glg log file before new start
			if QFile.exists(glg_file):
				QFile.remove(glg_file)

			self.update_message("Running autogrid")
			proc = psutil.Popen([autogrid, '-p', gpf_file, '-l', glg_file],
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				cwd = work_dir
			)

			read_start = 0

			while proc.poll() is None:
				QThread.sleep(1)

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

			if proc.returncode != 0:
				raise Exception(proc.stderr.read())

			params.maps.value = job.rn

		#set autodock vina parameters and make config file
		config_file = os.path.join(work_dir, "config.txt".format(job.rn))
		out_file = "out.pdbqt"
		log_file = "out.log"

		params.receptor.value = os.path.basename(rpdbqt)
		params.ligand.value = os.path.basename(lpdbqt)
		params.center_x.value = job.cx
		params.center_y.value = job.cy
		params.center_z.value = job.cz
		params.size_x.value = job.x
		params.size_y.value = job.y
		params.size_z.value = job.z
		params.spacing.value = job.spacing
		params.out.value = out_file
		params.make_config_file(config_file)

		#run autodock vina
		self.update_message("Running autodock vina")
		proc = psutil.Popen([vina, '--config', config_file],
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE,
			cwd = work_dir,
			encoding = 'utf8'
		)

		star = 0
		log_fw = open(log_file, 'w')

		while proc.poll() is None:
			char = proc.stdout.read(1)

			if char == '*':
				star += 1
				p = round(star/51*100, 2)
				self.update_progress(p)

			if char:
				log_fw.write(char)
			else:
				QThread.msleep(500)

		log_fw.write(proc.stdout.read())

		if proc.returncode != 0:
			raise Exception(proc.stderr.read())

		rows = []
		with open(log_file) as fh:
			for line in fh:
				if line.startswith('-----+'):
					break

			for line in fh:
				cols = line.strip().split()

				rows.append((None, job.id, int(cols[0]), float(cols[1]),
				 float(cols[2]), float(cols[3])))

		with open(out_file) as fh:
			for line in fh:
				if line.startswith('MODEL'):
					pass
				elif line.startswith('ENDMDL'):
					pass

					



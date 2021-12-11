import os
import time
import traceback
import subprocess

from PySide6.QtSql import *
from PySide6.QtCore import *

from param import *
from utils import *
from backend import *

__all__ = ['AutodockWorker']

class WorkerSignals(QObject):
	finished = Signal()
	error = Signal(str)
	progress = Signal(int)

class BaseWorker(QRunnable):
	def __init__(self, *args, **kwargs):
		super(BaseWorker, self).__init__()
		self.args = args
		self.kwargs = kwargs
		self.signals = WorkerSignals()
		self.settings = QSettings()

	@Slot()
	def run(self):
		try:
			self.pipline(*self.args, **self.kwargs)
		except:
			self.signals.error.emit(traceback.format_exc())
		else:
			pass #successfully
		finally:
			self.signals.finished.emit()

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
		self.log_size = 0
		self.log_start = 0
		self.log_file = None
		self.job_id = 0
		self.autogrid_stopped = False

	@Slot()
	def get_progress(self):
		if self.log_file is None:
			return

		if not QFile.exists(self.log_file):
			return

		s = QFileInfo(self.log_file).size()

		if s > self.log_size:
			self.log_size = s
		else:
			return

		with open(self.log_file) as fh:
			fh.seek(self.log_start)
			p = 0
			for line in fh:
				if '%' not in line:
					continue

				p = float(line.strip().split()[2].replace('%', ''))
					
			if p > 0:
				DB.query("UPDATE jobs SET progress=? WHERE id=?", (p, self.job_id))
				self.signals.progress.emit(self.job_id)
				print(p)

			self.log_start = fh.tell()

	def get_commands(self):
		autodock = self.settings.value('Tools/autodock_4')
		autogrid = self.settings.value('Tools/autogrid_4')

		return autodock, autogrid

	def pipline(self, job, params):
		self.job_id = job.id

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

		#proc = subprocess.Popen([autogrid, '-p', gpf_file, '-l', glg_file], cwd=work_dir)
		#read_start = 0
		parent = QObject()
		proc = QProcess(parent)
		proc.setWorkingDirectory(work_dir)
		proc.start(autogrid, ['-p', gpf_file, '-l', glg_file])
		proc.waitForStarted(-1)

		read_start = 0
		p = 0

		while 1:
			if QFileInfo(glg_file).size() > read_start:
				p = 0

				with open(glg_file) as log:
					log.seek(read_start)

					for line in log:
						if '%' in line:
							p = float(line.strip().split()[2].replace('%', ''))

					if p > 0:
						DB.query("UPDATE jobs SET progress=? WHERE id=?", (p, job.id))
						self.signals.progress.emit(job.id)

					read_start = log.tell()

			if p == 100:
				break

			QThread.sleep(1)

		proc.waitForFinished(-1)

		#set autodock4 parameters and make dpf parameter file
		#dpf_file = params.make_dpf_file(rpdbqt, lpdbqt)
		#dlg_file = dpf_file.replace('.dpf', '.dlg')

		#run autodock4
		#proc = QProcess()
		#proc.setWorkingDirectory(work_dir)
		#proc.start(autodock, ['-p', dpf_file, '-l', dlg_file])
		#proc.waitForFinished(-1)


class AutodockVinaWorker(BaseWorker):
	pass




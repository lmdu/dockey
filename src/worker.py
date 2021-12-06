import os
import traceback

from PySide6.QtSql import *
from PySide6.QtCore import *

from param import *
from utils import *
from backend import *

__all__ = ['AutodockWorker']

class WorkerSignals(QObject):
	finished = Signal()
	error = Signal(str)
	progress = Signal()



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

	def prepare_receptor(self, infile, outfile):
		py27, folder = self.get_mgltools()
		script = os.path.join(folder, 'prepare_receptor4.py')

		args = [script, '-r', infile, '-o', outfile, '-A', 'bonds_hydrogens']

		proc = QProcess()
		proc.setWorkingDirectory(os.path.dirname(infile))
		proc.start(py27, args)
		proc.waitForFinished(-1)

	def prepare_ligand(self, infile, outfile):
		py27, folder = self.get_mgltools()
		script = os.path.join(folder, 'prepare_ligand4.py')

		args = [script, '-l', infile, '-o', outfile, '-A', 'bonds_hydrogens']

		proc = QProcess()
		proc.setWorkingDirectory(os.path.dirname(infile))
		proc.start(py27, args)
		proc.waitForFinished(-1)

class AutodockWorker(BaseWorker):
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

		autodock, autogrid = self.get_commands()

		#run autogrid4
		proc = QProcess()
		proc.setWorkingDirectory(work_dir)
		proc.start(autogrid, ['-p', gpf_file, '-l', glg_file])
		proc.waitForFinished(-1)

		#set autodock4 parameters and make gdf parameter file

		#run autodock4
		#proc = QProcess()
		#proc.setWorkingDirectory(work_dir)
		#proc.start(autodock, ['-p', gdf_file, '-l', dlg_file])
		#proc.waitForFinished(-1)

class AutodockVinaWorker(BaseWorker):
	pass




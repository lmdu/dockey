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

class AutodockWorker(BaseWorker):
	def pipline(self, job):
		work_dir = os.path.join(job.root, 'jobs', "job{}{}vs{}".format(
			job.id, job.rn, job.ln
		))

		#create receptor and ligand dock dir
		os.mkdir(work_dir)

		#convert receptor and ligand to pdbqt format
		rpdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.rn))
		convert_other_to_pdbqt(job.rp, job.rf, rpdbqt)

		lpdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.ln))
		convert_other_to_pdbqt(job.lp, job.lf, lpdbqt)

		#set autogrid parameter and create gpf parameter file
		ag_param = AutogridParameter(rpdbqt, lpdbqt, (job.x, job.y, job.z),
			(job.cx, job.cy, job.cz), job.spacing
		)
		gpf_file = ag_param.make_gpf_file()
		glg_file = gpf_file.replace('.gpf', '.glg')

		#run autogrid4
		proc = QProcess()
		proc.start('autogrid4.exe', ['-p', gpf_file, '-l', glg_file])
		proc.waitForFinished(-1)

		print(proc.exitStatus())
		print(proc.readAllStandardError())

		print('finished yes')

		#set autogrid4 parameters and make gdf parameter file

		#run autogrid4

class AutodockVinaWorker(BaseWorker):
	pass




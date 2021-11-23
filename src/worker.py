import os
import traceback

from PySide6.QtSql import *
from PySide6.QtCore import *

from db import *
from param import *
from utils import *

__all__ = ['AutodockWorker']

class BaseWorker(QObject):
	finished = Signal()

	def __init__(self, parent):
		super(BaseWorker, self).__init__()
		self.parent = parent
		self.job_dir = os.path.join(self.parent.project_folder, 'jobs')

		#job query
		self.jq = QSqlQuery("SELECT * FROM jobs")

	@property
	def jobs(self):
		while self.jq.next():
			job = AttrDict({
				'id': self.jq.value('id'),
				'rid': self.jq.value('rid'),
				'lid': self.jq.value('lid')
			})

			job['receptor'] = DB.get_molecular_by_id(job.rid)
			job['ligand'] = DB.get_molecular_by_id(job.lid)
			job['grid'] = DB.get_grid_box(job.rid)

			yield job

	def __error(self):
		err_msg = traceback.format_exc()
		print(err_msg)

	def run(self):
		try:
			self.pipline()
		except:
			self.__error()

		self.finished.emit()


class AutodockWorker(BaseWorker):
	def pipline(self):
		for job in self.jobs:
			work_dir = os.path.join(self.job_dir, "j{}{}vs{}".format(
				job.id, job.receptor.name, job.ligand.name
			))

			#create receptor and ligand dock dir
			os.mkdir(work_dir)

			#convert receptor and ligand to pdbqt format
			receptor_pdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.receptor.name))
			convert_other_to_pdbqt(job.receptor.content, job.receptor.format, receptor_pdbqt)

			ligand_pdbqt = os.path.join(work_dir, "{}.pdbqt".format(job.ligand.name))
			convert_other_to_pdbqt(job.ligand.content, job.ligand.format, ligand_pdbqt)

			#set autogrid parameter and create gpf parameter file
			ag_param = AutogridParameter(receptor_pdbqt, ligand_pdbqt,
				(job.grid.x, job.grid.y, job.grid.z),
				(job.grid.cx, job.grid.cy, job.grid.cz),
				job.grid.spacing
			)
			gpf_file = ag_param.make_gpf_file()

			glg_file = gpf_file.replace('.gpf', '.glg')

			#run autogrid4
			proc = QProcess(self.parent)
			proc.start('autogrid4', ['-p', gpf_file, '-l', glg_file])
			proc.waitForFinished(-1)

			#set autogrid4 parameters and make gdf parameter file

			#run autogrid4

			





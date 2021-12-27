import os

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *

__all__ = ['BrowseInput', 'CreateProjectDialog', 'AutodockConfigDialog',
			'AutodockVinaConfigDialog', 'ExportImageDialog']

class BrowseInput(QWidget):
	def __init__(self, parent=None, is_file=True, is_save=False, _filter=""):
		super(BrowseInput, self).__init__(parent)
		self.filter = _filter
		self.input = QLineEdit(self)
		self.input.setReadOnly(True)
		self.browse = QPushButton(self)
		self.browse.setFlat(True)
		self.browse.setIcon(QIcon("icons/folder.svg"))

		if is_file:
			if is_save:
				self.browse.clicked.connect(self.select_save)
			else:
				self.browse.clicked.connect(self.select_file)
		else:
			self.browse.clicked.connect(self.select_folder)

		layout = QHBoxLayout()
		layout.setSpacing(0)
		layout.setContentsMargins(0,0,0,0)
		layout.addWidget(self.input)
		layout.addWidget(self.browse)

		self.setLayout(layout)

	@Slot()
	def select_save(self):
		exe, _ = QFileDialog.getSaveFileName(self, filter=self.filter)

		if exe:
			self.input.setText(exe)

	@Slot()
	def select_file(self):
		exe, _ = QFileDialog.getOpenFileName(self, filter=self.filter)

		if exe:
			self.input.setText(exe)

	@Slot()
	def select_folder(self):
		folder = QFileDialog.getExistingDirectory(self)

		if folder:
			self.input.setText(folder)

	def get_text(self):
		return self.input.text()

	def set_text(self, text):
		self.input.setText(text)

class CreateProjectDialog(QDialog):
	def __init__(self, parent=None):
		super(CreateProjectDialog, self).__init__(parent)
		self.setWindowTitle("Create a New Project")

		self.name_input = QLineEdit(self)
		self.location_input = BrowseInput(self, False)
		tips = (
			'<font color="gray">'
			'This will create a new folder with input name in the '
			'location you selected'
			'</font>'
		)

		self.tip_label = QLabel(tips, self)

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.check_input_name)
		action_box.rejected.connect(self.reject)

		layout = QFormLayout()
		layout.addRow("Name", self.name_input)
		layout.addRow("Location", self.location_input)
		layout.addRow(self.tip_label)
		layout.addRow(action_box)

		self.setLayout(layout)

	@Slot()
	def check_input_name(self):
		name = self.name_input.text()
		folder = self.location_input.get_text()
		if not name:
			return QMessageBox.warning(self, "Input Warning",
				"Please input a project name"
			)

		if not folder:
			return QMessageBox.warning(self, "Input Warning",
				"Please select a directory where to create project"
			)

		if ' ' in name:
			return QMessageBox.critical(self, "Input Error",
				"White-spaces are not allowed in project name"
			)
		
		if any(c in name for c in '\\/:*?"<>|'):
			return QMessageBox.critical(self, "Input Error",
				"Your input name has illegal character"
			)

		if os.path.exists(os.path.join(folder, name)):
			return QMessageBox.critical(self, "Input Error",
				"Your input project name <b>{}</b> exists in folder <b>{}</b>. Please input a new project name".format(
					name, folder
				)
			)

		self.accept()

	@classmethod
	def get_project_folder(cls, parent):
		dlg = cls(parent)
		if dlg.exec():
			name = dlg.name_input.text()
			location = dlg.location_input.get_text()
			return os.path.join(location, name)

class ProgramConfigDialog(QDialog):
	title = ""
	progs = []

	def __init__(self, parent=None):
		super(ProgramConfigDialog, self).__init__(parent)
		self.setWindowTitle(self.title)
		self.resize(QSize(500, 10))
		self.settings = QSettings()

		layout = QVBoxLayout()
		self.setLayout(layout)

		self.inputs = []

		for i, prog in enumerate(self.progs):
			k, l = prog
			v = self.settings.value(k, '')

			if v:
				self.inputs.append(None)
			else:
				self.inputs.append(BrowseInput(self))
				layout.addWidget(QLabel(l, self))
				layout.addWidget(self.inputs[i])

		action_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.save_config)
		action_box.rejected.connect(self.reject)

		layout.addWidget(action_box)

	@Slot()
	def save_config(self):
		for put in self.inputs:
			if put is None:
				continue

			v = put.get_text()

			if not v:
				return

		for i, put in enumerate(self.inputs):
			if put is None:
				continue

			v = put.get_text()
			k = self.progs[i][0]

			self.settings.setValue(k, v)

		self.accept()

	@classmethod
	def check_executable(cls, parent):
		settings = QSettings()
		configs = [settings.value(k) for k, l in cls.progs]

		if all(configs):
			return True

		dlg = cls(parent)
		return dlg.exec()

class AutodockConfigDialog(ProgramConfigDialog):
	title = "Autodock Config"
	progs = [
		('Tools/autodock_4', 'Autodock executable'),
		('Tools/autogrid_4', 'Autogrid executable'),
		('Tools/mgltools', 'MGLTools directory')
	]

class AutodockVinaConfigDialog(ProgramConfigDialog):
	title = "Autodock Vina Config"
	progs = [
		('Tools/autodock_vina', 'Autodock Vina executable'),
		('Tools/mgltools', 'MGLTools directory')
	]

class ExportImageDialog(QDialog):
	def __init__(self, parent=None):
		super(ExportImageDialog, self).__init__(parent)
		self.options = None

		self.setWindowTitle("Export as PNG Image")
		self.resize(QSize(450, 10))
		layout = QFormLayout()
		layout.setVerticalSpacing(10)
		self.setLayout(layout)

		self.file_input = BrowseInput(self, is_save=True, _filter="Image (*.png)")
		self.width_input = QSpinBox(self)
		self.width_input.setRange(0, 100000)
		self.width_input.setValue(1000)
		self.height_input = QSpinBox(self)
		self.height_input.setRange(0, 100000)
		self.height_input.setValue(1000)
		self.dpi_input = QComboBox(self)
		self.dpi_input.addItems(['72', '96', '150', '200', '300', '600'])
		self.dpi_input.setCurrentIndex(1)

		layout.addRow("Image file: ", self.file_input)

		sub_layout = QHBoxLayout()
		sub_layout.addWidget(QLabel("Width: ", self))
		sub_layout.addWidget(self.width_input, 1)
		sub_layout.addSpacing(20)
		sub_layout.addWidget(QLabel("Height: ", self))
		sub_layout.addWidget(self.height_input, 1)
		sub_layout.addSpacing(20)
		sub_layout.addWidget(QLabel("DPI: ", self))
		sub_layout.addWidget(self.dpi_input, 1)

		layout.addRow(sub_layout)

		action_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
		action_box.accepted.connect(self.save_image)
		action_box.rejected.connect(self.reject)

		layout.addRow(action_box)

	@Slot()
	def save_image(self):
		image = self.file_input.get_text()

		if not image:
			return QMessageBox.warning(self, "Input Warning",
				"Please select a save file"
			)

		self.options = AttrDict({
			'image': self.file_input.get_text(),
			'width': self.width_input.value(),
			'height': self.height_input.value(),
			'dpi': int(self.dpi_input.currentText()),
		})

		self.accept()

	@classmethod
	def save_to_png(cls, parent):
		dlg = cls(parent)
		if dlg.exec():
			return dlg.options

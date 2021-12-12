from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from config import *
from backend import *

__all__ = ['GridBoxSettingPanel']

class GridBoxParamter(QObject):
	updated = Signal()

	def __init__(self):
		super(GridBoxParamter, self).__init__()
		self.initialize()

	def __setitem__(self, k, v):
		setattr(self, k, v)

	def __getitem__(self, k):
		return self.property(k)

	def initialize(self):
		self.x = 0
		self.y = 0
		self.z = 0
		self.cx = 0
		self.cy = 0
		self.cz = 0
		self.show_face = True
		self.show_edge = True
		self.use_line = True
		self.use_cylinder = False
		self.edge_width = 1
		self.edge_color = (1.0, 1.0, 1.0)
		self.spacing = 0.375
		self.opacity = 50
		self.bg_x = (1.0, 0.0, 0.0)
		self.bg_y = (0.0, 1.0, 0.0)
		self.bg_z = (0.0, 0.0, 1.0)

	def update_dimension(self, points):
		x, y, z, cx, cy, cz = convert_coordinates_to_dimension(points, self.spacing)
		self.x = x
		self.y = y
		self.z = z
		self.cx = cx
		self.cy = cy
		self.cz = cz
		self.updated.emit()

	def update_grid(self, data):
		self.x = data.x
		self.y = data.y
		self.z = data.z
		self.cx = data.cx
		self.cy = data.cy
		self.cz = data.cz
		self.spacing = data.spacing
		self.updated.emit()

	def reset(self):
		self.initialize()
		self.updated.emit()

	def custom(self):
		self.initialize()
		self.x = 40
		self.y = 40
		self.z = 40
		self.updated.emit()

#https://www.pythonguis.com/widgets/qcolorbutton-a-color-selector-tool-for-pyqt/
class ColorButton(QPushButton):
	colorChanged = Signal(tuple)

	def __init__(self, *args, color=None, **kwargs):
		super(ColorButton, self).__init__(*args, **kwargs)

		self._color = None
		self._default = color
		self.pressed.connect(self.onColorPicker)
		self.setColor(self._default)

		self.setFixedSize(QSize(20, 15))

	def setColor(self, color):
		if color != self._color:
			self._color = color
			self.colorChanged.emit(QColor(color).getRgbF()[0:3])

		if self._color:
			self.setStyleSheet("background-color: {}".format(self._color))
		else:
			self.setStyleSheet("")

	def color(self):
		return self._color

	def onColorPicker(self):
		dlg = QColorDialog()

		if self._color:
			dlg.setCurrentColor(QColor(self._color))

		if dlg.exec():
			self.setColor(dlg.currentColor().name())

	def mousePressEvent(self, e):
		if e.button() == Qt.RightButton:
			self.setColor(self._default)

		return super(ColorButton, self).mousePressEvent(e)

class GridBoxSettingPanel(QWidget):
	updated = Signal()

	def __init__(self, parent=None):
		super(GridBoxSettingPanel, self).__init__(parent)
		self.parent = parent

		self.params = GridBoxParamter()

		self.layout = QFormLayout()
		self.setLayout(self.layout)

		self.create_controler()

		self.updated.connect(self.redraw_box)
		self.params.updated.connect(self.set_value)

	def create_widgets(self):
		#self.face_widget = QCheckBox("Show box face", self)
		#self.face_widget.setTristate(False)
		#self.face_widget.setCheckState(Qt.Checked)

		#self.edge_widget = QCheckBox("Show box edge", self)
		#self.edge_widget.setTristate(False)
		#self.edge_widget.setCheckState(Qt.Checked)
		pass

	def redraw_box(self):
		draw_gridbox(self.parent.cmd, self.params)

	def update_value(self, k, v):
		if self.params[k] == v:
			return

		self.params[k] = v
		self.updated.emit()

	@Slot()
	def set_value(self):
		self.size_x.setValue(self.params.x)
		self.size_y.setValue(self.params.y)
		self.size_z.setValue(self.params.z)
		self.pos_x.setValue(self.params.cx)
		self.pos_y.setValue(self.params.cy)
		self.pos_z.setValue(self.params.cz)
		self.wide.setValue(self.params.edge_width)

	def create_controler(self):
		#show box face
		self.face = QCheckBox("Show box face", self)
		self.face.setCheckState(Qt.Checked)
		self.face.setTristate(False)
		self.face.stateChanged.connect(lambda x: self.update_value('show_face', x))
		self.layout.addRow(self.face)

		#face background
		self.bg_x = ColorButton(self, color=QColor.fromRgbF(*self.params.bg_x).name())
		self.bg_x.colorChanged.connect(lambda x: self.update_value('bg_x', x))
		self.bg_y = ColorButton(self, color=QColor.fromRgbF(*self.params.bg_y).name())
		self.bg_y.colorChanged.connect(lambda x: self.update_value('bg_y', x))
		self.bg_z = ColorButton(self, color=QColor.fromRgbF(*self.params.bg_z).name())
		self.bg_z.colorChanged.connect(lambda x: self.update_value('bg_z', x))
		bg_layout = QHBoxLayout()
		bg_layout.addWidget(QLabel("Background", self))
		bg_layout.addWidget(self.bg_x)
		bg_layout.addWidget(self.bg_y)
		bg_layout.addWidget(self.bg_z)
		self.layout.addRow(bg_layout)

		self.opacity = QSlider(Qt.Horizontal, self)
		self.opacity.setRange(0, 99)
		self.opacity.setValue(self.params.opacity)
		self.opacity.valueChanged.connect(lambda x: self.update_value('opacity', x))
		self.layout.addRow(QLabel("Background opacity"))
		self.layout.addRow(self.opacity)

		#show box edge
		self.edge = QCheckBox("Show box line", self)
		self.edge.setCheckState(Qt.Checked)
		self.edge.setTristate(False)
		self.edge.stateChanged.connect(lambda x: self.update_value('show_edge', x))
		self.layout.addRow(self.edge)

		self.line = QRadioButton("Line", self)
		self.line.setChecked(True)
		self.line.toggled.connect(lambda x: self.update_value('use_line', x))
		self.cylinder = QRadioButton("Cylinder", self)
		self.cylinder.toggled.connect(lambda x: self.update_value('use_cylinder', x))
		type_group = QButtonGroup(self)
		type_group.setExclusive(True)
		type_group.addButton(self.line)
		type_group.addButton(self.cylinder)
		type_layout = QHBoxLayout()
		type_layout.addWidget(self.line)
		type_layout.addWidget(self.cylinder)
		self.layout.addRow(QLabel("Line type", self))
		self.layout.addRow(type_layout)

		self.wide = QSpinBox(self)
		self.wide.setRange(1, 100)
		self.wide.valueChanged.connect(lambda x: self.update_value('edge_width', x))
		self.layout.addRow("Line width", self.wide)
		self.color = ColorButton(self, color=QColor.fromRgbF(*self.params.edge_color).name())
		self.color.colorChanged.connect(lambda x: self.update_value('edge_color', x))
		color_layout = QHBoxLayout()
		color_layout.addWidget(QLabel("Line color", self))
		color_layout.addWidget(self.color)
		self.layout.addRow(color_layout)

		#spacing controler
		self.spacing = QDoubleSpinBox(self)
		self.spacing.setSingleStep(0.005)
		self.spacing.setRange(0, 1)
		self.spacing.setDecimals(3)
		self.spacing.setValue(self.params.spacing)
		self.spacing.valueChanged.connect(lambda x: self.update_value('spacing', x))
		#self.layout.addRow(QLabel('Padding space', self))
		self.layout.addRow("Spacing", self.spacing)

		self.size_x = QSpinBox(self)
		self.size_x.setSingleStep(2)
		self.size_x.setRange(0, 100000)
		self.size_x.valueChanged.connect(lambda x: self.update_value('x', x))
		self.size_y = QSpinBox(self)
		self.size_y.setSingleStep(2)
		self.size_y.setRange(0, 100000)
		self.size_y.valueChanged.connect(lambda x: self.update_value('y', x))
		self.size_z = QSpinBox(self)
		self.size_z.setSingleStep(2)
		self.size_z.setRange(0, 100000)
		self.size_z.valueChanged.connect(lambda x: self.update_value('z', x))

		self.layout.addRow(QLabel("Points in each dimension", self))
		size_layout = QGridLayout()
		size_layout.setColumnStretch(1, 1)
		size_layout.addWidget(QLabel("x", self), 0, 0)
		size_layout.addWidget(self.size_x, 0, 1)
		size_layout.addWidget(QLabel("y", self), 1, 0)
		size_layout.addWidget(self.size_y, 1, 1)
		size_layout.addWidget(QLabel("z", self), 2, 0)
		size_layout.addWidget(self.size_z, 2, 1)
		self.layout.addRow(size_layout)

		self.pos_x = QDoubleSpinBox(self)
		self.pos_x.setRange(-100000, 100000)
		self.pos_x.setDecimals(3)
		self.pos_x.setSingleStep(0.252)
		self.pos_x.valueChanged.connect(lambda x: self.update_value('cx', x))
		self.pos_y = QDoubleSpinBox(self)
		self.pos_y.setRange(-100000, 100000)
		self.pos_y.setDecimals(3)
		self.pos_y.setSingleStep(0.252)
		self.pos_y.valueChanged.connect(lambda x: self.update_value('cy', x))
		self.pos_z = QDoubleSpinBox(self)
		self.pos_z.setRange(-100000, 100000)
		self.pos_z.setDecimals(3)
		self.pos_z.setSingleStep(0.252)
		self.pos_z.valueChanged.connect(lambda x: self.update_value('cz', x))

		self.layout.addRow(QLabel("Centre coordinates", self))
		pos_layout = QGridLayout()
		pos_layout.setColumnStretch(1, 1)
		pos_layout.addWidget(QLabel("x", self), 0, 0)
		pos_layout.addWidget(self.pos_x, 0, 1)
		pos_layout.addWidget(QLabel("y", self), 1, 0)
		pos_layout.addWidget(self.pos_y, 1, 1)
		pos_layout.addWidget(QLabel("z", self), 2, 0)
		pos_layout.addWidget(self.pos_z, 2, 1)
		self.layout.addRow(pos_layout)

		self.save_btn = QPushButton("Save grid box", self)
		self.save_btn.clicked.connect(self.on_save)
		self.layout.addRow(self.save_btn)

	@Slot()
	def on_save(self):
		if self.parent.current_molecular:
			if self.parent.current_molecular.type == 1:
				sql = "INSERT INTO grid VALUES (?,?,?,?,?,?,?,?,?)"
				DB.query(sql, (
					None,
					self.parent.current_molecular.id,
					self.params.x,
					self.params.y,
					self.params.z,
					self.params.cx,
					self.params.cy,
					self.params.cz,
					self.params.spacing
				))

				QMessageBox.information(self.parent, "Grid box saved",
					"Successfully set grid box for receptor {}".format(
						self.parent.current_molecular.name)
				)
			else:
				QMessageBox.warning(self.parent, "Warning",
					"Could not save grid box for ligand"
				)
		else:
			QMessageBox.warning(self.parent, "Warning",
				"Please select a receptor to add grid box"
			)


from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from config import *

__all__ = ['GridBoxSettingPanel', 'GBD']

class GridBoxStorage(QObject):
	updated = Signal()

	def __init__(self):
		super(GridBoxStorage, self).__init__()
		self.initialize()

	def __setitem__(self, k, v):
		setattr(self, k, v)

	def __getitem__(self, k):
		return self.property(k)

	def initialize(self):
		self.x = 10
		self.y = 10
		self.z = 10
		self.cx = 0
		self.cy = 0
		self.cz = 0
		self.show_face = True
		self.show_edge = True
		self.use_line = True
		self.use_cylinder = False
		self.edge_width = 1
		self.edge_color = (1.0, 1.0, 1.0)
		self.padding = 0
		self.opacity = 50
		self.bg_x = (1.0, 0.0, 0.0)
		self.bg_y = (0.0, 1.0, 0.0)
		self.bg_z = (0.0, 0.0, 1.0)

	def update_dimension(self, points):
		x, y, z, cx, cy, cz = convert_coordinates_to_dimension(points)
		self.x = x
		self.y = y
		self.z = z
		self.cx = cx
		self.cy = cy
		self.cz = cz
		self.updated.emit()

	def reset(self):
		self.initialize()
		self.updated.emit()

GBD = GridBoxStorage()

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

		self.layout = QFormLayout()
		self.setLayout(self.layout)

		self.create_controler()

		self.updated.connect(self.redraw_box)

	def create_widgets(self):
		#self.face_widget = QCheckBox("Show box face", self)
		#self.face_widget.setTristate(False)
		#self.face_widget.setCheckState(Qt.Checked)

		#self.edge_widget = QCheckBox("Show box edge", self)
		#self.edge_widget.setTristate(False)
		#self.edge_widget.setCheckState(Qt.Checked)
		pass

	def redraw_box(self):
		draw_gridbox(self.parent.cmd, GBD)

	def update_value(self, k, v):
		if GBD[k] == v:
			return

		GBD[k] = v
		self.updated.emit()

	@Slot()
	def set_value(self):
		self.size_x.setValue(GBD.x)
		self.size_y.setValue(GBD.y)
		self.size_z.setValue(GBD.z)
		self.pos_x.setValue(GBD.cx)
		self.pos_y.setValue(GBD.cy)
		self.pos_z.setValue(GBD.cz)
		self.wide.setValue(GBD.edge_width)

	def create_controler(self):
		#show box face
		self.face = QCheckBox("Show box face", self)
		self.face.setCheckState(Qt.Checked)
		self.face.setTristate(False)
		self.face.stateChanged.connect(lambda x: self.update_value('show_face', x))
		self.layout.addRow(self.face)

		#face background
		self.bg_x = ColorButton(self, color=QColor.fromRgbF(*GBD.bg_x).name())
		self.bg_x.colorChanged.connect(lambda x: self.update_value('bg_x', x))
		self.bg_y = ColorButton(self, color=QColor.fromRgbF(*GBD.bg_y).name())
		self.bg_y.colorChanged.connect(lambda x: self.update_value('bg_y', x))
		self.bg_z = ColorButton(self, color=QColor.fromRgbF(*GBD.bg_z).name())
		self.bg_z.colorChanged.connect(lambda x: self.update_value('bg_z', x))
		bg_layout = QHBoxLayout()
		bg_layout.addWidget(QLabel("Background", self))
		bg_layout.addWidget(self.bg_x)
		bg_layout.addWidget(self.bg_y)
		bg_layout.addWidget(self.bg_z)
		self.layout.addRow(bg_layout)

		self.opacity = QSlider(Qt.Horizontal, self)
		self.opacity.setRange(0, 99)
		self.opacity.setValue(GBD.opacity)
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
		self.wide.valueChanged.connect(lambda x: self.update_value('edge_width', x))
		self.layout.addRow("Line width", self.wide)
		self.color = ColorButton(self, color=QColor.fromRgbF(*GBD.edge_color).name())
		self.color.colorChanged.connect(lambda x: self.update_value('edge_color', x))
		color_layout = QHBoxLayout()
		color_layout.addWidget(QLabel("Line color", self))
		color_layout.addWidget(self.color)
		self.layout.addRow(color_layout)

		#padding controler
		self.padding = QSlider(Qt.Horizontal, self)
		self.padding.valueChanged.connect(lambda x: self.update_value('padding', x))
		self.layout.addRow(QLabel('Padding space', self))
		self.layout.addRow(self.padding)

		self.size_x = QDoubleSpinBox(self)
		self.size_x.valueChanged.connect(lambda x: self.update_value('x', x))
		self.size_y = QDoubleSpinBox(self)
		self.size_y.valueChanged.connect(lambda x: self.update_value('y', x))
		self.size_z = QDoubleSpinBox(self)
		self.size_z.valueChanged.connect(lambda x: self.update_value('z', x))

		
		self.layout.addRow(QLabel("Box size", self))
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
		self.pos_x.valueChanged.connect(lambda x: self.update_value('cx', x))
		self.pos_y = QDoubleSpinBox(self)
		self.pos_y.valueChanged.connect(lambda x: self.update_value('cy', x))
		self.pos_z = QDoubleSpinBox(self)
		self.pos_z.valueChanged.connect(lambda x: self.update_value('cz', x))

		self.layout.addRow(QLabel("Centre position", self))
		pos_layout = QGridLayout()
		pos_layout.setColumnStretch(1, 1)
		pos_layout.addWidget(QLabel("x", self), 0, 0)
		pos_layout.addWidget(self.pos_x, 0, 1)
		pos_layout.addWidget(QLabel("y", self), 1, 0)
		pos_layout.addWidget(self.pos_y, 1, 1)
		pos_layout.addWidget(QLabel("z", self), 2, 0)
		pos_layout.addWidget(self.pos_z, 2, 1)
		self.layout.addRow(pos_layout)


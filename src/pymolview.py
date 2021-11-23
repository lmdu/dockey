import sys

import pymol
from pymol2 import SingletonPyMOL as PyMOL

from PySide6.QtCore import *
from PySide6.QtGui import *
from PySide6.QtWidgets import *
from PySide6.QtOpenGLWidgets import *

#install plugins
import pymolplugins

from utils import *
from config import *
from gridbox import *

__all__ = ['PymolGLWidget']

#Gesture = QtCore.QEvent.Gesture

#Module for translating Qt key codes to PyMOL key and "special" codes
keyMap = {
	Qt.Key_Escape: 27,
	Qt.Key_Tab: 9,
	Qt.Key_Backspace: 8,
	Qt.Key_Return: 13,
	Qt.Key_Enter: 13,
	Qt.Key_Delete: 127,
}

specialMap = {
	Qt.Key_Left: 100,
	Qt.Key_Up: 101,
	Qt.Key_Right: 102,
	Qt.Key_Down: 103,
	Qt.Key_PageUp: 104,
	Qt.Key_PageDown: 105,
	Qt.Key_Home: 106,
	Qt.Key_End: 107,
	Qt.Key_Insert: 108,
	Qt.Key_F1: 1,
	Qt.Key_F2: 2,
	Qt.Key_F3: 3,
	Qt.Key_F4: 4,
	Qt.Key_F5: 5,
	Qt.Key_F6: 6,
	Qt.Key_F7: 7,
	Qt.Key_F8: 8,
	Qt.Key_F9: 9,
	Qt.Key_F10: 10,
	Qt.Key_F11: 11,
	Qt.Key_F12: 12,
}


def get_modifiers(ev):
	'''Get modifers from event and translate into PyMOL modifier mask'''
	pymolmod = 0
	qtmodifiers = ev.modifiers()

	for mask, qtm in [
		(0x1, Qt.ShiftModifier),
		(0x2, Qt.MetaModifier),  # CTRL on Mac
		(0x2, Qt.ControlModifier),
		(0x4, Qt.AltModifier)
	]:
		if qtmodifiers & qtm:
			pymolmod |= mask

	return pymolmod


def keyPressEventToPyMOLButtonArgs(ev):
	# translate modifier mask
	pymolmod = get_modifiers(ev)

	# Qt::Key_*
	key = ev.key()

	if key in specialMap:
		k = specialMap[key]

		# PyMOL_Special
		state = -2
	else:
		# PyMOL_Key
		state = -1

		# doesn't work on Mac: k = keyMap.get(key, ev.nativeVirtualKey())
		k = keyMap.get(key, -1)
		if k == -1:
			text = ev.text()
			if text:
				k = ord(text)

		# CTRL-<key>
		if k == -1 and (pymolmod & 0x2):
			k = key - 64

		# ALT-<key>
		if k != -1 and (pymolmod & 0x4):
			k = key

		if k > 255 or k < 0:
			#if DEBUG:
			print('DEBUG: skipped: 0x%x 0x%x' % (key, k))
			return

	return (k, state, 0, 0, pymolmod)

def get_wheel_delta(ev):
	'''
	Get mouse wheel delta from event.
	Ignores horizontal scrolling (returns zero).
	'''
	try:
		# Qt4
		return ev.delta()
	except AttributeError:
		pass

	# Qt5
	angledelta = ev.angleDelta()
	delta_x = angledelta.x()
	delta_y = angledelta.y()

	if abs(delta_y) < abs(delta_x):
		# Shift+Wheel emulates horizontal scrolling
		if not (ev.modifiers() & Qt.ShiftModifier):
			return 0
		return delta_x

	return delta_y


def get_wheel_button(ev):
	'''
	Get mouse wheel button index (3 or 4) from event, or 0 if no vertial
	scrolling was detected.
	'''
	delta = get_wheel_delta(ev)
	if delta > 0:
		return 3
	if delta < 0:
		return 4
	return 0

# QOpenGLWidget is supposed to supersede QGLWidget, but has issues (e.g.
# no stereo support)
#USE_QOPENGLWIDGET = int(
#	os.getenv("PYMOL_USE_QOPENGLWIDGET") or
#	(pymol.IS_MACOS and QtCore.QT_VERSION >= 0x50400))

#if USE_QOPENGLWIDGET:
#	BaseGLWidget = QtWidgets.QOpenGLWidget
#	AUTO_DETECT_STEREO = False
#else:
#	from pymol.Qt import QtOpenGL
#	BaseGLWidget = QtOpenGL.QGLWidget
#	# only attempt stereo detection in Qt <= 5.6 (with 5.9+ on Linux I
#	# get GL_DOUBLEBUFFER=0 with flickering when requesting stereo)
#	AUTO_DETECT_STEREO = pymol.IS_WINDOWS or QtCore.QT_VERSION < 0x50700
#USE_QOPENGLWIDGET = True
#BaseGLWidget = QtOpenGLWidgets.
#AUTO_DETECT_STEREO = False

class PymolGLWidget(QOpenGLWidget):
	'''
	PyMOL OpenGL Widget
	Can be used as a context manager to make OpenGL context current.
	'''

	# mouse button map
	_buttonMap = {
		Qt.LeftButton: 0,
		Qt.MiddleButton: 1,
		Qt.RightButton: 2,
	}

	def __enter__(self):
		'''
		Context manager to make the OpenGL context "current".
		Fixes depth-buffer issue with QOpenGLWidget
		https://github.com/schrodinger/pymol-open-source/issues/25
		'''
		self.makeCurrent()
		pymol._cmd._pushValidContext(self.cmd._COb)

	def __exit__(self, exc_type, exc_value, traceback):
		pymol._cmd._popValidContext(self.cmd._COb)

	def __init__(self, parent):
		self.gui = parent
		self.fb_scale = 1.0

		# OpenGL context setup
		f = QSurfaceFormat()

		from pymol.invocation import options

		# logic equivalent to layer5/main.cpp:launch
		if options.multisample:
			f.setSamples(4)

		# disable splash information
		options.quiet = 0
		options.show_splash = 0
		options.external_gui = 0
		options.internal_gui = 0
		options.internal_feedback = 0
		options.no_quit = 1	
		
		super(PymolGLWidget, self).__init__(parent=parent)
		self.setFormat(f)
		self.setUpdateBehavior(QOpenGLWidget.PartialUpdate)

		# pymol instance
		self.pymol = PyMOL()
		self.pymol.start()
		self.cmd = self.pymol.cmd
		#self.cmd.set('internal_gui', True)
		#self.cmd.set('internal_feedback', 0)

		# capture python output for feedback
		#import pcatch
		#pcatch._install()

		# for passive move drag
		self.setMouseTracking(True)

		# for accepting keyboard input (command line, shortcuts)
		self.setFocusPolicy(Qt.ClickFocus)

		# for idle rendering
		self._timer = QTimer()
		self._timer.setSingleShot(True)
		self._timer.timeout.connect(self._pymol_process)

		# drag n drop
		self.setAcceptDrops(True)

		# pinch-zoom
		self.grabGesture(Qt.PinchGesture)

		self.receptor = None

	def sizeHint(self):
		# default 640 + internal_gui, 480 + internal_feedback
		return QSize(860, 498)

	##########################
	# Input Events
	##########################

	def event(self, ev):
		if ev.type() == QEvent.Gesture:
			return self.gestureEvent(ev)

		return super(PymolGLWidget, self).event(ev)

	def gestureEvent(self, ev):
		gesture = ev.gesture(Qt.PinchGesture)

		if gesture is None:
			return False

		if gesture.state() == Qt.GestureStarted:
			self.pinch_start_z = self.cmd.get_view()[11]

		changeFlags = gesture.changeFlags()

		if changeFlags & QPinchGesture.RotationAngleChanged:
			delta = gesture.lastRotationAngle() - gesture.rotationAngle()
			self.cmd.turn('z', delta)

		if changeFlags & QPinchGesture.ScaleFactorChanged:
			view = list(self.cmd.get_view())

			# best guess for https://bugreports.qt.io/browse/QTBUG-48138
			totalscalefactor = gesture.totalScaleFactor()
			if totalscalefactor == 1.0:
				totalscalefactor = gesture.scaleFactor()

			z = self.pinch_start_z / totalscalefactor
			delta = z - view[11]
			view[11] = z
			view[15] -= delta
			view[16] -= delta
			self.cmd.set_view(view)

		return True

	def _event_x_y_mod(self, ev):
		return (
			int(self.fb_scale * ev.position().x()),
			int(self.fb_scale * (self.height() - ev.position().y())),
			get_modifiers(ev),
		)

	def mouseMoveEvent(self, ev):
		self.pymol.drag(*self._event_x_y_mod(ev))

	def mousePressEvent(self, ev, state=0):
		if ev.button() not in self._buttonMap:
			return
		self.pymol.button(self._buttonMap[ev.button()], state,
						  *self._event_x_y_mod(ev))

	def mouseReleaseEvent(self, ev):
		self.mousePressEvent(ev, 1)

	def wheelEvent(self, ev):
		button = get_wheel_button(ev)
		if not button:
			return
		args = self._event_x_y_mod(ev)
		self.pymol.button(button, 0, *args)
		self.pymol.button(button, 1, *args)

	##########################
	# OpenGL
	##########################

	def paintGL(self):
		self.pymol.draw()
		self._timer.start(0)

	def resizeGL(self, w, h):
		w = int(w * self.fb_scale)
		h = int(h * self.fb_scale)

		self.pymol.reshape(w, h, True)

	def updateFbScale(self, context):
		'''Update PyMOL's display scale factor from the window or screen context
		@type context: QWindow or QScreen
		'''
		self.fb_scale = context.devicePixelRatio()
		try:
			self.cmd.set('display_scale_factor', int(self.fb_scale))
		except BaseException as e:
			# fails with modal draw (mpng ..., modal=1)
			print(e)

	def initializeGL(self):
		# Scale framebuffer for Retina displays
		try:
			window = self.windowHandle()

			# QOpenGLWidget workaround
			if window is None:
				window = self.parent().windowHandle()

			self.updateFbScale(window)
			window.screenChanged.connect(self.updateFbScale)
			window.screen().physicalDotsPerInchChanged.connect(
					lambda dpi: self.updateFbScale(window))

		except AttributeError:
			# Fallback for Qt4
			pass

	def _pymol_process(self):
		idle = self.pymol.idle()
		if idle or self.pymol.getRedisplay():
			self.update()

		self._timer.start(20)

	##########################
	# drag n drop
	##########################

	#def dragEnterEvent(self, event):
	#	if event.mimeData().hasUrls:
	#		event.accept()
	#	else:
	#		event.ignore()

	#def dropEvent(self, event):
	#	if event.mimeData().hasUrls:
	#		for url in event.mimeData().urls():
	#			if url.isLocalFile():
	#				url = url.toLocalFile()
	#			else:
	#				url = url.toString()
	#			self.gui.load_dialog(url)
	#		event.accept()

	def load_receptor(self, name, molecule):
		#view = self.cmd.get_view()
		self.receptor = name
		self.cmd.read_pdbstr(molecule, name)


		#self.cmd.draw_box(0, 0, 0, 10, 10, 10)

		#self.cmd.get_box()
		self.cmd.get_box(self.receptor)
		#xyz = self.cmd.get_coordset('las2', copy=0)
		#print(xyz)
		#self.cmd.set_view(view)

	def sidebar_controler(self, flag):
		self.cmd.set('internal_gui', flag)
		#w = int(self.fb_scale * self.width())
		#h = int(self.fb_scale * self.height())
		self.pymol.reshape(-1, -1, True)


import os
import sys

babel_path = os.path.join(sys._MEIPASS, 'openbabel')
#os.environ['BABEL_LIBDIR'] = babel_path
#os.environ['BABEL_DATADIR'] = os.path.join(babel_path, 'data')
os.environ['BABEL'] = babel_path

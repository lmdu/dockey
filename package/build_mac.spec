# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(['../src/main.py'],
             pathex=[],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=['hooks'],
             hooksconfig={},
             runtime_hooks=[],
             excludes=['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter', 'numpy',
                 'sqlite3', 'pkg_resources', 'multiprocessing', 'PySide6.QtWebEngineCore',
                 'PySide6.QtWebEngineWidgets', 'PySide6.QtWebEngineQuick', 'PySide6.QtQml',
                 'PySide6.QtDataVisualization', 'PySide6.QtCharts', 'PySide6.QtConcurrent',
                 'PySide6.QtDBus', 'PySide6.QtDesigner', 'PySide6.QtNetwork', 'PySide6.QtNetworkAuth',
                 'PySide6.QtRemoteObjects', 'PySide6.QtQml', 'PySide6.QtQuickControls2',
                 'PySide6.QtQuickWidgets', 'PySide6.QtScxml', 'PySide6.QtSql', 'PySide6.QtStateMachine',
                 'PySide6.QtSerialPort', 'PySide6.QtXml', 'PySide6.Qt3DAnimation', 'PySide6.Qt3DCore',
                 'PySide6.Qt3DExtras', 'PySide6.Qt3DInput', 'PySide6.Qt3DLogic', 'PySide6.Qt3DRender',
                 'PySide6.QtQuick'
             ],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts, 
          [],
          exclude_binaries=True,
          name='Dockey',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None,
          icon='../src/icons/logo.icns')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas, 
               strip=False,
               upx=True,
               upx_exclude=[],
               name='Dockey')
app = BUNDLE(coll,
             name='Dockey.app',
             icon='../src/icons/logo.icns',
             bundle_identifier=None)
# -*- mode: python ; coding: utf-8 -*-
import sys

block_cipher = None

if sys.version.startswith('3.8'):
    ob_libs = '/usr/lib/openbabel/3.0.0/*'
    ob_share = '/usr/share/openbabel/3.0.0/*'
else:
    ob_libs = '/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/*'
    ob_share = '/usr/share/openbabel/3.1.1/*'

a = Analysis(['../src/main.py'],
             pathex=[],
             binaries=[(ob_libs, 'openbabel/lib')],
             datas=[(ob_share, 'openbabel/data')],
             hiddenimports=[],
             hookspath=['hooks'],
             hooksconfig={},
             runtime_hooks=['hooks/pyi_rth_openbabel.py'],
             excludes=['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter'],
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
          icon='../src/icons/logo.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas, 
               strip=False,
               upx=True,
               upx_exclude=[],
               name='Dockey')

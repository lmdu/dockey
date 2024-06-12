from PyInstaller.compat import is_win, is_darwin

block_cipher = None

if is_win:
    icon_file = ['../src/icons/logo.ico', '../src/icons/dock.ico']
    datas = []
elif is_darwin:
    icon_file = '../src/icons/logo.icns'
    datas = [('../src/icons/dock.icns', '.')]
else:
    icon_file = '../src/icons/logo.ico'
    datas = []


a = Analysis(['../src/main.py'],
             pathex=[],
             binaries=[],
             datas=datas,
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
          upx=False,
          console=False,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None,
          icon=icon_file)
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas, 
               strip=False,
               upx=False,
               upx_exclude=[],
               name='Dockey')

if is_darwin:
    app = BUNDLE(coll,
                 name='Dockey.app',
                 icon=icon_file,
                 bundle_identifier=None,
                 info_plist={
                    'CFBundleDocumentTypes': [
                        {
                            'CFBundleTypeName': 'Dockey Project File',
                            'CFBundleTypeIconFile': 'dock.icns',
                            'CFBundleTypeRole': 'Editor',
                            'LSHandlerRank': 'Owner',
                            'LSItemContentTypes': ['app.Dockey.dock']
                        }
                    ],
                    'UTExportedTypeDeclarations': [{
                        'UTTypeIdentifier': 'app.Dockey.dock',
                        'UTTypeTagSpecification': {
                            'public.filename-extension': ['dock']
                        },
                        'UTTypeConformsTo': ['public.data'],
                        'UTTypeDescription': 'Dockey Project File',
                        'UTTypeIconFile': 'dock.icns'
                    }]
                })

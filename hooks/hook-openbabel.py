from PyInstaller.compat import is_win, is_linux, is_darwin
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs

if is_win:
	datas = collect_data_files('openbabel', subdir='data')
	binaries = collect_data_files('openbabel', excludes=['data/*', '*.exe'])

if is_linux:
	datas = collect_data_files('openbabel', subdir='data')
	binaries = collect_data_files('openbabel', subdir='plugin', include_py_files=True)
	binaries += collect_data_files('openbabel', subdir='../openbabel.libs', include_py_files=True)

if is_darwin:
	datas = collect_data_files('openbabel', subdir='data')
	binaries = collect_data_files('openbabel', subdir='plugin', include_py_files=True)
	binaries += collect_dynamic_libs('openbabel')

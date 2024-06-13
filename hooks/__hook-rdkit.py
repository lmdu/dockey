from PyInstaller.compat import is_win, is_linux, is_darwin
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs

#if is_win or is_linux:
#	binaries = collect_data_files('rdkit', subdir='../rdkit.libs', include_py_files=True)

#if is_darwin:
#	binaries = collect_dynamic_libs('rdkit')

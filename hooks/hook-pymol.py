from PyInstaller.compat import is_linux, is_darwin
from PyInstaller.utils.hooks import collect_submodules, collect_data_files, collect_dynamic_libs

hiddenimports = collect_submodules('pymol')

if is_linux:
	binaries = collect_data_files('pymol',
		subdir='../pymol.libs',
		include_py_files=True
	)

if is_darwin:
	binaries = collect_dynamic_libs('pymol')

datas = collect_data_files('pymol',
	subdir='pymol_path',
	include_py_files = True,
	excludes = ['**/__pycache__']
)

from PyInstaller.utils.hooks import collect_submodules, collect_data_files

hiddenimports = collect_submodules('pymol')

datas = collect_data_files('pymol', subdir='pymol_path',
	include_py_files = True,
	excludes = ['**/__pycache__']
)

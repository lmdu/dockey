import os
from PyInstaller.compat import is_win
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs

if is_win:
	#this can collect dll files
	files = collect_data_files('openbabel', subdir='data')
	datas = [(src, os.path.basename(des)) for src, des in files]
	files = collect_data_files('openbabel', excludes=['*.exe', '*.dll'])
	datas += [(src, '.') for src, des in files]
	files = collect_dynamic_libs('openbabel')
	binaries = [(src, '.') for src, des in files]

else:
	datas = collect_data_files('openbabel', subdir='data')
	binaries = collect_dynamic_libs('openbabel', search_patterns=['plugin/*.so'])

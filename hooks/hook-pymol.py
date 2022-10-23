from PyInstaller.compat import is_win
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

hiddenimports = collect_submodules('pymol')

if not is_win: 
	datas = collect_data_files('pymol', subdir='../pymol.libs')

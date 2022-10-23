from PyInstaller.compat import is_win
from PyInstaller.utils.hooks import collect_all, collect_data_files

#if is_win:
datas, binaries, _ = collect_all('openbabel')

if not is_win:
	datas += collect_data_files('openbabel', subdir='../openbabel.libs')

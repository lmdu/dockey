from PyInstaller.compat import is_win
from PyInstaller.utils.hooks import collect_all

if is_win:
	datas, binaries, _ = collect_all('openbabel')

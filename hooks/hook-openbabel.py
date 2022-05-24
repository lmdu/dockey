from PyInstaller.compat import is_linux
from PyInstaller.utils.hooks import collect_all
if not is_linux:
	datas, binaries, _ = collect_all('openbabel')
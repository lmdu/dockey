from PyInstaller.utils.hooks import collect_data_files
datas = collect_data_files('pdbfixer', includes=['soft.xml'])
datas += collect_data_files('pdbfixer', subdir='templates')

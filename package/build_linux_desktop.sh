#!/bin/sh

version=$1
packager=$2

wget https://github.com/goreleaser/nfpm/releases/download/v2.11.3/nfpm_2.11.3_Linux_x86_64.tar.gz
tar xzvf nfpm_2.11.3_Linux_x86_64.tar.gz

desktop="[Desktop Entry]
Version=${version}
Name=Dockey
Comment=a morden tool for molecular docking
GenericName=Molecular Docking
Keywords=Molecular;Docking;Drug
Exec=/usr/lib/Dockey/Dockey
Icon=/usr/share/icons/dockey/logo.svg
Terminal=false
Type=Application
Categories=Application
StartupNotify=true
"
echo "$desktop" > dockey.desktop

nfpmconfig="name: Dockey
arch: amd64
platform: linux
version: v${version}
section: default
priority: extra
maintainer: lmdu <adullb@qq.com>
description: a morden tool for molecular docking
vendor: Bioinformatics and Integrative Genomics
homepage: https://github.com/lmdu/dockey
license: MIT
contents:
  - src: ./Dockey
    dst: /usr/lib/Dockey
  - src: ./dockey.desktop
    dst: /usr/share/applications/dockey.desktop
  - src: ./logo.svg
    dst: /usr/share/icons/dockey/logo.svg
rpm:
  compression: lzma
deb:
  compression: xz
"
echo "$nfpmconfig" > nfpm.yaml

# copy logo file
cp ../src/icons/logo.svg ./logo.svg

if [ $packager == "deb" ]
then
  ./nfpm pkg -t Dockey-v${version}-amd64.deb
  tar -czvf Dockey-v${version}-ubuntu.tar.gz Dockey
elif [ $packager == "rpm" ]
  ./nfpm pkg -t Dockey-v${version}-amd64.rpm
  tar -czvf Dockey-v${version}-fedora.tar.gz Dockey
else
  echo $version
fi
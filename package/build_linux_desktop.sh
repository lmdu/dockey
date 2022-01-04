#!/bin/sh

version=$1

wget https://github.com/goreleaser/nfpm/releases/download/v2.11.3/nfpm_2.11.3_Linux_x86_64.tar.gz
tar xzvf nfpm_2.11.3_Linux_x86_64.tar.gz

desktop="[Desktop Entry]
Version=${version}
Name=Dockey
Comment=a morden tool for molecular docking
Exec=/usr/lib/Dockey/Dockey
Icon=/usr/share/icons/dockey_logo.svg
Terminal=false
Type=Application
Categories=Application
"
echo "$desktop" > dockey.desktop

nfpmconfig="name: Dockey
arch: amd64
platform: linux
version: v${version}
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
  - src: ./dockey_logo.svg
    dst: /usr/share/icons/dockey_logo.svg
rpm:
  compression: lzma
deb:
  compression: xz
"
echo "$nfpmconfig" > nfpm.yaml

# copy logo file
cp ../src/icons/logo.svg ./dockey_logo.svg

./nfpm pkg -p deb -t Dockey-v${version}-amd64.deb
./nfpm pkg -p rpm -t Dockey-v${version}-amd64.rpm

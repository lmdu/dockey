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
arch: x86_64
platform: linux
version: ${version}
maintainer: lmdu <adullb@qq.com>
description: a morden tool for molecular docking
vendor: Bioinformatics and Integrative Genomics
homepage: https://github.com/lmdu/dockey
license: MIT
files:
  ./Dockey/**/*: /usr/lib/Dockey
  ./dockey.desktop: /usr/share/applications/dockey.desktop
  ./dockey_logo.svg: /usr/share/icons/dockey_logo.svg
"
echo "$nfpmconfig" > nfpm.yaml

# copy logo file
cp ../src/icons/logo.svg ./dockey_logo.svg

./nfpm pkg -t Dockey-v${version}-x86_64.deb
./nfpm pkg -t Dockey-v${version}-x86_64.rpm

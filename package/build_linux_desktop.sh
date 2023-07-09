#!/bin/bash

version=$1
#packager=$2

wget https://github.com/goreleaser/nfpm/releases/download/v2.20.0/nfpm_2.20.0_Linux_x86_64.tar.gz
tar xzvf nfpm_2.20.0_Linux_x86_64.tar.gz
rm nfpm_2.20.0_Linux_x86_64.tar.gz

cat > dockey.desktop <<EOF
[Desktop Entry]
Version=${version}
Name=Dockey
Comment=a modern tool for molecular docking
GenericName=Molecular Docking
Keywords=Molecular;Docking;Drug
Exec=/usr/lib/Dockey/Dockey %f
Icon=dockey.svg
Terminal=false
Type=Application
Categories=Education
StartupNotify=true
MimeType=application/x-dock
EOF

cat > application-x-dock.xml <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-dock">
    <comment>Dockey Project File</comment>
    <glob pattern="*.dock"/>
  </mime-type>
</mime-info>
EOF

cat > nfpm.yaml <<EOF
name: Dockey
arch: amd64
platform: linux
version: v${version}
section: default
priority: extra
maintainer: lmdu <adullb@qq.com>
description: a modern tool for molecular docking
vendor: Bioinformatics and Integrative Genomics
homepage: https://github.com/lmdu/dockey
license: MIT
contents:
  - src: ./Dockey
    dst: /usr/lib/Dockey
  - src: ./dockey.desktop
    dst: /usr/share/applications/dockey.desktop
  - src: ./application-x-dock.xml
    dst: /usr/share/mime/packages/application-x-dock.xml
  - src: ./logo.svg
    dst: /usr/share/icons/hicolor/scalable/apps/dockey.svg
  - src: ./dock.svg
    dst: /usr/share/icons/hicolor/scalable/mimetypes/application-x-dock.svg
rpm:
  compression: lzma
deb:
  compression: xz
EOF

# copy logo file
cp ../src/icons/logo.svg ./logo.svg
cp ../src/icons/dock.svg ./dock.svg

uver=$(cat /etc/issue | cut -d " " -f2 | cut -d "." -f1)

if [ "$uver" = "20" ]
then
  linux="linux"
else
  linux="linux-modern"
fi

./nfpm pkg -t Dockey-v$version-$linux.deb
tar -czvf Dockey-v$version-$linux.tar.gz Dockey

#build appimage
wget --no-check-certificate --quiet https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage
chmod +x appimagetool-x86_64.AppImage

cp dockey.desktop Dockey
cp logo.svg Dockey/dockey.svg

mkdir -p Dockey/usr/share/icons/hicolor/scalable/apps
cp logo.svg Dockey/usr/share/icons/hicolor/scalable/apps

cat > Dockey/AppRun <<EOF
#!/bin/bash
appdir=$(dirname $0)
exec "$appdir/Dockey" "$@"
EOF
chmod 755 Dockey/AppRun

cat > Dockey/dockey.desktop <<EOF
[Desktop Entry]
Name=Dockey
Comment=a modern tool for molecular docking
GenericName=Molecular Docking
Keywords=Molecular;Docking;Drug
Exec=Dockey %F
Icon=dockey
Terminal=false
Type=Application
Categories=Education
MimeType=application/x-dock
X-AppImage-Version=${version}
EOF

mkdir -p Dockey/usr/share/metainfo
cat > Dockey/usr/share/metainfo/dockey.appdata.xml <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<component type="desktop-application">
<id>com.big.dockey</id>
<metadata_license>CC0-1.0</metadata_license>
<project_license>MIT</project_license>
<name>Dockey</name>
<summary>molecular docking tool</summary>
<description>
  <p>Dockey is an integrated tool for molecular docking and virtual screening</p>
</description>
<screenshots>
  <screenshot type="default">
    <caption>Dockey</caption>
    <image>https://raw.githubusercontent.com/lmdu/dockey/main/src/icons/logo.svg</image>
  </screenshot>
</screenshots>
<url type="homepage">https://github.com/lmdu/dockey</url>
</component>
EOF

./appimagetool-x86_64.AppImage --appimage-extract-and-run Dockey Dockey-v$version-$linux.AppImage
rm appimagetool-x86_64.AppImage

#./nfpm pkg -t Dockey-v${version}-amd64.rpm

#if [ "$packager" = "deb" ]
#then
#  ./nfpm pkg -t Dockey-v${version}-amd64.deb
#  #tar -czvf Dockey-v${version}-ubuntu.tar.gz Dockey
#elif [ "$packager" = "rpm" ]
#then
#  ./nfpm pkg -t Dockey-v${version}-amd64.rpm
#  #tar -czvf Dockey-v${version}-centos.tar.gz Dockey
#else
#  echo $version
#fi
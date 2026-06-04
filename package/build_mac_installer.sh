#!/bin/bash

version=$1
arch=$2

brew install create-dmg
cd dist
create-dmg \
	--volname "Dockey Installer" \
	--volicon "../src/icons/logo.icns" \
	--hdiutil-quiet \
	"Dockey-v${version}-macos-${arch}.dmg" \
	"Dockey.app"

pkgbuild \
	--identifier dulab.big.dockey \
	--component \
	Dockey.app \
	Dockey-Component.pkg \
	--install-location /Applications

productbuild \
	--synthesize \
	--package Dockey-Component.pkg \
	Dockey-distribution.xml

productbuild \
	--distribution Dockey-distribution.xml \
	--package-path . \
	Dockey-v${version}-macos-${arch}.pkg
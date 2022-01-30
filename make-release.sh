#TARGET=SailfishOS-4.3.0.12-armv7hl
TARGET=SailfishOS-4.3.0.12-aarch64

LANGS_RPM=https://openrepos.net/sites/default/files/packages/990/okboard-full-0.6.34-1.armv7hl.rpm

mkdir okb && cd okb
git clone https://git.tuxfamily.org/okboard/okb-engine.git
git clone https://git.tuxfamily.org/okboard/okb-keyboard.git

cd okb-engine/db
curl "$LANGS_RPM" |sfdk engine exec rpm2cpio |cpio -i -d './usr/share/okboard/*.*'
mv usr/share/okboard/* .
gunzip -d *.gz

cd ../../okb-keyboard
~/SailfishOS/bin/sfdk config target="$TARGET"
~/SailfishOS/bin/sfdk tools package-install "$TARGET.default" python3-devel \
    qt5-qtdeclarative-qtquick-devel meego-rpm-config git fakeroot \
    libsailfishapp-devel
~/SailfishOS/bin/sfdk engine exec mkdir -p /home/mersdk/rpmbuild/{BUILD,BUILDROOT,RPMS,SOURCES,SPECS,SRPMS}
~/SailfishOS/bin/sfdk build-shell ./release.sh
~/SailfishOS/bin/sfdk engine exec cp -r /home/mersdk/rpmbuild/RPMS "$PWD/"

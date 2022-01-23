# OkBoard
hints to build okboard

# Build OkBoard
Okboard-Engine contains a readme.md file with How to build and deploy
The section is a bit outdated
You can connect to build vm using this correct string:
`ssh -i <SDK install dir>/vmshare/ssh/private_keys/sdk -p 2222 mersdk@localhost`

which you can actually get from:
QtCreator->Options->SailfishOS
- on build engine tab
- in connection

after successfully connected, i did skip this: sb2-config -d SailfishOS-armv7hl`
and did that:

list all available targets:
sb2-config -l
set propper target using:
sb2-config -d SailfishOS-4.3.0.12-aarch64

then followed article to use sb2

then exit sb2

then created all subfolders as required

then created a folder okboard-root

then copy extracted stuff to 

scp -r -i ~/SailfishOS/vmshare/ssh/private_keys/sdk -P 2222  okb-engine-0.6.20 mersdk@localhost:/home/mersdk/okboard-root

scp -r -i ~/SailfishOS/vmshare/ssh/private_keys/sdk -P 2222  OkBoard-Engine-Git mersdk@localhost:/home/mersdk/okboard-root

rem: check what is in what .. folder
copy engine as okb-engine into ..
rename ok-board-6.. to ok-board
cp okb-engine-0.6.20 okb-engine -r

#cp okb-engine-0.6.20 OkBoard-Engine-Git/okb-engine -r
then run release.sh from okboard-engine-git

then i had an error in build: fatal error: Python.h: No such file or directory

sudo zypper in python3-devel

you have really to copy the language files to db folder and you have to unzip all the .gz files there
using gzip -d ....gz

then build runs till 

fatal: not a git repository (or any of the parent directories): .git

scp -r -i ~/SailfishOS/vmshare/ssh/private_keys/sdk -P 2222  /home/pawel/Downloads/okboard-0.6.34.tar.gz   mersdk@localhost:/home/mersdk/rpmbuild/SOURCES
scp -r -i ~/SailfishOS/vmshare/ssh/private_keys/sdk -P 2222 /home/pawel/Downloads/okb-engine-0.6.20.tar.gz  mersdk@localhost:/home/mersdk/rpmbuild/SOURCES

both versions need to be same, so i had to rename okb-engine to 6.34

finally a real build error:
rror: Failed build dependencies:
	pkgconfig(Qt5Core) is needed by okboard-full-0.6.34-1.i386
	pkgconfig(Qt5Gui) is needed by okboard-full-0.6.34-1.i386
	pkgconfig(Qt5Qml) is needed by okboard-full-0.6.34-1.i386
	pkgconfig(Qt5Quick) is needed by okboard-full-0.6.34-1.i386
	pkgconfig(sailfishapp) >= 0.0.10 is needed by okboard-full-0.6.34-1.i386


for some reasons i do need in build root board and engine twice there, once with, once without version
then there was an error maybe introduced by me in curve_match.cpp on line 908 or so TH should be MAX_DEVICE_WIDTH or something, go and check in real files

then finaly i get an error:
in okboard-full.spec ln 55:
mv cfslm*.so cfslm.so

this returns target cfslm.so is not a directory 
possible to something like here

https://github.com/triton-inference-server/server/issues/1888








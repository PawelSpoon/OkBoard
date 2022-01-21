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









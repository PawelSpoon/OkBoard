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



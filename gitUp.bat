echo off

set major=0
set minor=5
set sub=2
set build=2

git status
pause

git add *
git status
pause

echo Enter commit message
set /P "Message="

git commit -m "v%major%.%minor%.%sub% #%build% - %Message%"

git push https://github.com/dilara255/Phase-Field-Monolayers.git

pause
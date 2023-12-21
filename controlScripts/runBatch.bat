@echo off

REM You can run a sequence of runners like this, no input required
REM You can also create copies of this with different runners, for parallel processing : p

echo.|python runner.py > logRunner.txt

timeout /t -1
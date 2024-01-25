@echo off

echo.|python runner.py > logRunner.txt
REM You can create multiple runners and run them in a sequence like this

REM You can also create copies of this with different runners, for parallel processing : p

timeout /t -1
#!/usr/bin/env python3
import subprocess
import hashlib
import time
import os
import shutil

makeGifFilename = "makeGif.py"

#all fields can be set to 'default'

simToRun = 'default' #only 1 (single layer) is currently implemented
lambda_ = 'default' #def = 3, 7.824813 used on some tests
gamma = 'default' #def = 0.06, 0.043986 used on some tests
dt = '2' #def = 1, may be brought down depending on other parameters
cells = 'default' # Must be > 0, and if initial conditions ask for cells, is set to at least 1
width = 'default' #default depends on the sim chosen. Must be > 0
height = 'default' #default depends on the sim chosen Must be > 0
initialCond = 'default' #default depends on the sim chosen
bias = 'default' #bias only applies to randomized initial conditions
seed = 'default' #some known-good seeds are available on depend/fAux/include/fAux/API/prng.hpp
method = 'default' #Numerical method. 0-3, each slower then the previous : )
startPaused = 'default' #starting paused here will just hang forever, so leave default or 0
changePerElmPerStepToStop = 'default' #default waits for fairly minimal changes
maxSteps = '15000' #default gives enough time for most tested networks to converge
stepsPerCheck = 'default' #default is 5k on debug and 25k on release (~100s)
changePerCheck = 'default' #default is 5k on debug and 25k on release (~100s)
useAdaptativeDt = '0' #default is 0 (false), use 0 or 1

time = str(time.time())

keyInt = 0

while keyInt == 0:
    key = hashlib.blake2b(time.encode(), digest_size=4, 
                      key=lambda_.encode(), salt=gamma.encode(), person=width.encode())

    keyInt = int(key.hexdigest(), 16)

argumentsString = (r'..\bin\Release-windows-x86_64\controlCL\controlCL.exe '
+ simToRun + ' '
+ lambda_ + ' '
+ gamma + ' '
+ dt + ' '
+ cells + ' '
+ width + ' '
+ height + ' '
+ initialCond + ' '
+ bias + ' '
+ seed + ' '
+ method + ' '
+ startPaused + ' '
+ changePerElmPerStepToStop + ' '
+ maxSteps + ' '
+ stepsPerCheck + ' '
+ changePerCheck + ' '
+ str(keyInt) + ' '
+ useAdaptativeDt
)

print('\n' + argumentsString + '\n')

subprocess.run(argumentsString)

print('\nWill try to acquire the path to the saved data for processing...\n')

keyFileName = str(keyInt) + ".txt"

keyFile = open(keyFileName)
dirPath = keyFile.readline()

keyFile.close()
os.remove(keyFileName)

os.system('python ' + makeGifFilename + ' ' + dirPath)
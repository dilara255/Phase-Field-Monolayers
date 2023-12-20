#!/usr/bin/env python3
import subprocess

#all fields can be set to 'default'

simToRun = '1' #only 1 (single layer) is currently implemented
lambda_ = '7.824813' #def = 3, 7.824813 used on some tests
gamma = '0.043986' #def = 0.06, 0.043986 used on some tests
dt = '1.4' #def = 1, may be brought down depending on other parameters
cells = '5' # Must be > 0, and if initial conditions ask for cells, is set to at least 1
width = '64' #default depends on the sim chosen. Must be > 0
height = '96' #default depends on the sim chosen Must be > 0
initialCond = '0' #default depends on the sim chosen
bias = '0.5' #bias only applies to randomized initial conditions
seed = '2' #some known-good seeds are available on depend/fAux/include/fAux/API/prng.hpp
method = '1' #Numerical method. 0-3, each slower then the previous : )
startPaused = 'default' #starting paused here will just hang forever
changePerElmPerStepToStop = '0.00001' #default waits for fairly minimal changes
maxSteps = '100' #default gives enough time for most tested networks to converge

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
+ maxSteps)

print(argumentsString)

subprocess.run(argumentsString)
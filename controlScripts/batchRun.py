#!/usr/bin/env python3
import subprocess

simToRun = 'default'
lambda_ = 'default'
gamma = 'default'
dt = 'default'
cells = 'default'
width = 'default'
height = 'default'
initialCond = 'default'
bias = 'default'
seed = 'default'
method = 'default'
startPaused = 'default'
changePerElmPerStepToStop = '0.00001'
maxSteps = 'default'

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
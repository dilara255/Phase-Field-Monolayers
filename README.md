# Phase-Field-Monolayers

A simulation of cell monolayers using the **phase-field model**.

Main contents:
- **simulation** holds the actual simulation. It's compiled as a static library.
- **controlCL** is thin, command-line only control application.
- **controlGUI** is a control application with a GUI that allows some real-time control and tinkering with the simulation. Mainly meant for exploration and testing.
- **controlScripts** are scripts which call controlCL, to automate batches of simulations and possibly some post-processing.

## Dependencies and requirements:

- Uses fAux as a utility library (included).
- The **GUI application** is built on top of fViz2D (also included), and **requires openGL 3.3** or better.
- The build system is controlled by **premake5**.
- **Expects 64 bit Windows or Linux**.
- The scripts use **Python 3**.

## Compilation:

- premake5.lua describes the solution's configuration. It's set up to use **clang**.
- run **"./premake5 gmake"** to use **make**, or run **generateProjects_VS2019.bat** (on windows) for **VS**.
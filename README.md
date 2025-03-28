# CH Phase-Field-Monolayers (wip)

A simulation of cell monolayers using the **phase-field model** and **Cahn-Hilliard potential**.

Main contents:
- **simulation** holds the actual simulation. It's compiled as a static library.
- **controlCL** is a thin, command-line only control application. Ideal to be called from scripts.
- **controlGUI** is a control application with a GUI that allows some real-time control and tinkering with the simulation. Mainly meant for exploration and testing.
- **controlScripts** are scripts which call controlCL, to automate batches of simulations and possibly some post-processing.

## Current state:

Simulation of a single Cahn-Hilliard field, with support for FTCS, Heun, RK2 and RK4.
Periodic density and absolute change checks and saving of these and of the field data (in .pgm, .jpg and/or .bin).
Parameters and configurations can be defined via command line or GUI (supports restart with new configurations and realtime change of parameters).
- Next goal: support for multiple fields (using subdomains) and intractions between them.

## Dependencies and requirements:

- Uses fAux as a utility library (included). The **simulation** project doesn't need to link to it (but uses some headers).
- The **GUI application** is built on top of fViz2D (also included), and **requires openGL 3.3** or better.
- The build system is controlled by **premake5** (tested for gmake and Visual Studio). Windows binary included in the depend folder.
- **Expects 64 bit Windows or Linux**.
- The scripts use **Python 3**.

## Compilation:

- premake5.lua describes the solution's configuration. It's set up to use **clang**.
- run **"premake5 gmake"** to use **make**, run **generateProjects_VS2022.bat** (on windows) for **VS**, or run premake5 with any other supported option.

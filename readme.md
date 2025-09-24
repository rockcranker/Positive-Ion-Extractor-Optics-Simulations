## Beam optics simulations for MINION positive ion extractor diagnostic
Igor Morozov (igoryau@gmail.com), 26.09.2025

### Overview
Provides tools for optic simulations and related analysis using IBSIMU simulation package. Can be used to perform 3D beam optic simulations for the triode extractor geometry, modeling extraction of positive hydrogen ions from plasma. Based on these simulations analysis of extracted particles and beam divergence and a simulation of Beam Emission Spectroscopy data can be done.  


### Requirements
IBSIMU: Follow installation instructions (https://ibsimu.sourceforge.net/installation.html), and update makefile based on install location. See examples on the website for reference.

ROOT: Installed by default on RAT2 and RAT4 servers. If installed correctly, should work with no modifications to makefile. Used for histograms, plots, curve fitting. Code could be modified to remove dependency if needed or if external tools are used.

### Workflow
The suggested workflow is:

Optics simulations &rarr; postprocessing &rarr; BES simulation &rarr; BES data analysis

This is done by moving files through the corresponding folders and running the programs in the order:

Triode3D &rarr; Postp &rarr; BES_sim


### Triode3D
This is the main simulation file. Relevant simulation settings can be adjusted in the triode3d.cpp to specify parameters, geometry, iterations, etc. Then, save the file and compile it using the terminal command
```console
[foo@bar]$ make
```
The simulation can then run directly through terminal using the commands:
```console
[foo@bar]$ ./triode3d
```
Alternatively, for large simulations the job can be submitted to easyPBS to run on the server in the background. First, set the correct settings and allocate appropriate resources in the triode.pbs file. Then, submit with the command:
```console
[foo@bar]$ qsub triode.pbs
```

Once the simulation finishes, the output data files epot.dat, geom.dat, pdb.dat, and parameters.txt will be created, alongside two beam trajectory plots. These files are used in the next analysis steps

### Postp
This is the main file used for all analysis of simulation results. Existing analysis includes divergence calculations, beam power deposition, and beam profile visualizations. The particle database is also cleaned up by culling irrelevant particles, reducing filesize. ROOT is required to create and fit the histograms, but this dependency can be removed if different analysis tools are used. Additional analysis and output can be built into this program as needed.

To use, the output files of a simulation (epot.dat, geom.dat, pdb.dat, and parameters.txt) should be moved to this folder. Postp.cpp can then be compiled and run with 
```console
[foo@bar]$ make
[foo@bar]$ ./postp
```
By default, outputs are moved into a new folder with a name determined by the string on the first line of the parameters.txt file. This string can be generated automatically in the simulation code or set manually as needed.

### BES_sim
Monte-carlo simulation of expected BES results based on particle trajectories. Requires pdb.dat and parameters.txt from an optics simulation as input. ROOT is once again needed to create and fit histograms, though results can be output to a file instead if needed. As before, the files should be moved into the folder, settings are adjusted in the .cpp files, and the code can be compiled with:
```console
[foo@bar]$ make
```
The BES simulation can be run with 
```console
[foo@bar]$ ./BES_sim
```
or
```console
[foo@bar]$ qsub BES_sim.pbs
```

After the program finishes, the resulting histogram can be fitted quickly using analysis.cpp. Make sure the angle setting is the same as in BES_sim, then run
```console
[foo@bar]$ ./analysis
```
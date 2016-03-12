PeldorFit2015
=========
The program PeldorFit performs analysis of the orientation-selective Pulsed ELectron-electron Double Resonance (PELDOR or DEER) signals. The orientation-selective PELDOR signals contain information about the distance between two spin centers in an ensemble of identical molecules and about the relative orientations of these centers. To extract this information, the program uses a simplified model of a spin pair for which several assumptions are taken:

1) spins are considered as single-point objects;

2) geometric parameters of the model have either normal or uniform distribution; 

3) correlation between individual geometric parameters of the model is neglected.

The geometry of the model is optimized via genetic algorithm provide the best match between experimental PELDOR signals and the corresponding PELDOR signals simulated for the model.

***

PeldorFit2015 is the second version of the program PeldorFit. The changes to the first version are:

1) The ~2 times faster performance (still in order of hours on usual PC).

2) The design of the configuration file is reconsidered to make it easier and more intuitive.

3) All fitting parameters can have eigher uniform or normal distribution.

***

Further description of the program can be found in the manual and in the paper (see below).

General Information
=========
The source code of the program PeldorFit is written using C++. The program uses two external open-access libraries, Intel TBB (https://www.threadingbuildingblocks.org/) and libconfig (http://www.hyperrealm.com/libconfig/). Both of them are required for compiling the program.

The compiled executables of the program for the Linux and Windows operating systems are gathered in the archive PeldorFit2015.zip.

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite:
D. Abdullin, G. Hagelueken, R. I. Hunter, G. M. Smith, O. Schiemann, Geometric model-based fitting algorithm for orientation-selective PELDOR data, Mol. Phys. 2015, 113, 544-560.

PeldorFit2015
=========
The program PeldorFit is developed for analysis of orientation-selective Pulsed ELectron-electron Double Resonance (PELDOR or DEER) data.

***

PeldorFit2015 is the second version of the program PeldorFit. The changes made to the first version are:

1) A two times faster performance (still in order of hours on usual PC).

2) The structure of a configuration file is reconsidered to make it easier and more intuitive.

3) All fitting parameters can now have eigher uniform or normal distribution.

***

Further description of the program can be found in the manual and in the paper (see below).

General Information
=========
The source code of the PeldorFit program is written in C++. The program uses two external open-access libraries, Intel TBB (https://www.threadingbuildingblocks.org/) and libconfig (http://www.hyperrealm.com/libconfig/). 

The program was already compiled for Linux and Windows, the corresponding executables are stored at https://github.com/dinarabdullin/PeldorFit2015/releases.

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite:
D. Abdullin, G. Hagelueken, R. I. Hunter, G. M. Smith, O. Schiemann, Geometric model-based fitting algorithm for orientation-selective PELDOR data, Mol. Phys. 2015, 113, 544-560.

PeldorFit2015
=========
A program PeldorFit is developed for the analysis of orientation-selective PELDOR (or alterntively DEER) data. This data encodes information about the inter-spin distance distribution and the relative orientation of spin centers in a certain spin pair. To extract this information, PeldorFit uses a simplified geometric model of the spin pair. This model is optimized via genetic algorithm until the PELDOR signals that are simulated for the model provide the best fit to the experimental PELDOR signals.

***

PeldorFit2015 is the second version of the program PeldorFit. The changes made to the first version are:

1) The ~2 times faster performance (still in order of hours on usual PC).

2) The design of the configuration file is reconsidered to make it easier and more intuitive.

3) All fitting parameters can now have eigher uniform or normal distribution.

***

Further description of the program can be found in the manual and in the paper (see below).

General Information
=========
The source code of the PeldorFit program is written in C++. The program uses two external open-access libraries, Intel TBB (https://www.threadingbuildingblocks.org/) and libconfig (http://www.hyperrealm.com/libconfig/). Both libraries are required for compiling the program. The program was already compiled for Linux and Windows, the corresponding executables are provided.

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite:
D. Abdullin, G. Hagelueken, R. I. Hunter, G. M. Smith, O. Schiemann, Geometric model-based fitting algorithm for orientation-selective PELDOR data, Mol. Phys. 2015, 113, 544-560.

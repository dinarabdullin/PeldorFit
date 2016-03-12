PeldorFit
=========
The program PeldorFit performs analysis of the orientation-selective Pulsed ELectron-electron Double Resonance (PELDOR or DEER) signals. These signals usually encode information about the distance between two spin centers and about their relative orientation. To extract this information, the program employs a simplified model of a spin pair for which several important assumptions are taken:

1) spins are considered as single-point objects, 

2) geometric parameters of the model have either normal or uniform distribution, 

3) correlation between individual geometric parameters of the model is neglected.

The geometry of the model is optimized via genetic algprithm provide a match between experimental PELDOR signals and the corresponding PELDOR signals simulated for the model.  

Further description of the optimization procedure can be found in the manual and in the paper (see below).

General Information
=========
The source code of the program PeldorFit is written using C++. The program uses two external open-access libraries, Intel TBB (https://www.threadingbuildingblocks.org/) and libconfig (http://www.hyperrealm.com/libconfig/). Both of them are required for compiling the program.

The compiled executables of the program for the Linux and Windows operating systems are gathered in the archive PeldorFit.zip.

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite:
D. Abdullin, G. Hagelueken, R. I. Hunter, G. M. Smith, O. Schiemann, Geometric model-based fitting algorithm for orientation-selective PELDOR data, Mol. Phys. 2015, 113, 544-560.

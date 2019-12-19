Listed below are the additional data to run the SWAT-MRMT-R model apart from those for the native SWAT model:

1) Reaction input file for single hyporheic zone in vertical and lateral direction is located in PFLOTRAN_INPUT/SingleHZ/.

2) Reaction input file for 10 sub-storage zones in vertical and lateral direction is located in PFLOTRAN_INPUT/MultipleHZ/.

3) Reaction database is located in reaction_database/.

4) The native NEXSS output is located in NEXSS_OUTPUT/.

5) NEXSS mapped to the simulated watershed is located in NEXSStoSWAT/SingleHZ/ for mean residence time and single hyporheic zone, and NEXSStoSWAT/MulitpleHZ/ for distributed residence time in 10 sub-storage zones.
   The source code and input files for the mapping are in NEXSStoSWAT/mapping/. 


The following notes are for running the SWAT-MRMT-R model on PNNL's constance supercomputer.

Notes for source codes:
1) SWAT originated from revision 664
2) PFLOTRAN originated from hash bc0bb2d 
3) petsc hash 1a9d3c3 is needed to compile PFLOTRAN and SWAT

To generate an executable:
1) Install petsc
   git clone https://gitlab.com/petsc/petsc petsc
   cd petsc
   git checkout 1a9d3c3

   module purge
   module load precision/i4
   module load intel
   module load openmpi
   module load mkl
   module load gcc/4.8.2
   module load cmake

   export PETSC_DIR=`pwd`
   export PETSC_ARCH=contance_intel15_O

   ./config/configure.py      \
   --with-fc=mpif90           \
   --with-cxx=mpicxx          \
   --with-clanguage=c         \
   --with-blas-lapack-dir=/share/apps/intel/2015u1/mkl \
   --with-shared-libraries=1 \
   --download-hdf5=yes \
   --download-parmetis=yes \
   --download-metis=yes \
   --with-debugging=0
 
   make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all

2) Download swatmrmtr
   
3) Compile PFLOTRAN to create library libpflotranchem.a
   cd swatmrmtr/pflotran/src/pflotran
   make pflotran_rxn

3) Compile SWAT-MRMT-R
   cd swatmrmtr/SWAT-MRMT-R
   make swatmrmtr.e

To run the model:
   swatmrmtr.e

Apart from the input files for the native swat model, four more files are required to run swatmrmtr.e:
1) Input file for the reaction parameters.  Default file name: pflotran.in.  
2) Input file for user defined database.  
3) Input file for hyporheic exchange flow rates and residence times.  Default file name: nxs2swat.txt.  
4) Initial chemical composition in the hyporheic zones.  Default file name: hz_init_conc.dat

The definition of parameters related to multirate mass transfer and hyporheic biological reactions can be found in function CyberRead in pflotran/src/pflotran/reaction_sandbox_pnnl_cyber_mrmt.F90.
Note: reactions in reaction_sandbox_pnnl_cyber_mrmt.F90 only assume two-step reactions for denitrification and an aerobic respiration reaction.  Otherwise, the code has to be modified following steps documented in https://www.pflotran.org/documentation/theory_guide/reaction_sandbox.html

Disclaimer
SWAT-MRMT-R is an open-source software distributed under the terms of the GNU Lesser General Public License as published by the Free Software Foundation either version 2.1 of the License, or any later version.



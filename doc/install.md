Installation       {#pageInstall}
============

- @subpage pageSource
- @subpage pageDependency
- @subpage pageBuild
- @subpage pageTrouble

<!-- ************************************************************ -->
@page pageSource Obtaining the source code

**TBD**

<!-- ************************************************************ -->
@page pageDependency Dependencies


%PEXSI is built on top of SelInv (for sequential LDLT factorization
and selected inversion), and SuperLU\_DIST (for distributed memory
parallel LU factorization).  Furthermore, the graph partitioning
module of SuperLU\_DIST requires ParMETIS.

Build ParMETIS
--------------

Download ParMETIS (latest version 4.0.2) from

http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz

\attention After untar the parmetis package, in Install.txt

    Edit the file metis/include/metis.h and specify the width (32 or
    64 bits) of the elementary data type used in ParMetis (and
    METIS). This is controled by the IDXTYPEWIDTH constant.

    For now, on a 32 bit architecture you can only specify a width
    of 32, whereas for a 64 bit architecture you can specify a width
    of either 32 or 64 bits.

In most cases,

    #define IDXTYPEWIDTH 32

works fine.

Build SuperLU_DIST
------------------

Download SuperLU_DIST (latest version 3.2) from

http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_3.2.tar.gz

\attention SuperLU_DIST is **NOT** guaranteed to work on OSX system.


<!-- ************************************************************ -->
@page pageBuild Build %PEXSI

Edit make.inc
-------------

The building sequenence of %PEXSI is not yet fully automatic (will
be updated in later versions).  Currently the first step is to edit
make.inc file.  There are two examples: make.hopper (NERSC machine)
and make.colosseum (MAC OSX system).  

For porting to other machines, copy make.inc.hopper to another
file (<hostid> can be the name of your machine)

    cp make.inc.hopper make.inc.<hostid> 

Set the environment variable (this line can be put into your
.bashrc)

    export HOSTID=<hostid> 
    export PEXSI_DIR = <directory>

If you use csh, use

    setenv HOSTID <hostid>
    setenv PEXSI_DIR = <directory>

Edit the variables in make.inc.<hostid> to your directory

    PEXSI_DIR     = 
    DSUPERLU_DIR  = 
    PARMETIS_DIR  =

Build the sequential SelInv
---------------------------

    cd ${PEXSI_DIR}
    cd external/SelInv
    make

Build the %PEXSI library
------------------------
    
    cd ${PEXSI_DIR}
    cd src
    make

Make examples
-------------

    cd ${PEXSI_DIR}
    cd examples
    make

For running the test examples, see \ref pag_example.

\attention The `tests` folder are for developing purpose, and not
all test subroutines are readily to be compiled.


<!-- ************************************************************ -->
@page pageTrouble Troubleshooting

**TBD**

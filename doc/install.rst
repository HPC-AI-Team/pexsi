Installation
----------------

Dependencies
============

PEXSI requires an external parallel :math:`LU` factorization or
:math:`LDL^T` factorization routine, and an external parallel matrix
reordering routine to reduce the fill-in of the factorization routine.

Starting from v1.0, PEXSI requires both symPACK and SuperLU_DIST.
symPACK is the default option for the :math:`LDL^T` factorization of
symmetric matrices, and use SuperLU_DIST as the default option for the
:math:`LU` factorization of unsymmetric matrices.  SuperLU_DIST can
also be used for symmetric matrices, by means of treating the matrix as
a general matrix but use symmetric reordering.

Starting from v1.0, PEXSI uses the PT-Scotch as the default package
for matrix reordering.  The ParMETIS package can also be used.

The installation procedure and dependencies of every version of the PEXSI
package may be slightly different. Please follow the documentation of the version
of the PEXSI package you are working with.
(provided in the :ref:`Download Page <pageDownload>` )


Build PT-Scotch
=============================

PT-Scotch can be downloaded from (latest version 6.0.0)
https://gforge.inria.fr/frs/download.php/31831/scotch_6.0.0.tar.gz

**PT-Scotch 6.0.5 seems to be incompatible with PEXSI. For the moment
please use 6.0.0 (contributed by Victor Yu, 6/20/2018) **

Follow the installation step to install PT-Scotch.
**In INSTALL.TXT, pay special attention to the following
sections in order to compile PT-Scotch correctly.**

    2.3) Integer size issues

    2.5) Threads issues


PT-Scotch is also METIS-Compatible.  See the following section in
INSTALL.TXT for more information.

    2.9) MeTiS compatibility library

In `src/` directory, you need
:: 
    make ptscotch 
    
to compile PT-Scotch.


.. note::  

  Just typing ``make`` will generate the Scotch library but not PT-Scotch.  
  Then all libraries will be given in ``lib/`` directory.**


Build symPACK
=============================


symPACK is a sparse symmetric matrix direct linear solver.
More information can be found at http://www.sympack.org/.

To use symPACK, first, download the package as follows
::
    git clone https://github.com/symPACK/symPACK.git  /path/to/sympack

Several environment variables can be set before configuring the build:

``SCOTCH_DIR`` = Installation directory for SCOTCH and PT-SCOTCH

Then, create a build directory, enter that directory and type::

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/sympack ...OPTIONS... /path/to/sympack


The ``...OPTIONS...`` can be one of the following:

- ``-DENABLE_METIS=ON|OFF``   to make METIS ordering available in symPACK (``METIS_DIR`` must be set in the environment)
- ``-DENABLE_PARMETIS=ON|OFF``   to make ParMETIS ordering available in symPACK (``PARMETIS_DIR`` must be set in the environment, ``METIS_DIR`` is required as well)
- ``-DENABLE_SCOTCH=ON|OFF``   to make SCOTCH / PT-SCOTCH orderings available in symPACK (``SCOTCH_DIR`` must be set in the environment)



Some platforms have preconfigured toolchain files which can be used by
adding the following option to the `cmake` command (To build on NERSC
Edison machine for instance)::

    -DCMAKE_TOOLCHAIN_FILE=/path/to/sympack/toolchains/edison.cmake     
    


A sample toolchain file can be found in `/path/to/sympack/toolchains/build_config.cmake` and customized for the target platform.


The `cmake` command will configure the build process, which can now start by typing::

    make
    make install

Additionally, a standalone driver for symPACK can be built by typing `make examples`

.. note:: 

  Since cmake also compiles UPCxx and GASNET, the compilation
  time may be long especially on certain clusters.


Build SuperLU_DIST
==================


Download SuperLU_DIST (latest version 6.1.0) from

http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_6.1.0.tar.gz

Follow the installation step to install SuperLU_DIST.

Our experience shows that on some machines it may be better
to build SuperLU_DIST with -O2 option than the more aggresive
optimization options provided by vendors.

 - In SuperLU_DIST, some functions conflict when both real
   and complex arithmetic factorization is needed. This can be temporarily
   solved by adding  `-Wl,--allow-multiple-definition` in the linking
   option.

 - In SuperLU_DIST, there could be some excessive outputs.
   This can be removed by going to the SRC/ directory of superlu, and
   comment out the line starting with `printf(".. dQuery_Space` in
   dmemory_dist.c. Do the same thing for the line starting with
   `printf(".. zQuery_Space..)` in zmemory_dist.c.

 - Please note that the number of processors for symbolic
   factorization cannot be too large when PARMETIS is used together with
   SuperLU. The exact number of processors for symbolic factorization is
   unfortunately a **magic parameter**. See :ref:`FAQ page <pageFAQ>`.



(Optional) Build ParMETIS
=========================

Download ParMETIS (latest version 4.0.3) from

http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz

Follow the installation step to install ParMETIS.

**After untar the ParMETIS package, in Install.txt**

    Edit the file metis/include/metis.h and specify the width (32 or
    64 bits) of the elementary data type used in ParMetis (and
    METIS). This is controled by the IDXTYPEWIDTH constant.

    For now, on a 32 bit architecture you can only specify a width
    of 32, whereas for a 64 bit architecture you can specify a width
    of either 32 or 64 bits.

**In our experience for most cases, the following setup work
fine.**::

    #define IDXTYPEWIDTH 32


Build PEXSI
===========

There are two ways to build PEXSI: 1) Using CMake 2) Using the standard
makefile system.

Build option 1: Use CMake
^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: 

  PEXSI requires CMake version 3.17+** (latest CMake can be
  downloaded at https://cmake.org/download/)


CMake is a meta-build system provided by Kitware. In essence, the purpose of
the CMake build system is to generate Makefiles which are customized to the
user's particular build environment. Generally, CMake operates by taking
information provided by the user in the form of CMake variables to notify
the build generator of things such as the location of dependency installations,
the enablement/disablement of software features, etc. In practice, this process
generally takes the form ::

    cmake -H<TOP SOURCE DIR> -B<BINARY DIR> -D<VAR1>=<VAL1> -D<VAR2>=<VAL2> ...

The project may then be compiled via ::

    make -C <BINARY DIR>

The following is a table of CMake variables which are influencial to the
PEXSI project

.. list-table:: PEXSI CMake Variables 
   :widths: 25 50 50
   :header-rows: 1

   * - Variable Name
     - Description
     - Possible Values (Default) 
   * - CMAKE_<Lang>_COMPILER
     - The <Lang> (C/CXX/Fortran) Compiler
     - (System Default)
   * - CMAKE_BUILD_TYPE
     - Build type
     - Release/Debug/RelWithDebInfo (Release)
   * - PEXSI_DEBUG_LEVEL
     - Level of PEXSI Debug print (Debug only)
     - 1-3 (1)
   * - PEXSI_ENABLE_OPENMP
     - Enable PEXSI OpenMP bindings
     - ON/OFF (ON)
   * - PEXSI_ENABLE_FORTRAN
     - Enable PEXSI Fortran bindings
     - ON/OFF (ON)
   * - PEXSI_ENABLE_PROFILE
     - Enable PEXSI Profiling (GNUProf)
     - ON/OFF (OFF)
   * - CMAKE_INSTALL_PREFIX
     - PEXSI Installation path
     - (None)
   * - CMAKE_PREFIX_PATH
     - Common installation prefix of dependencies 
     - (None)

PEXSI Also allows for the manualy specification of dependency locations
as either a prefix path or as a full linker

.. list-table:: PEXSI Dependency Variables 
   :widths: 25 50 50
   :header-rows: 1

   * - Variable Name
     - Description
     - Possible Values (Default) 

   * - <DEP>_PREFIX
     - Installation prefix of <DEP>
     - (None)
   * - <DEP>_LIBRARIES
     - A full linker for <DEP>
     - (None)
   * - <DEP>_INCLUDE_DIR
     - Location of <DEP> header files
     - (None)

Here, ``<DEP>`` is one of ``SuperLU_DIST``, ``METIS``, ``ParMETIS``,
``BLAS``, or ``LAPACK``. Note that the ``(PT-)SCOTCH`` and ``symPACK``
build paths are not supported through the build system at this time.

.. note:: 

  When specifying ``<DEP>_LIBRARIES``, the value must be a full linker,
  i.e. all of the libraries required to link to said dependency. e.g. ::

    SuperLU_LIBRARIES="-lsuperlu_dist -lparmetis -lmetis -lblas"

  We generally suggest that users specify ``<DEP>_PREFIX`` in preference
  over ``<DEP>_LIBRARIES`` whenever possible to avoid explicit specification
  of dependency trees such as these.
  


CMake also offers a mechanism to combine configuration parameters into
a single "toolchain" file, e.g. ::

  # my_toolchain.cmake
  set( CMAKE_C_COMPILER       gcc      )
  set( CMAKE_CXX_COMPILER     g++      )
  set( CMAKE_Fortran_COMPILER gfortran )
  set( SuperLU_DIST_PREFIX    "/home/linlin/Software/SuperLU_DIST_install/v6.1.0" )
  set( ParMETIS_PREFIX        "/home/linlin/Software/parmetis-4.0.3_install" )

Toolchains may be specified by ``CMAKE_TOOLCHAIN_FILE`` as a full path::

  cmake -H<TOP DIR> -B<BINARY DIR> -DCMAKE_TOOLCHAIN_FILE=$PWD/my_toolchain.cmake






..
  A few examples of the configuration options are given in the
  ``config/`` directory.
  
  
  Find ``build.sh`` with the most similar architecture, and copy to the main
  PEXSI directory (using Cori for example at NERSC, a CRAY X40 machine).
  ``${PEXSI_DIR}`` stands for the main directory of PEXSI. ::
  
      cd ${PEXSI_DIR}
      cp config/build.sh.CRAY_XC40.intel ./build.sh
      mkdir build; cd build;
  
  Edit the variables in ``build.sh``  ::
     
      PEXSI_INSTALL_DIR=Directory to install PEXSI
      DSUPERLU_DIR=Directory for SuperLU_DIST
      PARMETIS_DIR=Directory for ParMETIS 
      PTSCOTCH_DIR=Directory for PT-Scotch
  
  Edit the compiler options, for instance ::
  
      CC=cc
      CXX=CC
      FC=ftn
  
  Modify locations for other libraries if needed.  Then ::
      
      ../build.sh
  
  should prepare the ``build/`` directory.  If the configuration does not
  generate error messages, then ::
      
      make 
      make install
  
  should install PEXSI in ``PEXSI_INSTALL_DIR``. The examples files 
  are also compiled in ``build/examples/``. 



Tests
"""""

In the ``examples/`` folder::

    examples$ mpirun -n 1 ./driver_pselinv_complex_(suffix)

should return the diagonal of the matrix
:math:`(A + i I)^{-1}`
saved on the 0-th processor, where :math:`A` is the five-point
discretization of a Laplacian operator on a 2D domain.  The result can
be compared with `examples/driver_pselinv_complex.out` to check the
correctness of the result. 


The FORTRAN examples are given in ``build/fortran/``.  For more
information on the examples, see :ref:`Tutorial Page <pageTutorial>`.


.. note:: 

  If error messages occur, after debugging the compilation file,
  it is recommended to remove all files under ``build/`` first and then
  rerun ``build.sh``.





Build option 2: Use standard Makefile system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Configuration of PEXSI is controlled by a single ``make.inc`` file.
Examples of the ``make.inc`` file are given under the ``config/`` directory.

Find ``make.inc`` with the most similar architecture, and copy to the main
PEXSI directory (using Edison at NERSC for example, a CRAY X30 machine).
``${PEXSI_DIR}`` stands for the main
directory of PEXSI. ::

    cd ${PEXSI_DIR}
    cp config/make.inc.CRAY_XC30.intel make.inc

Edit the variables in make.inc.  ::
   
    PEXSI_DIR     = Main directory for PEXSI
    DSUPERLU_DIR  = Main directory for SuperLU_DIST
    PARMETIS_DIR  = Main directory for ParMETIS 
    PTSCOTCH_DIR  = Main directory for PT-Scotch

Edit the compiler options, for instance ::

    CC           = cc
    CXX          = CC
    FC           = ftn
    LOADER       = CC


The ``USE_SYMPACK`` option can be set to use the symPACK solver in
PEXSI. It is set to 0 by default. When set to 1, the ``SYMPACK_DIR`` variable
must be pointing to symPACK's installation directory.


.. note::

  - Starting from PEXSI v0.8.0, ``-std=c++11`` is required in ``CXXFLAGS``. 
  
  - Starting from PEXSI v0.9.2, ``-std=c99`` is required in ``CFLAGS`` to be
    compatible with SuperLU_DIST starting from v4.3.
  
  - For **FORTRAN** users, ``CPP_LIB=-lstdc++ -lmpi -lmpi_cxx`` is often needed.
    Check this if there is link error.
  
  - PEXSI can be compiled using ``debug`` or ``release`` mode in
    by the variable ``COMPILE_MODE`` in ``make.inc``.  This variable mainly controls the
    compiling flag ``-DRELEASE``.  The ``debug`` mode introduces tracing of call
    stacks at all levels of functions, and may significantly slow down the
    code.  For production runs, use ``release`` mode.
  
  - The ``USE_PROFILE`` option is for internal test purpose. Usually set this to 0.


The installation procedure and dependencies of every version of the PEXSI
package may be different. Please follow the documentation of the version
of the PEXSI package you are working with 
(provided in the :ref:`Download Page <pageDownload>` )

If make.inc is configured correctly,::
    
    make 
    make install

Should build the PEXSI library under the `build` directory ready to be
used in an external package.  If the FORTRAN interface is needed, type::

    make finstall

If examples are needed (not necessary if you use PEXSI in an external
package), type ::

    make examples

which will generate C examples in `examples/` directory and FORTRAN examples in
`fortran/` directory, respectively.::

    make all

will make the library and the examples. 



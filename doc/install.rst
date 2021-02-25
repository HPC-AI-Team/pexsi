Installation
----------------

.. warning::
  Often a software package has a "quick installation guide" which allows
  the software to be installed with a handful of generic lines (at least
  for the serial mode).  Unfortunately due to dependencies, we have not
  been able to work out a simple solution for PEXSI. We are very
  conscious about this and would like to find a better solution in the
  future. For now please be patient and go through the following steps..
  

Dependencies
============

PEXSI requires an external parallel :math:`LU` factorization or
:math:`LDL^T` factorization routine (default is SuperLU_DIST), and an
external parallel matrix reordering routine (default is ParMETIS) to
reduce the fill-in of the factorization routine. These packages need to
be downloaded and installed separately, before building the PEXSI
package.

..
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

Although the standard makefile system is supported, starting from v2.0,
we recommend the usage of the CMake system.


Build ParMETIS
==============

Download ParMETIS (latest version 4.0.3) from

http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz

..
  After untar the ParMETIS package, in Install.txt
  Edit the file metis/include/metis.h and specify the width (32 or
  64 bits) of the elementary data type used in ParMetis (and
  METIS). This is controled by the IDXTYPEWIDTH constant.
  For now, on a 32 bit architecture you can only specify a width
  of 32, whereas for a 64 bit architecture you can specify a width
  of either 32 or 64 bits.
  In our experience for most cases, the following setup work
  well.::
    #define IDXTYPEWIDTH 32

.. note::

  ParMETIS has its own (slightly non-standard) CMake system. Be sure to
  follow `Install.txt` and `BUILD.txt` after untaring the file.

  We provide an example bash script for building ParMETIS (see
  `build_parmetis.sh <https://bitbucket.org/berkeleylab/pexsi/src/cmake-refactor/config/build_parmetis.sh>`_)


Build SuperLU_DIST
==================


Download SuperLU_DIST (latest version 6.4.0) from

https://github.com/xiaoyeli/superlu_dist/archive/v6.4.0.tar.gz

.. note::

  We provide an example bash script for building SuperLU_DIST (see
  `build_superlu_dist.sh <https://bitbucket.org/berkeleylab/pexsi/src/cmake-refactor/config/build_superlu_dist.sh>`_)


..
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
  

.. warning::
   When the number of processors for symbolic factorization cannot be
   too large when PARMETIS is used together with SuperLU. The exact
   number of processors for symbolic factorization is unfortunately
   sometimes a **magic parameter**. See :ref:`FAQ page <pageFAQ>`.


Build PEXSI
===========

The recommended method to build PEXSI is to use CMake. You may also use
the :ref:`old makefile system <Makefile>`.

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
     - ON/OFF (OFF)
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
  set( SuperLU_DIST_PREFIX    "/yourdirectory/SuperLU_DIST_install/v6.1.0" )
  set( ParMETIS_PREFIX        "/yourdirectory/parmetis-4.0.3_install" )

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






Optional
========

.. toctree::
   :maxdepth: 1

   optpackage


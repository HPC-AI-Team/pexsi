.. _pageoptpackage:

PT-Scotch
^^^^^^^^^

PT-Scotch can be used to replace ParMETIS. (We prefer ParMETIS since
this is the default for SuperLU_DIST)


PT-Scotch can be downloaded from (latest version 6.0.0)

https://gforge.inria.fr/frs/download.php/31831/scotch_6.0.0.tar.gz

.. note::
  PT-Scotch 6.0.5 seems to be incompatible with PEXSI. For the moment
  please use 6.0.0 (contributed by Victor Yu, 6/20/2018) 

Follow the installation step to install PT-Scotch.
**In INSTALL.TXT, pay special attention to the following
sections in order to compile PT-Scotch correctly.**

    2.3) Integer size issues

    2.5) Threads issues

PT-Scotch is also METIS-Compatible.  See the following section in
INSTALL.TXT for more information.

    2.9) MeTiS compatibility library

In `src/` directory, you need :: 
    make ptscotch 
    
to compile PT-Scotch.


.. note::  
  Just typing ``make`` will generate the Scotch library but not PT-Scotch.  
  Then all libraries will be given in ``lib/`` directory.**




symPACK
^^^^^^^

symPACK can be used to replace SuperLU_DIST (for :math:`LDL^T` factorization but not :math:`LU` factorization)

symPACK is a sparse symmetric matrix direct linear solver.
More information can be found at http://www.sympack.org/.

To use symPACK, first, download the package as follows ::
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

The `cmake` command will configure the build process, which can now start by typing ::
    make
    make install

Additionally, a standalone driver for symPACK can be built by typing `make examples`

.. note:: 

  Since cmake also compiles UPCxx and GASNET, the compilation
  time may be long especially on certain clusters.

.. _Makefile:

Build PEXSI with the old Makefile system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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


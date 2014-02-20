Frequently asked questions {#pageFAQ}
==========================

General questions
-----------------

**Q: How do I know whether %PEXSI works for my application?**

A: %PEXSI may not necessarily be faster than the diagonalization method or
other competitive methods.  The simplest way to see whether %PEXSI
brings acceleration for your applications is to use %PEXSI to compute
the selected elements of the inverse for a typical matrix from your
applications.  See @ref pagePselinvComplex and @ref
driver_pselinv_complex.c for
how to do this.

**Q: Can I just use the selected inversion routine?**

A: The parallel selected inversion (PSelInv) is a standalone routine.
See @ref pagePselinvComplex for an example.


**Q: Does %PEXSI accelerate dense matrix computation?**

A: No.  The acceleration is based on the sparsity of the LU factor or
the Cholesky factor.  %PEXSI should not be fast if the matrix is dense
or nearly dense.

**Q: Does %PEXSI work for asymmetric matrices?**

A: Currently %PEXSI only works for (real or complex) symmetric matrices.
We plan to support asymmetric matrices in the future.  If you have
applications in mind and %PEXSI may bring significant acceleration to
your application, please contact us as in @ref pageTrouble, with some
description of your matrix and application.

**Q: I only found the driver routine for selected inversion of complex
symmetric matrices.  Why %PEXSI does not provide driver routines for
real arithmetic selected inversion?**

A: Real arithmetic selected inversion *is available* in %PEXSI.
However, due to legacy reason it is not available in an intuitive way,
and therefore we do not provide the driver routine in the first release.
We plan to make major change of the structure of the code in the next
version.  If you need the real version of selected inversion, the
easiest workaround is to use the complex arithmetic interface.  If speed
is crucial for your application, please contact us as in @ref
pageTrouble, with some description of your matrix and application.


Installation
------------

**To be added**


Performance
-----------

**Q: What if %PEXSI is numerically unstable or inaccurate for my
application?**

A: If you are using the pole expansion, the expansion converges
exponentially with respect to the number of poles.  So first increase
the number of poles until the error saturates.  After this step, error
only comes from the selected inversion phase. The selected inversion is
in principle an exact method, but may possibly suffer from numerical
instability issue due to the round off erros.  This is possible due to
the lack of dynamic pivoting strategies.  If this problem persists,
please contact us as in @ref pageTrouble, with some description of your
matrix and application.

**Q: When using ParMETIS/PT-Scotch, I got segmentaiton fault in the
factorization phase.**

A: We have observed that when ParMETIS/PT-Scotch is used, the number of
processors for symbolic factorization (`npSymbFact`) cannot be larger
than a magic number depending on the matrix and the machine.  This may
to be related to the parallel symbolic factorization routine in
SuperLU\_DIST.  If this problem happens, try to reduce `npSymbFact` to
a small number (such as 4 or 16), or even try to use the sequential
symbolic factorization if feasible.


**Q: Why the PPEXSI driver routines do not specify nprow/npcol?**

A: Due to legacy reason the number of processors used in the interface
routines for each parallel selected inversion process has to be a
**square number** (1,4,9,16,64 etc).  This will be changed in the next
version of the user interface.

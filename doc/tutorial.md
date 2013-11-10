@page pageTutorial Tutorial 

- @subpage pagePselinvComplex
- @subpage pagePEXSISolve


@page pagePselinvComplex Parallel selected inversion for a complex matrix
\tableofcontents

The computation of selected elements of an inverse matrix is a
standalone functionality of %PEXSI. If you are a C/C++ programmer, the
parallel selected inversion routine can be used as follows.

~~~~~~~~~~{.cpp}
    #include <fftw3.h>
    ...
    {
        fftw_complex *in, *out;
        fftw_plan p;
        ...
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        ...
        fftw_execute(p); /* repeat as needed */
        ...
        fftw_destroy_plan(p);
        fftw_free(in); fftw_free(out);
    }
~~~~~~~~~~

void @ref PPEXSISelInvInterface "PPEXSISelInvInterface" (




@page pagePEXSISolve Electronic structure calculation using PEXSI
 

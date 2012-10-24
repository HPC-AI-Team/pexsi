#ifndef __GETPOLE
#define __GETPOLE


// Avoid conflict with doublecomplex defined elsewhere (such as in
// SuperLU)
#ifndef _DOUBLECOMPLEX_
typedef struct { double r, i; } doublecomplex;
#endif

#ifdef __cplusplus
extern "C"{
#endif
int getpole_rho(doublecomplex* zshift, doublecomplex* zweight, 
	     int* Npole, double* temp, double* gap, double* deltaE,
	     double* mu);

int getpole_hmz(doublecomplex* zshift, doublecomplex* zweight, 
	     int* Npole, double* temp, double* gap, double* deltaE,
	     double* mu);

int getpole_egy(doublecomplex* zshift, doublecomplex* zweight, 
	     int* Npole, double* temp, double* gap, double* deltaE,
	     double* mu);
#ifdef __cplusplus
}
#endif

#endif

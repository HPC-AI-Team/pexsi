#ifndef __GETPOLE
#define __GETPOLE
typedef struct { double r, i; } doublecomplex;

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

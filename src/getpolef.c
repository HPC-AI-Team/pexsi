#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "getpolef.h"

double complex fd(double complex z, double beta, double mu);
double complex egy(double complex z, double beta, double mu);
double complex hmz(double complex z, double beta, double mu);
void ellipkkp(double* K, double* Kp, double *pL);
void ellipjc(double complex* psn, double complex* pcn,
	     double complex* pdn, double complex* pu, 
	     double *pL, int *flag) ;
int getpolef(double complex (*func)(double complex, double, double),
	      double complex* zshift, double complex* zweight, int* Npole,
	      double* temp, double* gap, double* deltaE, double* mu);


// Fermi-Dirac distribution fd(z) = 2/(1+exp(beta*z))
double complex fd(double complex z, double beta, double mu){
  double complex val, ez;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( creal(z) >= 0 ){
    ez = cexp(-beta*z);
    val = 2.0 * ez / (1.0 + ez);
  }
  else{
    ez = cexp(beta* z);
    val = 2.0 / (1.0 + ez);
  }
  return val;
}

// Energy function egy(z) = (z+mu) * fd(z) 
double complex egy(double complex z, double beta, double mu){
  double complex val, ez;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( creal(z) >= 0 ){
    ez = cexp(-beta*z);
    val = (z + mu) * 2.0 * ez / (1.0 + ez);
  }
  else{
    ez = cexp(beta* z);
    val = (z + mu) * 2.0 / (1.0 + ez);
  }
  return val;
}

// Helmholtz free energy function hmz(z) = -2/beta*log(1+exp(-beta*z))
double complex hmz(double complex z, double beta, double mu){
  double complex val, ez;
  /* LL: VERY IMPORTANT TO AVOID OVERFLOW/UNDERFLOW! */
  if( creal(z) >= 0 ){
    ez = cexp(-beta*z);
    val = -2.0/beta*clog(1.0+ez);
  }
  else{
    ez = cexp(beta* z);
    val = 2.0*z - 2.0/beta*clog(1+ez);
  }
  return val;
}


/*********************************************************************
   ELLIPKKP Computes the complete elliptic integral of the first kind,
   with complement.

   Input parameters:
     L             

   Output parameters
     K, Kp

   K is the value of the complete elliptic integral of the first kind,
   evaluated at M=exp(-2*pi*L), 0 < L < Inf. Kp is the result for
   complementary parameter, which is useful when M < EPS.  Even when M
   < 1e-6, the built-in ELLIPKE of MATLAB can lose digits of accuracy
   for KP.

   Recall that the elliptic modulus k is related to the parameter
   M by M = k^2.

   ELLIPKKP uses the method of the arithmetic-geometric mean described
   in 17.6 of M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
   Functions," Dover, 1965.  Same method as in ELLIPKE, only
   interchanging 1 and 1-m to find KP.

   When m=exp(-2*pi*L) is extremely small, use O(m) approximations.  

   Originally written by Toby Driscoll in 1999.

   Rewritten in C by 
   Lin Lin
   Computer Research Division, Lawrence Berkeley National Lab
   Last modified:  09-15-2011
*********************************************************************/
void ellipkkp(double* K, double* Kp, double *pL){
  double pi, m, a0, b0, s0, a1, b1, c1, w1;
  double mm, eps = 1e-15;
  int i1; 
  double L = (*pL);

  pi = atan(1.0)*4.0;
  if( L == 0 ){
    fprintf(stderr, "L == 0. STOP!\n");
    return;
  }
  if( L > 10.0 ){
    *K = pi / 2.0;
    *Kp = pi * L + log(4.0);
    return;
  }

  m = exp(-2*pi*L);
  
  /* Calculate K */
 
  a0 = 1.0;
  b0 = sqrt(1.0 - m);
  s0 = m;
  i1 = 0; mm = 1.0;
  while( mm > eps ){
    a1 = (a0+b0) * 0.5;
    b1 = sqrt(a0*b0);
    c1 = (a0-b0)/2;
    i1++;
    w1 = pow(2.0, i1) * c1 * c1;
    mm = w1;
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  }
  
  *K = pi / (2.0 * a1);

  /* Calculate Kp */
  a0 = 1.0;
  b0 = sqrt(m);
  s0 = 1.0-m;
  i1 = 0;  mm = 1.0;
  while (mm > eps ){
    a1 = (a0+b0)/2.0;
    b1 = sqrt(a0 * b0);
    c1 = (a0-b0)/2.0;
    i1++;
    w1 = pow(2.0, i1) * c1 * c1;
    mm = w1;
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  }
  *Kp = pi / (2.0 * a1);
  return;
}

/*********************************************************************
  ELLIPJC Calculates Jacobi elliptic functions for complex argument.

  Input parameters:
    u, L             

  Output parameters
    sn, cn, dn

  [sn,cn,dn] are the values of the Jacobi elliptic functions
  evaluated at complex argument U and parameter M=exp(-2*pi*L), 0 < L
  < Inf.  Recall that M = k^2, where k is the elliptic modulus.

  The entries of U are expected to lie within the rectangle |Re U| <
  K, 0 < Im U < Kp, where K,Kp are evaluated from ELLIPKKP.

  The built-in ELLIPJ of MATLAB can't handle compelx arguments, and
  standard transformations to handle this would require ELLIPJ
  called with parameter 1-M. When M < eps (or is even close),
  this can't be done accurately.

  The algorithm is the descending Landen transformation,
  described in L. Howell's PhD thesis from MIT. Additional
  formulas from Gradshteyn & Ryzhik, 5th ed., and Abramowitz
  & Stegun.  

  NOTE:  When calling this function outside *flag = 0
         *flag = 1 is only for recursive use.

  Originally written by Toby Driscoll in 1999.

  Rewritten in C by 
  Lin Lin
  Computer Research Division, Lawrence Berkeley National Lab
  Last modified:  10-28-2011
*********************************************************************/
void ellipjc(double complex* psn, double complex* pcn,
	     double complex* pdn, double complex* pu, 
	     double *pL, int *flag) {

  double K, Kp; 
  double L = (*pL);
  double m;
  double pi = atan(1.0)*4.0;
  double eps = 1e-15;
  double x, kappa, mu;
  int high;
  int ione = 1;
 
  /* C complex numbers */
  double complex sinu, cosu;
  double complex u;
  double complex snh, cnh, dnh, v;
  double complex sn1, cn1, dn1, denom;

  u = (*pu);

  /* Check and transform u in the upper half of the rectangle */
  high = 0;
  if( *flag == 0 ){
    ellipkkp(&K, &Kp, &L);
    if( cimag(u) > Kp * 0.5 ){
      high = 1;
      u = Kp - u;
    }
    m = exp(-2.0*pi*L);
  }
  else{
    m = L;
  }

  
  /* Case 1 : m is already small */
  if( m < 4.0 * eps ){
    sinu = csin(u);
    cosu = ccos(u);

    *psn = sinu + m/4.0 * (sinu * cosu - u) * cosu;

    *pcn = cosu + m/4.0 * (-sinu * cosu + u) * sinu;
    
    *pdn = 1.0 +  m/4.0 * (cosu * cosu - sinu * sinu - 1.0); 

  }
  /* Case 2 : m is big, call recursive formula */
  else{
    if( m > 1e-3 )
      kappa = (1.0-sqrt(1.0-m))/(1.0+sqrt(1.0-m));
    else{
      x = m/4;
      kappa = x*(1.0+x*(2.0+x*(5.0+x*(14.0+x*(42.0+x*132)))));
    }
    mu = kappa * kappa;
    v  = u / (1.0 + kappa);
    /* Call ellipjc recursively */
    ellipjc(&sn1, &cn1, &dn1, &v, &mu, &ione);
    denom = 1.0 + kappa * sn1 * sn1;
    
    *psn = (1.0+kappa) * sn1 / denom;
    *pcn = cn1 * dn1 / denom; 
    *pdn = (1.0-kappa*sn1*sn1) / denom;
  }


  if( high ){
    snh = (*psn);
    cnh = (*pcn);
    dnh = (*pdn);

    (*psn)  = -1.0/ (sqrt(m) * snh);
    (*pcn) = I * dnh / (sqrt(m)*snh);
    (*pdn) = I * cnh / snh;
  }
  return;
}





/*********************************************************************
  GETPOLEF generates the poles and weights for any function f that
  shares the same analytic structure with the Fermi-Dirac distribution
  with chemical potential mu and inverse temperature beta.

  Input:
		func      :    input function to be expanded by pole
		               expansion 
                Npole     :    the number of poles to be used.
		temp      :    temperature, unit(K)
		Gap       :    Energy gap defined to be min(abs(EV-mu)).
			       EV is the eigenvalue set of Hamiltonian,
			       unit(hatree) 
		deltaE    :    Spectrum width defined to be
			       max(EV)-min(EV). EV is the eigenvalue set
			       of Hamiltonian, unit(hartree) 
                mu        :    Chemical potential, unit(hartree)

  Output:
		
		zshift    :    Complex shift of poles.
		zweight   :    Weight of poles.
  
 
  Example:
    Pseudocode (MATLAB notation) using pole expansion to reconstruct the
    electron density.  NOTE: mu is included in zshift.

    Rho = zeros(N, 1);  
    for i = 1 : Npoles
      Rho = Rho + diag(imag( zweight(i) * inv(H - zshift(i))));
    end 

  Reference:

    L. Lin, J. Lu, L. Ying and W. E, Pole-based approximation of the
    Fermi-Dirac function, Chin. Ann. Math.  30B, 729, 2009 

 Author: 
 Lin Lin
 Computer Research Division, Lawrence Berkeley National Lab
 Last modified:  10-28-2011
*********************************************************************/
int getpolef(double complex (*func)(double complex, double, double),
	      double complex* zshift, double complex* zweight, int* Npole,
	      double* temp, double* gap, double* deltaE, double* mu){
  double K2au = 3.166815e-6, beta;
  double M, mshift, m2, M2, kr, L, K, Kp, coef;
  double complex t, sn, cn, dn, z, dzdt, zsq, funczsq;
  double pi = atan(1.0)*4.0;
  int i, j, Npolehalf, flag=0;
  beta = 1.0 / ((*temp) * K2au);

  if( (*Npole) % 2 != 0 ){
    fprintf(stderr, "Npole has to be an even number!\n");
    return 1;
  } 
  
  Npolehalf   = (*Npole) / 2;
  M           = *deltaE;
  mshift      = pow((pi/beta), 2.0);
  m2          = mshift + (*gap)*(*gap);
  M2          = M*M;
  kr          = (sqrt(M2/m2)-1.0)/(sqrt(M2/m2)+1.0);
  L           = -log(kr)/pi;
  ellipkkp(&K, &Kp, &L);

  for( j = 0; j < Npolehalf; j++){
    t   = (-K + (0.5 + j) / Npolehalf * 2.0 * K) + I * 0.5 * Kp;
    ellipjc(&sn, &cn, &dn, &t, &L, &flag);
   
    z    = sqrt(m2*M2) * (1.0/kr + sn) / (1.0/kr-sn) - mshift;

    /* General formula for zshift and zweight
     *
       zshift(j)  = zsqrt(j) + mu;
       zweight(j) = 2*K*sqrt(m2*M2)/(kr*pi*Npolehalf) / 
                    zsqrt(j) * dzdt(j) * func(zsqrt(j)); 
    */

    dzdt = cn * dn / ((1.0/kr-sn) * (1.0/kr-sn));   
    
    coef = 2.0 * K * sqrt(m2*M2) / (kr*pi*Npolehalf);
    
    
    /* The first Npolehalf poles */
    zsq     = csqrt(z);
    
    funczsq = (*func)(zsq, beta, *mu);
    
    zshift[j]  = (*mu) + zsq;
    zweight[j] = funczsq * dzdt / zsq * coef;

    /* The second Npolehalf poles */
    zsq  = -csqrt(z);
    
    funczsq = (*func)(zsq, beta, *mu);
    
    zshift[j+Npolehalf]  = (*mu) + zsq;
    zweight[j+Npolehalf] = funczsq * dzdt / zsq * coef;

  }

  return 0;
}


/* Wrapper function */

int getpole_rho(doublecomplex* zshift, doublecomplex* zweight, 
	     int* Npole, double* temp, double* gap, double* deltaE,
	     double* mu){
  return getpolef(&fd, (double complex*)zshift, (double complex*) zweight,
	   Npole, temp, gap, deltaE, mu);
}

int getpole_hmz(doublecomplex* zshift, doublecomplex* zweight, 
	     int* Npole, double* temp, double* gap, double* deltaE,
	     double* mu){
  return getpolef(&hmz, (double complex*)zshift, (double complex*) zweight,
	   Npole, temp, gap, deltaE, mu);
}

int getpole_egy(doublecomplex* zshift, doublecomplex* zweight, 
	     int* Npole, double* temp, double* gap, double* deltaE,
	     double* mu){
  return getpolef(&egy, (double complex*)zshift, (double complex*) zweight,
	   Npole, temp, gap, deltaE, mu);
}



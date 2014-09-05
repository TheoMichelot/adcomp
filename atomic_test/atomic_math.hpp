/* Namespace with double versions of R special math library */
namespace Rmath {
  #include <Rmath.h>
  // Macros do not respect the namespace limits.
  // Some of them will conflict with other TMB functions.
  #undef dnorm
  #undef pnorm
  #undef qnorm

  #include <R_ext/Applic.h>
  void integrand_D_incpl_gamma_shape(double *x, int n, void *ex){
    double* parms=(double*)ex;
    for(int i=0;i<n;i++) x[i] = exp(-x[i]) * pow(x[i],parms[0]-1.0) * pow(log(x[i]),parms[1]);
  }
  /* n'th order derivative of incomplete gamma wrt. shape parameter */
  double D_incpl_gamma_shape(double x, double shape, double n){
    double a=0;
    double b=x;
    double epsabs=1e-8;
    double epsrel=1e-8;
    double result=0;
    double abserr=10000;
    int neval=10000;
    int ier=0;
    int limit=100;
    int lenw = 4 * limit;
    int last=0;
    int* iwork =  Calloc(limit, int);
    double* work = Calloc(lenw, double);
    double ex[2];
    ex[0]=shape;
    ex[1]=n;
    Rdqags(integrand_D_incpl_gamma_shape, ex, &a, &b,
	   &epsabs, &epsrel,
	   &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    Free(iwork);
    Free(work);
    if(ier!=0)warning("Integrate incomplete gamma function unreliable");
    return result;
  }


}

#include "ugly_macro.hpp"

TMB_ATOMIC_UNARY_FUNCTION(
			  pnorm,
			  Rmath::pnorm5(x,0,1,1,0);,
			  ty[0] = pnorm(tx[0]);,
			  px[0] = dnorm(tx[0],Type(0),Type(1),false) * py[0];
			  )

TMB_ATOMIC_UNARY_FUNCTION(
			  qnorm,
			  Rmath::qnorm5(x,0,1,1,0);,
			  ty[0] = qnorm(tx[0]);,
			  px[0] = Type(1) / dnorm(ty[0],Type(0),Type(1),false) * py[0];
			  )

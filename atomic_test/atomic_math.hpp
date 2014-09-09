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
			  // ATOMIC_NAME
			  pnorm,
			  // ATOMIC_DOUBLE
			  Rmath::pnorm5(x,0,1,1,0);,
			  // ATOMIC_FORWARD
			  ty[0] = pnorm(tx[0]);,
			  // ATOMIC_REVERSE
			  px[0] = dnorm(tx[0],Type(0),Type(1),false) * py[0];
			  )

TMB_ATOMIC_UNARY_FUNCTION(
			  // ATOMIC_NAME
			  qnorm,
			  // ATOMIC_DOUBLE
			  Rmath::qnorm5(x,0,1,1,0);,
			  // ATOMIC_FORWARD
			  ty[0] = qnorm(tx[0]);,
			  // ATOMIC_REVERSE
			  px[0] = Type(1) / dnorm(ty[0],Type(0),Type(1),false) * py[0];
			  )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   D_incpl_gamma_shape
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   vy[0]=Rmath::D_incpl_gamma_shape(vx[0],vx[1],vx[2]);
			   ,
			   // ATOMIC_REVERSE
			   px[0] = exp(-tx[0])*pow(tx[0],tx[1]-Type(1.0))*pow(log(tx[0]),tx[2]) * py[0];
			   CppAD::vector<Type> tx_(tx);
			   tx_[2] = tx_[2] + Type(1.0);  // Add one to get partial wrt. tx[1]
			   px[1] = D_incpl_gamma_shape(tx_)[0] * py[0];
			   px[2] = Type(0);
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   matmul
			   ,
			   // OUTPUT_DIM
			   vx.size()/2
			   ,
			   // ATOMIC_DOUBLE
			   int n=sqrt(vx.size()/2);
			   matrix<double> left(n,n);
			   matrix<double> right(n,n);
			   for(int i=0;i<n*n;i++){left(i)=vx[i];right(i)=vx[i+n*n];}
			   matrix<double> res=left*right; // Use Eigen matrix multiply
			   for(int i=0;i<n*n;i++)vy[i]=res(i);
			   ,
			   // ATOMIC_REVERSE
			   int n=sqrt(ty.size());
			   CppAD::vector<Type> arg1(2*n*n);  // For mat-mult W*Y^T
			   CppAD::vector<Type> arg2(2*n*n);  // For mat-mult X^T*W
			   int k=0;
			   for(int i=0;i<n;i++){
			     for(int j=0;j<n;j++){
			       arg1[k] = py[k];
			       arg1[k+n*n] = tx[n*n + i + j*n];
			       arg2[k] = tx[i + j*n];
			       arg2[k+n*n] = py[k];
			       k++;
			     }
			   }
			   CppAD::vector<Type> res1=matmul(arg1);
			   CppAD::vector<Type> res2=matmul(arg2);
			   for(int i=0;i<n*n;i++){px[i]=res1[i];px[i+n*n]=res2[i];}
			   )


/* Utilities for easier use of matmul symbol */
template<class Type>
matrix<Type> vec2mat(CppAD::vector<Type> x){
  int n=sqrt(x.size());
  matrix<Type> res(n,n);
  for(int i=0;i<n*n;i++)res(i)=x[i];
  return res;
}
template<class Type>
CppAD::vector<Type> mat2vec(matrix<Type> x){
  int n=x.size();
  CppAD::vector<Type> res(n);
  for(int i=0;i<n;i++)res[i]=x(i);
  return res;
}
template<class Type>
matrix<Type> matmul(matrix<Type> x, matrix<Type> y){
  CppAD::vector<Type> arg(x.size()+y.size());
  for(int i=0;i<x.size();i++){arg[i]=x(i);}
  for(int i=0;i<y.size();i++){arg[i+x.size()]=y(i);}
  return vec2mat(matmul(arg));
}

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   matinv
			   ,
			   // OUTPUT_DIM
			   vx.size()
			   ,
			   // ATOMIC_DOUBLE
			   int n=sqrt(vx.size());
			   matrix<double> X(n,n);
			   for(int i=0;i<n*n;i++){X(i)=vx[i];}
			   matrix<double> res=X.inverse();   // Use Eigen matrix inverse (LU)
			   for(int i=0;i<n*n;i++)vy[i]=res(i);
			   ,
			   // ATOMIC_REVERSE  (-f(X)^T*W*f(X)^T)
			   int n=sqrt(ty.size());
			   matrix<Type> W=vec2mat(py);       // Range direction
			   matrix<Type> Y=vec2mat(ty);       // f(X)
			   matrix<Type> Yt=Y.transpose();    // f(X)^T
			   matrix<Type> tmp=matmul(W,Yt);    // W*f(X)^T
			   matrix<Type> res=-matmul(Yt,tmp); // -f(X)^T*W*f(X)^T
			   px=mat2vec(res);
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   logdet
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   matrix<double> X=vec2mat(vx);
			   matrix<double> LU=X.lu().matrixLU();    // Use Eigen LU decomposition
			   vector<double> LUdiag = LU.diagonal();
			   double res=LUdiag.abs().log().sum();    // TODO: currently PD only - take care of sign.
			   vy[0] = res;
			   ,
			   // ATOMIC_REVERSE  (X^-1*W[0])
			   CppAD::vector<Type> invX=matinv(tx);
			   for(int i=0;i<tx.size();i++)px[i]=invX[i]*py[0];
			   )


/* Temporary test of dmvnorm implementation based on atomic symbols.
   Should reduce tape size from O(n^3) to O(n^2).
*/
template<class Type>
Type nldmvnorm(vector<Type> x, matrix<Type> Sigma){
  matrix<Type> Q=vec2mat(matinv(mat2vec(Sigma)));
  Type logdetQ = -logdet(mat2vec(Sigma))[0];
  Type quadform = (x*(Q*x)).sum();
  return -Type(.5)*logdetQ + Type(.5)*quadform + x.size()*Type(log(sqrt(2.0*M_PI)));
}

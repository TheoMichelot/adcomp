#include <TMB.hpp>
#include "atomic_math.hpp"
#include "incpl_gamma.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(a);
  PARAMETER_VECTOR(x);
  Type res=0;
  if(a==0){
    for(int i=0;i<x.size();i++)res+=pnorm(x[i]);
  } else if(a==1){
    for(int i=0;i<x.size();i++)res+=qnorm(x[i]);
  } else if(a==2){
    CppAD::vector<Type> arg(3);
    arg[0] = x[0];
    arg[1] = x[1];
    arg[2] = Type(0);
    res += D_incpl_gamma_shape(arg)[0];
  } else {
    error("Invalid a");
  }
  return res;
}


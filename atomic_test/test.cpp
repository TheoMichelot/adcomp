#include <TMB.hpp>
#include "atomic_math.hpp"

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
  } else {
    error("Invalid a");
  }
  return res;
}


library(TMB)
##compile("test.cpp","-O0 -g -Woverloaded-virtual")
compile("test.cpp")
dyn.load(dynlib("test"))
##config(optimize.instantly=0)
x <- (2:9) / 10

## =============== Test pnorm
obj <- MakeADFun(data=list(a=0),parameters=list(x=x),DLL="test")
obj$fn(x)
sum(pnorm(x)) ## Check
obj$gr(x)
dnorm(x) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=0),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
library(numDeriv)
diag(hessian(obj$env$f,x))

## =============== Test qnorm
obj <- MakeADFun(data=list(a=1),parameters=list(x=x),DLL="test")
obj$fn(x)
sum(qnorm(x)) ## Check
obj$gr(x)
1/(dnorm(qnorm(x))) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=1),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
library(numDeriv)
diag(hessian(obj$env$f,x))

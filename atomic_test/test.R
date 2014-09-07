library(TMB)
library(numDeriv)
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
diag(hessian(obj$env$f,x))

## =============== Test incomplete gamma
incpl_gamma <- function(x)pgamma(x[1],x[2])*gamma(x[2])
obj <- MakeADFun(data=list(a=2),parameters=list(x=x),DLL="test")
obj$fn(x)
incpl_gamma(x[1:2]) ## Check
obj$gr(x)
grad(incpl_gamma,x[1:2]) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=2),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
hessian(incpl_gamma,x[1:2]) ## Check

## =============== Test matrix multiply
f <- function(x){
    n <- length(x)/2
    m1 <- matrix(x[1:n],sqrt(n))
    m2 <- matrix(x[-(1:n)],sqrt(n))
    sum(m1%*%m2)
}
obj <- MakeADFun(data=list(a=3),parameters=list(x=x),DLL="test")
obj$fn(x)
f(x) ## Check
obj$gr(x)
grad(f,x) ## Check


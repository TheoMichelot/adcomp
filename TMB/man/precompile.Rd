% Generated by roxygen2 (4.0.1): do not edit by hand
\name{precompile}
\alias{precompile}
\title{Precompile the TMB library in order to speed up compilation of templates.}
\usage{
precompile(...)
}
\arguments{
\item{...}{Passed to \code{compile}.}
}
\description{
Precompile the TMB library
}
\details{
The precompilation should only be run once, typically right after installation of TMB.
Note that the precompilation requires write access to the TMB package folder.
Two versions of the library - with/without the openmp flag - will be generated. After this,
compilation times of templates should be reduced.
\itemize{
\item To precompile on Linux run \code{precompile()}.
\item To precompile on OS X run \code{precompile(PKG_LIBS = "-install_name `pwd`/$@")}.
}
Note that precompilation has side effects: It is not possible to work with more than one
model at a time for a single R instance.
}


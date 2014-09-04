/* Namespace with double versions of R special math library */
namespace Rmath {
  #include <Rmath.h>
  // Macros do not respect the namespace limits.
  // Some of them will conflict with other TMB functions.
  #undef dnorm
  #undef pnorm
  #undef qnorm
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

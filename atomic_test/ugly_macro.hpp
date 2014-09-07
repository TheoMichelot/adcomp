// Example:
// ========
// ATOMIC_NAME:      pnorm
// ATOMIC_DOUBLE:    Rmath::pnorm5(x,0,1,1,0);
// ATOMIC_FORWARD:   ty[0] = pnorm(tx[0]);
// ATOMIC_REVERSE:   px[0] = dnorm(tx[0],Type(0),Type(1),false) * py[0];

#define TMB_ATOMIC_UNARY_FUNCTION(ATOMIC_NAME,ATOMIC_DOUBLE,ATOMIC_FORWARD,ATOMIC_REVERSE) \
double ATOMIC_NAME(double x);						\
AD<double> ATOMIC_NAME(AD<double> x);					\
AD<AD<double> > ATOMIC_NAME(AD<AD<double> > x);				\
AD<AD<AD<double> > > ATOMIC_NAME(AD<AD<AD<double> > > x);		\
template <class Type>							\
class atomic_##ATOMIC_NAME : public CppAD::atomic_base<Type> {		\
public:									\
  atomic_##ATOMIC_NAME(const char* name) : CppAD::atomic_base<Type>(name){ \
    std::cout << "Constructing atomic " << #ATOMIC_NAME << "\n" ;	\
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);		\
  }									\
private:								\
  virtual bool forward(size_t p,					\
		       size_t q,					\
		       const CppAD::vector<bool>& vx,			\
		       CppAD::vector<bool>& vy,				\
		       const CppAD::vector<Type>& tx,			\
		       CppAD::vector<Type>& ty				\
		       )						\
  {									\
    if( vx.size() > 0 )vy[0] = vx[0];					\
    ATOMIC_FORWARD;							\
    return true;							\
  }									\
  virtual bool reverse(size_t q,					\
		       const CppAD::vector<Type>& tx,			\
		       const CppAD::vector<Type>& ty,			\
		       CppAD::vector<Type>& px,				\
		       const CppAD::vector<Type>& py			\
		       )						\
  {									\
    ATOMIC_REVERSE;							\
    return true;							\
  }									\
  virtual bool rev_sparse_jac(size_t q,					\
  			      const CppAD::vector<bool>& rt,		\
  			      CppAD::vector<bool>& st)			\
  {									\
    for(size_t i=0;i<st.size();i++)st[i]=true;				\
    return true;							\
  }									\
  virtual bool rev_sparse_jac(size_t q,					\
			      const CppAD::vector< std::set<size_t> >& rt, \
			      CppAD::vector< std::set<size_t> >& st)	\
  {									\
    error("Should not be called");					\
  }									\
};									\
double ATOMIC_NAME(double x){						\
  return ATOMIC_DOUBLE;							\
}									\
atomic_##ATOMIC_NAME<double> afun1##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
AD<double> ATOMIC_NAME(AD<double> x){					\
  CppAD::vector<AD<double> > vx(1);					\
  CppAD::vector<AD<double> > vy(1);					\
  vx[0]=x;								\
  afun1##ATOMIC_NAME(vx,vy);						\
  return vy[0];								\
}									\
atomic_##ATOMIC_NAME<AD<double> > afun2##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
AD<AD<double> > ATOMIC_NAME(AD<AD<double> > x){				\
  CppAD::vector<AD<AD<double> > > vx(1);				\
  CppAD::vector<AD<AD<double> > > vy(1);				\
  vx[0]=x;								\
  afun2##ATOMIC_NAME(vx,vy);						\
  return vy[0];								\
}									\
atomic_##ATOMIC_NAME<AD<AD<double> > > afun3##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
AD<AD<AD<double> > > ATOMIC_NAME(AD<AD<AD<double> > > x){		\
  CppAD::vector<AD<AD<AD<double> > > > vx(1);				\
  CppAD::vector<AD<AD<AD<double> > > > vy(1);				\
  vx[0]=x;								\
  afun3##ATOMIC_NAME(vx,vy);						\
  return vy[0];								\
}


#define TMB_ATOMIC_VECTOR_FUNCTION(ATOMIC_NAME,OUTPUT_DIM,ATOMIC_DOUBLE,ATOMIC_REVERSE) \
CppAD::vector<double> ATOMIC_NAME(CppAD::vector<double> x);		\
CppAD::vector<AD<double > > ATOMIC_NAME(CppAD::vector<AD<double> > x);	\
CppAD::vector<AD<AD<double> > > ATOMIC_NAME(CppAD::vector<AD<AD<double> > > x);	\
CppAD::vector<AD<AD<AD<double> > > > ATOMIC_NAME(CppAD::vector<AD<AD<AD<double> > > > x); \
template <class Type>							\
class atomic##ATOMIC_NAME : public CppAD::atomic_base<Type> {		\
public:									\
  atomic##ATOMIC_NAME(const char* name) : CppAD::atomic_base<Type>(name){ \
    std::cout << "Constructing atomic " << #ATOMIC_NAME << "\n" ;	\
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);		\
  }									\
private:								\
  virtual bool forward(size_t p,					\
		       size_t q,					\
		       const CppAD::vector<bool>& vx,			\
		       CppAD::vector<bool>& vy,				\
		       const CppAD::vector<Type>& tx,			\
		       CppAD::vector<Type>& ty				\
		       )						\
  {									\
    if( vx.size() > 0 ){						\
      bool anyvx = false;						\
      for(int i=0;i<vx.size();i++)anyvx |= vx[i];			\
      for(int i=0;i<vy.size();i++)vy[i] = anyvx;			\
    }									\
    ty = ATOMIC_NAME(tx);	       					\
    return true;							\
  }									\
  virtual bool reverse(size_t q,					\
		       const CppAD::vector<Type>& tx,			\
		       const CppAD::vector<Type>& ty,			\
		       CppAD::vector<Type>& px,				\
		       const CppAD::vector<Type>& py			\
		       )						\
  {									\
    ATOMIC_REVERSE;							\
    return true;							\
  }									\
  virtual bool rev_sparse_jac(size_t q,					\
			      const CppAD::vector<bool>& rt,		\
			      CppAD::vector<bool>& st)			\
  {									\
    for(size_t i=0;i<st.size();i++)st[i]=true;				\
    return true;							\
  }									\
  virtual bool rev_sparse_jac(size_t q,					\
			      const CppAD::vector< std::set<size_t> >& rt, \
			      CppAD::vector< std::set<size_t> >& st)	\
  {									\
    error("Should not be called");					\
  }									\
};									\
CppAD::vector<double> ATOMIC_NAME(CppAD::vector<double> vx){		\
  CppAD::vector<double> vy(OUTPUT_DIM);					\
  ATOMIC_DOUBLE;							\
  return vy;								\
}									\
atomic##ATOMIC_NAME<double> afun1##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
CppAD::vector<AD<double > > ATOMIC_NAME(CppAD::vector<AD<double> > vx){	\
  CppAD::vector<AD<double> > vy(OUTPUT_DIM);				\
  afun1##ATOMIC_NAME(vx,vy);						\
  return vy;								\
}									\
atomic##ATOMIC_NAME<AD<double> > afun2##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
CppAD::vector<AD<AD<double> > > ATOMIC_NAME(CppAD::vector<AD<AD<double> > > vx){ \
  CppAD::vector<AD<AD<double> > > vy(OUTPUT_DIM);			\
  afun2##ATOMIC_NAME(vx,vy);						\
  return vy;								\
}									\
atomic##ATOMIC_NAME<AD<AD<double> > > afun3##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
CppAD::vector<AD<AD<AD<double> > > > ATOMIC_NAME(CppAD::vector<AD<AD<AD<double> > > > vx){ \
  CppAD::vector<AD<AD<AD<double> > > > vy(OUTPUT_DIM);			\
  afun3##ATOMIC_NAME(vx,vy);						\
  return vy;								\
}



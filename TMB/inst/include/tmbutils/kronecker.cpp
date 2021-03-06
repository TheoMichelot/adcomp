/** \file
    \brief Kronecker product of two matrices
*/

/** \brief Kronecker product of two matrices */
template <class scalartype, int n1, int n2, int n3, int n4>
Matrix<scalartype,n1*n3,n2*n4> kronecker(Matrix<scalartype,n1,n2> x, Matrix<scalartype,n3,n4> y){
  Matrix<scalartype,n1*n3,n2*n4> ans;
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      for(int k=0;k<n3;k++)
	for(int l=0;l<n4;l++)
	  ans(i*n3+k,j*n4+l)=x(i,j)*y(k,l);
  return ans;
}

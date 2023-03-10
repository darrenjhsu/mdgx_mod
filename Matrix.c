#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include "mdgxVector.h"
#include "Macros.h"

//-----------------------------------------------------------------------------
// CreateImat: create an M x N integer matrix, initialized to zero.     
//-----------------------------------------------------------------------------
imat CreateImat(int M, int N)
{
  int i;
  imat A;

  A.row = M;
  A.col = N;
  A.map = (int**)malloc(M*sizeof(int*));
  A.data = (int*)calloc(M*N, sizeof(int));
  for (i = 0; i < M; i++) {
    A.map[i] = &A.data[N*i];
  }

  return A;
}

//-----------------------------------------------------------------------------
// DestroyImat: destroy an integer matrix.                              
//-----------------------------------------------------------------------------
void DestroyImat(imat *A)
{
  free(A->data);
  free(A->map);
}

//-----------------------------------------------------------------------------
// ReallocImat: re-allocate an integer matrix.  Data in the original matrix
//              will remain "in place" in that indices of the new matrix,
//              so long as they existed in the original matrix, will contain 
//              the same data as before.  New indices will contain zeros.
//
// Arguments:                                                           
//   A:     the original matrix                                         
//   M:     the new number of rows                                      
//   N:     the new number of columns                                   
//-----------------------------------------------------------------------------
imat ReallocImat(imat *A, int M, int N)
{
  int i, j, minM, minN;
  imat Ap;

  Ap = CreateImat(M, N);
  minM = (M < A->row) ? M : A->row;
  minN = (N < A->col) ? N : A->col;
  for (i = 0; i < minM; i++) {
    for (j = 0; j < minN; j++) {
      Ap.map[i][j] = A->map[i][j];
    }
  }
  DestroyImat(A);

  return Ap;
}

//-----------------------------------------------------------------------------
// CreateDmat: create an M x N double-precision real matrix, initialized to
//             zero.
//
// Arguments:                                                           
//   M:        the number of rows                                       
//   N:        the number of columns                                    
//   prepFFT:  fset to 1 to activate preparation for FFTs, 0 otherwise  
//-----------------------------------------------------------------------------
dmat CreateDmat(int M, int N, int prepFFT)
{
  int i, Np;
  dmat A;

  // Prepare for FFTs
  A.pfft = prepFFT;
  Np = (prepFFT == 1) ? 2*(N/2 + 1) : N;

  // Allocate memory
  A.data = (double*)calloc(M*Np, sizeof(double));

  // Make the map
  A.map = (double**)malloc(M*sizeof(double*));
  for (i = 0; i < M; i++) {
    A.map[i] = &A.data[i*Np];
  }
  
  // Set up complex maps if FFT preparation is required
#if 0
  if (prepFFT == 1) {
    A.fdata = (fftw_complex*)A.data;
    A.fmap = (fftw_complex**)malloc(M*sizeof(fftw_complex*));
    for (i = 0; i < M; i++) {
      A.fmap[i] = &A.fdata[(i*Np)/2];
    }
  }
#endif
  
  // Store the dimensionality
  A.row = M;
  A.col = N;

  return A;
}

//-----------------------------------------------------------------------------
// DestroyDmat: destroy a double-precision real matrix.                 
//-----------------------------------------------------------------------------
void DestroyDmat(dmat *A)
{
  free(A->data);
  free(A->map);
}

//-----------------------------------------------------------------------------
// ReallocDmat: re-allocate a double-precision real matrix.  Data in the
//              original matrix will remain "in place" in that indices  
//              of the new matrix, so long as they existed in the original
//              matrix, will contain the same data as before.  New indices
//              will contain zeros.
//
// Arguments:                                                           
//   A:     the original matrix                                         
//   M:     the new number of rows                                      
//   N:     the new number of columns                                   
//-----------------------------------------------------------------------------
dmat ReallocDmat(dmat *A, int M, int N)
{
  int i, j, minM, minN;
  double *odtmp, *ndtmp;
  dmat Ap;

  Ap = CreateDmat(M, N, A->pfft);
  minM = (M < A->row) ? M : A->row;
  minN = (N < A->col) ? N : A->col;
  for (i = 0; i < minM; i++) {
    odtmp = A->map[i];
    ndtmp = Ap.map[i];
    for (j = 0; j < minN; j++) {
      ndtmp[j] = odtmp[j];
    }
  }
  DestroyDmat(A);

  return Ap;
}

//-----------------------------------------------------------------------------
// CopyDmat: make a copy of a double-precision matrix structure.  The matrix
//           may be a Fourier-space representation of a real matrix, in which
//           case the Fourier-space representation will be copied completely.
//
// Arguments:                                                           
//   Ac:    the matrix copy                                             
//   A:     the original matrix                                         
//   Acex:  flag to specify that memory for Ac is already allocated     
//-----------------------------------------------------------------------------
void CopyDmat(dmat *Ac, dmat *A, int Acex)
{
  int i, j;
  double *dA, *dAc;

  if (Acex == 0) {
    *Ac = CreateDmat(A->row, A->col, A->pfft);
  }
  const int jlim = (A->pfft == 1) ? 2*(A->col/2 + 1) : A->col;
  for (i = 0; i < A->row; i++) {
    dA = A->map[i];
    dAc = Ac->map[i];
    for (j = 0; j < jlim; j++) {
      dAc[j] = dA[j];
    }
  }
}

//-----------------------------------------------------------------------------
// AplusBmat: add double precision matrix B to matrix A, store the result in
//            A.  If the pom variable is set to any value other than 1, B will
//            be subtracted from A instead.        
//-----------------------------------------------------------------------------
void AplusBmat(dmat *A, dmat *B, int pom)
{
  int i, j, ilim, jlim;
  double *atmp, *btmp;

  // Check
  if (A->row != B->row) {
    printf("AplusBmat >> Error.  Rows of A and B do not match (%d, %d).\n",
	   A->row, B->row);
  }
  if (A->col != B->col) {
    printf("AplusBmat >> Error.  Columns of A and B do not match (%d, %d).\n",
	   A->col, B->col);
  }
  ilim = A->row;
  jlim = (A->pfft == 1 && B->pfft == 1) ? 2*(A->col/2+1) : A->col;
  for (i = 0; i < ilim; i++) {
    atmp = A->map[i];
    btmp = B->map[i];
    if (pom == 1) {
      for (j = 0; j < jlim; j++) {
	atmp[j] += btmp[j];
      }
    }
    else {
      for (j = 0; j < jlim; j++) {
	atmp[j] -= btmp[j];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// AddDmat: add two double-precision real matrices, with scaling factors in
//          front of each of them.                                   
//                                                                      
// Arguments:                                                           
//   A:     the first matrix                                            
//   B:     the second matrix                                           
//   Ascl:  the first matrix's scaling factor                           
//   Bscl:  the second matrix's scaling factor                          
//-----------------------------------------------------------------------------
void AddDmat(dmat *A, dmat *B, double Ascl, double Bscl)
{
  int i, j, ilim, jlim;
  double *atmp, *btmp;

  // Check
  if (A->row != B->row) {
    printf("AplusBmat >> Error.  Rows of A and B do not match (%d, %d).\n",
	   A->row, B->row);
  }
  if (A->col != B->col) {
    printf("AplusBmat >> Error.  Columns of A and B do not match (%d, %d).\n",
	   A->col, B->col);
  }
  ilim = A->row;
  jlim = (A->pfft == 1 && B->pfft == 1) ? 2*(A->col/2+1) : A->col;

  for (i = 0; i < ilim; i++) {
    atmp = A->map[i];
    btmp = B->map[i];
    for (j = 0; j < jlim; j++) {
      atmp[j] = Ascl*atmp[j] + Bscl*btmp[j];
    }
  }
}

//-----------------------------------------------------------------------------
// CreateCmat: create an M x N character matrix, initialized to null.   
//-----------------------------------------------------------------------------
cmat CreateCmat(int M, int N)
{
  int i;
  cmat A;

  A.row = M;
  A.col = N;
  A.map = (char**)malloc(M*sizeof(char*));
  A.data = (char*)calloc(M*N, sizeof(char));
  for (i = 0; i < M; i++) {
    A.map[i] = &A.data[N*i];
  }

  return A;
}

//-----------------------------------------------------------------------------
// DestroyCmat: destroy a character matrix.                             
//-----------------------------------------------------------------------------
void DestroyCmat(cmat *A)
{
  free(A->data);
  free(A->map);
}

//-----------------------------------------------------------------------------
// ReallocCmat: extend a character matrix to hold more rows or columns. 
//-----------------------------------------------------------------------------
cmat ReallocCmat(cmat *A, int M, int N)
{
  int i, j, minM, minN;
  char *rA, *rnA;
  cmat nA;

  nA = CreateCmat(M, N);
  minM = (M < A->row) ? M : A->row;
  minN = (N < A->col) ? N : A->col;
  for (i = 0; i < minM; i++) {
    rA = A->map[i];
    rnA = nA.map[i];
    for (j = 0; j < minN; j++) {
      rnA[j] = rA[j];
    }
  }
  DestroyCmat(A);

  return nA; 
}

//-----------------------------------------------------------------------------
// DMatMult: multiply two double-precision real matrices A and B to form
//           another matrix C.                                          
//-----------------------------------------------------------------------------
void DMatMult(dmat A, dmat B, dmat C)
{
  int i, j, k;
  double tempval;
  double* tmp;
  double* tm2p;

  if (A.col != B.row) {
    printf("DMatMult >> Error: internal dimensions do not agree.\n"
           "DMatMult >> Taking minimum value.\n\n");
    if (A.col > B.row) {
      A.col = B.row;
    }
    else {
      B.row = A.col;
    }
  }

  for (i = 0; i < A.row; i++) {
    tmp = A.map[i];
    tm2p = C.map[i];
    for (j = 0; j < B.col; j++) {
      tempval = 0.0;
      for (k = 0; k < A.col; k++) {
	tempval += tmp[k]*B.map[k][j];
      }
      tm2p[j] = tempval;
    }
  }
}

//-----------------------------------------------------------------------------
// DMatVecMult: multiply matrix A times vector x to produce result b.  Assumes
//              that x has at least as many elements as A has columns, and
//              that b has at least as many elements as A has rows.
//-----------------------------------------------------------------------------
void DMatVecMult(dmat *A, double* x, double* b)
{
  int i, j;
  double res;
  double *dtmp;

  for (i = 0; i < A->row; i++) {
    res = 0.0;
    dtmp = A->map[i];
    for (j = 0; j < A->col; j++) {
      res += dtmp[j]*x[j];
    }
    b[i] = res;
  }
}

//-----------------------------------------------------------------------------
// IMatMult: multiply two integer matrices A and B to form matrix C.    
//-----------------------------------------------------------------------------
void IMatMult(imat A, imat B, imat C)
{
  int i, j, k;
  int tempval;
  int* tmp;
  int* tm2p;

  if (A.col != B.row) {
    printf("IMatMult >> Error: internal dimensions do not agree.\n"
           "IMatMult >> Taking minimum value.\n\n");
    if (A.col > B.row) {
      A.col = B.row;
    }
    else {
      B.row = A.col;
    }
  }

  for (i = 0; i < A.row; i++) {
    tmp = A.map[i];
    tm2p = C.map[i];
    for (j = 0; j < B.col; j++) {
      tempval = 0.0;
      for (k = 0; k < A.col; k++) {
	tempval += tmp[k]*B.map[k][j];
      }
      tm2p[j] = tempval;
    }
  }
}

//-----------------------------------------------------------------------------
// ttInvRow: accessory to ttInv below, makes the rows of the inverse.   
//-----------------------------------------------------------------------------
static void ttInvRow(dmat invA, double* pq, double* rs, int j, double deta)
{
  invA.map[0][j] = deta*(pq[1]*rs[2]-rs[1]*pq[2]);
  invA.map[1][j] = deta*(pq[0]*rs[2]-rs[0]*pq[2]);
  invA.map[2][j] = deta*(pq[0]*rs[1]-rs[0]*pq[1]);
}

//-----------------------------------------------------------------------------
// ttInv: inverts a 3 x 3 matrix A and returns the result as invA.  The 
//        algorithm is HIGHLY INEFFICIENT (scales as N!) but for the    
//        3 x 3 case (3! = 6) we don't care.                            
//-----------------------------------------------------------------------------
void ttInv(dmat A, dmat invA)
{
  double deta;
  double* pq;
  double* rs;

  // Compute the determinant of a
  deta = A.map[0][0]*(A.map[1][1]*A.map[2][2] - A.map[1][2]*A.map[2][1]) +
         A.map[0][1]*(A.map[1][2]*A.map[2][0] - A.map[1][0]*A.map[2][2]) +
         A.map[0][2]*(A.map[1][0]*A.map[2][1] - A.map[1][1]*A.map[2][0]);
  deta = 1.0/deta;

  pq = A.map[1];
  rs = A.map[2];
  ttInvRow(invA, pq, rs, 0, deta);
  pq = A.map[0];
  ttInvRow(invA, pq, rs, 1, deta);
  rs = A.map[1];
  ttInvRow(invA, pq, rs, 2, deta);
  invA.map[0][1] *= -1.0;
  invA.map[1][0] *= -1.0;
  invA.map[1][2] *= -1.0;
  invA.map[2][1] *= -1.0;
}

//-----------------------------------------------------------------------------
// DMatAdd: computes the expression da*A + db*B = C, where da and db are
//          scalars and A, B, and C are matrices.  Note that C may point to 
//          addresses of either A or B so that the result will overwrite A
//          or B.                                           
//-----------------------------------------------------------------------------
void DMatAdd(dmat A, dmat B, double da, double db, dmat C)
{
  int i;

  // Check
  if (A.row != B.row || A.col != B.col || A.row != C.row || A.col != C.col) {
    printf("DMatAdd >> Error.  Matrix dimensions do not agree.\n"
	   "DMatAdd >> A = [ %6d x %6d ]\nDMatAdd >> B = [ %6d x %6d ]\n"
	   "DMatAdd >> C = [ %6d x %6d ]\n", A.row, A.col, B.row, B.col,
	   C.row, C.col);
    exit(1);
  }

  for (i = 0; i < A.row*A.col; i++) {
    C.data[i] = da*A.data[i] + db*B.data[i];
  }
}

//-----------------------------------------------------------------------------
// RotationMatrix: constructs a rotation matrix for rotating about a    
//                 particular vector.                                   
//-----------------------------------------------------------------------------
void RotationMatrix(dmat *mat, double* vec, double angle)
{
  double cT = cos(angle), cT1 = 1-cT;
  double sT = sin(angle);
  double x = vec[0], y = vec[1], z = vec[2];

  // Generate the rotation matrix
  mat->map[0][0] = 1 + cT1 * (x*x - 1);
  mat->map[0][1] = -z*sT + cT1*x*y;
  mat->map[0][2] = y*sT + cT1*x*z;
  mat->map[1][0] = z*sT + cT1*x*y;
  mat->map[1][1] = 1 + cT1*(y*y - 1);
  mat->map[1][2] = -x*sT + cT1*y*z;
  mat->map[2][0] = -y*sT + cT1*x*z;
  mat->map[2][1] = x*sT + cT1*y*z;
  mat->map[2][2] = 1 + cT1*(z*z - 1);
}

//-----------------------------------------------------------------------------
// AxbQRRxc: function for solving linear least-squares problems.  This performs
//           the first stage of the solution by taking a matrix problem Ax = b,
//           where A { m by n, m >= n, and decomposes A by the "modern
//           classical" QR algorithm to recast the problem as Rx = c, where R
//           is the the upper-triangular component of A and the vector c is
//           implicitly computed as (Q*b).  This algorithm may be found in
//           [REF]:  
//                                                                      
//           Trefethen, Lloyd N. and Bau, David III. "Numerical Linear  
//           Algebra." pp.73.  Society for Industrial and Applied Mathematics,
//           Philadelphia.  1997.                          
//-----------------------------------------------------------------------------
void AxbQRRxc(dmat A, double* b, int update_user)
{
  int i, j, k;
  double tnm_v, tnm_v2, tempval, sign_v;
  double* v;
  double* vprime;
  double* tmp;

  if (A.row >= A.col) {
    v = (double*)malloc(A.row*sizeof(double));
    vprime = (double*)malloc(A.row*sizeof(double));
    for (k = 0; k < A.col; k++) {

      // Update the user
      if (update_user == 1) {
        fprintf(stderr, "\rmdgx >> Solving LLSP for column %5d of %5d", k+1,
		A.col);
        fflush(stderr);
      }

      // Compute the kth column of Q*
      tnm_v2 = 0.0;
      for (i = 0; i < A.row-k; i++) {
        v[i] = A.map[i+k][k];
        tnm_v2 += v[i]*v[i];
      }
      sign_v = SIGN(v[0]);
      tnm_v = sqrt(tnm_v2);
      tnm_v2 -= v[0]*v[0];
      v[0] += sign_v*tnm_v;
      tnm_v = 1.0/sqrt(tnm_v2 + v[0]*v[0]);
      for (i = 0; i < A.row-k; i++) {
        v[i] = v[i]*tnm_v;
      }

      // Update A as R evolves
      for (i = 0; i < A.col-k; i++) {
	vprime[i] = 0.0;
      }
      for (i = 0; i < A.row-k; i++) {
	tmp = &A.map[i+k][k];
        tempval = v[i];
	for (j = 0; j < A.col-k; j++) {
	  vprime[j] += tempval*tmp[j];
	}
      }
      for (i = 0; i < A.row-k; i++) {
        tmp = &A.map[i+k][k];
        tempval = 2.0*v[i];
        for (j = 0; j < A.col-k; j++) {
          tmp[j] -= tempval*vprime[j];
        }
      }

      // Update b as Q* evolves
      tmp = &b[k];
      tempval = 2.0*DotP(v, tmp, A.row-k);
      for (i = 0; i < A.row-k; i++) {
        tmp[i] -= tempval*v[i];
      }
    }

    // Free Allocated Memory
    free(v);
    free(vprime);

    if (update_user == 1) {
      printf("\n");
    }
  }
  else {
    printf("mdgx >> Error: Matrix A is rank deficient.\n");
  }
}

//-----------------------------------------------------------------------------
// BackSub: solve the equation Rx = b, where R is an upper triangular matrix of
//          dimension n, and b is a vector of dimension n.  Results are
//          returned in the vector b.                       
//-----------------------------------------------------------------------------
void BackSub(dmat R, double* b)
{
  int i, j;
  double multval, temp_b, pivot;

  for (i = R.col-1; i > 0; i--) {
    pivot = 1.0/R.map[i][i];
    temp_b = b[i];
    for (j = i-1; j >= 0; j--) {
      multval = R.map[j][i]*pivot;
      b[j] -= multval*temp_b;
    }
    b[i] *= pivot;
  }
  b[0] /= R.map[0][0];
}

//-----------------------------------------------------------------------------
// TRED2: reduces a symmetric matrix to tridiagonal form.  It is a modification
//        more suited for C programs from [REF]:           
//                                                                      
//        William H. Press, Saul A. Teukolsky, William T. Vetterling,   
//        and Brian P. Flannery.  Numerical Recipes in C, Second Edition. 
//        Cambridge University Press, 1992.                   
//-----------------------------------------------------------------------------
void TRED2(double** A, int n, double* d, double* e)
{
  int i, j, k, l;
  double scale, hh, h, g, f;
  double* tmp;
  double* tm2p;

  for (i = n-1; i >= 1; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    tmp = A[i];
    if (l > 0) {
      for (k = 0; k <= l; k++) {
	scale += fabs(tmp[k]);
      }
      if (scale == 0.0) {
	e[i] = tmp[l];
      }
      else {
	for (k = 0; k <= l; k++) {
	  tmp[k] /= scale;
	  h += tmp[k]*tmp[k];
	}
	f = tmp[l];
	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i] = scale*g;
	h -= f*g;
	tmp[l] = f - g;
	f = 0.0;
	for (j = 0; j <= l; j++) {
	  tm2p = A[j];
	  tm2p[i] = tmp[j]/h;
	  g = 0.0;
	  for (k = 0; k <= j; k++) {
	    g += tm2p[k]*tmp[k];
	  }
	  for (k = j+1; k <=l; k++) {
	    g += A[k][j]*tmp[k];
	  }
	  e[j] = g/h;
	  f += e[j]*tmp[j];
	}
	hh = f/(h + h);
	for (j = 0; j <= l; j++) {
	  f = tmp[j];
	  e[j] = g = e[j] - hh*f;
	  tm2p = A[j];
	  for (k = 0; k <= j; k++) {
	    tm2p[k] -= (f*e[k] + g*tmp[k]);
	  }
	}
      }
    }
    else {
      e[i] = tmp[l];
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    tmp = A[i];
    l = i - 1;
    if (d[i]) {
      for (j = 0; j <= l; j++) {
        g = 0.0;
        for (k = 0; k <= l; k++) {
          g += tmp[k]*A[k][j];
        }
        for (k = 0; k <= l; k++) {
          A[k][j] -= g*A[k][i];
        }
      }
    }
    d[i] = tmp[i];
    tmp[i] = 1.0;
    for (j = 0; j <= l; j++) {
      A[j][i] = tmp[j] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// TQLI: this function is also modified from [REF]:                     
//                                                                      
//        William H. Press, Saul A. Teukolsky, William T. Vetterling,   
//        and Brian P. Flannery.  Numerical Recipes in C, Second        
//        Edition.  Cambridge University Press, 1992.                   
//                                                                      
//        The input is the super- and sub- diagonal vectors d and e.    
//        The eigenvalues are returned in d and the eigenvectors are    
//        returned in Z.                                                
//-----------------------------------------------------------------------------
void TQLI(double* d, double* e, int n, double** z)
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  double* tmp;

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;

    // Added to initialize m
    m = l - 1;

    while (m != l) {
      for (m = l; m < n - 1; m++) {
	dd = fabs(d[m]) + fabs(d[m+1]);
	if (fabs(e[m]+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  printf("TQLI >> Error: Too many iterations in tqli\n");
	}
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = pythag(g,1.0);
	g = d[m]-d[l]+e[l]/(g+SIGN2(r,g));
	c = 1.0;
	s = 1.0;
	p = 0.0;
	for (i = m-1; i >= l; i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1]-p;
	  r = (d[i] - g)*s + 2.0*c*b;
	  d[i+1] = g + (p = s*r);
	  g = c*r - b;
	  for (k = 0; k < n; k++) {
	    tmp = z[k];
	    f = tmp[i+1];
	    tmp[i+1] = s*tmp[i] + c*f;
	    tmp[i] = c*tmp[i] - s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// DetectInCmat: detect a character string from within a character matrix.
//               Returns 1 if found, 0 if not.
//
// Arguments:
//   C:       the character matrix to search
//   s:       the string to search for
//   row:     the row of C at which the string is found
//   col:     the column of C at which the string ENDS
//   initl:   the first line at which to start reading (this can greatly
//            accelerate the search) 
//   direc:   direction to search starting at initl, 0 = forward, 1 = reverse
//-----------------------------------------------------------------------------
int DetectInCmat(cmat *C, char* s, int *row, int *col, int initl, int direc)
{
  int i, j, k, jlim, slen, pass;
  char *ctmp;

  slen = strlen(s);
  jlim = C->col-slen;
  if (direc == 0) {
    for (i = initl; i < C->row; i++) {
      ctmp = C->map[i];
      for (j = 0; j < jlim; j++) {
        pass = 1;
        for (k = 0; k < slen; k++) {
	  if (ctmp[j+k] != s[k]) {
	    pass = 0;
	    break;
	  }
        }
        if (pass == 1) {
	  *row = i;
	  *col = j + slen;
          return 1;
	}
      }
    }
  }
  else {
    for (i = initl; i >= 0; i--) {
      ctmp = C->map[i];
      for (j = 0; j < jlim; j++) {
        pass = 1;
        for (k = 0; k < slen; k++) {
	  if (ctmp[j+k] != s[k]) {
	    pass = 0;
	    break;
	  }
        }
        if (pass == 1) {
	  *row = i;
	  *col = j + slen;
          return 1;
	}
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// CompareCmat: compare two character matrices, to the extent that each line
//              contains a word (two matrices will be considered identical if
//              each of their rows contains identical characters up to the 
//              first instance of '\0' on that row, regardless of what
//              characters may come after the '\0' and how may columns there
//              are.  Returns 1 if the matrices are not identical, 0 if they
//              are.
// 
// Arguments:
//   C[1,2]     : the character matrices to compare
//-----------------------------------------------------------------------------
int CompareCmat(cmat *C1, cmat *C2)
{
  int i, j, mincol;
  char *ctmp1, *ctmp2;

  if (C1->row != C2->row) {
    return 1;
  }
  mincol = (C1->col < C2->col) ? C1->col : C2->col;
  for (i = 0; i < C1->row; i++) {
    ctmp1 = C1->map[i];
    ctmp2 = C2->map[i];
    j = 0;
    while (j < mincol && ctmp1[j] != '\0' && ctmp2[j] != '\0') {
      if (ctmp1[j] != ctmp2[j]) {
	return 1;
      }
      j++;
    }
  }

  // Both matrices are identical in the ways that count
  return 0;
}

//-----------------------------------------------------------------------------
// DistinguishImat: compare all rows of two integer matrices, returning 0 if
//                  they are identical, 1 if each contains unique rows, 2 if
//                  only the first contains unique rows (the second is a subset
//                  of the first), and 3 if only the second contains unique
//                  rows not found in the first.
//
// Arguments:
//   I[1,2]     : the integer matrices to compare
//-----------------------------------------------------------------------------
int DistinguishImat(imat *I1, imat *I2)
{
  int i, j, k, unique1, unique2, match;
  int *itmp, *jtmp;

  if (I1->col != I2->col) {
    return 1;
  }
  unique1 = I1->row;
  unique2 = I2->row;
  for (i = 0; i < I1->row; i++) {
    itmp = I1->map[i];
    for (j = 0; j < I2->row; j++) {
      match = 1;
      jtmp = I2->map[j];
      for (k = 0; k < I1->col; k++) {
	if (itmp[k] != jtmp[k]) {
	  match = 0;
	  break;
	}
      }
      if (match == 1) {
	unique1--;
      }
    }
  }
  for (i = 0; i < I2->row; i++) {
    itmp = I2->map[i];
    for (j = 0; j < I1->row; j++) {
      match = 1;
      jtmp = I1->map[j];
      for (k = 0; k < I1->col; k++) {
        if (itmp[k] != jtmp[k]) {
          match = 0;
          break;
        }
      }
      if (match == 1) {
        unique2--;
      }
    }
  }

  // Report the result
  if (unique1 > 0 && unique2 > 0) {
    return 1;
  }
  else if (unique1 > 0) {
    return 2;
  }
  else if (unique2 > 0) {
    return 3;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// RotateA2B: produce the rotation matrix to rotate vector A into vector B.
//            This routine will make normalized copies of the two vectors, so
//            no special preparation is needed.
//
// Arguments:
//   A, B:    the vectors to consider--three elements each
//-----------------------------------------------------------------------------
dmat RotateA2B(double* A, double* B)
{
  int i;
  double ma, mb, dotab, mcrab, mbarej;
  double va[3], vb[3], crab[3], barej[3];
  dmat G, F, invF, U;
  
  ma = 1.0/sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
  mb = 1.0/sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  dotab = 0.0;
  for (i = 0; i < 3; i++) {
    va[i] = A[i]*ma;
    vb[i] = B[i]*mb;
    dotab += va[i]*vb[i];
  }
  CrossP(va, vb, crab);
  mcrab = sqrt(crab[0]*crab[0] + crab[1]*crab[1] + crab[2]*crab[2]);
  G = CreateDmat(3, 3, 0);
  invF = CreateDmat(3, 3, 0);
  F = CreateDmat(3, 3, 0);
  U = CreateDmat(3, 3, 0);
  G.map[0][0] = dotab;
  G.map[1][1] = dotab;
  G.map[0][1] = -mcrab;
  G.map[1][0] = mcrab;
  G.map[2][2] = 1.0;
  mbarej = 0.0;
  for (i = 0; i < 3; i++) {
    invF.map[i][0] = va[i];
    barej[i] = vb[i] - dotab*va[i];
    mbarej += barej[i]*barej[i];
  }
  mbarej = 1.0/sqrt(mbarej);
  for (i = 0; i < 3; i++) {
    invF.map[i][1] = barej[i]*mbarej;
    invF.map[i][2] = -crab[i];
  }
  ttInv(invF, F);
  DMatMult(G, F, U);
  DMatMult(invF, U, G);
  ReflectDVec(U.data, G.data, 9);

  // Free allocated memory
  DestroyDmat(&G);
  DestroyDmat(&invF);
  DestroyDmat(&F);
  
  return U;
}

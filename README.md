# Math
Mainly linear algebra function by various computer languages
---

## Complex Eigenvalue Problem
### Single QR method

1) Extraction of Eigenvalus by <a href="HouseholderTransformation.vba">HouseholderTransformation.vba</a> and <a href="SingleShiftQRforComplexEigenProblem.vba">SingleShiftQRforComplexEigenProblem.vba</a>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You can use <a href="DoubleShiftQR.vba">DoubleShiftQR.vba</a> instead of SingleShiftQRforComplexEigenProblem.vba, but it is NOT stable.<br>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <i> [Extraction procedure] </i>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Asymmetric matrix -->(Householder transformation)-->

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Hessenberg matrix -->(Single QR shift)--> <u>Extraction of Complex Eigen value</u>


2) Eigenvectors are extracted from original matrix (asymmetric matrix) and Eigenvalues by <a href="EigenVectorForComplex.vba">EigenVectorForComplex.vba</a>

&nbsp;&nbsp;![#008c15](https://placehold.it/15/008c15/000000?text=+)
<em> How to use?</em><br>

0) Turn visual basic window (ALT + F11 or develop on ribon -> Visual Basic) on ,and insert standard module. Then copy and paste user functions.

1) Select range on Excel worksheet. (Excel extracts calculation values on this range)

2) Type user function.

3) Ctrl + Shift + Enter (It makes those cells as a chunk)

&nbsp;&nbsp; <em> For example </em>

&nbsp;&nbsp;&nbsp;&nbsp;=Householder(A1:D4)

&nbsp;&nbsp;&nbsp;&nbsp;In this case, range A1:D4 are reference matrix (asymmetric matrix you want to extract Eigenvalues)

---
## Real Eigenvalue Problem
### <a href="EigenValuesForSymmetricMatirx.vba">Jacobi Eigenvalue algorithm</a>

&nbsp;&nbsp;A symmetric matrix is transformed to diagonal matrix by Givens rotation.

&nbsp;&nbsp;Each diagonal elements are Eigenvalues.

### <a href="CholeskyDecomposition.vba">Cholesky decompostion</a>

&nbsp;&nbsp;<a href="CholeskyDecomposition.vba">Cholesky decompostion</a> is for generalized Eigenvalue problem.<br>

&nbsp;&nbsp;For example noise and vibration (MK type) ,buckling (Eular buckling) and etc..

&nbsp;&nbsp;Cholesky decomposition makes lower triangle matrix L from A square matrix (or upper triangle matrix U).

&nbsp;&nbsp;A=L\*L^t or A=U^t\*U ( L^t is transpose(L) )

&nbsp;&nbsp;Concept of Cholesky decomposition is Root of matrix. (I think so)

&nbsp;&nbsp;LU decomposition is NOT same, but decomposition for Gauss elimination.

&nbsp;&nbsp;

---

# Under Construction

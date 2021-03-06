# Math
Mainly linear algebra function by various computer languages
---

## Complex Eigenvalue Problem
### Single QR method

1) Extraction of Eigenvalus by <a href="HouseholderTransformation.vba">HouseholderTransformation.vba</a> and <a href="SingleShiftQRforComplexEigenProblem.vba">SingleShiftQRforComplexEigenProblem.vba</a>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You can use <a href="DoubleShiftQR.vba">DoubleShiftQR.vba</a> instead of SingleShiftQRforComplexEigenProblem.vba, but it is NOT stable.<br>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <i> [Extraction procedure] </i>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Asymmetric matrix -->(Householder transformation)-->

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Hessenberg matrix -->(Single QR shift)--> <u>Extraction of Complex Eigenvalue</u>


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

### Applications

&nbsp;&nbsp;Standard Eigenvalue problems are mainly for tensors probably.

&nbsp;&nbsp;Tensors are a kind of number, are similar to complex value. And those notation are same as matrix.

&nbsp;&nbsp;In real number and non-direction dependency, tensor is symmetric matrix.

&nbsp;&nbsp;As N * N matrix, Tensor includes N pairs of principle value and principle axis.

&nbsp;&nbsp;And principle values are as Eigenvalues, princple axes are as Eigenvector.

&nbsp;&nbsp;So, you can use real Eigenvalue problem for solving tensor to principle values and axes.

&nbsp;&nbsp;<a href="TensorConverter.htm">TensorConverter.htm</a> is able to extract principles from tensor by browser easily.

&nbsp;&nbsp;This Eingenvalue solver is translation from vba (Jacobi Eigenvalue algorithm) to javascript.

---

## <a href="SVD.vba">Singular Value Decomposition (SVD)</a>

My understanding is that SVD is a type of Eigenvalue problem for non-square matrix.

SVD has not only similarity of Eigenvalue problem but also similarity of inverse matrix too.

In major use , it is for pseudoinverse in multiple regression equation. 

<a href="http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm">Proper example of SVD</a> (It is much easier to understand than my explanation.)

Another way for pseudoinverse is for example <a href="http://help.matheass.eu/en/Pseudoinverse.html">here</a>.

---

## Library for matrix calculation in Visual Studio C++

<a href="Matrix.h">Matrix.h</a> is library for matrix.

I wrote by Visual Studio 2013 C++. But almost of C++ compilers are able to compile by small modification.

See prototype to each function you want to use.

Sorry for comments are in Japanese, I will add explanations if I felt like to do. XD


# Under Construction

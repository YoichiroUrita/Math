//Header file :matrix00
/*******************************************
 * 行列計算用ライブラリ
 *     ver.2018/01/05
 *******************************************/
#define _USE_MATH_DEFINES
#include <math.h>

//#include <amp.h>
//#include "stdafx.h"
//#include <amp_math.h> //for M_PI
#ifndef Matrix00_H
#define Matrix00_H

typedef double **Mata, **Matb, **Matc, **Matd,**Mate,*Veca,*Vecb,*Vecc;

int madd(Mata, Matb, Matc, int, int);//行列の和
int msub(Mata, Matb, Matc, int, int);//行列の差
int mtimes(Mata, double, Matc, int, int);//行列のｎ倍
int mmult(Mata, Matb, Matc, int, int, int);//行列の内積
int mtrans(Mata, Matb, int, int);//転置
int minv(Mata, Matb, int);//逆行列
int mcomp(Mata, Matb, Matc, int, int);//実行列による複素表現
int mreal(Mata, Matb,int, int);//複素表現から実部のみ取りだす
int msimi(Mata, Matb, Matc, int, int);//相似変換
int mchol(Mata, Matb, int);//コレスキー分解
int meigj(Mata, Veca, Matb, int);//標準実固有値分解　Jacobi法
int mskew(Veca, Mata);//歪対称行列（３＊３）
int moutpro(Veca, Vecb,Vecc);//外積（３＊３）
int mdiag(Veca, Mata, int);//対角行列化
int mcompdiag(Mata, Matb, int);//複素対角行列化
int mhouse(Mata, Matb, int);//Householder変換（Hessenberg行列化）
int mqr(Mata, Matb, int);//QR分解（Hessenberg用）
int mevqr(Mata, Matb, int);//複素固有値（QR分解用）
int mevecqr(Mata, Matb, Matc, int);//複素固有ベクトル（ガウスの消去法）
int mceig(Mata, Matb, Matc, int);//標準複素固有値・固有ベクトル
int mmerge_row(Mata, Matb, Matc, int, int, int);//行列の連結　縦方向
int mmerge_col(Mata, Matb, Matc, int, int, int);//行列の連結　横方向
int mgevecqr(Mata, Matb, Matc, Matd, Mate, int);//一般複素固有ベクトル
int GCEIGQR(Mata, Matb, Matc, Matd, Mate, int);//一般複素固有値解析　固有値・固有ベクトル　QR法
int GREIGJ(Mata, Matb, Veca, Matc, int);//一般実固有値解析　固有値・固有ベクトル　Jacobi法
int msplit(Mata, Matb, int, int, int, int);//行列の抜き取り
int mgmat(Mata, Matb, Matc, Matd, double, int);//Gマトリクス作成
int mpinv(Mata, Matb, int , int );//疑似逆行列
int mrotX(Mata, double);//オイラーの回転行列　X軸廻り
int mrotY(Mata, double);//オイラーの回転行列　Y軸廻り
int mrotZ(Mata, double);//オイラーの回転行列　Z軸廻り
double radians(double);//rad <- deg
double degrees(double);//deg <- rad
int mrotyxyzyx(Mata, double, double, double, double, double, double);//オイラー回転行列の積
int mDmat(Mata, Matb, double, double, double);//DistanceMatrix作成
int mX0byX1(Mata, Matb, int, int);//Gain用計算(X1/X0)
int mGuyan(Mata, Matb, int, int);//Guyan縮約
int mTrv(Mata, Matb, int , int );//振動伝達率(F1/F0)
int mSEC(Mata, Matb, Matc, int);//歪エネルギー分布計算
int mKEC(Mata, Matb, Matc, int);//運動エネルギー分布計算
int GCEIGQR2(Mata, Matb, Matc, Matd, Mate, int);//一般複素固有値解析　固有値・固有ベクトル　複素ばね用（減衰なし）
int mhouse2(Mata, Matb,Matc, int);//Householder変換（Hessenberg行列化） Qベクトル付
int mqr2(Mata, Matb,Matc, int);//QR分解（Hessenberg用） Qベクトル付
int mceig2(Mata, Matb, Matc, int);//標準複素固有値・固有ベクトル 固有ベクトルは、元の行列と同サイズ
int mqr3(Mata, Matb,int, int);//QR分解（Hessenberg用） イタレーション回数指定
int mqr4(Mata, Matb, Matc, int,int);//QR分解（Hessenberg用） Qベクトル付　イタレーション回数指定
int mdqr(Mata, Matb, int, int);//シフト付きダブルQR分解(Hessenberg用）
double mysqrt(double);//sqrt関数
int momit3(double **mass, double**damp, double**stiff, double**stiff_im, int *index, int row);
int mred(double Mata, double Matb, int, int);//縮約

#endif

//static const double PI = 3.14159265358979323846264338327950288;

int madd(double **a, double **b, double **c, int row, int col){
	//6/19/2014
	int  i, j;


	//行列の各項の足し算
	{

		if (row < 1 || col < 1) return (999);
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				c[i][j] = a[i][j] + b[i][j];
			}
		}
	}
	return 0;
}
int msub(double **a, double **b, double **c, int row, int col){
	//6/19/2014
	int  i, j;

	//行列の各項の引き算
	{
		if (row < 1 || col < 1)return(999);
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				c[i][j] = a[i][j] - b[i][j];
			}
		}
	}
	return 0;
}
int mtimes(double **a, double b, double **c, int row, int col){
	//6/19/2014
	int i, j;
	//行列のｂ倍
	{
		if (row < 1 || col < 1)return(999);
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				c[i][j] = b*a[i][j];
			}
		}
	}
	return 0;
}
int mmult(double **a, double **b, double **ab, int row_a, int col_a, int col_b){
	//6/19/2014
	int i, j, k;
	//行列の内積
	//a[row][col]*b[col][col2]


	{
		if (row_a < 1 || col_a < 1)return(999);

		for (i = 0; i < row_a; i++)
		{
			for (j = 0; j < col_b; j++)
			{
				ab[i][j] = 0.0;
				for (k = 0; k < col_a; k++)
				{
					ab[i][j] += a[i][k] * b[k][j];
				}
			}
		}
	}

	return 0;
}
int mtrans(double **a, double **aT, int row, int col){
	//6/20/2014
	int i, j;
	//転置

	{
		if (row < 1 || col < 1)return(999);
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				aT[j][i] = a[i][j];
			}
		}
	}
	return 0;
}

int minv(double** inMat, double** outMat, int n)
//6/20/2014
//inverse matrix
//http://seesaawiki.jp/w/pafuhana1213/d/n%BC%A1%A4%CE%B5%D5%B9%D4%CE%F3 参照
{
	double tmp;
	int i, j, k;

	// 配列動的確保
	double **a = new double*[n];
	double **inv_a = new double*[n];
	for (i = 0; i < n; i++)
	{
		a[i] = new double[n];
		inv_a[i] = new double[n];
	}

	// 入力行列のコピー
	for (j = 0; j<n; j++){
		for (i = 0; i<n; i++){
			a[j][i] = inMat[i][j];
		}
	}
	// n次の単位行列の生成
	for (j = 0; j<n; j++){
		for (i = 0; i<n; i++){
			inv_a[j][i] = (i == j) ? 1.0 : 0.0;
		}
	}
	// ピボット選択を行ったGauss-Jordan法
	for (k = 0; k<n; k++){
		//ピボット選択 k行k列目の要素の絶対値が最大に
		int max = k;
		for (j = k + 1; j<n; j++){
			if (fabs(a[j][k]) > fabs(a[max][k])){
				max = j;
			}
		}
		// 行の入れ替え
		if (max != k){
			for (i = 0; i<n; i++){
				// 入力行列側
				tmp = a[max][i];
				a[max][i] = a[k][i];
				a[k][i] = tmp;
				// 単位行列側
				tmp = inv_a[max][i];
				inv_a[max][i] = inv_a[k][i];
				inv_a[k][i] = tmp;
			}
		}
		// k行k列目の要素が1になるように
		tmp = a[k][k];
		for (i = 0; i<n; i++){
			a[k][i] /= tmp;
			inv_a[k][i] /= tmp;
		}
		// k行目のk列目以外の要素が0になるように
		for (j = 0; j<n; j++){
			if (j != k){
				tmp = a[j][k] / a[k][k];
				for (i = 0; i<n; i++){
					a[j][i] = a[j][i] - a[k][i] * tmp;
					inv_a[j][i] = inv_a[j][i] - inv_a[k][i] * tmp;
				}
			}
		}

	}
	//copy
	for (j = 0; j<n; j++){
		for (i = 0; i<n; i++){
			outMat[i][j] = inv_a[j][i];
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < n; i++)
	{
		delete[] a[i], inv_a[i];
	}
	delete[] a, inv_a;
	return 0;
}
int mcomp(double **Re, double **Im, double **cc, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			cc[2 * i][2 * j] = Re[i][j];
			cc[2 * i + 1][2 * j + 1] = Re[i][j];
			cc[2 * i][2 * j + 1] = -Im[i][j];
			cc[2 * i + 1][2 * j] = Im[i][j];
		}
	}
	return 0;
}
int mreal(double **comp, double **Re, int m, int n)//m,nはcompの行列サイズ
{
	int i, j;
	for (i = 0; i < m; i+=2)
	{
		for (j = 0; j < n; j+=2)
		{
			Re[i/2][j/2]=comp[i][j];
		}
	}
	return 0;
}
int mimag(double **comp, double **Im, int m, int n)//m,nはcompの行列サイズ
{
	int i, j;
	for (i = 0; i < m; i += 2)
	{
		for (j = 0; j < n; j += 2)
		{
			Im[i / 2][j / 2+1] = -comp[i][j+1];
		}
	}
	return 0;
}

//c=b^t*a*b
//a:[row][row]
//b:[row][col]
int msimi(double **a, double **b, double **c, int row, int col){
	//6/21/2014
	int i, j, k;
	//行列の相似変換(similarity transformation)
	//b^t[col][row]*a[row][row]*b[row][col]
	// 配列動的確保
	double **ab = new double*[row];
	for (i = 0; i < row; i++)
	{
		ab[i] = new double[col];
	}

	{
		if (row < 1 || col < 1)return(999);
		//ab=a*b
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				ab[i][j] = 0.0;
				for (k = 0; k < row; k++)
				{
					ab[i][j] += a[i][k] * b[k][j];
				}
			}
		}
	}
	//c=b^t*ab
	for (i = 0; i < col; i++)
	{
		for (j = 0; j < col; j++)
		{
			c[i][j] = 0.0;
			for (k = 0; k < row; k++)
			{
				c[i][j] += b[k][i] * ab[k][j];
			}
		}
	}

	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] ab[i];
	}
	delete[] ab;
	return 0;
}
int mchol(double **a, double **u,int row){
	//6/21/2014
	//Cholesky decomposition
	int i, j, k;
	double b1, b2,b4;

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			u[i][j] = 0.;
		}
	}

	if (a[0][0]== 0)return(999);
	//upper right triangle
	//u11
	 b1= sqrt(a[0][0]);
	 u[0][0] = b1;
	//u1j
	for (j = 1; j < row; j++)
	{
		u[0][j] = a[0][j] / b1;
		//System::Console::WriteLine(u[0][j]+L" "+b1);
	}

	//uii & uij
	for (i = 1; i < row; i++)
	{
		b1 = a[i][i];
		b4 = 0.0;
		//uii
		for (k = 0; k < i ; k++)
		{
			b1 -= u[k][i] * u[k][i];
		}
		b2= sqrt(b1);
		u[i][i] = b2;
		//uij
		if (i < row - 1)
		{
			for (j = i + 1; j < row; j++)
			{
				b1 = a[i][j];
				for (k = 0; k < i ; k++)
				{
					b1 -= u[k][i] * u[k][j];
				}
				u[i][j] = b1 /b2;
			}
		}
	}
	return 0;
}
int meigj(double **aa, double *eig, double **vec, int row){
	//6/21/2014
	//Eigen value with Jacobian
	int i, j,  k_new,  i_rot, j_rot, kmax;
	double b1, b2, b3,eps,k_iter,ax,v_cos,v_sin,v_tan;
	// 配列動的確保
	double **a = new double*[row];
	for (i = 0; i < row; i++)
	{
		a[i] = new double[row];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			a[i][j] = aa[i][j];
		}
	}

	eps = 10e-16;//define acceptable error
	kmax = 5000;//maximum iteration number
	//vec=unit matrix
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			if (i == j)
			{
				vec[i][i] = 1.;
			}
			else
			{
				vec[i][j] = 0.;
			}
		}
	}
	for (k_iter = 0; k_iter < kmax;k_iter++)
	{
		k_new = 0;//flag:if nondiagonal elementhas unacceptable error,flag was 1
		ax = 0.;
		//looking for maximum element and replace ax
		for (i = 0; i < row; i++)
		{
			for (j = i + 1; j < row; j++)
			{
				b1 = fabs(a[i][j]);
				if (b1>ax)
				{
					i_rot = i;
					j_rot = j;
					ax = b1;
				}
			}
		}
		if (ax > eps)
		{
			k_new = 1;


			//calculate cosine and sine
			b1 = a[i_rot][i_rot] - a[j_rot][j_rot];
			b2 = -b1 - sqrt(b1*b1 + 4 * ax*ax);
			b3 = 2 * a[i_rot][j_rot];
			v_tan = b2 / b3;
			v_cos = 1 / sqrt(1 + v_tan*v_tan);
			v_sin = v_cos*v_tan;

			//rotate vec
			for (i = 0; i < row; i++)
			{
				b1 = vec[i][i_rot];
				vec[i][i_rot] = v_cos*b1 + v_sin*vec[i][j_rot];
				vec[i][j_rot] = v_cos*vec[i][j_rot] - v_sin*b1;
			}

			//rotate A matrix
			for (i = 0; i < i_rot; i++)
			{
				b1 = a[i][i_rot];
				a[i][i_rot] = v_cos*b1 + v_sin*a[i][j_rot];
				a[i][j_rot] = v_cos*a[i][j_rot] - v_sin*b1;
			}
			for (i = i_rot + 1; i < j_rot; i++)
			{
				b1 = a[i_rot][i];
				a[i_rot][i] = v_cos*b1 + v_sin*a[i][j_rot];
				a[i][j_rot] = v_cos*a[i][j_rot] - v_sin*b1;
			}
			for (i = j_rot + 1; i < row; i++)
			{
				b1 = a[i_rot][i];
				a[i_rot][i] = v_cos*b1 + v_sin*a[j_rot][i];
				a[j_rot][i] = v_cos*a[j_rot][i] - v_sin*b1;
			}

			//erase nondiagonal element
			b1 = a[i_rot][i_rot];
			b2 = 2 * v_cos*v_sin*a[i_rot][j_rot];
			b3 = v_cos*v_cos*b1 + b2;

			a[i_rot][i_rot] = b3 + v_sin*v_sin*a[j_rot][j_rot];
			b3 = v_cos*v_cos*a[j_rot][j_rot] + v_sin*v_sin*b1;
			a[j_rot][j_rot] = b3 - b2;
			a[i_rot][j_rot] = 0;
			a[j_rot][i_rot] = a[i_rot][j_rot];

		}
		if (k_new = 0)
		{
			break;
		}
	}

	//eig <- diagonal elements
	for (i = 0; i < row; i++)
	{
		eig[i] = a[i][i];
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] a[i];
	}
	delete[] a;
	return 0;
}
int mskew(double* vec, double** outMat){
	//6/22/2014
	//skew symmetric matrix(multiple operator)
	int i;
	for (i = 0; i < 3; i++)
	{
		outMat[i][i] = 0.;
	}
	outMat[0][1] = -vec[2];
	outMat[0][2] = vec[1];
	outMat[1][2] = -vec[0];

	outMat[1][0] = vec[2];
	outMat[2][0] = -vec[1];
	outMat[2][1] = vec[0];

	return 0;
}

int moutpro(double*vec1, double*vec2, double*outvec)
{
	//1/12/2018
	//outerproduct
	int i;
	double**skewmat = new double*[3];
	double**vecmat1 = new double*[3];
	double**vecmat2 = new double*[3];
	for (i = 0; i < 3; i++)
	{
		skewmat[i] = new double[3];
		vecmat1[i] = new double[1];
		vecmat2[i] = new double[1];
		vecmat1[i][0] = vec2[i];
	}
	mskew(vec1, skewmat);
	mmult(skewmat, vecmat1, vecmat2,3,3,1);
	for (i = 0; i < 3; i++)
	{
		outvec[i] = vecmat2[i][0];
	}
	return 0;
}

int mdiag(double* vec, double** outMat, int row){
	//6/22/2014
	//diagonal matrix
	int i, j;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			if (i == j)
			{
				outMat[i][i] = vec[i];
			}
			else
			{
				outMat[i][j] = 0.;
			}
		}
	}
	return 0;
}
int mcompdiag(double** vec, double** outMat, int row){
	//6/22/2014
	//diagonal matrix
	int i, j;
	for (i = 0; i < row; i+=2)
	{
		for (j = 0; j < row; j+=2)
		{
			if (i == j)
			{
				outMat[i][i] = vec[i][0];
				outMat[i + 1][i + 1] = vec[i + 1][1];
				outMat[i][i + 1] = vec[i][1];
				outMat[i + 1][i] = vec[i + 1][0];
			}
			else
			{
				outMat[i][j] = 0.;
				outMat[i + 1][j] = 0.;
				outMat[i][j + 1] = 0.;
				outMat[i + 1][i] = 0.;
			}
		}
	}
	return 0;
}
int mhouse(double** inmat, double** outmat, int row){
	//6/22/2014
	//Householder transformation
	int i, j, k,m;
	double b1, b2, b3, sigma, c1;
	// 配列動的確保
	double **qmat = new double*[row];
	double **hessen = new double*[row];
	double **buf = new double*[row];
	double *uvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		hessen[i] = new double[row];
		buf[i] = new double[row];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			hessen[i][j] = inmat[i][j];
			buf[i][j] = 0.;
		}
	}
	//
	for (k = 0; k < row - 2; k++)
	{
		//vector u
		b1 = 0.;
		for (i = k + 1; i < row; i++)
		{
			b1 += hessen[i][k] * hessen[i][k];
		}
		b2 = sqrt(b1);
		b3 = hessen[k + 1][k];
		if (b3 >= 0)
		{
			sigma = b2;
		}
		else
		{
			sigma = -b2;
		}

		uvec[k + 1] = b3 + sigma;
		c1 = sigma*uvec[k + 1];
		uvec[k] = 0;

		for (i = k + 2; i < row; i++)////
		{
			uvec[i] = hessen[i][k];
		}

		//Matrix Q
		//initialize
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0.;
			}
		}
		for (i = 0; i < row; i++)
		{	
			qmat[i][i] = 1.;
		}

		for (i = k; i < row; i++)
		{
			for (j = k; j < row; j++)
			{
				if (c1 != 0)
				{
					qmat[i][j] -= uvec[i] * uvec[j] / c1;
				}
			}
		}

		//QA calculation
		for (i = 0; i <= k; i++)////
		{
			for (j = k; j < row; j++)
			{
				buf[i][j] = hessen[i][j];
				//double tmp = outmat[i][j];
			}
		}
		for (i = k + 1; i < row; i++)
		{
			for (j = k; j < row; j++)
			{
				b1 = 0.;
				for (m = k; m < row; m++)
				{
					b1 += qmat[i][m] * hessen[m][j];
					//double tmp = qmat[i][m];
				}
				buf[i][j] = b1;
			}
		}
		//QAQ^t
		for (j = 0; j < row; j++)
		{
			for (i = k; i < row; i++)
			{
				hessen[i][j] = buf[i][j];
			}
		}
		for (j = k + 1; j < row; j++)
		{
			for (i = 0; i < row; i++)
			{
				b1 = 0.;
				for (m = k; m < row; m++)
				{
					b1 += buf[i][m] * qmat[m][j];
				}
				hessen[i][j] = b1;
			}
		}
	}
	//lower left =0
	for (i = 2; i < row; i++)
	{
		for (j = 0; j <= i - 2; j++)
		{
			hessen[i][j] = 0.;
		}
	}
	for (i = 0; i < row  ; i++)
	{
		for (j = 0; j < row; j++)
		{
			outmat[i][j] = hessen[i][j];
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] buf[i],qmat[i],hessen[i];
	}
	delete[] buf,qmat,uvec,hessen;

	return 0;
}
//QR factorization
int mqr(double** inmat, double** outmat, int row){
	//6/22/2014
	//pukiwiki for PBCG Lab 参照
	int i, j, k, m,cnt,itermax;
	double c, s, b1, rq, eps,er;
	// 配列動的確保
	double **qmat = new double*[row];
	double **rmat = new double*[row];
	double **tmat = new double*[row];
	double *uvec = new double[row];
	double *vvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		rmat[i] = new double[row];
		tmat[i] = new double[row];
	}

	//copy Hessenberg matrix to rmat & outmat
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			rmat[i][j] = inmat[i][j];
			outmat[i][j] = inmat[i][j];
			tmat[i][j] = 0.;
		}
	}

	itermax = 5000;//iteration maximum number
	eps = 1.0e-8;

	//
	for (cnt = 0; cnt < itermax; cnt++)
	{
		//double tmp = rmat[0][1];
		//Q <- unit matrix
		for (i = 0; i < row;i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0.;
			}
		}
		for (i = 0; i < row; i++)
		{
			qmat[i][i] = 1.;
		}
		for (k = 0; k < row - 1; k++)
		{
			//sine and cosine
			//double tmp1 = rmat[k][k];
			//double tmp2 = rmat[k + 1][k];
			b1 = sqrt(rmat[k][k] * rmat[k][k] + rmat[k + 1][k] * rmat[k + 1][k]);
			c = rmat[k][k] / b1;
			s = -rmat[k + 1][k] / b1;

			//R calculation
			for (j = k + 1; j < row; j++)
			{
				uvec[j] = c*rmat[k][j] - s*rmat[k + 1][j];
				vvec[j] = s*rmat[k][j] + c*rmat[k + 1][j];
			}

			rmat[k][k] = b1;
			rmat[k + 1][k] = 0.;
			
			//double tmp0 =rmat[0][1];

			for (j = k + 1; j < row; j++)
			{
				rmat[k][j] = uvec[j];
				rmat[k + 1][j] = vvec[j];
			}

			//Q calculation
			for (j = 0; j <= k; j++)//
			{
				uvec[j] = c*qmat[k][j];
				vvec[j] = s*qmat[k][j];
			}
			qmat[k][k + 1] = -s;
			qmat[k + 1][k + 1] = c;
			for (j = 0; j <= k; j++)
			{
				qmat[k][j] = uvec[j];
				qmat[k + 1][j] = vvec[j];
			}
		}

		//A=RQ (outmat=RQ)
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rq = 0.;
				for (m = 0; m < row; m++)
				{
					rq += rmat[i][m] * qmat[j][m];
					//double tmp1 = rmat[i][m];
					//double tmp2 = qmat[j][m];
				}
				tmat[i][j] = rq;
			}
		}

		//convergence check
		er = 0;
		for (i = 0; i < row; i++)
		{
			er += fabs(tmat[i][i] - inmat[i][i]);
		}
		if (er < eps)
		{
			break;
		}
		//R,A <-RQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rmat[i][j] = tmat[i][j];
				outmat[i][j] = tmat[i][j];
			}
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] rmat[i], qmat[i],tmat[i];
	}
	delete[] rmat,qmat,tmat,uvec,vvec;

	return cnt;
}
//Eigen value from QR decomposition
int mevqr(double** inmat, double** eig ,int row){
	//6/22/2014
	//Eigen value from QR decomposition
	int i,j,k;
	int tmp_size;
	double eps, bb, cc;

	eps = 1.e-8;

	i = row-1;
	//
	while (i>=0)
	{
		//search 0
		//single Eigen value (not pair)
		if (i == 0)
		{
			//real part
			eig[0][row-1] = inmat[0][0];
			//imag part
			eig[1][row-1] = 0;
			break;
		}
		if (fabs(inmat[i][i - 1]) < eps)
		{
			//real part
			eig[0][row - i-1] = inmat[i][i];
			double tmp = eig[0][row - i - 1];
			//imag part
			eig[1][row - i-1] = 0;
			i--;
		}
		else
			//Eigen value pair
		{
			if (i == 1)
			{
				bb = inmat[0][0] + inmat[1][1];
				cc = inmat[0][0] * inmat[1][1] - inmat[0][1] * inmat[1][0];
				if ((bb*bb - 4 * cc) > 0)
				{
					//real pair
					//real part
					eig[0][row-2] = (bb - sqrt(bb*bb - 4 * cc)) / 2;
					eig[0][row-1] = (bb + sqrt(bb*bb - 4 * cc)) / 2;
					//imag part
					eig[1][row-2] = 0;
					eig[1][row-1] = 0;
				}
				else
				{
					//complex pair
					//real part
					eig[0][row-2] = bb / 2;
					eig[0][row-1] = bb / 2;
					//imag part
					eig[1][row-2] = -(sqrt(fabs(bb*bb - 4 * cc))) / 2;
					eig[1][row-1] = (sqrt(fabs(bb*bb - 4 * cc))) / 2;
				}
				break;
			}
			else
			{
				//in convegence
				if (fabs(inmat[i - 1][i - 2]) < eps)
				{
					bb = inmat[i - 1][i - 1]+inmat[i][i];
					cc = inmat[i - 1][i - 1] * inmat[i][i] - inmat[i - 1][i] * inmat[i][i - 1];
					if ((bb*bb - 4 * cc)>0)
					{
						//real pair
						//real part
						eig[0][row - i-1] = (bb - sqrt(bb*bb - 4 * cc)) / 2;
						eig[0][row - i] = (bb + sqrt(bb*bb - 4 * cc)) / 2;
						//imag part
						eig[1][row - i-1] = 0;
						eig[1][row - i] = 0;
					}
					else
					{
						//complex pair
						//real part
						eig[0][row - i-1] = bb / 2;
						eig[0][row - i] = bb / 2;
						//imag part
						eig[1][row - i-1] = -(sqrt(fabs(bb*bb - 4 * cc))) / 2;
						eig[1][row - i] = (sqrt(fabs(bb*bb - 4 * cc))) / 2;
					}
					i -= 2;
				}
				else
				{
					return 1000+(row-i-1);
					break;
				}
			}
		}

	}
	//transformation error check
	for (j = 2; j < row; j++)
	{
		for (k = 0; k < j - 2; k++)
		{
			if (fabs(inmat[j][k])> eps)
			{
				tmp_size = j + 1;
			}
		}
	}


	return tmp_size;
}
//complex standard Eigen vector by Gaussian elimination
//eig[2][row]
//ceigvec[2*row][row]
int mevecqr(double** inmat, double** eig,double** ceigvec, int row){
	//6/23/2014

	int i, j, k, kk,kset,icc,nn,ib,ic,flag;
	double b1, pivotvalue, valueb,alp,bet,mag;
	//double check[8][8];
	// 配列動的確保
	double **a = new double*[2 * row];
	double **eigvec = new double*[row];
	for (i = 0; i < 2*row; i++)
	{
		a[i] = new double[2 * row];
	}
	for (i = 0; i < row; i++)
	{
		eigvec[i] = new double[2*row];
	}

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < 2 * row; j++)
		{
			eigvec[i][j] = 0.;
		}
	}


	flag = 0;
	//
	icc = 0;
	for (kset = 0; kset < row; kset++)
	{
		for (i = 0; i < (2 * row); i++)
		{
			for (j = 0; j < (2 * row); j++)
			{
				a[i][j] = 0.;
			}
		}
		if (icc == 0)//case:not 2nd complex pair
		{
			alp = eig[0][kset];
			bet = eig[1][kset];
			if (bet == 0)//case:real only
			{
				icc = 0;
				for (i = 0; i < row; i++)
				{
					for (j = 0; j < row; j++)
					{
						a[i][j] = inmat[i][j];
					}
					a[i][i] -= alp;
				}
				nn = row;
			}
			else//case:1st complex pair
			{
				icc = 1;
				for (i = 0; i < row; i++)
				{
					ib = row + i;
					for (j = 0; j < row; j++)
					{
						a[i][j] = inmat[i][j];
						a[ib][row + j] = inmat[i][j];
					}
					a[i][i] = a[i][i] - alp;
					a[ib][ib] = a[ib][ib] - alp;
					a[i][ib] = bet;
					a[ib][i] = -bet;
				}
				//check array
				//for (i = 0; i < row * 2; i++)
				//{
				//	for (j = 0; j < row * 2; j++)
				//	{
				//		check[i][j] = a[i][j];
				//	}
				//}
				nn = 2 * row;
			}
		}
		else
		{
			icc = 2;//case:2nd complex pair
		}

		if (icc < 2)
		{
			eigvec[kset][0] = 1.;//dammy number
			if ((a[0][0] != 0) && (a[1][0] == 0))//pivot
			{
				for (j = 0; j < nn; j++)
				{
					a[1][j] = a[0][j];
				}
			}
			for (i = 1; i < nn; i++)
			{
				eigvec[kset][i] = -a[i][0];
				a[i][0] = 0;
			}

			//Gaussian elimination
			for (k = 1; k < nn; k++)
			{
				ic = 0;//flag for pivoting
				for (i = k; i < nn; i++)
				{
					if (a[i][k] != 0)//search pivoting point
					{
						ic = 1;
						kk = i;
						break;
					}
				}
				if (ic == 0)//case:error
				{
					return 999;
					flag = 999;
				}
				pivotvalue = a[kk][k];
				for (j = k; j < nn; j++) //pivot k rows to kk rows each other
				{
					b1 = a[kk][j];
					a[kk][j] = a[k][j];
					a[k][j] = b1 / pivotvalue;
				}
				b1 = eigvec[kset][kk];
				eigvec[kset][k] = b1 / pivotvalue;
				for (i = 1; i < nn; i++)//elimination without pivot row
				{
					if (i != k)
					{
						valueb = a[i][k];
						for (j = k; j < nn; j++)
						{
							a[i][j] -= valueb*a[k][j];
						}
						eigvec[kset][i] -= valueb*eigvec[kset][k];
						//double tmp1 = eigvec[kset][i];
					}
				}
			}//Gaussian elimination end
			//for (i = 0; i < 2 * row; i++)
			//{
			//	for (j = 0; j < 2 * row; j++)
			//	{
			//		check[i][j] = a[i][j];
			//	}
			//}
			int x = 0;
		}
		if (icc == 2)//2nd complex pair
		{
			for (i = 0; i < nn; i++)
			{
				if (i < row )
				{
					eigvec[kset][i] = eigvec[kset - 1][i];
				}
				else
				{
					eigvec[kset][i] = -eigvec[kset - 1][i];
				}
			}
			icc = 0;
		}

		//normalize magunitude = 1
		mag = 0.;
		for (i = 0; i < nn; i++)
		{
			mag += eigvec[kset][i] * eigvec[kset][i];
			//double tmp0 = eigvec[kset][i];
		}
		mag = sqrt(mag);
		for (i = 0; i < nn; i++)
		{
			eigvec[kset][i] /= mag;
			//double tmp1 = eigvec[kset][i];
		}
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < (2 * row); j++)
		{
			ceigvec[j][i] = eigvec[i][j];//change vector direction (row to column)
			//double tmp = ceigvec[j][i];
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
	delete[] eigvec[i];
	}
	for (i = 0; i < 2 * row; i++)
	{
		delete[] a[i];
	}
	delete[] eigvec, a;

	return flag;
}
//Standard complex eigenvalu problem
int mceig(double** inmat, double** eig, double** eigvec, int row){
	//6/23/2014
	//Standard complex eigenvalu problem
	int flag;
	// 配列動的確保
	double **hessenberg = new double*[row];
	double **rightupper = new double*[row];
	for (int i = 0; i < row; i++)
	{
		hessenberg[i] = new double[row];
		rightupper[i] = new double[row];
	}
	mhouse(inmat, hessenberg, row);//ヘッセンベルグ化
	mqr(hessenberg, rightupper, row);//右上三角行列化（部分的に副対角項有り）
	mevqr(rightupper, eig, row);//複素固有値算出
	flag=mevecqr(inmat, eig, eigvec, row);//複素固有ベクトル算出
	// 確保したメモリ領域の解放
	for (int i = 0; i < row; i++)
	{
		delete[] hessenberg[i],rightupper[i];
	}
	delete[] hessenberg, rightupper;
	return flag;
}
int mmerge_row(double**inmat1, double**inmat2, double**outmat, int row1, int row2, int col){
	//6/23/2014
	//Merge matrix along row direction
	//inmat1[row1][col] ,inmat2[row2][col] => outmat[row1+row2][col]
	int i, j;
	for (i = 0; i < row1; i++)
	{
		for (j = 0; j < col; j++)
		{
			outmat[i][j] = inmat1[i][j];
		}
	}
	for (i = 0; i < row2; i++)
	{
		for (j = 0; j < col; j++)
		{
			outmat[row1 + i][j] = inmat2[i][j];
		}
	}
	return 0;
}
int mmerge_col(double**inmat1, double**inmat2, double**outmat, int row, int col1, int col2){
	//6/23/2014
	//Merge matrix along col direction
	//inmat1[row][col1],inmat2[row][col2] => outmat[row][col1+col2]
	int i, j;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col1; j++)
		{
			outmat[i][j] = inmat1[i][j];
		}
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col2; j++)
		{
			outmat[i][col1 + j] = inmat2[i][j];
		}
	}
	return 0;
}
int mgevecqr(double** mass,double**damp,double**stiff, double** eig, double** ceigvec, int row){
	//6/23/2014
	//generalized Eigen vector by Gaussian elimination
	//eig[2][2*row]
	//ceigvec[2*row][row]  [2*row][2*row]?
	int i, j, k, kk, kset, icc, nn, ib, ic;
	double b1, pivotvalue, valueb, alp, bet, mag;
	//double check[8][8];
	// 配列動的確保
	double **a = new double*[2 * row];
	double **eigvec = new double*[2 * row];
	for (i = 0; i < 2 * row; i++)
	{
		a[i] = new double[2 * row];
		eigvec[i] = new double[2*row];
	}
	for (i = 0; i < 2*row; i++)
	{
		for (j = 0; j < 2 * row; j++)
		{
			eigvec[i][j] = 0.;
		}
	}

	//
	icc = 0;
	for (kset = 0; kset < row*2; kset++)
	{
		for (i = 0; i < (2 * row); i++)
		{
			for (j = 0; j < (2 * row); j++)
			{
				a[i][j] = 0.;
			}
		}
		if (icc == 0)//case:not 2nd complex pair
		{
			alp = eig[0][kset];
			bet = eig[1][kset];
			//case:1st complex pair
			{
				icc = 1;
				for (i = 0; i < row; i++)
				{
					ib = row + i;
					for (j = 0; j < row; j++)
					{
						a[i][j] = (alp*alp-bet*bet)*mass[i][j]+alp*damp[i][j]+stiff[i][j];
						a[ib][row + j] = a[i][j];
						a[i][row + j] = -2 * alp*bet*mass[i][j] - bet*damp[i][j];
						a[ib][j] = 2 * alp*bet*mass[i][j] + bet*damp[i][j];
					}

				}
				//check array
				//for (i = 0; i < row * 2; i++)
				//{
				//	for (j = 0; j < row * 2; j++)
				//	{
				//		check[i][j] = a[i][j];
				//	}
				//}
				nn = 2 * row;
			}
		}
		else
		{
			icc = 2;//case:2nd complex pair
		}

		if (icc < 2)
		{
			eigvec[kset][0] = 1.;//dammy number
			if ((a[0][0] != 0) && (a[1][0] == 0))//pivot
			{
				for (j = 0; j < nn; j++)
				{
					a[1][j] = a[0][j];
				}
			}
			for (i = 1; i < nn; i++)
			{
				eigvec[kset][i] = -a[i][0];
				a[i][0] = 0;
			}

			//Gaussian elimination
			for (k = 1; k < nn; k++)
			{
				ic = 0;//flag for pivoting
				for (i = k; i < nn; i++)
				{
					if (a[i][k] != 0)//search pivoting point
					{
						ic = 1;
						kk = i;
						break;
					}
				}
				if (ic == 0)//case:error
				{
					return 999;
					break;
				}
				pivotvalue = a[kk][k];
				for (j = k; j < nn; j++) //pivot k rows to kk rows each other
				{
					b1 = a[kk][j];
					a[kk][j] = a[k][j];
					a[k][j] = b1 / pivotvalue;
				}
				b1 = eigvec[kset][kk];
				eigvec[kset][k] = b1 / pivotvalue;
				for (i = 1; i < nn; i++)//elimination without pivot row
				{
					if (i != k)
					{
						valueb = a[i][k];
						for (j = k; j < nn; j++)
						{
							a[i][j] -= valueb*a[k][j];
						}
						eigvec[kset][i] -= valueb*eigvec[kset][k];
						//double tmp1 = eigvec[kset][i];
					}
				}
			}//Gaussian elimination end
			for (i = 0; i < 2 * row; i++)
			{
				for (j = 0; j < 2 * row; j++)
				{
				//	check[i][j] = a[i][j];
				}
			}
			int x = 0;
		}
		if (icc == 2)//2nd complex pair
		{
			for (i = 0; i < nn; i++)
			{
				if (i < row)
				{
					eigvec[kset][i] = eigvec[kset - 1][i];
				}
				else
				{
					eigvec[kset][i] = -eigvec[kset - 1][i];
				}
			}
			icc = 0;
		}

		//normalize magunitude = 1
		mag = 0.;
		for (i = 0; i < nn; i++)
		{
			mag += eigvec[kset][i] * eigvec[kset][i];
			//double tmp0 = eigvec[kset][i];
		}
		mag = sqrt(mag);
		for (i = 0; i < nn; i++)
		{
			eigvec[kset][i] /= mag;
			//double tmp1 = eigvec[kset][i];
		}
	}
	for (i = 0; i < 2*row; i++)
	{
		for (j = 0; j < (2 * row); j++)
		{
			ceigvec[j][i] = eigvec[i][j];
			//double tmp = ceigvec[j][i];
		}
	}
	// 確保したメモリ領域の解放

		for (i = 0; i < 2 * row; i++)
		{
			delete[] a[i], eigvec[i];;
		}
		delete[] eigvec, a;

	return 0;
}
//μ^2*M*z+μ*C*z+K*z=0 
//M,C,K are real number.
int GCEIGQR(double**mass,double**damp,double**stiff,double**eigenvalue,double**eigenvector,int row){
	//6/23/2014
	//Generalized Complex Eigen Value
	int i, j;

	double **invM = new double *[row];
	double **invMK = new double *[row];
	double **invMC = new double *[row];

	double **O = new double *[2 * row];
	double **KC = new double *[2 * row];

	double **a = new double *[2 * row];//generalized matirx a
	double **h = new double *[2 * row];//Hessenberg
	double **r = new double *[2 * row];//right upper

	for (i = 0; i < row; i++)//列を作る
	{
		invM[i] = new double[row];
		invMK[i] = new double[row];
		invMC[i] = new double[row];
	}
	for (i = 0; i < 2 * row; i++)
	{
		O[i] = new double[row];
		KC[i] = new double[row];
		a[i] = new double[2 * row];
		h[i] = new double[2 * row];
		r[i] = new double[2 * row];
	}
	//0E
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row * 2; j++)
		{
			O[i][j] = 0.;
		}
	}
	for (i = 0; i < row; i++)
	{
		O[i][row + i] = 1.;
	}

	minv(mass, invM, row);	//M^-1
	mmult(invM, stiff, invMK, row, row, row);//M^-1*K
	mmult(invM, damp, invMC, row, row, row);//M^-1*C
	mtimes(invMK, -1, invMK, row, row);//-M^-1*K
	mtimes(invMC, -1, invMC, row, row);//-M^-1*C
	mmerge_col(invMK, invMC, KC, row, row, row);//-M^-1*K,-M^-1*C
	mmerge_row(O, KC, a, row, row, 2 * row);//a matrix
	mhouse(a, h, 2 * row);//Hessenberg
	mqr(h, r, 2 * row);//QR
	mevqr(r, eigenvalue, 2 * row);//complex Eigen value
	mgevecqr(mass, damp, stiff, eigenvalue, eigenvector, row);//complex Eigen vector

	//メモリの解放
	for (i = 0; i < row; i++)
	{
		delete[] invM[i], invMK[i], invMC[i];
	}
	delete[] invM, invMC, invMK;
	for (i = 0; i < row * 2; i++)
	{
		delete[] O[i], KC[i], a[i], h[i], r[i];
	}
	delete[] O, KC, a, h, r;
	return 0;
}
int GREIGJ(double**mass, double**stiff, double*eigenvalue, double**eigenvector, int row){
	//6/23/2014
	//Generalized Real Eigen Value
	int i, j;
	double mag;
	double **U = new double *[row];
	double **invU = new double *[row];
	double **a = new double *[row];
	double **y = new double *[row];

	for (i = 0; i < row; i++)//列を作る
	{
		U[i] = new double[row];
		invU[i] = new double[row];
		a[i] = new double[row];
		y[i] = new double[row];
	}

	mchol(mass, U, row);
	minv(U, invU, row);
	msimi(stiff, invU, a,row,row);
	meigj(a, eigenvalue, y, row);
	mmult(invU, y, eigenvector, row,row,row);

	//omega^ => omega
	/*for (i = 0; i < row; i++)
	{
		eigenvalue[i] = sqrt(eigenvalue[i]);
	}*/

	//normalized vector(magnitude=1)
	for (i = 0; i < row; i++)
	{
		mag = 0.;
		for (j = 0; j < row; j++)
		{
			mag += eigenvector[j][i] * eigenvector[j][i];
		}
		mag = sqrt(mag);
		for (j = 0; j < row; j++)
		{
			eigenvector[j][i] /= mag;
		}
	}

	for (i = 0; i < row; i++)
	{
		delete[] U[i], invU[i],a[i],y[i];
	}
	delete[] U, invU,a,y;
	return 0;
}
int msplit(double**inmat, double**outmat,int start_row, int start_col, int row_length, int col_length){
	//6/23/2014
	//left upper (row,col) 
	int i, j;
	for (i = 0; i <  row_length; i++)
	{
		for (j =0; j <  col_length; j++)
		{
			outmat[i][j] = inmat[i + start_row][j + start_col];
		}
	}
	return 0;
}
int mgmat(double**mass, double**damp, double**stiff, double**outmat,double omega,int row){
	//6/23/2014
	//Dynamic stiffness(complex)
	int i;
	double **M = new double *[row];
	double **C = new double *[row];
	for (i = 0; i < row; i++)//列を作る
	{
		M[i] = new double[row];
		C[i] = new double[row];
	}
	mtimes(mass, (-1*omega*omega), M, row, row);
	mtimes(damp, omega, C, row, row);
	madd(M, stiff, M, row, row);
	mcomp(M, C, outmat, row, row);
	for (i = 0; i < row; i++)
	{
		delete[] M[i], C[i];
	}
	delete[] M, C;
	return 0;
}
//Pseudo inverse matrix
//NEED row_a > col_a 
int mpinv(double**a, double**pinv_a, int row_a, int col_a){
	//6/25/2014
	//Pseudo inverse matrix
	int i;
	//dynamic array
	//rows
	double**aT = new double*[col_a];
	double**inv_aTa = new double*[col_a];
	double**aTa = new double*[col_a];
	//columns
	for (i = 0; i < col_a; i++)
	{
		aT[i] = new double[row_a];
		inv_aTa[i] = new double[col_a];
		aTa[i] = new double[col_a];
	}

	mtrans(a, aT, row_a, col_a);
	mmult(aT, a, inv_aTa, col_a, row_a, col_a);
	minv(inv_aTa, inv_aTa, col_a);
	mmult(inv_aTa, aT, pinv_a, col_a, col_a, row_a);

	//empty memory
	//rows
	for (i = 0; i < col_a; i++)
	{
		delete[] aT[i], inv_aTa[i],aTa[i];
	}
	//columns
	delete[] aT, inv_aTa,aTa[i];
	return 0;
}
//Make out Euler's rotation matrix around X axis 
int mrotX(double**rot, double theta_rad){
	//6/27/2014
	//Euler rotation around X axis
	rot[0][0] = 1.0; rot[0][1] = 0.0; rot[0][2] = 0.0;
	rot[1][0] = 0.0; rot[1][1] = cos(theta_rad); rot[1][2] = -sin(theta_rad);
	rot[2][0] = 0.0; rot[2][1] = sin(theta_rad); rot[2][2] = cos(theta_rad);
	return 0;
}
//Make out Euler's rotation matrix around Y axis 
int mrotY(double**rot, double theta_rad){
	//6/27/2014
	//Euler rotation around Y axis
	rot[0][0] = cos(theta_rad); rot[0][1] = 0.0; rot[0][2] = sin(theta_rad);
	rot[1][0] = 0.0; rot[1][1] = 1.0; rot[1][2] =0.0 ;
	rot[2][0] = -sin(theta_rad); rot[2][1] = 0.0; rot[2][2] = cos(theta_rad);
	return 0;
}
//Make out Euler's rotation matrix around Z axis 
int mrotZ(double**rot, double theta_rad){
	//6/27/2014
	//Euler rotation around Y axis
	rot[0][0] = cos(theta_rad); rot[0][1] = -sin(theta_rad); rot[0][2] = 0.0;
	rot[1][0] = sin(theta_rad); rot[1][1] = cos(theta_rad); rot[1][2] = 0.0;
	rot[2][0] = 0.0; rot[2][1] = 0.0; rot[2][2] = 1.0;
	return 0;
}
//convert degrees to radians
double radians(double deg){
	//6/27/2014
	//convert degrees to radians
	return deg*M_PI / 180.0;
}
//convert radians to degrees
double degrees(double rad){
	//6/27/2014
	//convert radians to degrees
	return rad * 180.0/ M_PI;
}
// [rot_y]*[rot_x]*[rot_y]*[rot_z]*[rot_y]*[rot_x]
//example1. rot_Z -> rot_X -> rot_Y:mrotyxyzyx(OutMatrix,radY,radX,0,radZ,0,0);
//example2. rot_X -> rot_Y -> rot_Z:mrotyxyzyx(OutMatrix,0,0,0,radZ,radY,radX);
int mrotyxyzyx(double **outmat,
	double radY3, double radX2, double radY2, double radZ1, double radY1, double radX1){
	//6/27/2014
	int i;
	//dynamic array
	//rows
	double**tmp1 = new double*[3];
	double**tmp2 = new double*[3];
	double**tmp3 = new double*[3];
	//columns
	for (i = 0; i < 3; i++)
	{
		tmp1[i] = new double[3];
		tmp2[i] = new double[3];
		tmp3[i] = new double[3];
	}
	mrotX(tmp1, radX1);
	mrotY(tmp2, radY1);
	mmult(tmp2, tmp1, tmp3,3, 3, 3);
	mrotZ(tmp1, radZ1);
	mmult(tmp1, tmp3, tmp2, 3, 3, 3);
	mrotY(tmp3, radY2);
	mmult(tmp3, tmp2, tmp1, 3, 3, 3);
	mrotX(tmp2, radX2);
	mmult(tmp2, tmp1, tmp3, 3, 3, 3);
	mrotY(tmp1, radY3);
	mmult(tmp1, tmp3, outmat, 3, 3, 3);

	//rows
	for (i = 0; i < 3; i++)
	{
		delete[] tmp1[i], tmp2[i],tmp3[i];
	}
	//columns
	delete[] tmp1, tmp2,tmp3;

	return 0;
}
//make Distance matrix(6*6)
int mDmat(double **Rotation_matrix, double **Distanse_Matrix, double X, double Y, double Z){
	//6/28/2014
	int i, j;
	//dynamic array
	//rows
	double **d_mat = new double*[3];
	double **rTd_mat = new double*[3];
	double **rT_mat = new double*[3];
	//columns
	for (i = 0; i < 3; i++)
	{
		d_mat[i] = new double[3];
		rTd_mat[i] = new double[3];
		rT_mat[i] = new double[3];
	}
	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			Distanse_Matrix[i][j] = 0.;
		}
	}
	for (i = 0; i < 3; i++)d_mat[i][i]=0.;
	d_mat[0][1] = Z;
	d_mat[0][2] = -Y;
	d_mat[1][2] = X;

	d_mat[1][0] = -Z;
	d_mat[2][0] = Y;
	d_mat[2][1] = -X;
/*	d_mat[0][1] = -Z;
	d_mat[0][2] = Y;
	d_mat[1][2] = -X;

	d_mat[1][0] = Z;
	d_mat[2][0] = -Y;
	d_mat[2][1] = X;
*/
	mtrans(Rotation_matrix, rT_mat, 3, 3);
	mmult(rT_mat, d_mat, rTd_mat, 3, 3, 3);
	
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Distanse_Matrix[i][j] = rT_mat[i][j];
			Distanse_Matrix[3+i][3+j] = rT_mat[i][j];
			//Distanse_Matrix[i][j] = Rotation_matrix[i][j];
			//Distanse_Matrix[3 + i][3 + j] = Rotation_matrix[i][j];
			Distanse_Matrix[i][3 + j] = rTd_mat[i][j];
		}
	}
	//rows
	for (i = 0; i < 3; i++)
	{
		delete[] rT_mat[i], d_mat[i], rTd_mat[i];
	}
	//columns
	delete[] rT_mat, d_mat, rTd_mat;

	return 0;

}
//num=G_matrix size(complex rows or complex columns)
//function for Gain calculation (Gain=x0/x1)
//when
//G_matrix * X_vector = Force_vector
//X_vector=(x0,x1)^t
//Gain_matrix:[num-G00num][G00num]
int mX0byX1(double **Gmat,double **Gain_mat,int G00num,int num){
	//6/28/2014
	int i,G11num;
	G11num = num - G00num;

	//行を作る
	double **G00 = new double *[G00num];
	double **invG00 = new double *[G00num];
	double **G01 = new double*[G00num];

	//列を作る
	for (i = 0; i < G00num; i++)
	{
		G00[i] = new double[G00num];
		invG00[i] = new double[G00num];
		G01[i] = new double[G11num];
	}

	//
	msplit(Gmat, G00, 0, 0, G00num,G00num);
	msplit(Gmat, G01, 0, G00num, G00num, G11num);
	minv(G00, invG00, G00num);
	mmult(invG00, G01, Gain_mat, G00num, G00num, G11num);
	mtimes(Gain_mat, -1, Gain_mat, G00num, G00num);

	for (i = 0; i < G00num; i++)
	{
		delete[] G00[i], invG00[i], G01[i];
	}
	delete[] G00, invG00, G01;

	return 0;

}
//Guyan's Reduction
//Gmat:				base matrix
//Greduction_mat:	matrix after Guyan's reduction
//Gred_num:			size of reduced matrix(rows or columns)
//Gnum:				size of base matrix 
int mGuyan(double **Gmat, double **Greduction_mat, int Gred_num, int Gnum){
	//6/29/2014
	int i, G11num;
	G11num = Gnum - Gred_num;

	//行を作る
	double **G00 = new double *[Gred_num];
	double **G11 = new double*[G11num];
	double **invG11 = new double *[G11num];
	double **G10 = new double*[G11num];
	double **tmat = new double *[Gred_num];
	
	//列を作る
	for (i = 0; i < Gred_num; i++)
	{
		G00[i] = new double[Gred_num];
		tmat[i] = new double[Gred_num];
	}
	for (i = 0; i < G11num; i++)
	{
		G11[i] = new double[G11num];
		invG11[i] = new double[G11num];
		G10[i] = new double[G11num];
	}

	//
	msplit(Gmat, G00, 0, 0, Gred_num, Gred_num);
	msplit(Gmat, G10, Gred_num, 0, G11num, Gred_num);
	msplit(Gmat, G11, Gred_num, Gred_num, G11num, G11num);
	minv(G11, invG11, G11num);
	msimi(invG11, G10, tmat, G11num, Gred_num);
	msub(G00, tmat, Greduction_mat, Gred_num, Gred_num);
	//mmult(invG00, G01, Gain_mat, G00num, G00num, G11num);

	for (i = 0; i < Gred_num; i++)
	{
		delete[] G00[i], tmat[i];
	}
	for (i = 0; i < G11num; i++)
	{
		delete[] G11[i], invG11[i], G10[i];
	}
	delete[] G00, tmat,G11,invG11, G10;

	return 0;

}
// transmissibility of vibration (F1/F0)
// notice::F0 is INPUT force.
//when
//G_matrix * X_vector = Force_vector
//F_vector=(F0,F1)^t
//trv_matrix:[num-G00num][G00num]
//num=G_matrix size(complex rows or complex columns)
int mTrv(double **Gmat, double **trv_mat, int G00num, int num){
	//6/28/2014
	int i, G11num;
	G11num = num - G00num;

	//行を作る
	double **G00 = new double *[G00num];
	double **invG00 = new double *[G00num];
	double **G10 = new double*[G11num];

	//列を作る
	for (i = 0; i < G00num; i++)
	{
		G00[i] = new double[G00num];
		invG00[i] = new double[G00num];
		
	}
	for (i = 0; i < G11num; i++)
	{
		G10[i] = new double[G00num];
	}

	//
	msplit(Gmat, G00, 0, 0, G00num, G00num);
	msplit(Gmat, G10, G00num, 0, G11num, G00num);
	minv(G00, invG00, G00num);
	mmult(G10, invG00, trv_mat, G00num, G00num, G11num);

	for (i = 0; i < G00num; i++)
	{
		delete[] G00[i], invG00[i], G10[i];
	}
	delete[] G00, invG00, G10;

	return 0;

}
//c=b*a*b^t
//a:[col][col]
//b:[row][col]
int msimi2(double **a, double **b, double **c, int row, int col){
	//6/29/2014
	int i, j, k;
	//行列の相似変換(similarity transformation)
	//b[row][col]*a[col][col]*b^t[col][row]
	// 配列動的確保
	double **abT = new double*[col];
	for (i = 0; i < col; i++)
	{
		abT[i] = new double[row];
	}

	{
		if (row < 1 || col < 1)return(999);
		//ab=a*b^t
		for (i = 0; i < col; i++)
		{
			for (j = 0; j < row; j++)
			{
				abT[i][j] = 0.0;
				for (k = 0; k < col; k++)
				{
					abT[i][j] += a[i][k] * b[j][k];
				}
			}
		}
	}
	//c=b*ab
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			c[i][j] = 0.0;
			for (k = 0; k < col; k++)
			{
				c[i][j] += b[i][k] * abT[k][j];
			}
		}
	}

	// 確保したメモリ領域の解放
	for (i = 0; i < col; i++)
	{
		delete[] abT[i];
	}
	delete[] abT;
	return 0;
}
//calculate elastic axis
//positon:[3][6] position of each axis center
//vector:[3][6] direction of each axis
//eig_k[6] principle value of each axis
int mElas_old(double **stiff, double **potision, double **vector,double *eig_k){
	//6/29/2014
	int i, j;
	double **K00 = new double*[3];//stiffness matrix for translation
	double **K11 = new double*[3];//stiffness matrix for rotation
	double **R00 = new double*[3];//rotation matrix for translation
	double **R11 = new double*[3];//rotation matrix for rotation
	double **R = new double*[6];
	double **tmp_pos = new double*[3];
	double **tmp_pos2 = new double*[3];
	double **RtKR = new double*[6];
	double *Ktr = new double[3];//translation principle stiffness 
	double *Krot = new double[3];//rotation principle stiffness
	for (i = 0; i < 3; i++)
	{
		K00[i] = new double[3];
		K11[i] = new double[3];
		R00[i] = new double[3];
		R11[i] = new double[3];
		tmp_pos[i] = new double[3];
		tmp_pos2[i] = new double[3];
	}
	for (i = 0; i < 6; i++)
	{
		R[i] = new double[6];
		RtKR[i] = new double[6];
	}
	///////translation
	msplit(stiff, K00, 0, 0, 3, 3);
	meigj(K00, Ktr, R00, 3 );
	for (i = 0; i < 3; i++)//make R mat
	{
		for (j = 0; j < 3; j++)
		{
			R[i][j] = R00[i][j];
			R[i + 3][j + 3] = R00[i][j];
			R[i][j + 3] = 0.;
			R[i + 3][j] = 0.;
		}
	}
	msimi(stiff, R, RtKR, 6, 6);
	tmp_pos[0][0] = RtKR[0][3] / RtKR[0][0];
	tmp_pos[1][0] = -1 * RtKR[0][5] / RtKR[0][0];
	tmp_pos[2][0] = RtKR[0][4] / RtKR[0][0];
	
	tmp_pos[0][1] = RtKR[1][5] / RtKR[1][1];
	tmp_pos[1][1] = RtKR[1][4] / RtKR[1][1];
	tmp_pos[2][1] = -1 * RtKR[1][3] / RtKR[1][1];

	tmp_pos[0][2] = -1 * RtKR[2][4] / RtKR[2][2];
	tmp_pos[1][2] = RtKR[2][3] / RtKR[2][2];
	tmp_pos[2][2] = RtKR[2][5] / RtKR[2][2];

	mmult(R00, tmp_pos, tmp_pos2, 3, 3, 3);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			potision[i][j] = tmp_pos2[i][j];
		}
	}

	///////rotation
	msplit(stiff, K11, 3, 3, 3, 3);
	meigj(K11, Krot, R11, 3);
	for (i = 0; i < 3; i++)//make R mat
	{
		for (j = 0; j < 3; j++)
		{
			R[i][j] = R11[i][j];
			R[i + 3][j + 3] = R11[i][j];
			R[i][j + 3] = 0.;
			R[i + 3][j] = 0.;
		}
	}
	msimi(stiff, R, RtKR, 6, 6);
	tmp_pos[0][0] = RtKR[0][3] / RtKR[3][3];
	tmp_pos[1][0] = -1 * RtKR[2][3] / RtKR[3][3];
	tmp_pos[2][0] = RtKR[1][3] / RtKR[3][3];

	tmp_pos[0][1] = RtKR[2][4] / RtKR[4][4];
	tmp_pos[1][1] = RtKR[1][4] / RtKR[4][4];
	tmp_pos[2][1] = -1 * RtKR[0][4] / RtKR[4][4];

	tmp_pos[0][2] = -1 * RtKR[1][5] / RtKR[5][5];
	tmp_pos[1][2] = RtKR[0][5] / RtKR[5][5];
	tmp_pos[2][2] = RtKR[2][5] / RtKR[5][5];
	mmult(R11, tmp_pos, tmp_pos2, 3, 3, 3);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			potision[i][j+3] = tmp_pos2[i][j];
		}
	}

	////copy
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			vector[i][j] = R00[i][j];
			vector[i][j + 3] = R11[i][j];
		}
	}
	for (i = 0; i < 3; i++)
	{
		eig_k[i] = Ktr[i];
		eig_k[i + 3] = Krot[i];
	}


	for (i = 0; i < 3; i++)
	{
		delete[] K00[i],K11[i], R00[i], R11[i],tmp_pos[i],tmp_pos2[i];
	}
	for (i = 0; i < 6; i++)
	{
		delete[] R[i], RtKR[i];
	}
	delete[] K00,K11, R00,R11,tmp_pos,tmp_pos2, R, RtKR,Ktr,Krot;

	return 0;
}
//calculate strain enegy contribution
//U=1/2*k*x^2
int mSEC(double **stiff, double **eigen_vector, double **contribution_matrix, int size_of_number){
	//7/2/2014
	int i, j;
	double **vecPow2 = new double *[size_of_number];
	double **tmp_mat = new double*[size_of_number];
	for (i = 0; i < size_of_number; i++)
	{
		vecPow2[i] = new double[size_of_number];
		tmp_mat[i] = new double[size_of_number];
	}

	//make eigen_vector^2
	for (i = 0; i < size_of_number; i++)
	{
		for (j = 0; j < size_of_number; j++)
		{
			vecPow2[i][j] = eigen_vector[i][j] * eigen_vector[i][j];
		}
	}
	mmult(stiff, vecPow2, tmp_mat, size_of_number,size_of_number,size_of_number);
	
	//normalize
	for (i = 0; i < size_of_number; i++)
	{
		double tmp = 0;
		for (j = 0; j < size_of_number; j++)
		{
			tmp += tmp_mat[j][i] * tmp_mat[j][i];
		}
		tmp = sqrt(tmp);
		for (j = 0; j < size_of_number; j++)
		{
			contribution_matrix[j][i] = fabs(tmp_mat[j][i] / tmp);
		}
	}

	//clear memory
	for (i = 0; i < size_of_number; i++)
	{
		delete[] vecPow2[i],tmp_mat[i];
	}
	delete[] vecPow2,tmp_mat;
	return 0;
}
//calculate kinetic enegy contribution
//U=1/2*m*x^2
int mKEC(double **mass, double **eigen_vector, double **contribution_matrix, int size_of_number){
	//7/2/2014
	int i, j;
	double **vecPow2 = new double *[size_of_number];
	double **tmp_mat = new double*[size_of_number];
	for (i = 0; i < size_of_number; i++)
	{
		vecPow2[i] = new double[size_of_number];
		tmp_mat[i] = new double[size_of_number];
	}

	//make eigen_vector^2
	for (i = 0; i < size_of_number; i++)
	{
		for (j = 0; j < size_of_number; j++)
		{
			vecPow2[i][j] = eigen_vector[i][j] * eigen_vector[i][j];
		}
	}
	mmult(mass, vecPow2, tmp_mat, size_of_number, size_of_number, size_of_number);

	//normalize
	for (i = 0; i < size_of_number; i++)
	{
		double tmp = 0;
		for (j = 0; j < size_of_number; j++)
		{
			tmp += tmp_mat[j][i] * tmp_mat[j][i];
		}
		tmp = sqrt(tmp);
		for (j = 0; j < size_of_number; j++)
		{
			contribution_matrix[j][i] = fabs(tmp_mat[j][i] / tmp);
		}
	}

	//clear memory
	for (i = 0; i < size_of_number; i++)
	{
		delete[] vecPow2[i], tmp_mat[i];
	}
	delete[] vecPow2, tmp_mat;
	return 0;
}
//Kx=λMx
//K is complex number.
int GCEIGQR2(double**Mass, double**Stiff_Re, double**Stiff_Im,
	double**eig_value, double**eig_vector, int Mass_size){
	int i, j,flag;
	double **comp_mass = new double*[Mass_size * 2];
	double **comp_stiff = new double*[Mass_size * 2];
	double **comp_U = new double*[Mass_size * 2];
	double **invU = new double*[Mass_size * 2];
	double **a = new double*[Mass_size * 2];
	double **y = new double*[Mass_size * 2];
	double **tmpVec = new double*[Mass_size * 2];
	double **dammy_mass = new double*[Mass_size];
	for (i = 0; i < Mass_size * 2; i++)
	{
		comp_mass[i] = new double[Mass_size * 2];
		comp_stiff[i] = new double[Mass_size * 2];
		comp_U[i] = new double[Mass_size * 2];
		invU[i] = new double[Mass_size * 2];
		a[i] = new double[Mass_size * 2];
		y[i] = new double[Mass_size * 2];
		tmpVec[i] = new double[Mass_size * 2];
	}
	for (i = 0; i < Mass_size; i++)
	{
		dammy_mass[i] = new double[Mass_size];
	}
	//0
	for (i = 0; i < Mass_size; i++)
	{
		for (j = 0; j < Mass_size; j++)
		{
			dammy_mass[i][j] = 0.0;
		}
	}
	double chk[12][12];
	mcomp(Mass, dammy_mass, comp_mass, Mass_size, Mass_size);
	mcomp(Stiff_Re, Stiff_Im, comp_stiff, Mass_size, Mass_size);
	mchol(comp_mass, comp_U, Mass_size * 2);

	minv(comp_U, invU, Mass_size*2);
	msimi(comp_stiff, invU, a, Mass_size * 2, Mass_size * 2);
	flag=mceig2(a, eig_value, y, Mass_size * 2);////y -> eig_vector
	for (i = 0; i < Mass_size * 2; i++)
	{
		double r = sqrt(eig_value[0][i] * eig_value[0][i]
			+ eig_value[1][i] * eig_value[1][i]);
		double theta = 0.5*atan2(eig_value[1][i] , eig_value[0][i]);
		
		eig_value[0][i] = sqrt(r)*cos(theta)/2/M_PI;
		eig_value[1][i] =sqrt(r)*sin(theta);
	}

	mmult(invU, y, tmpVec, Mass_size * 2, Mass_size * 2, Mass_size*2);
//	//normalize
	for (i = 0; i < Mass_size*2; i++)
	{
		double tmp = 0.;
		for (j = 0; j < Mass_size * 2; j++)
		{
			tmp += tmpVec[j][i] * tmpVec[j][i];
		}
		tmp = sqrt(tmp);
		for (j = 0; j < Mass_size * 2; j++)
		{
			tmpVec[j][i] /= tmp;
		}
	}
	//change order
	for (i = Mass_size - 1; i >= 0; i--)
	{
		int k = Mass_size - i - 1;
		for (j = 0; j<Mass_size*2 ; j++)
		{
			eig_vector[j][i * 2] = tmpVec[j][k * 2];
			eig_vector[j][i * 2 + 1] = tmpVec[j][k * 2 + 1];
		}
	}
	/*
	for (i = 0; i < Mass_size * 2;i++)
	{
		for (j = 0; j < Mass_size * 2; j++)
		{
			eig_vector[i][j] = a[i][j];
		}
	}
	*/

	//clear memory
	for (i = 0; i < 12; i++)
	{
		for (j = 0; j < 12; j++)
		{
			chk[i][j] = comp_U[i][j];
		}
	}

	for (i = 0; i < Mass_size * 2; i++)
	{
		delete[] comp_mass[i], comp_stiff[i], comp_U[i],invU[i],a[i],y[i],tmpVec[i];
	}
	for (i = 0; i < Mass_size; i++)
	{
		delete[] dammy_mass[i];
	}

	delete[] comp_mass, comp_stiff, comp_U,invU,a,y,tmpVec, dammy_mass;

	return flag;
}
//Householder transformation
//with Q vector
int mhouse2(double** inmat, double** outmat,double **Qmat, int row){
	//7/5/2014
	//Householder transformation
	int i, j, k, m;
	double b1, b2, b3, sigma, c1;
	// 配列動的確保
	double **qmat = new double*[row];
	double **buf = new double*[row];
	double **tempQ = new double*[row];
	double *uvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		buf[i] = new double[row];
		tempQ[i] = new double[row];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			outmat[i][j] = inmat[i][j];
			buf[i][j] = 0.;
			tempQ[i][j] = 0.;
		}
		tempQ[i][i] = 1.;
	}
	//
	for (k = 0; k < row - 2; k++)
	{
		//vector u
		b1 = 0.;
		for (i = k + 1; i < row; i++)
		{
			b1 += outmat[i][k] * outmat[i][k];
		}
		b2 = sqrt(b1);
		b3 = outmat[k + 1][k];
		if (b3 >= 0)
		{
			sigma = b2;
		}
		else
		{
			sigma = -b2;
		}

		uvec[k + 1] = b3 + sigma;
		c1 = sigma*uvec[k + 1];
		uvec[k] = 0;

		for (i = k + 2; i < row; i++)////
		{
			uvec[i] = outmat[i][k];
		}

		//Matrix Q
		//initialize
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0.;
			}
		}
		for (i = 0; i < row; i++)
		{
			qmat[i][i] = 1.;
		}

		for (i = k; i < row; i++)
		{
			for (j = k; j < row; j++)
			{
				if (c1 != 0)
				{
					qmat[i][j] -= uvec[i] * uvec[j] / c1;
				}
			}
		}

		//QA calculation
		for (i = 0; i <= k; i++)////
		{
			for (j = k; j < row; j++)
			{
				buf[i][j] = outmat[i][j];
				//double tmp = outmat[i][j];
			}
		}
		for (i = k + 1; i < row; i++)
		{
			for (j = k; j < row; j++)
			{
				b1 = 0.;
				for (m = k; m < row; m++)
				{
					b1 += qmat[i][m] * outmat[m][j];
					//double tmp = qmat[i][m];
				}
				buf[i][j] = b1;
			}
		}
		//QAQ^t
		for (j = 0; j < row; j++)
		{
			for (i = k; i < row; i++)
			{
				outmat[i][j] = buf[i][j];
			}
		}
		for (j = k + 1; j < row; j++)
		{
			for (i = 0; i < row; i++)
			{
				b1 = 0.;
				for (m = k; m < row; m++)
				{
					b1 += buf[i][m] * qmat[m][j];
				}
				outmat[i][j] = b1;
			}
		}
		//Qmat(k)=Qmat(k-1)*qmat
		//Qmat(k-1)=tempQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				Qmat[i][j] = 0.;
				for (m = 0; m < row; m++)
				{
					Qmat[i][j] += tempQ[i][m] * qmat[m][j];
				}
			}
		}
		//tempQ=Q
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				tempQ[i][j] = Qmat[i][j];
			}
		}

	}

	//lower left =0
	for (i = 2; i < row; i++)
	{
		for (j = 0; j <= i - 2; j++)
		{
			outmat[i][j] = 0.;
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] buf[i], qmat[i],tempQ[i];
	}
	delete[] buf, qmat, uvec,tempQ;

	return 0;
}
//QR factorization
//with Q vector
int mqr2(double** inmat, double** outmat, double**Qmat,int row){
	//7/5/2014
	//pukiwiki for PBCG Lab 参照
	int i, j, k, m, cnt, itermax;
	double c, s, b1, rq, eps, er;
	// 配列動的確保
	double **qmat = new double*[row];
	double **rmat = new double*[row];
	double **tmat = new double*[row];
	double **tempQ = new double*[row];
	double *uvec = new double[row];
	double *vvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		rmat[i] = new double[row];
		tmat[i] = new double[row];
		tempQ[i] = new double[row];
	}

	//copy Hessenberg matrix to rmat & outmat
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			rmat[i][j] = inmat[i][j];
			outmat[i][j] = inmat[i][j];
			tmat[i][j] = 0.;
			tempQ[i][j] = 0.;
		}
		tempQ[i][i] = 1;
	}

	itermax = 5000;//iteration maximum number
	eps = 1.0e-8;

	//
	for (cnt = 0; cnt < itermax; cnt++)
	{
		//double tmp = rmat[0][1];
		//Q <- unit matrix
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0.;
			}
		}
		for (i = 0; i < row; i++)
		{
			qmat[i][i] = 1.;
		}
		for (k = 0; k < row - 1; k++)
		{
			//sine and cosine
			//double tmp1 = rmat[k][k];
			//double tmp2 = rmat[k + 1][k];
			b1 = sqrt(rmat[k][k] * rmat[k][k] + rmat[k + 1][k] * rmat[k + 1][k]);
			c = rmat[k][k] / b1;
			s = -rmat[k + 1][k] / b1;

			//R calculation
			for (j = k + 1; j < row; j++)
			{
				uvec[j] = c*rmat[k][j] - s*rmat[k + 1][j];
				vvec[j] = s*rmat[k][j] + c*rmat[k + 1][j];
			}

			rmat[k][k] = b1;
			rmat[k + 1][k] = 0.;

			//double tmp0 =rmat[0][1];

			for (j = k + 1; j < row; j++)
			{
				rmat[k][j] = uvec[j];
				rmat[k + 1][j] = vvec[j];
			}

			//Q calculation
			for (j = 0; j <= k; j++)//
			{
				uvec[j] = c*qmat[k][j];
				vvec[j] = s*qmat[k][j];
			}
			qmat[k][k + 1] = -s;
			qmat[k + 1][k + 1] = c;
			for (j = 0; j <= k; j++)
			{
				qmat[k][j] = uvec[j];
				qmat[k + 1][j] = vvec[j];
			}
		}

		//A=RQ (outmat=RQ)
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rq = 0.;
				for (m = 0; m < row; m++)
				{
					rq += rmat[i][m] * qmat[j][m];
					//double tmp1 = rmat[i][m];
					//double tmp2 = qmat[j][m];
				}
				tmat[i][j] = rq;
			}
		}

		//convergence check
		er = 0;
		for (i = 0; i < row; i++)
		{
			er += fabs(tmat[i][i] - inmat[i][i]);
		}
		if (er < eps)
		{
			break;
		}
		//R,A <-RQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rmat[i][j] = tmat[i][j];
				outmat[i][j] = tmat[i][j];
			}
		}

		//Qmat(cnt)=Qmat(cnt-1)*qmat
		//Qmat(cnt-1)=tempQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				Qmat[i][j] = 0;
				for (k = 0; k < row; k++)
				{
					Qmat[i][j] += tempQ[i][k] * qmat[j][k];
				}
			}
		}
		//tempQ=Qmat
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				tempQ[i][j] = Qmat[i][j];
			}
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] rmat[i], qmat[i], tmat[i],tempQ[i];
	}
	delete[] rmat, qmat, tmat, uvec, vvec,tempQ;

	return cnt;
}
//Standard complex eigenvalu problem
//when solving complex stiffness
//eigvec=[row][row];using Householder's Q vector and QR's Q vector
int mceig2(double** inmat, double** eig, double** eigvec, int row){
	//7/5/2014
	//Standard complex eigenvalu problem
	//int flag;
	// 配列動的確保
	double **hessenberg = new double*[row];
	double **rightupper = new double*[row];
	double **Qmat1 = new double*[row];
	double **Qmat2 = new double*[row];
	for (int i = 0; i < row; i++)
	{
		hessenberg[i] = new double[row];
		rightupper[i] = new double[row];
		Qmat1[i] = new double[row];
		Qmat2[i] = new double[row];
	}
	mhouse2(inmat, hessenberg, Qmat1,row);//ヘッセンベルグ化
	int cnt=mqr4(hessenberg, rightupper,Qmat2, row,5000);//右上三角行列化（部分的に副対角項有り）
	mevqr(rightupper, eig, row);//複素固有値算出
	//flag = mevecqr(inmat, eig, eigvec, row);//複素固有ベクトル算出
	mmult(Qmat1, Qmat2, eigvec, row, row, row);//複素固有ベクトル算出
	// 確保したメモリ領域の解放
	for (int i = 0; i < row; i++)
	{
		delete[] hessenberg[i], rightupper[i],Qmat1[i],Qmat2[i];
	}
	delete[] hessenberg, rightupper,Qmat1,Qmat2;
	return cnt;
}
//how many numbers to omit rows(columns)
//and make matrix which is omited rows and columns
//for square matrix
int momit(double **inmat,int row){
	//7/10/2014
	//
	int i, j,k;
/*	double **tmp = new double*[row];
	for (i = 0; i < row; i++)
	{
		tmp[i] = new double[row];
	}*/
	int omit_times = 0;
	for (i = 0; i < row; i++)
	{
		if (inmat[i][i] == 0)
		{
			omit_times += 1;

		}//endif
	}


	double tmp;
	for (k = 0; k < row; k++)
	{


		for (int loop = 0; loop < omit_times; loop++)
		{
			if (inmat[k][k] == 0)
			{
				for (i = k; i < row - 1; i++)
				{
					//exchange rows 
					for (j = 0; j < row; j++)
					{
						tmp = inmat[i][j];
						inmat[i][j] = inmat[i + 1][j];
						inmat[i + 1][j] = tmp;
					}
					//exchange columns
					for (j = 0; j < row; j++)
					{
						tmp = inmat[j][i];
						inmat[j][i] = inmat[j][i + 1];
						inmat[j][i + 1] = tmp;
					}
				}
			}
		} 

	}
	return omit_times;
}
//how many numbers to omit rows(columns)
//and make matrix which is omited rows and columns
//for square matrix
//omit condition:mass[i][i]&&stiff[i][i]==0
int momit2(double **mass,double**damp,double**stiff,double**stiff_im, int row){
	//7/10/2014
	//
	int i, j, k;
	int omit_times = 0;
	for (i = 0; i < row; i++)
	{
		if (mass[i][i] == 0 || stiff[i][i]==0)
		{
			omit_times += 1;

		}//endif
	}


	double tmp;
	for (k = 0; k < row; k++)
	{
		for (int loop = 0; loop < omit_times; loop++)
		{
			if (mass[k][k] == 0 || stiff[k][k]==0)
			{
				for (i = k; i < row - 1; i++)
				{
					//exchange rows 
					for (j = 0; j < row; j++)
					{
						tmp = mass[i][j];
						mass[i][j] = mass[i + 1][j];
						mass[i + 1][j] = tmp;
						tmp = damp[i][j];
						damp[i][j] = damp[i + 1][j];
						damp[i + 1][j] = tmp;
						tmp = stiff[i][j];
						stiff[i][j] = stiff[i + 1][j];
						stiff[i + 1][j] = tmp;
						tmp = stiff_im[i][j];
						stiff_im[i][j] = stiff_im[i + 1][j];
						stiff_im[i + 1][j] = tmp;
					}
					//exchange columns
					for (j = 0; j < row; j++)
					{
						tmp = mass[j][i];
						mass[j][i] = mass[j][i + 1];
						mass[j][i + 1] = tmp;
						tmp = damp[j][i];
						damp[j][i] = damp[j][i + 1];
						damp[j][i + 1] = tmp;
						tmp = stiff[j][i];
						stiff[j][i] = stiff[j][i + 1];
						stiff[j][i + 1] = tmp;
						tmp = stiff_im[j][i];
						stiff_im[j][i] = stiff_im[j][i + 1];
						stiff_im[j][i + 1] = tmp;
					}
				}
			}
		}

	}
	return omit_times;
}
//QR factorization
int mqr3(double** inmat, double** outmat, int row,int iteration_num){
	//6/22/2014
	//pukiwiki for PBCG Lab 参照
	int i, j, k, m, cnt, itermax,elim_i;
	double c, s, b1, rq, eps;
	// 配列動的確保
	double **qmat = new double*[row];
	double **rmat = new double*[row];
	double **tmat = new double*[row];
	double *uvec = new double[row];
	double *vvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		rmat[i] = new double[row];
		tmat[i] = new double[row];
	}

	//copy Hessenberg matrix to rmat & outmat
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			rmat[i][j] = inmat[i][j];
			outmat[i][j] = inmat[i][j];
			tmat[i][j] = 0.;
		}
	}

	itermax = iteration_num;//iteration maximum number
	eps = 1.0e-16;
	elim_i = row;
	//
	for (cnt = 0; cnt < itermax; cnt++)
	{
		//double tmp = rmat[0][1];
		//Q <- unit matrix
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0.;
			}
		}
		for (i = 0; i < row; i++)
		{
			qmat[i][i] = 1.;
		}
		for (k = 0; k < elim_i - 1; k++)
		{
			//sine and cosine
			//double tmp1 = rmat[k][k];
			//double tmp2 = rmat[k + 1][k];
			b1 = sqrt(rmat[k][k] * rmat[k][k] + rmat[k + 1][k] * rmat[k + 1][k]);
			c = rmat[k][k] / b1;
			s = -rmat[k + 1][k] / b1;

			//R calculation
			for (j = k + 1; j < elim_i; j++)
			{
				uvec[j] = c*rmat[k][j] - s*rmat[k + 1][j];
				vvec[j] = s*rmat[k][j] + c*rmat[k + 1][j];
			}

			rmat[k][k] = b1;
			rmat[k + 1][k] = 0.;

			//double tmp0 =rmat[0][1];

			for (j = k + 1; j < elim_i; j++)
			{
				rmat[k][j] = uvec[j];
				rmat[k + 1][j] = vvec[j];
			}

			//Q calculation
			for (j = 0; j<elim_i; j++)//
			{
				uvec[j] = c*qmat[k][j];
				vvec[j] = s*qmat[k][j];
			}
			qmat[k][k + 1] = -s;
			qmat[k + 1][k + 1] = c;
			for (j = 0; j <= k; j++)
			{
				qmat[k][j] = uvec[j];
				qmat[k + 1][j] = vvec[j];
			}
		}

		//A=RQ (outmat=RQ)
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rq = 0.;
				for (m = 0; m < row; m++)
				{
					rq += rmat[i][m] * qmat[j][m];
					//double tmp1 = rmat[i][m];
					//double tmp2 = qmat[j][m];
				}
				tmat[i][j] = rq;
			}
		}


		//R,A <-RQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rmat[i][j] = tmat[i][j];
				outmat[i][j] = tmat[i][j];
			}
		}
		//elimination
		if (elim_i>1)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{
				elim_i -= 1;
			}
		}
		else if (elim_i>2)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{
				elim_i -= 2;
			}
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] rmat[i], qmat[i], tmat[i];
	}
	delete[] rmat, qmat, tmat, uvec, vvec;

	return cnt;
}
//QR factorization
//with Q vector
int mqr4(double** inmat, double** outmat, double**Qmat, int row,int iteration_num){
	//7/5/2014
	//pukiwiki for PBCG Lab 参照
	int i, j, k, m, cnt, itermax,elim_i;
	double c, s, b1, rq, eps;
	// 配列動的確保
	double **qmat = new double*[row];
	double **rmat = new double*[row];
	double **tmat = new double*[row];
	double **tempQ = new double*[row];
	double *uvec = new double[row];
	double *vvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		rmat[i] = new double[row];
		tmat[i] = new double[row];
		tempQ[i] = new double[row];
	}

	//copy Hessenberg matrix to rmat & outmat
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			rmat[i][j] = inmat[i][j];
			outmat[i][j] = inmat[i][j];
			tmat[i][j] = 0.;
			tempQ[i][j] = 0.;
		}
		tempQ[i][i] = 1;
	}

	itermax = iteration_num;//iteration maximum number
	eps = 1.0e-8;
	elim_i = row;
	//
	for (cnt = 0; cnt < itermax; cnt++)
	{
		//double tmp = rmat[0][1];
		//Q <- unit matrix
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0.;
			}
		}
		for (i = 0; i < row; i++)
		{
			qmat[i][i] = 1.;
		}
		for (k = 0; k < elim_i - 1; k++)
		{
			//sine and cosine
			//double tmp1 = rmat[k][k];
			//double tmp2 = rmat[k + 1][k];
			b1 = sqrt(rmat[k][k] * rmat[k][k] + rmat[k + 1][k] * rmat[k + 1][k]);
			c = rmat[k][k] / b1;
			s = -rmat[k + 1][k] / b1;

			//R calculation
			for (j = k + 1; j < elim_i; j++)
			{
				uvec[j] = c*rmat[k][j] - s*rmat[k + 1][j];
				vvec[j] = s*rmat[k][j] + c*rmat[k + 1][j];
			}

			rmat[k][k] = b1;
			rmat[k + 1][k] = 0.;

			//double tmp0 =rmat[0][1];

			for (j = k + 1; j < elim_i; j++)
			{
				rmat[k][j] = uvec[j];
				rmat[k + 1][j] = vvec[j];
			}

			//Q calculation
			for (j = 0; j <elim_i; j++)//
			{
				uvec[j] = c*qmat[k][j];
				vvec[j] = s*qmat[k][j];
			}
			qmat[k][k + 1] = -s;
			qmat[k + 1][k + 1] = c;
			for (j = 0; j <= k; j++)
			{
				qmat[k][j] = uvec[j];
				qmat[k + 1][j] = vvec[j];
			}
		}

		//A=RQ (outmat=RQ)
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rq = 0.;
				for (m = 0; m < row; m++)
				{
					rq += rmat[i][m] * qmat[j][m];
					//double tmp1 = rmat[i][m];
					//double tmp2 = qmat[j][m];
				}
				tmat[i][j] = rq;
			}
		}


		//R,A <-RQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				rmat[i][j] = tmat[i][j];
				outmat[i][j] = tmat[i][j];
			}
		}

		//Qmat(cnt)=Qmat(cnt-1)*qmat
		//Qmat(cnt-1)=tempQ
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				Qmat[i][j] = 0;
				for (k = 0; k < row; k++)
				{
					Qmat[i][j] += tempQ[i][k] * qmat[j][k];
				}
			}
		}
		//tempQ=Qmat
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				tempQ[i][j] = Qmat[i][j];
			}
		}
		//elimination
		if (elim_i>1)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{
				elim_i -= 1;
			}
		}
		else if (elim_i>2)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{
				elim_i -= 2;
			}
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] rmat[i], qmat[i], tmat[i], tempQ[i];
	}
	delete[] rmat, qmat, tmat, uvec, vvec, tempQ;

	return cnt;
}
//doubleQR with shift factorization
int mdqr(double** inmat, double** outmat, int row, int iteration_num){
	//7/16/2014

	int i, j, k, m, cnt, itermax,elim_i;
	double c, s, b1,b2, rq, eps,bsum,bpro;
	// 配列動的確保
	double **qmat = new double*[row];
	double **rmat = new double*[row];
	double **tmat = new double*[row];
	double *uvec = new double[row];
	double *vvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		rmat[i] = new double[row];
		tmat[i] = new double[row];
	}
	//double a[5][5],b[5][5];
	//copy Hessenberg matrix to rmat & outmat
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			//a[i][j] = inmat[i][j];
			outmat[i][j] = inmat[i][j];
			tmat[i][j] = 0.;
			rmat[i][j] = 0.;

		}
	}
	mmult(inmat, inmat, rmat, row, row, row);
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			//a[i][j] = rmat[i][j];
		}
	}


	itermax = iteration_num;//iteration maximum number
	eps = 1.0e-16;
	elim_i = row;
	bsum = 0;
	bpro = 0;
	//
	for (cnt = 0; cnt < itermax; cnt++)
	{
		//double tmp = rmat[0][1];
		//Q <- unit matrix
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < row; j++)
			{
				qmat[i][j] = 0;
				//b[i][j] = 0.;
			}
			qmat[i][i] = 1;
			//b[i][i] = 1.;
		}

		for (k = 0; k < elim_i - 1; k++)
		{
			if (k < elim_i - 2)//a[k+2][k]=0
			{
				//sine and cosine
				//double tmp1 = rmat[k][k];
				//double tmp2 = rmat[k + 1][k];
				b1 = sqrt(rmat[k][k] * rmat[k][k] + rmat[k + 2][k] * rmat[k + 2][k]);
				c = rmat[k][k] / b1;
				s = -rmat[k + 2][k] / b1;

				//R calculation
				for (j = k; j < elim_i; j++)
				{
					uvec[j] = c*rmat[k][j] - s*rmat[k + 2][j];
					vvec[j] = s*rmat[k][j] + c*rmat[k + 2][j];
				}

				rmat[k][k] = b1;
				rmat[k + 2][k] = 0.;
				//a[k][k] = b1;
				//a[k + 2][k] = 0.;

				//double tmp0 =rmat[0][1];

				for (j = k ; j < elim_i; j++)
				{
					rmat[k][j] = uvec[j];
					rmat[k + 2][j] = vvec[j];
					//a[k][j] = uvec[j];
					//a[k + 2][j] = vvec[j];
				}

				//Q^t calculation
				for (j = 0; j <elim_i; j++)//
				{
					b1 = qmat[k][j];
					b2 = qmat[k + 2][j];
					qmat[k][j] = c*b1-s*b2;
					qmat[k + 2][j] = s*b1+c*b2;
					//b[k][j] = c*b1 - s*b2;
					//b[k + 2][j] = s*b1 + c*b2;
				}
			}
			//a[k+1][k]=0

			//sine and cosine
			//double tmp1 = rmat[k][k];
			//double tmp2 = rmat[k + 1][k];
			b1 = sqrt(rmat[k][k] * rmat[k][k] + rmat[k + 1][k] * rmat[k + 1][k]);
			c = rmat[k][k] / b1;
			s = -rmat[k + 1][k] / b1;

			//R calculation
			for (j = k + 1; j < elim_i; j++)
			{
				uvec[j] = c*rmat[k][j] - s*rmat[k + 1][j];
				vvec[j] = s*rmat[k][j] + c*rmat[k + 1][j];
			}

			rmat[k][k] = b1;
			rmat[k + 1][k] = 0.;
			//a[k][k] = b1;
			//a[k + 1][k] = 0.;
			//double tmp0 =rmat[0][1];

			for (j = k + 1; j < elim_i; j++)
			{
				rmat[k][j] = uvec[j];
				rmat[k + 1][j] = vvec[j];
				//a[k][j] = uvec[j];
				//a[k + 1][j] = vvec[j];
			}

			//Q^2t calculation
			for (j = 0; j <elim_i; j++)//
			{
				b1 = qmat[k][j];
				b2 = qmat[k + 1][j];
				qmat[k][j] = c*b1-s*b2;
				qmat[k + 1][j] = s*b1+c*b2;
				//b[k][j] = c*b1 - s*b2;
				//b[k + 1][j] = s*b1 + c*b2;
			}
		}

		//A_cnt+1=Q^t*A_cnt*Q 
		for (i = 0; i < elim_i; i++)
		{
			for (j = 0; j < elim_i; j++)
			{
				rq = 0.;
				for (m = 0; m < elim_i; m++)
				{
					rq += outmat[i][m] * qmat[j][m];
					//double tmp1 = rmat[i][m];
					//double tmp2 = qmat[j][m];
				}
				tmat[i][j] = rq;
			}
		}
		for (i = 0; i < elim_i; i++)
		{
			for (j = 0; j < elim_i; j++)
			{
				rq = 0.;
				for (m = 0; m < elim_i; m++)
				{
					rq += qmat[i][m] * tmat[m][j];
				}
				outmat[i][j] = rq;
				//a[i][j] = outmat[i][j];
			}
		}

		//elimination
		if (elim_i>3)
		{
			if (fabs(outmat[elim_i-1][elim_i - 2]) < eps)
			{
				//bsum = 0;
				//bpro = outmat[elim_i-1][elim_i-1];
				bsum = outmat[elim_i - 3][elim_i - 3] + outmat[elim_i - 2][elim_i - 2];
				bpro = outmat[elim_i -3][elim_i -3] * outmat[elim_i - 2][elim_i - 2]
					- outmat[elim_i - 3][elim_i - 2] * outmat[elim_i - 2][elim_i - 3];
				elim_i -= 1;
			}
		}
		else if (elim_i>4)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{
				//bsum = outmat[elim_i - 2][elim_i - 2] + outmat[elim_i-1][elim_i-1];
				//bpro = outmat[elim_i - 2][elim_i - 2] * outmat[elim_i-1][elim_i-1]
				//	- outmat[elim_i - 2][elim_i-1] * outmat[elim_i-1][elim_i - 2];
				bsum = outmat[elim_i - 4][elim_i - 4] + outmat[elim_i - 3][elim_i - 3];
				bpro = outmat[elim_i - 4][elim_i - 4] * outmat[elim_i - 3][elim_i - 3]
					- outmat[elim_i - 4][elim_i - 3] * outmat[elim_i - 3][elim_i - 4];
				elim_i -= 2;
			}
		}
		/*else
		{
			bsum = 0;
			bpro = 0;
		}*/
		//R<-A^2-(μ+μ^*)A+(μμ^*)E
		for (i = 0; i < elim_i; i++)
		{
			for (j = 0; j < elim_i; j++)
			{
				rmat[i][j] = 0;
				for (m = 0; m < elim_i; m++)
				{
					rmat[i][j] += outmat[i][m]*outmat[m][j];
				}
				rmat[i][j] -= bsum*outmat[i][j];
			}
			rmat[i][i] += bpro;
		}
		for (i = 0; i < elim_i; i++)
		{
			for (j = 0; j < elim_i; j++)
			{
				//a[i][j] = rmat[i][j];
			}
		}

	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] rmat[i], qmat[i], tmat[i];
	}
	delete[] rmat, qmat, tmat, uvec, vvec;

	return cnt;
}
//sqrt
double mysqrt(double x){
	double s, last;
	if (x > 0)
	{
		if (x > 1)
		{
			s = x;
		}
		else
		{
			s = 1;
		}
		do
		{
			last = s;
			s = (x / s + s) / 2;
		} while (s < last);
		return last;
	}
	return 0;
}
int mhqr(double** inmat, double** outmat, int row,int iteration){
	//6/22/2014
	//Householder transformation
	int i, j, k, m,cnt,elim_i;
	double b1, b2, b3, sigma, c1,eps;
	// 配列動的確保
	double **qmat = new double*[row];
	double **h = new double*[row];
	double **buf = new double*[row];
	double *uvec = new double[row];
	for (i = 0; i < row; i++)
	{
		qmat[i] = new double[row];
		h[i] = new double[row];
		buf[i] = new double[row];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			h[i][j] = inmat[i][j];
			buf[i][j] = 0.;
		}
	}
	//
	eps = 1.0e-16;
	elim_i = row;
	for (cnt = 0; cnt < iteration; cnt++)
	{
		for (k = 0; k < elim_i - 2; k++)
		{
			//vector u
			b1 = 0.;
			for (i = k; i < elim_i; i++)
			{
				b1 += h[i][k] * h[i][k];
			}
			b2 = mysqrt(b1);
			b3 = h[k][k];
			if (b3 >= 0)
			{
				sigma = -b2;
			}
			else
			{
				sigma = b2;
			}

			uvec[k] = b3 - sigma;

			for (i = k + 1; i < elim_i; i++)////
			{
				uvec[i] = h[i][k];
			}
			c1 = 0;
			for (i = k; i < elim_i; i++)
			{
				c1 += uvec[i] * uvec[i];
			}
			if (c1 != 0)
			{
				c1 = mysqrt(c1);
				for (i = k; i < elim_i; i++)
				{
					uvec[i] /= c1;
				}
			}
			//Matrix Q
			//initialize
			for (i = 0; i < elim_i; i++)
			{
				for (j = 0; j < elim_i; j++)
				{
					qmat[i][j] = 0.;
				}
				qmat[i][i] = 1.;
			}

			for (i = k; i < elim_i; i++)
			{
				for (j = k; j < elim_i; j++)
				{
					qmat[i][j] -= 2*uvec[i] * uvec[j];

				}
			}

			//QAQ calculation
			//msimi(h, qmat, h, elim_i, elim_i);
			//QA calculation
			for (i = 0; i <= k; i++)////
			{
				for (j = k; j < elim_i; j++)
				{
					buf[i][j] = h[i][j];
					//double tmp = outmat[i][j];
				}
			}
			for (i = k ; i < elim_i; i++)
			{
				for (j = 0; j < elim_i; j++)
				{
					b1 = 0.;
					for (m = k; m < elim_i; m++)
					{
						b1 += qmat[i][m] * h[m][j];
						//double tmp = qmat[i][m];
					}
					buf[i][j] = b1;
				}
			}
			//QAQ^t
			for (j = 0; j < k; j++)
			{
				for (i = k; i < elim_i; i++)
				{
					h[i][j] = buf[i][j];
				}
			}
			for (j = k ; j < elim_i; j++)
			{
				for (i = 0; i < elim_i; i++)
				{
					b1 = 0.;
					for (m = k; m < elim_i; m++)
					{
						b1 += buf[i][m] * qmat[m][j];
					}
					h[i][j] = b1;
				}
			}

		}
		//elimination
		if (elim_i>1)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{

				elim_i -= 1;
			}
		}
		else if (elim_i>2)
		{
			if (fabs(outmat[elim_i - 2][elim_i - 3]) < eps)
			{

				elim_i -= 2;
			}
		}
		cnt++;
	}
	

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < row; j++)
		{
			outmat[i][j] = h[i][j];
		}
	}
	// 確保したメモリ領域の解放
	for (i = 0; i < row; i++)
	{
		delete[] buf[i], qmat[i], h[i];
	}
	delete[] buf, qmat, uvec, h;

	return 0;
}
//how many numbers to omit rows(columns)
//and make matrix which is omited rows and columns
//for square matrix
//omit condition:mass[i][i]&&stiff[i][i]==0
//with index
int momit3(double **mass, double**damp, double**stiff, double**stiff_im, int *index,int row){
	//7/18/2014
	//
	int i, j, k;
	int omit_times = 0;
	for (i = 0; i < row; i++)
	{
		if (mass[i][i] == 0 || stiff[i][i] == 0)
		{
			omit_times += 1;

		}//endif
	}
	for (i = 0; i < row; i++)
	{
		index[i] = i;
	}


	double tmp;
	for (k = 0; k < row; k++)
	{
		for (int loop = 0; loop < omit_times; loop++)
		{
			if (mass[k][k] == 0 || stiff[k][k] == 0)
			{
				for (i = k; i < row - 1; i++)
				{
					//exchange rows 
					for (j = 0; j < row; j++)
					{
						tmp = mass[i][j];
						mass[i][j] = mass[i + 1][j];
						mass[i + 1][j] = tmp;
						tmp = damp[i][j];
						damp[i][j] = damp[i + 1][j];
						damp[i + 1][j] = tmp;
						tmp = stiff[i][j];
						stiff[i][j] = stiff[i + 1][j];
						stiff[i + 1][j] = tmp;
						tmp = stiff_im[i][j];
						stiff_im[i][j] = stiff_im[i + 1][j];
						stiff_im[i + 1][j] = tmp;

					}
					//exchange columns
					for (j = 0; j < row; j++)
					{
						tmp = mass[j][i];
						mass[j][i] = mass[j][i + 1];
						mass[j][i + 1] = tmp;
						tmp = damp[j][i];
						damp[j][i] = damp[j][i + 1];
						damp[j][i + 1] = tmp;
						tmp = stiff[j][i];
						stiff[j][i] = stiff[j][i + 1];
						stiff[j][i + 1] = tmp;
						tmp = stiff_im[j][i];
						stiff_im[j][i] = stiff_im[j][i + 1];
						stiff_im[j][i + 1] = tmp;
					}
					int tmp_i = index[i];
					index[i] = index[i + 1];
					index[i + 1] = tmp_i;
				}
			}
		}

	}
	return omit_times;
}
//calculate elastic axis
//positon:[3][6] position of each axis center
//vector:[3][6] direction of each axis
//eig_k[6] principle value of each axis
int mElas(double **stiff, double **potision, double **vector, double *eig_k){
	//6/29/2014
	int i, j,k;
	double **K00 = new double*[3];//stiffness matrix for translation
	double **K11 = new double*[3];//stiffness matrix for rotation
	double **R00 = new double*[3];//rotation matrix for translation
	double **R11 = new double*[3];//rotation matrix for rotation
	double **D00 = new double*[3];//distance matrix for translation
	double **D11 = new double*[3];//distance matrix for rotation
	double **Kt_tr = new double*[3];//Kt for translation
	double **Kt_rot = new double*[3];//Kt for rotation
	double **invKt_tr = new double*[3];//inverse Kt for translation
	double **invKt_rot = new double*[3];//inverse Kt for translation
	double **R = new double*[6];
	double **tmp_pos = new double*[3];
	double **tmp_pos2 = new double*[3];
	double **RtKR = new double*[6];
	double *K_tr = new double[3];//translation principle stiffness 
	double *K_rot = new double[3];//rotation principle stiffness
	double **DtKTD = new double*[3];
	double **tmp_R = new double*[3];
	for (i = 0; i < 3; i++)
	{
		K00[i] = new double[3];
		K11[i] = new double[3];
		R00[i] = new double[3];
		R11[i] = new double[3];
		D00[i] = new double[3];
		D11[i] = new double[3];
		Kt_tr[i] = new double[3];
		Kt_rot[i] = new double[3];
		invKt_tr[i] = new double[3];
		invKt_rot[i] = new double[3];
		tmp_pos[i] = new double[3];
		tmp_pos2[i] = new double[3];
		DtKTD[i] = new double[3];
		tmp_R[i] = new double[3];

	}
	for (i = 0; i < 6; i++)
	{
		R[i] = new double[6];
		RtKR[i] = new double[6];
	}
	///////translation
	msplit(stiff, K00, 0, 0, 3, 3);
	meigj(K00, K_tr, R00, 3);//並進方向固有ベクトル
	for (i = 0; i < 3; i++)//make R mat
	{
		for (j = 0; j < 3; j++)
		{
			R[i][j] = R00[i][j];
			R[i + 3][j + 3] = R00[i][j];
			R[i][j + 3] = 0.;
			R[i + 3][j] = 0.;
		}
	}
	msimi(stiff, R, RtKR, 6, 6);//並進軸上にＫを座標変換
	tmp_pos[0][0] = RtKR[0][3] / RtKR[0][0];//距離マトリクス
	tmp_pos[1][0] = -1 * RtKR[0][5] / RtKR[0][0];
	tmp_pos[2][0] = RtKR[0][4] / RtKR[0][0];

	tmp_pos[0][1] = RtKR[1][5] / RtKR[1][1];
	tmp_pos[1][1] = RtKR[1][4] / RtKR[1][1];
	tmp_pos[2][1] = -1 * RtKR[1][3] / RtKR[1][1];

	tmp_pos[0][2] = -1 * RtKR[2][4] / RtKR[2][2];
	tmp_pos[1][2] = RtKR[2][3] / RtKR[2][2];
	tmp_pos[2][2] = RtKR[2][5] / RtKR[2][2];

	mmult(R00, tmp_pos, tmp_pos2, 3, 3, 3);//グローバル座標変換
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			potision[i][j] = tmp_pos2[i][j]*1000;
		}
	}

	///////rotation
	mdiag(K_tr, Kt_tr, 3);//固有値を対角行列化
	minv(Kt_tr, invKt_tr, 3);//並進方向固有値マトリクスの逆行列
	for (i = 0; i < 3; i++)//距離行列を求める
	{
		for (j = 0; j < 3; j++)
		{
			D00[i][j] = 0;
			for (k = 0; k < 3; k++)
			{
				D00[i][j] += invKt_tr[i][k] * RtKR[k][j + 3];
			}
		}
	}
	msimi(Kt_tr, D00, DtKTD, 3, 3);//並行軸の定理　並進成分の影響分
	for (i = 0; i < 3; i++)//並進成分の影響を除いた回転成分
	{
		for (j = 0; j < 3; j++)
		{
			K11[i][j] = RtKR[i + 3][j + 3] - DtKTD[i][j];
		}
	}
	
	meigj(K11, K_rot, R11, 3);//回転成分の固有値
	for (i = 0; i < 3; i++)//make R mat
	{
		for (j = 0; j < 3; j++)
		{
			R[i][j] = R11[i][j];
			R[i + 3][j + 3] = R11[i][j];
			R[i][j + 3] = 0.;
			R[i + 3][j] = 0.;
		}
	}
	msimi(stiff, R, RtKR, 6, 6);
	tmp_pos[0][0] = RtKR[0][3] / RtKR[3][3];
	tmp_pos[1][0] = -1 * RtKR[2][3] / RtKR[3][3];
	tmp_pos[2][0] = RtKR[1][3] / RtKR[3][3];

	tmp_pos[0][1] = RtKR[2][4] / RtKR[4][4];
	tmp_pos[1][1] = RtKR[1][4] / RtKR[4][4];
	tmp_pos[2][1] = -1 * RtKR[0][4] / RtKR[4][4];

	tmp_pos[0][2] = -1 * RtKR[1][5] / RtKR[5][5];
	tmp_pos[1][2] = RtKR[0][5] / RtKR[5][5];
	tmp_pos[2][2] = RtKR[2][5] / RtKR[5][5];
	mmult(R11, tmp_pos, tmp_pos2, 3, 3, 3);//グローバル座標変換
	mmult(R00, tmp_pos2, tmp_pos, 3, 3, 3);//グローバル座標変換
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			potision[i][j + 3] = tmp_pos[i][j]*1000;
		}
	}

	////copy
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			vector[i][j] = R00[i][j];
			vector[i][j + 3] = R11[i][j];
		}
	}
	for (i = 0; i < 3; i++)
	{
		eig_k[i] = K_tr[i]/1000;
		eig_k[i + 3] = K_rot[i]*M_PI/180;
	}


	for (i = 0; i < 3; i++)
	{
		delete[] K00[i], K11[i], R00[i], R11[i], tmp_pos[i], tmp_pos2[i],DtKTD[i];
		delete[] D00[i] ,D11[i],Kt_tr[i],Kt_rot[i],invKt_tr[i],invKt_rot[i],tmp_R[i];
	}
	for (i = 0; i < 6; i++)
	{
		delete[] R[i], RtKR[i];
	}
	delete[] K00, K11, R00, R11, tmp_pos, tmp_pos2, R, RtKR, K_tr, K_rot,DtKTD;
	delete[] D00, D11, Kt_tr, Kt_rot, invKt_tr, invKt_rot,tmp_R;
	return 0;
}
//Reduction
//Gmat:				base matrix
//Greduction_mat:	matrix after Guyan's reduction
//Gred_num:			size of reduced matrix(rows or columns)
//Gnum:				size of base matrix 
int mred(double **Gmat, double **Greduction_mat, int Gred_num, int Gnum){
	//6/29/2014
	int i, G11num;
	G11num = Gnum - Gred_num;

	//行を作る
	double **invG = new double *[Gnum];

	//列を作る
	for (i = 0; i < Gnum; i++)
	{
		invG[i] = new double[Gnum];
	}

	//
	minv(Gmat, invG, Gnum);
	minv(invG, Greduction_mat, G11num);

	for (i = 0; i < Gnum; i++)
	{
		delete[] invG[i];
	}
	delete[] invG;

	return 0;

}
int vproduct(double*a, double*b, double*ab, int num){
	int i;
	double *tmp_a = new double[num * 2];
	double *tmp_b = new double[num * 2];

	for (i = 0; i < num; i++)
	{
		//a
		tmp_a[i] = a[i];
		tmp_a[i + num] = a[i];
		//b
		tmp_b[i] = b[i];
		tmp_b[i + num] = b[i];
	}

	for (i = 0; i < num; i++)
	{
		ab[i] = tmp_a[i + 1] * tmp_b[i + 2] - tmp_b[i + 1] * tmp_a[i + 2];
	}
	return  0;

	delete[] tmp_a,tmp_b;
}

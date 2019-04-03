#ifndef _MY_MATH_H_
#define _MY_MATH_H_



int maxIndex(const  double *x,int width);


double Norm( const double *x,int width);


void copyVec(double *desc,int dInc, const double *src,int sInc,int size);


void subVec(const double *op1,const double *op2,double *result,int size);


void MatMultiVec(double *A, double* x, double* y, int row, int col);


void Submat(double *x, double *x_sub, int m, int n);


void MatColCopyMat(double *A,double *B,int *index,int num,int m,int n);


void ppp(double* a,double* e,double* s,double* v,int m,int n);
void sss(double fg[2],double cs[2]);

/*SVD*/
int bmuav(double* a,int m,int n,double* u,double* v,double eps,int ka);


int bginv(double* a,int m,int n,double* aa,double eps,double* u,double* v,int ka);


void Sort(double *a,int *index,int n);


void BubbleSort(int *a,int num);


int Check(int *a,int num);

#endif
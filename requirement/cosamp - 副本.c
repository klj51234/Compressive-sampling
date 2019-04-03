#include"cosamp.h"
#include"mymath.h"
#include"stdlib.h"
#include"stdio.h" /*chen*/
/*measure:观测向量，recSignal:待恢复信号，m:观测数，n:原始信号长度，error:控制迭代停止的相对残差，k:稀疏度 */
void cosamp(double *measure, double *recSignal, double *A, int m, int n, int k)/* double err,*/
{
	double *rError,*h,*temp1,*temp2,*temp3,*temp4,*temp5,*temp6,*subMat_A,eps, temp3_temp;
	int num,ka,*support,alpha,*h_index,*a_index,*a_index2,i,support_num,result;
	int j=0; /* chen */
	
	result = 0;
	alpha = 2;
	support_num = k + alpha*k;

	support = (int *)malloc(sizeof(int) * m);
	h_index = (int *)malloc(sizeof(int) * n);
	a_index = (int *)malloc(sizeof(int) * n);
	a_index2 = (int *)malloc(sizeof(int) * k);

	rError = (double *)malloc(m * sizeof(double));
	subMat_A = (double *)malloc(n * m * sizeof(double) );
	h = (double *)malloc(n * sizeof(double));
	temp1 = (double *)malloc(m * sizeof(double));
	temp2 = (double *)malloc(m * sizeof(double));

	for (i=0;i<n;i++)
		recSignal[i] = 0.0;

	eps = 0.000001;		//求伪逆的误差上限
	ka= m+1;

	if (n>m)
		ka = n+1;

	copyVec(rError,1,measure,1,m);
	num = 0;

	Submat(A,subMat_A,m,n);
	while (j<40)
	{
		MatMultiVec(subMat_A,rError,h,n,m);
		
		Sort(h,h_index,n);
		
		if ( num == 0 )
		{
			for (i=0;i<alpha*k;i++)
				support[i] = h_index[i];
			support_num = alpha*k;
			BubbleSort(support,alpha*k);
		}
		else
		{
			for (i=0;i<k;i++)
				support[i] = a_index2[i];
			for (i=0;i<alpha*k;i++)
				support[k+i] = h_index[i];
			support_num = k + alpha*k;
			BubbleSort(support,(alpha+1)*k);
		}

		result = Check(support,support_num);
		support_num = support_num - result;
		num++;

		temp3 = (double *)malloc(support_num * m * sizeof(double));	
		temp4 = (double *)malloc(support_num * m * sizeof(double));
		temp5 = (double *)malloc(m * m * sizeof(double));
		temp6 = (double *)malloc(support_num * support_num * sizeof(double));
		
		MatColCopyMat(A,temp3,support,support_num,m,n);
		bginv(temp3, m, support_num,temp4, eps, temp5, temp6, ka);
		MatMultiVec(temp4,measure , temp1, support_num, m);
		
		if (num==1)
		{
			for (i=0;i<support_num;i++)
				recSignal[support[i]] = temp1[i];
		}
		else
		{
			for (i=0;i<support_num;i++)
				recSignal[support[i]] = temp1[i];
		}
		Sort(recSignal,a_index,n);
		for (i=0;i<n-k;i++)
			recSignal[a_index[i+k]] = 0;

		for (i=0;i<k;i++)
			a_index2[i] = a_index[i];
		
		MatMultiVec(A,recSignal ,temp2, m, n);	
		
		subVec(measure,temp2,rError,m);	

		free(temp3);
		free(temp4);
		free(temp5);
		free(temp6);
		j++; /*Chen*/
	}
}



#ifndef _COSAMP_H_
#define _COSAMP_H_

/*measure:观测向量，recSignal:待恢复信号，m:观测数，n:原始信号长度，error:控制迭代停止的相对残差,k:稀疏度*/
void cosamp(double *measure, double *recSignal, double *A, int width, int height, int k);/*double error,*/


#endif
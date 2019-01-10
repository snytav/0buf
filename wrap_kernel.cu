/*
 * wrap_kernel.cu
 *
 *  Created on: Jul 23, 2018
 *      Author: snytav
 */







__global__
void GPU_getCellEnergy(
		GPUCell **cells,double *d_ee,
		double *d_Ex,double *d_Ey,double *d_Ez)
{
	unsigned int i = blockIdx.x;
	unsigned int l= blockIdx.y;
	unsigned int k = blockIdx.z;
	//int i,l,k;
	Cell *c0 = cells[0],nc;
	double t,ex,ey,ez;
//	__shared__ extern CellDouble fd[9];
	//double *src;//,*dst;
	int n  = c0->getGlobalCellNumber(i,l,k);

	ex = d_Ex[n];
	ey = d_Ey[n];
	ez = d_Ez[n];

	t = ex*ex+ey*ey+ez*ez;

	cuda_atomicAdd(d_ee,t);
}



/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include "archAPI.h"

int getRank()
{
    int rank;


    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    return rank;
}

int getSize()
{
    int rank;


    MPI_Comm_size(MPI_COMM_WORLD,&rank);

    return rank;
}

int InitMPI(int argc,char *argv[])
{
	MPI_Init(&argc,&argv);
}

int sumMPI(int size,double *d_jx,double *d_jy,double *d_jz)
{
    double *snd,*rcv,*jx,*jy,*jz;
    int i;
    
    if(getSize() == 1) return 0;

    jx = (double *)malloc(size*sizeof(double));
    jy = (double *)malloc(size*sizeof(double));
    jz = (double *)malloc(size*sizeof(double));

    int err  = MemoryCopy(jx,d_jx,size*sizeof(double),DEVICE_TO_HOST);
    int err1 = MemoryCopy(jy,d_jy,size*sizeof(double),DEVICE_TO_HOST);
    int err2 = MemoryCopy(jz,d_jz,size*sizeof(double),DEVICE_TO_HOST);
    printf("sumMPI: err %d err1 %d err2 %d \n",err,err1,err2);

    snd = (double *)malloc(3*size*sizeof(double));
    rcv = (double *)malloc(3*size*sizeof(double));
    
    for(i = 0;i < size;i++)
    {   
        snd[i]          = jx[i];
        snd[i +   size] = jy[i];
        snd[i + 2*size] = jz[i];
    }

    MPI_Allreduce(snd,rcv,size,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
    
    for(i = 0;i < size;i++)
    {   
        jx[i] = rcv[i];
        jy[i] = rcv[i +   size];
        jz[i] = rcv[i + 2*size];
    }

    err  = MemoryCopy(d_jx,jx,size*sizeof(double),HOST_TO_DEVICE);
    err1 = MemoryCopy(d_jy,jy,size*sizeof(double),HOST_TO_DEVICE);
    err2 = MemoryCopy(d_jz,jz,size*sizeof(double),HOST_TO_DEVICE);
    printf("sumMPI-after: err %d err1 %d err2 %d \n",err,err1,err2);


    return 0;
}

int sumMPIenergy(double *e)
{
    double snd,rcv;
    
    snd = *e;

    MPI_Allreduce(&snd,&rcv,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
    
    *e = rcv;
    
    return 0;
}



int  CloseMPI()
{
	MPI_Finalize();
}

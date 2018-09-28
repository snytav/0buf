/*
 * gpu_plasma.h
 *
 *  Created on: Aug 21, 2013
 *      Author: snytav
 */


#ifndef GPU_PLASMA_H_
#define GPU_PLASMA_H_

#include<stdlib.h>
#include<stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <string>

#include "load_data.h"

//#include <unistd.h>
//#include <stdio.h>
#include <errno.h>

#ifdef __CUDACC__
#include <nvToolsExtCuda.h>
#include <nvToolsExtCudaRt.h>
#endif

#include "maxwell.h"

#include <time.h>

//#ifdef __OMP__
#include <omp.h>
//#endif

#ifdef __CUDACC__
#include <cuda.h>
#endif

#include "archAPI.h"
#include "maxwell.h"
#include "plasma.h"
#include "gpucell.h"
#include "mpi_shortcut.h"

#include "service_functions.h"

#include <sys/resource.h>
#include <stdint.h>

#include <sys/sysinfo.h>
#include <sys/time.h>


#include "init.h"
#include "diagnose.h"

#include<string>
#include <iostream>

#include "particle_target.h"

#include "params.h"

#include "memory_control.h"

#include "add.cu"
#include "wrap_kernel.cu"
//#include "run_control.h"
#include "kernels.cu"

#include <vector>



using namespace std;






#define FORTRAN_ORDER







class Plasma
{
public:
   GPUCell **h_CellArray,**d_CellArray;
   GPUCell **cp;
   thrust::device_vector<GPUCell> *d_AllCells;
   double *d_Ex,*d_Ey,*d_Ez,*d_Hx,*d_Hy,*d_Hz,*d_Jx,*d_Jy,*d_Jz,*d_Rho,*d_npJx,*d_npJy,*d_npJz;
   double *d_Qx,*d_Qy,*d_Qz;
   double *dbg_x,*dbg_y,*dbg_z,*dbg_px,*dbg_py,*dbg_pz;
   int total_particles;

   int h_controlParticleNumberArray[4000];

   int  jx_wrong_points_number;
   int3 *jx_wrong_points,*d_jx_wrong_points;

//#ifdef ATTRIBUTES_CHECK
   double *ctrlParticles,*d_ctrlParticles,*check_ctrlParticles;
//#endif

   int jmp,size_ctrlParticles;
   double ami,amb,amf;
   int real_number_of_particle[3];
   FILE *f_prec_report;



   int CPU_field;




void copyCells(std::string where,int nt)
{
	static int first = 1;
	size_t m_free,m_total;
	int size = (*AllCells).size();
	struct sysinfo info;
//	unsigned long c1,c2;

    if(first == 1)
    {
    	cp = (GPUCell **)malloc(size*sizeof(GPUCell *));
    }

	unsigned long m1,m2,delta,accum;
	memory_monitor("beforeCopyCells",nt);

	for(int i = 0;i < size;i++)
	{
//		if(i == 141)
//		{
//			int qq = 0;
//		}
		cudaError_t err = cudaMemGetInfo(&m_free,&m_total);
		sysinfo(&info);
		m1 = info.freeram;
	 	GPUCell c,*d_c,*z0;
	 	z0 = h_CellArray[i];
	 	if(first == 1)
	 	{
	       d_c = c.allocateCopyCellFromDevice();
     	   cp[i] = d_c;
	 	}
	 	else
	 	{
	 	   d_c = cp[i];
	 	}
	    c.copyCellFromDevice(z0,d_c,where,nt);
		m2 = info.freeram;

	    delta = m1-m2;
        accum += delta;

	}

	if(first == 1)
	{
		first = 0;
	}

	memory_monitor("afterCopyCells",nt);
}


double checkGPUArray(double *a,double *d_a,std::string name,std::string where,int nt)
{
	 static double *t;
	 static int f1 = 1;
	 char fname[1000];
	 double res;
//	 FILE *f;
#ifndef CHECK_ARRAY_OUTPUT
   return 0.0;
#endif

	 sprintf(fname,"diff_%s_at_%s_nt%08d.dat",name.c_str(),where.c_str(),nt);


	 if(f1 == 1)
	 {
		 t = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
		 f1 = 0;
	 }
	 int err;
	 err = MemoryCopy(t,d_a,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),DEVICE_TO_HOST);
	 if(err != cudaSuccess)
	         {
	          	printf("bCheckArray err %d %s \n",err,getErrorString(err));
	        	exit(0);
	         }

#ifdef CHECK_ARRAY_DETAIL_PRINTS
	 if((f = fopen(fname,"wt")) != NULL)
	 {
		 res = CheckArray(a,t,f);
		 fclose(f);
	 }
#else
//	 res = CheckArray(a,t,f);
	 int size = (Nx+2)*(Ny+2)*(Nz+2);
	 res = CheckArraySilent(a,t,size);
#endif

	 return res;
}

double checkGPUArray(double *a,double *d_a)
{
	 static double *t;
	 static int f = 1;

	 if(f == 1)
	 {
		 t = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
		 f = 0;
	 }
	 int err;
	 err = MemoryCopy(t,d_a,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),DEVICE_TO_HOST);
	 if(err != cudaSuccess)
	         {
	          	printf("bCheckArray err %d %s \n",err,getErrorString(err));
	        	exit(0);
	         }

	 return CheckArray(a,t);

}


public:

void virtual emeGPUIterate(int3 s,int3 f,
			double *E,double *H1, double *H2,
			double *J,double c1,double c2, double tau,
			int3 d1,int3 d2)
{
	dim3 dimGrid(f.x-s.x+1,1,1),dimBlock(1,f.y-s.y+1,f.z-s.z+1);

//    GPU_eme//<<<dimGrid,dimBlock>>>(d_CellArray,s,
//    		                            E,H1,H2,
//    					    	  		J,c1,c2,tau,
//    					    	  		d1,d2
//    		);
    void* args[] = { (void* )&d_CellArray,
            		           (void *)&s,
            		           (void *)&E,
            		           (void *)&H1,
            		           (void *)&H2,
            		           (void *)&J,
                               (void *)&c1,
                               (void *)&c2,
                               (void *)&tau,
                               (void *)&d1,
                               (void *)&d2,
            		           0};

              cudaError_t cudaStatus = cudaLaunchKernel(
              	            (const void*)GPU_eme, // pointer to kernel func.
              	            dimGrid,                       // grid
              	            dimBlock,                      // block
              	            args,                          // arguments
              	            0,
              	            0
              	            );


}

void GetElectricFieldStartsDirs(
		int3 *start,
		int3 *d1,
		int3 *d2,
		int dir
		)
{
      start->x = (dir == 0)*0 + (dir == 1)*1 + (dir == 2)*1;
      start->y  = (dir == 0)*1 + (dir == 1)*0 + (dir == 2)*1;
      start->z = (dir == 0)*1 + (dir == 1)*1 + (dir == 2)*0;

      d1->x = (dir == 0)*0    + (dir == 1)*0    + (dir == 2)*(-1);
      d1->y = (dir == 0)*(-1) + (dir == 1)*0    + (dir == 2)*0;
      d1->z = (dir == 0)*0    + (dir == 1)*(-1) + (dir == 2)*0;

      d2->x = (dir == 0)*0    + (dir == 1)*(-1) + (dir == 2)*0;
      d2->y = (dir == 0)*0    + (dir == 1)*0    + (dir == 2)*(-1);
      d2->z = (dir == 0)*(-1) + (dir == 1)*0    + (dir == 2)*0;
}

int virtual ElectricFieldTrace(
  double *E,double *H1,double *H2,double *J,
  int dir,double c1,double c2,double tau)
  {
      int3 start,d1,d2,finish = make_int3(Nx,Ny,Nz);


      GetElectricFieldStartsDirs(&start,&d1,&d2,dir);

         emeGPUIterate(start,finish,E,H1,H2,
    	        		    	  		J,c1,c2,tau,
    	        		    	  		d1,d2);

    return 0;
  }

int checkFields_beforeMagneticStageOne(double *t_Ex,double *t_Ey,double *t_Ez,
		                               double *t_Hx,double *t_Hy,double *t_Hz,
		                               double *t_Qx,double *t_Qy,double *t_Qz,
		                               double *t_check,int nt)
{

	 memory_monitor("beforeComputeField_FirstHalfStep",nt);

//	         t_check[0] = checkControlMatrix("emh1",nt,"qx",t_Qx);
//			 t_check[1] = checkControlMatrix("emh1",nt,"qy",t_Qy);
//			 t_check[2] = checkControlMatrix("emh1",nt,"qz",t_Qz);
//
//			 t_check[3] = checkControlMatrix("emh1",nt,"ex",t_Ex);
//			 t_check[4] = checkControlMatrix("emh1",nt,"ey",t_Ey);
//			 t_check[5] = checkControlMatrix("emh1",nt,"ez",t_Ez);
//
//			 t_check[6] = checkControlMatrix("emh1",nt,"hx",t_Hx);
//			 t_check[7] = checkControlMatrix("emh1",nt,"hy",t_Hy);
//			 t_check[8] = checkControlMatrix("emh1",nt,"hz",t_Hz);
	return 0;
}

int checkFields_afterMagneticStageOne(double *t_Hx,double *t_Hy,double *t_Hz,
		                              double *t_Qx,double *t_Qy,double *t_Qz,
		                              double *t_check,int nt)
{
//	         t_check[9] = checkControlMatrix("emh1",nt,"qx",t_Qx);
//			 t_check[10] = checkControlMatrix("emh1",nt,"qy",t_Qy);
//			 t_check[11] = checkControlMatrix("emh1",nt,"qz",t_Qz);
//
//			 t_check[12] = checkControlMatrix("emh1",nt,"hx",t_Hx);
//			 t_check[13] = checkControlMatrix("emh1",nt,"hy",t_Hy);
//			 t_check[14] = checkControlMatrix("emh1",nt,"hz",t_Hz);


			 CPU_field = 1;



			 checkControlPoint(50,nt,0);
			 memory_monitor("afterComputeField_FirstHalfStep",nt);

	return 0;
}

void  ComputeField_FirstHalfStep(
		   int nt
		   )
{
	 double t_check[15];
    cudaError_t err;
	err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

     checkFields_beforeMagneticStageOne(d_Ex,d_Ey,d_Ez,
		 		                               d_Hx,d_Hy,d_Hz,
		 		                               d_Qx,d_Qy,d_Qz,
		 		                               t_check,nt);


	err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	 MagneticStageOne(d_Qx,d_Qy,d_Qz,d_Hx,d_Hy,d_Hz,d_Ex,d_Ey,d_Ez);

	err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	 checkFields_afterMagneticStageOne(d_Hx,d_Hy,d_Hz,
		 		                           d_Qx,d_Qy,d_Qz,
		 		                           t_check,nt);
	err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	 AssignCellsToArraysGPU();
	err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

}

virtual void ComputeField_SecondHalfStep(int nt)
{

     SetPeriodicCurrents(nt);



     MagneticFieldStageTwo(d_Hx,d_Hy,d_Hz,nt,d_Qx,d_Qy,d_Qz);




     ElectricFieldEvaluate(d_Ex,d_Ey,d_Ez,nt,d_Hx,d_Hy,d_Hz,d_Jx,d_Jy,d_Jz);


}

void ElectricFieldComponentEvaluateTrace(
		  double *E,double *H1,double *H2,double *J,
		  int dir,
		  double c1,double c2,double tau,
		  int dir_1,int start1_1,int end1_1,int start2_1,int end2_1,int N_1,
		  int dir_2,int start1_2,int end1_2,int start2_2,int end2_2,int N_2
		  )

{
     ElectricFieldTrace(E,H1,H2,J,dir,c1,c2,tau);

//     PeriodicBoundaries(E, dir_1,start1_1,end1_1,start2_1,end2_1,N_1);
//     PeriodicBoundaries(E, dir_2,start1_2,end1_2,start2_2,end2_2,N_2);
}


void ElectricFieldComponentEvaluatePeriodic(
		  double *E,double *H1,double *H2,double *J,
		  int dir,
		  double c1,double c2,double tau,
		  int dir_1,int start1_1,int end1_1,int start2_1,int end2_1,int N_1,
		  int dir_2,int start1_2,int end1_2,int start2_2,int end2_2,int N_2
		  )

{
//     ElectricFieldTrace(E,H1,H2,J,dir,c1,c2,tau);

     if (dir != 0)
     {
    	 PeriodicBoundaries(E, dir_1,start1_1,end1_1,start2_1,end2_1,N_1);
     }

     if (dir != 2)
     {
    	 PeriodicBoundaries(E, dir_2,start1_2,end1_2,start2_2,end2_2,N_2);
     }
}


void ElectricFieldEvaluate(double *locEx,double *locEy,double *locEz,
		   int nt,
		   double *locHx,double *locHy,double *locHz,
		   double *loc_npJx,double *loc_npJy,double *loc_npJz)
{
	 CPU_field = 0;
      double3 c1 = getMagneticFieldTimeMeshFactors();

      ElectricFieldComponentEvaluateTrace(
    		  locEx,locHz,locHy,loc_npJx,
    		  0,c1.y,c1.z,tau,
              1,0,Nx,1,Nz,Ny,
              2,0,Nx,0,Ny+1,Nz);

      ElectricFieldComponentEvaluateTrace(
          		  locEy,locHx,locHz,loc_npJy,
          		  1,c1.z,c1.x,tau,
                    0,0,Ny,1,Nz,Nx,
                    2,0,Nx+1,0,Ny,Nz);



      ElectricFieldComponentEvaluateTrace(
         		  locEz,locHy,locHx,loc_npJz,
         		  2,c1.x,c1.y,tau,
                   0,1,Ny,0,Nz,Nx,
                   1,0,Nx+1,0,Nz,Ny);

      checkControlPoint(550,nt,0);


      ElectricFieldComponentEvaluatePeriodic(
     		  locEx,locHz,locHy,loc_npJx,
     		  0,c1.y,c1.z,tau,
               1,0,Nx,1,Nz,Ny,
               2,0,Nx,0,Ny+1,Nz);

      ElectricFieldComponentEvaluatePeriodic(
    		  locEy,locHx,locHz,loc_npJy,
    		  1,c1.z,c1.x,tau,
              0,0,Ny,1,Nz,Nx,
              2,0,Nx+1,0,Ny,Nz);

//      SinglePeriodicBoundary(locEy,1,0,Nx+1,0,Nz+1,Ny);


      ElectricFieldComponentEvaluatePeriodic(
        		  locEz,locHy,locHx,loc_npJz,
        		  2,c1.x,c1.y,tau,
                  0,1,Ny,0,Nz,Nx,
                  1,0,Nx+1,0,Nz,Ny);

         checkControlPoint(600,nt,0);

         memory_monitor("after_ComputeField_SecondHalfStep",nt);
}

double3 getMagneticFieldTimeMeshFactors()
{
    Cell c = (*AllCells)[0];
	double hx = c.get_hx(),hy = c.get_hy(),hz = c.get_hz();
	double3 d;
    d.x = tau/(hx);
    d.y = tau/(hy);
    d.z = tau/hz;

	return d;
}

virtual void MagneticStageOne(
                  double *Qx,double *Qy,double *Qz,
                  double *Hx,double *Hy,double *Hz,
   	              double *Ex,double *Ey,double *Ez
		           )
{
	double3 c1 = getMagneticFieldTimeMeshFactors();

    MagneticFieldTrace(Qx,Hx,Ey,Ez,Nx+1,Ny,Nz,c1.z,c1.y,0);

    MagneticFieldTrace(Qy,Hy,Ez,Ex,Nx,Ny+1,Nz,c1.x,c1.z,1);

    MagneticFieldTrace(Qz,Hz,Ex,Ey,Nx,Ny,Nz+1,c1.y,c1.x,2);

}

virtual void MagneticFieldStageTwo(double *Hx,double *Hy,double *Hz,
		            int nt,
		            double *Qx,double *Qy,double *Qz)
{
    Cell c = (*AllCells)[0];

    SimpleMagneticFieldTrace(c,Qx,Hx,Nx+1,Ny,Nz);
    SimpleMagneticFieldTrace(c,Qy,Hy,Nx,Ny+1,Nz);
    SimpleMagneticFieldTrace(c,Qz,Hz,Nx,Ny,Nz+1);

    checkControlPoint(500,nt,0);
}


int PushParticles(int nt)
{
	double mass = -1.0/1836.0,q_mass = -1.0;

	memory_monitor("before_CellOrder_StepAllCells",nt);

    CellOrder_StepAllCells(nt,	mass,q_mass,1);
    puts("cell_order");

	memory_monitor("after_CellOrder_StepAllCells",nt);

	//checkParticleAttributes(nt);

	checkControlPoint(270,nt,1);

	return 0;
}


int readStartPoint(int nt)
{
	char fn[100];

	 if(nt == START_STEP_NUMBER)
	 {
		 readControlPoint(NULL,fn,0,nt,0,1,Ex,Ey,Ez,Hx,Hy,Hz,Jx,Jy,Jz,Qx,Qy,Qz,
							   dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz
			 );

		copyFieldsToGPU(
		    		                        d_Ex,d_Ey,d_Ez,
		    								d_Hx,d_Hy,d_Hz,
		    								d_Jx,d_Jy,d_Jz,
		    								d_npJx,d_npJy,d_npJz,
		    								d_Qx,d_Qy,d_Qz,
		    								Ex,Ey,Ez,
		    				        		Hx,Hy,Hz,
		    				        		Jx,Jy,Jz,
		    				        		npJx,npJy,npJz,
		    				                Qx,Qy,Qz,
		    				                Nx,Ny,Nz
		    		);
	 }

     checkControlPoint(0,nt,1);

	 return 0;
}



	void Step(int nt)
	 {
		ComputeField_FirstHalfStep(nt);

		PushParticles(nt);

		puts("push ended");

		ComputeField_SecondHalfStep(nt);
		puts("field computed-2");

		sumMPI((Nx+2)*(Ny+2)*(Nz+2),d_Jx,d_Jy,d_Jz);

		 Diagnose(nt);

	 }
	virtual double getElectricEnergy()
	{
		dim3 dimGrid(Nx+2,Ny+2,Nz+2),dimGridOne(1,1,1),dimBlock(MAX_particles_per_cell/2,1,1),
    		 dimBlockOne(1,1,1),dimBlockGrow(1,1,1),dimBlockExt(CellExtent,CellExtent,CellExtent);
		static int first = 1;
		static double *d_ee;
		double ee;

		if(first == 1)
		{
			cudaMalloc(&d_ee,sizeof(double));
			first = 0;
		}
		cudaMemset(d_ee,0,sizeof(double));

//		GPU_getCellEnergy//<<<dimGrid, dimBlockOne,16000>>>(d_CellArray,d_ee,d_Ex,d_Ey,d_Ez);
		void* args[] = { (void* )&d_CellArray,
		        		           (void *)&d_ee,
		        		           (void *)&d_Ex,
		        		           (void *)&d_Ey,
		        		           (void *)&d_Ez,
		        		           0};
		          cudaError_t cudaStatus = cudaLaunchKernel(
		          	            (const void*)GPU_getCellEnergy, // pointer to kernel func.
		          	            dimGrid,                       // grid
		          	            dimBlockOne,                      // block
		          	            args,                          // arguments
		          	            0,
		          	            0
		          	            );


        MemoryCopy(&ee,d_ee,sizeof(double),DEVICE_TO_HOST);

        return ee;

	}
	void Diagnose(int nt)
	{
		double ee;
		static FILE *f;
		static int first = 1;

		if(first == 1)
		{
			f = fopen("energy.dat","wt");
			first = 0;
		}
		else
		{
			f = fopen("energy.dat","at");

		}



        ee = getElectricEnergy();
       // sumMPIenergy(&ee);

        //if(getRank() == 0)
        	fprintf(f,"%10d %25.15e \n",nt,ee);

        fclose(f);
		//puts("GPU-Plasma");

	}
	virtual ~Plasma(){
		//~Plasma<Cell>();
		}

	  int Nx,Ny,Nz;

	  int n_per_cell;

	  int meh;

	  int magf;

	  double ion_q_m,tau;

	  double Lx,Ly,Lz;

	  double ni;

	  double *Qx,*Qy,*Qz,*dbg_Qx,*dbg_Qy,*dbg_Qz;

	  double *Ex,*Ey,*Ez,*Hx,*Hy,*Hz,*Jx,*Jy,*Jz,*Rho,*npJx,*npJy,*npJz;
	  double *dbgEx,*dbgEy,*dbgEz,*dbgHx,*dbgHy,*dbgHz,*dbgJx,*dbgJy,*dbgJz;
	  double *dbgEx0,*dbgEy0,*dbgEz0;
	  double *npEx,*npEy,*npEz;



	  std::vector<GPUCell> *AllCells;

	  int getBoundaryLimit(int dir){return ((dir == 0)*Nx  + (dir == 1)*Ny + (dir == 2)*Nz + 2);}

#include "init.cu"



	int getMagneticFieldTraceShifts(int dir,int3 &d1,int3 &d2)
	{
	      d1.x = (dir == 0)*0 + (dir == 1)*1 + (dir == 2)*0;
	      d1.y = (dir == 0)*0 + (dir == 1)*0 + (dir == 2)*1;
	      d1.z = (dir == 0)*1 + (dir == 1)*0 + (dir == 2)*0;

	      d2.x = (dir == 0)*0 + (dir == 1)*0 + (dir == 2)*1;
	      d2.y = (dir == 0)*1 + (dir == 1)*0 + (dir == 2)*0;
	      d2.z = (dir == 0)*0 + (dir == 1)*1 + (dir == 2)*0;

	      return 0;
	}

	int MagneticFieldTrace(double *Q,double *H,double *E1,double *E2,int i_end,int l_end,int k_end,double c1,double c2,int dir)
	{
	      int3 d1,d2;

	      getMagneticFieldTraceShifts(dir,d1,d2);

   		dim3 dimGrid(i_end+1,l_end+1,k_end+1),dimBlock(1,1,1);

//	    GPU_emh1//<<<dimGrid,dimBlock>>>(d_CellArray,Q,H,E1,E2,c1,c2,
//	    		d1,d2);
	    void* args[] = { (void* )&d_CellArray,
	            		           (void *)&Q,
	            		           (void *)&H,
	            		           (void *)&E1,
	            		           (void *)&E2,
	            		           (void *)&c1,
	                               (void *)&c2,
	            		           (void *)&d1,
	                               (void *)&d2,
	            		           0};
	              cudaError_t cudaStatus = cudaLaunchKernel(
	              	            (const void*)GPU_emh1, // pointer to kernel func.
	              	            dimGrid,                       // grid
	              	            dimBlock,                      // block
	              	            args,                          // arguments
	              	            0,
	              	            0
	              	            );



	      return 0;
	  }

	int SimpleMagneticFieldTrace(Cell &c,double *Q,double *H,int i_end,int l_end,int k_end)
	{


		   		dim3 dimGrid(i_end+1,l_end+1,k_end+1),dimBlock(1,1,1);

//			    GPU_emh2//<<<dimGrid,dimBlock>>>(d_CellArray,0,0,0,Q,H);
			    int i_s = 0;
			    int l_s = 0;
			    int k_s = 0;

			    void* args[] = { (void* )&d_CellArray,
			    		         (void *)&i_s,
			    		         (void *)&l_s,
			    		         (void *)&k_s,
			    		         (void *)&Q,
			    		         (void *)&H,
			   	            		     0};
			   	              cudaError_t cudaStatus = cudaLaunchKernel(
			   	              	            (const void*)GPU_emh2, // pointer to kernel func.
			   	              	            dimGrid,                       // grid
			   	              	            dimBlock,                      // block
			   	              	            args,                          // arguments
			   	              	            0,
			   	              	            0
			   	              	            );



	      return 0;
	  }


	  int PeriodicBoundaries(double *E,int dir,int start1,int end1,int start2,int end2,int N)
	  {
	      Cell  c = (*AllCells)[0];

//	      if(CPU_field == 0)
//	      {
	    		dim3 dimGrid(end1-start1+1,1,end2-start2+1),dimBlock(1,1,1);

//	    	    GPU_periodic//<<<dimGrid,dimBlock>>>(d_CellArray,start1,start2,E,dir,0,N);
	    	    int zero = 0;
	    	    void* args[] = { (void* )&d_CellArray,
	    	   			         (void *)&start1,
	    	   			    	 (void *)&start2,
	    	   			         (void *)&E,
	    	   			    	 (void *)&dir,
	    	   			    	 (void *)&zero,
	    	   			    	 (void *)&N,
	    	   			   	     0};
	    	     cudaError_t cudaStatus = cudaLaunchKernel(
	    	   	   	            (const void*)GPU_periodic, // pointer to kernel func.
	    	   			   	    dimGrid,                       // grid
	    	   			   	    dimBlock,                      // block
	    	   			   	    args,                          // arguments
	    	   			   	    0,
	    	   			   	    0
	    	   			   	    );

//	    	    GPU_periodic//<<<dimGrid,dimBlock>>>(d_CellArray,start1,start2,E,dir,N+1,1);
	    	    int one = 1;
	    	    int N1 = N+1;

				void *args1[] = { (void* )&d_CellArray,
								 (void *)&start1,
								 (void *)&start2,
								 (void *)&E,
								 (void *)&dir,
								 (void *)&N1,
								 (void *)&one,
								 0};
				 cudaStatus = cudaLaunchKernel(
								(const void*)GPU_periodic, // pointer to kernel func.
								dimGrid,                       // grid
								dimBlock,                      // block
								args1,                          // arguments
								0,
								0
								);
//
//	      }
//	      else
//	      {
//
//	      for(int k = start2;k <= end2;k++)
//	      {
//		  for(int i = start1;i <= end1;i++)
//		  {
//			  periodicElement(&c,i,k,E,dir,0,N);
//		  }
//	      }
//	      for(int k = start2;k <= end2;k++)
//	      {
//	         for(int i = start1;i <= end1;i++)
//	      	 {
//	        	 periodicElement(&c,i,k,E,dir,N+1,1);
//		  }
//	      }
//	      }
	      return 0;
	}

int SinglePeriodicBoundary(double *E,int dir,int start1,int end1,int start2,int end2,int N)
{
    Cell  c = (*AllCells)[0];

    if(CPU_field == 0)
    {
    	dim3 dimGrid(end1-start1+1,1,end2-start2+1),dimBlock(1,1,1);

//   	    GPU_periodic//<<<dimGrid,dimBlock>>>(d_CellArray,start1,start2,E,dir,N+1,1);
   	    int N1 = N+1;
   	    int one = 1;
     	void* args[] = { (void* )&d_CellArray,
   	         		           (void *)&start1,
   	         		           (void *)&start2,
   	         		           (void *)&E,
   	         		           (void *)&dir,
   	         		           (void *)&N1,
   	                            (void *)&one,
   	         		           0};
   	           cudaError_t cudaStatus = cudaLaunchKernel(
   	           	            (const void*)GPU_periodic, // pointer to kernel func.
   	           	            dimGrid,                       // grid
   	           	            dimBlock,                      // block
   	           	            args,                          // arguments
   	           	            16000,
   	           	            0
   	           	            );


    }
    else
    {
       for(int k = start2;k <= end2;k++)
       {
	  	  for(int i = start1;i <= end1;i++)
	  	  {
	  		  int3 i0,i1;

	                    int n   = c.getGlobalBoundaryCellNumber(i,k,dir,N+1);
	  		            int n1  = c.getGlobalBoundaryCellNumber(i,k,dir,1);
	  		            E[n]    = E[n1];
	  		            i0= c.getCellTripletNumber(n);
	  		            i1= c.getCellTripletNumber(n1);
	  		            std::cout << "ex1 "<< i0.x+1 << " "<< i0.y+1 << " " << i0.z+1  <<" " <<  i1.x+1 << " " << i1.y+1 << " " << i1.z+1  << " " << E[n]  << " " << E[n1] << std::endl;
	  		   	  }
	        }
    }
    return 0;
}





	  int SetPeriodicCurrentComponent(GPUCell **cells,double *J,int dir,int Nx,int Ny,int Nz)
	  {
		  dim3 dimGridX(Ny+2,1,Nz+2),dimGridY(Nx+2,1,Nz+2),dimGridZ(Nx+2,1,Ny+2),dimBlock(1,1,1);


//          GPU_CurrentPeriodic//<<<dimGridX,dimBlock>>>(cells,J,dir,0,0,0,Nx+2);
          int dir2 = 0;
          int i_s  = 0;
          int k_s  = 0;
          int N    = Nx+2;
          void* args[] = { (void* )&cells,
        		           (void *)&J,
        		           (void *)&dir,
        		           (void *)&dir2,
        		           (void *)&i_s,
        		           (void *)&k_s,
                           (void *)&N,
        		           0};
          cudaError_t cudaStatus = cudaLaunchKernel(
          	            (const void*)GPU_CurrentPeriodic, // pointer to kernel func.
          	            dimGridX,                       // grid
          	            dimBlock,                      // block
          	            args,                          // arguments
          	            16000,
          	            0
          	            );

//	      GPU_CurrentPeriodic//<<<dimGridY,dimBlock>>>(cells,J,dir,1,0,0,Ny+2);
	               dir2 = 1;
	               N    = Ny+2;
          cudaStatus = cudaLaunchKernel(
	               	            (const void*)GPU_CurrentPeriodic, // pointer to kernel func.
	               	            dimGridY,                       // grid
	               	            dimBlock,                      // block
	               	            args,                          // arguments
	               	            16000,
	               	            0
	               	            );

//	      GPU_CurrentPeriodic//<<<dimGridZ,dimBlock>>>(cells,J,dir,2,0,0,Nz+2);
	      dir2 = 2;
	      N    = Nz+2;
	      cudaStatus = cudaLaunchKernel(
	     	               	            (const void*)GPU_CurrentPeriodic, // pointer to kernel func.
	     	               	            dimGridZ,                       // grid
	     	               	            dimBlock,                      // block
	     	               	            args,                          // arguments
	     	               	            16000,
	     	               	            0
	     	               	            );

		  return 0;
	  }

	  void SetPeriodicCurrents(int nt)
	  {

		  memory_monitor("before275",nt);

		  checkControlPoint(275,nt,0);
		  SetPeriodicCurrentComponent(d_CellArray,d_Jx,0,Nx,Ny,Nz);
		  SetPeriodicCurrentComponent(d_CellArray,d_Jy,1,Nx,Ny,Nz);
		  SetPeriodicCurrentComponent(d_CellArray,d_Jz,2,Nx,Ny,Nz);


	     checkControlPoint(400,nt,0);

	   }

	  void InitQdebug(std::string fnjx,std::string fnjy,std::string fnjz)
	  {


	     read3Darray(fnjx, dbg_Qx);
	     read3Darray(fnjy, dbg_Qy);
	     read3Darray(fnjz, dbg_Qz);
	  }



void AssignCellsToArraysGPU()
{
	dim3 dimGrid(Nx,Ny,Nz),dimBlockExt(CellExtent,CellExtent,CellExtent);
	cudaError_t err = cudaGetLastError();
    printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); 

	size_t sz;
	err = cudaDeviceGetLimit(&sz,cudaLimitStackSize);
	printf("%s:%d - stack limit %d err = %d\n",__FILE__,__LINE__,((int)sz),err);
	err = cudaDeviceSetLimit(cudaLimitStackSize, 64*1024);
	printf("%s:%d - set stack limit %d \n",__FILE__,__LINE__,err);
	err= cudaDeviceGetLimit(&sz,cudaLimitStackSize);
	printf("%s:%d - stack limit %d err %d\n",__FILE__,__LINE__,((int)sz),err);

//	GPU_SetFieldsCel..ls//<<<dimGrid, dimBlockExt>>>(d_CellArray,d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz);
	void* args[] = { (void* )&d_CellArray,&d_Ex,&d_Ey,&d_Ez,&d_Hx,&d_Hy,&d_Hz,0};
	cudaError_t cudaStatus = cudaLaunchKernel(
	                        (const void*)GPU_SetFieldsToCells, // pointer to kernel func.
	                        dimGrid,                       // grid
	                        dimBlockExt,                      // block
	                        args,                          // arguments
	                        16000,
	                        0
	                        );

    
	cudaDeviceSynchronize();
	err = cudaGetLastError();
    printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); 

}


	  void AssignCellsToArrays()
	{
	     for(int n = 0;n < (*AllCells).size();n++)
	     {
	         Cell  c = (*AllCells)[n];
		 c.writeAllToArrays(Jx,Jy,Jz,Rho,0);

	     }
	     CheckArray(Jx, dbgJx);
	     SetPeriodicCurrents(0);
	     CheckArray(Jx, dbgJx);
	}



	  void write3Darray(char *name,double *d)
	  {
	    char fname[100];
	    GPUCell c = (*AllCells)[0];
	    FILE *f;

	    sprintf(fname,"%s_fiel3d.dat",name);

	    if((f = fopen(fname,"wt")) == NULL) return;

	    for(int i = 1;i < Nx+1;i++)
	    {
	        for(int l = 1;l < Ny+1;l++)
	        {
	            for(int k = 1;k < Nz+1;k++)
		    {
		        int n = c.getGlobalCellNumber(i,l,k);

			fprintf(f,"%15.5e %15.5e %15.5e %25.15e \n",c.getNodeX(i),c.getNodeY(l),c.getNodeZ(k),d[n]);
		    }
		}
	    }

	    fclose(f);
	}

void write3D_GPUArray(char *name,double *d_d)
{
//	double *d;

	return;

}

void readControlPoint(FILE **f1,char *fncpy,int num,int nt,int part_read,int field_assign,
		double *ex,double *ey,double *ez,
		double *hx,double *hy,double *hz,
		double *jx,double *jy,double *jz,
		double *qx,double *qy,double *qz,
		double *x,double *y,double *z,
		double *px,double *py,double *pz
		)
{
	char fn[100],fn_next[100];
	FILE *f;

	sprintf(fn,"mumu%03d%08d.dat",num,nt);
	strcpy(fncpy,fn);
	sprintf(fn_next,"mumu%03d%05d.dat",num,nt+1);
	if((f = fopen(fn,"rb")) == NULL) return;
	if(part_read)
	{
	   *f1 = f;
	}

	readFortranBinaryArray(f,ex);
	readFortranBinaryArray(f,ey);
	readFortranBinaryArray(f,ez);
	readFortranBinaryArray(f,hx);
	readFortranBinaryArray(f,hy);
	readFortranBinaryArray(f,hz);
	readFortranBinaryArray(f,jx);
	readFortranBinaryArray(f,jy);
	readFortranBinaryArray(f,jz);

	readFortranBinaryArray(f,qx);
	readFortranBinaryArray(f,qy);
	readFortranBinaryArray(f,qz);

	if(field_assign == 1) AssignArraysToCells();

}

double checkControlMatrix(char *wh,int nt,char *name, double *d_m)
{
	double t_jx;//,t_jy;//,t_jz;
	char fn[100];
	FILE *f;

#ifndef CHECK_CONTROL_MATRIX
	return 0.0;
#endif

	sprintf(fn,"wcmx_%4s_%08d_%2s.dat",wh,nt,name);
	if((f = fopen(fn,"rb")) == NULL) return -1.0;

	readFortranBinaryArray(f,dbgJx);

	t_jx = checkGPUArray(dbgJx,d_m);

    return t_jx;
}


void checkCurrentControlPoint(int j,int nt)
{
	 double /*t_ex,t_ey,t_ez,t_hx,t_hy,t_hz,*/t_jx,t_jy,t_jz;
		char fn[100];//,fn_next[100];
		FILE *f;

		sprintf(fn,"curr%03d%05d.dat",nt,j);
		if((f = fopen(fn,"rb")) == NULL) return;

		readFortranBinaryArray(f,dbgJx);
		readFortranBinaryArray(f,dbgJy);
		readFortranBinaryArray(f,dbgJz);

		int size = (Nx+2)*(Ny+2)*(Nz+2);

	 t_jx = CheckArraySilent(Jx,dbgJx,size);
	 t_jy = CheckArraySilent(Jy,dbgJy,size);
	 t_jz = CheckArraySilent(Jz,dbgJz,size);

     printf("Jx %15.5e Jy %15.5e Jz %15.5e \n",t_jx,t_jy,t_jz);
}

void checkControlPoint(int num,int nt,int check_part)
{
	 double t_ex,t_ey,t_ez,t_hx,t_hy,t_hz,t_jx,t_jy,t_jz;
	 double t_qx,t_qy,t_qz;//,t_njx,t_njy,t_njz;

	 if((nt != TOTAL_STEPS) || (num != 600))
	 {
#ifndef CONTROL_POINT_CHECK
	     return;
#endif
	 }

	 FILE *f;
	 char fn_copy[100];
//	 struct sysinfo info;

	 memory_monitor("checkControlPoint1",nt);

	 if(nt % FORTRAN_NUMBER_OF_SMALL_STEPS != 0) return;

	 memory_monitor("checkControlPoint2",nt);

	 readControlPoint(&f,fn_copy,num,nt,1,0,dbgEx,dbgEy,dbgEz,dbgHx,dbgHy,dbgHz,dbgJx,dbgJy,dbgJz,dbg_Qx,dbg_Qy,dbg_Qz,
                     dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz);

	 memory_monitor("checkControlPoint3",nt);

	 int size = (Nx+2)*(Ny+2)*(Nz+2);


	 t_ex = CheckArraySilent(Ex,dbgEx,size);
	 t_ey = CheckArraySilent(Ey,dbgEy,size);
	 t_ez = CheckArraySilent(Ez,dbgEz,size);
	 t_hx = CheckArraySilent(Hx,dbgHx,size);
	 t_hy = CheckArraySilent(Hy,dbgHy,size);
	 t_hz = CheckArraySilent(Hz,dbgHz,size);
	 t_jx = CheckArraySilent(Jx,dbgJx,size);
	 t_jy = CheckArraySilent(Jy,dbgJy,size);
	 t_jz = CheckArraySilent(Jz,dbgJz,size);

	 memory_monitor("checkControlPoint4",nt);

	 t_ex = CheckGPUArraySilent(dbgEx,d_Ex);
	 t_ey = CheckGPUArraySilent(dbgEy,d_Ey);
	 t_ez = CheckGPUArraySilent(dbgEz,d_Ez);
	 t_hx = CheckGPUArraySilent(dbgHx,d_Hx);
	 t_hy = CheckGPUArraySilent(dbgHy,d_Hy);
	 t_hz = CheckGPUArraySilent(dbgHz,d_Hz);

	 t_qx = CheckGPUArraySilent(dbg_Qx,d_Qx);
	 t_qy = CheckGPUArraySilent(dbg_Qy,d_Qy);
	 t_qz = CheckGPUArraySilent(dbg_Qz,d_Qz);

	 if(num >= 500)
	 {
		 char wh[100];

		 sprintf(wh,"%d",num);

		 t_jx = checkGPUArray(dbgJx,d_Jx,"Jx",wh,nt); //checkGPUArrayСomponent(dbgEx,d_Ex,"Ex",num);
		 t_jy = checkGPUArray(dbgJy,d_Jy,"Jy",wh,nt);
		 t_jz = checkGPUArray(dbgJz,d_Jz,"Jz",wh,nt);

		 t_ex = checkGPUArray(dbgEx,d_Ex,"Ex",wh,nt); //checkGPUArrayСomponent(dbgEx,d_Ex,"Ex",num);
		 t_ey = checkGPUArray(dbgEy,d_Ey,"Ey",wh,nt);
		 t_ez = checkGPUArray(dbgEz,d_Ez,"Ez",wh,nt);

	 }
	 else
	 {
		 t_ex = CheckGPUArraySilent(dbgEx,d_Ex);
		 t_ey = CheckGPUArraySilent(dbgEy,d_Ey);
		 t_ez = CheckGPUArraySilent(dbgEz,d_Ez);
	 }

	 t_jx = CheckGPUArraySilent(dbgJx,d_Jx);
	 t_jy = CheckGPUArraySilent(dbgJy,d_Jy);
	 t_jz = CheckGPUArraySilent(dbgJz,d_Jz);

//	 t_njx = CheckGPUArraySilent(dbgJx,d_Jx);
//	 t_njy = CheckGPUArraySilent(dbgJy,d_Jy);
//	 t_njz = CheckGPUArraySilent(dbgJz,d_Jz);

	 memory_monitor("checkControlPoint5",nt);

	 double t_cmp_jx = checkGPUArray(dbgJx,d_Jx,"Jx","step",nt);
	 double t_cmp_jy = checkGPUArray(dbgJy,d_Jy,"Jy","step",nt);
	 double t_cmp_jz = checkGPUArray(dbgJz,d_Jz,"Jz","step",nt);

#ifdef CONTROL_DIFF_GPU_PRINTS
     printf("GPU: Ex %15.5e Ey %15.5e Ez %15.5e \n",t_ex,t_ey,t_ez);
     printf("GPU: Hx %15.5e Hy %15.5e Ez %15.5e \n",t_hx,t_hy,t_hz);
     printf("GPU: Jx %15.5e Jy %15.5e Jz %15.5e \n",t_jx,t_jy,t_jz);
     printf("GPU compare : Jx %15.5e Jy %15.5e Jz %15.5e \n",t_cmp_jx,t_cmp_jy,t_cmp_jz);
#endif

     memory_monitor("checkControlPoint6",nt);

     double cp = checkControlPointParticles(num,f,fn_copy,nt);

     f_prec_report = fopen("control_points.dat","at");
     fprintf(f_prec_report,"nt %5d num %3d Ex %15.5e Ey %15.5e Ez %15.5e Hx %15.5e Hy %15.5e Hz %15.5e Jx %15.5e Jy %15.5e Jz %15.5e Qx %15.5e Qy %15.5e Qz %15.5e particles %15.5e\n",
    		 nt,num,
    		 t_ex,t_ey,t_ez,
    		 t_hx,t_hy,t_hz,
    		 t_jx,t_jy,t_jz,
    		 t_qx,t_qy,t_qz,
    		 cp
    		 );
     fclose(f_prec_report);

     memory_monitor("checkControlPoint7",nt);

     fclose(f);
}




void copyCellCurrentsToDevice(CellDouble *d_jx,CellDouble *d_jy,CellDouble *d_jz,
		                      CellDouble *h_jx,CellDouble *h_jy,CellDouble *h_jz)
{
	int err;

 	err = MemoryCopy(d_jx,h_jx,sizeof(CellDouble),HOST_TO_DEVICE);
 	if(err != cudaSuccess)
 	        {
 	         	printf("1copyCellCurrentsToDevice err %d %s \n",err,getErrorString(err));
 	       	exit(0);
 	        }
 	err = MemoryCopy(d_jy,h_jy,sizeof(CellDouble),HOST_TO_DEVICE);
 	if(err != cudaSuccess)
 	        {
 	         	printf("2copyCellCurrentsToDevice err %d %s \n",err,getErrorString(err));
 	       	exit(0);
 	        }
 	err = MemoryCopy(d_jz,h_jz,sizeof(CellDouble),HOST_TO_DEVICE);
 	if(err != cudaSuccess)
 	        {
 	         	printf("3copyCellCurrentsToDevice err %d %s \n",err,getErrorString(err));
 	       	exit(0);
 	        }

}


double CheckArray	(double* a, double* dbg_a,FILE *f)
	{
	    Cell c = (*AllCells)[0];
//	    int wrong = 0;
	    double diff = 0.0;



#ifdef CHECK_ARRAY_DETAIL_PRINTS
	    fprintf(f,"begin array checking=============================\n");
#endif
	    for(int n = 0;n < (Nx + 2)*(Ny + 2)*(Nz + 2);n++)
	    {
//	        double t  = a[n];
//		    double dt = dbg_a[n];
            diff += pow(a[n] - dbg_a[n],2.0);

	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
		    {

		       int3 i = c.getCellTripletNumber(n);
#ifdef CHECK_ARRAY_DETAIL_PRINTS
		       fprintf(f,"n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
#endif
     		}
	    }
#ifdef CHECK_ARRAY_DETAIL_PRINTS
	    fprintf(f,"  end array checking============================= %.4f less than %15.5e diff %15.5e \n",
	    		(1.0-((double)wrong/((Nx + 2)*(Ny + 2)*(Nz + 2)))),TOLERANCE,
	    		pow(diff/((Nx + 2)*(Ny + 2)*(Nz + 2)),0.5)
	    	  );
#endif

	    return pow(diff,0.5);
	}

double CheckArray	(double* a, double* dbg_a)
	{
	    Cell c = (*AllCells)[0];
	    int wrong = 0;
	    double diff = 0.0;
#ifdef CHECK_ARRAY_DETAIL_PRINTS
	    puts("begin array checking2=============================");
#endif
	    for(int n = 0;n < (Nx + 2)*(Ny + 2)*(Nz + 2);n++)
	    {
//	        double t  = a[n];
//		    double dt = dbg_a[n];
            diff += pow(a[n] - dbg_a[n],2.0);

	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
		    {

		       int3 i = c.getCellTripletNumber(n);
#ifdef CHECK_ARRAY_DETAIL_PRINTS
		       printf("n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
#endif
     		}
	    }
#ifdef CHECK_ARRAY_DETAIL_PRINTS
	    printf("  end array checking============================= %.4f less than %15.5e diff %15.5e \n",
	    		(1.0-((double)wrong/((Nx + 2)*(Ny + 2)*(Nz + 2)))),TOLERANCE,
	    		pow(diff/((Nx + 2)*(Ny + 2)*(Nz + 2)),0.5)
	    	  );
#endif

	    return (1.0-((double)wrong/((Nx + 2)*(Ny + 2)*(Nz + 2))));
	}




double CheckGPUArraySilent	(double* a, double* d_a)
	{
	    static double *t;
	    static int f = 1;
	    cudaError_t err;


	    if(f == 1)
	    {
	    	 t = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
	    	 f = 0;
	    }
	    MemoryCopy(t,d_a,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),DEVICE_TO_HOST);
	    err = cudaGetLastError();
	    if(err != cudaSuccess)
	            {
	             	printf("CheckArraySilent err %d %s \n",err,cudaGetErrorString(err));
	            	exit(0);
	            }


	   return CheckArraySilent(a,t,(Nx+2)*(Ny+2)*(Nz+2));
	}




	int CheckValue(double *a, double *dbg_a, int n)
	{
	    Cell  c = (*AllCells)[0];
//	    double t  = a[n];
//	    double dt = dbg_a[n];

	    if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
	    {

	       int3 i = c.getCellTripletNumber(n);
#ifdef CHECK_VALUE_DEBUG_PRINTS
	       printf("value n %5d i %3d l %3d k %3d %15.5e dbg %1.5e diff %15.5e \n",n,i.x,i.y,i.z,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]));
#endif

	       return 0;

	    }

	    return 1;
	}


	void read3DarrayLog(char* name, double* d, int offset, int col)
	{
	    char str[1000];
	    Cell c = (*AllCells)[0];
	    FILE *f;

	    //sprintf(fname,"%s_fiel3d.dat",name);

	    if((f = fopen(name,"rt")) == NULL) return;

	    while(fgets(str,1000,f) != NULL)
	    {
	          //str += offset;

	          int i = atoi(str + offset)      - 1;
	          int l = atoi(str + offset + 5)  - 1;
	          int k = atoi(str + offset + 10) - 1;
//	          if(i == 1 && l == 1 && k == 1 && col == 0)
//	          {
//	        	  int qq = 0;
//	          }
		  double t = atof(str + offset + 15 + col*25);
		  int n = c.getGlobalCellNumber(i,l,k);
		  d[n] = t;
#ifdef READ_ARRAY_LOG_PRINTS
		  printf("%d %d %5d %5d %15.5e \n",i,l,k,n,t);
#endif
	    }

	    fclose(f);

	}

	void read3Darray(char* name, double* d)
	{
	  char str[100];
	  Cell  c = (*AllCells)[0];
	  FILE *f;

	  //sprintf(fname,"%s_fiel3d.dat",name);

	  if((f = fopen(name,"rt")) == NULL) return;

	  while(fgets(str,100,f) != NULL)
	  {
	        int i = atoi(str);
	        int l = atoi(str + 10);
	        int k = atoi(str + 20);
		  double t = atof(str + 30);
		  //int n = c.getGlobalCellNumber(i,l,k);
		  int i1,l1,k1,n = c.getFortranCellNumber(i,l,k);
		  c.getFortranCellTriplet(n,&i1,&l1,&k1);
		  d[n] = t;
	  }

	  fclose(f);

	}


	void read3Darray(string name, double* d)
	{
	  char str[100];

	  strcpy(str,name.c_str());

	  read3Darray(str,d);

	}

	int PeriodicCurrentBoundaries(double* E, int dirE,int dir, int start1, int end1, int start2, int end2)
	{
	      Cell  c = (*AllCells)[0];

	      int N = getBoundaryLimit(dir);

	      for(int k = start2;k <= end2;k++)
	      {
	    	  for(int i = start1;i <= end1;i++)
		  {
		      int n1    = c.getGlobalBoundaryCellNumber(i,k,dir,1);
		      int n_Nm1 = c.getGlobalBoundaryCellNumber(i,k,dir,N-1);
	#ifdef DEBUG_PLASMA
		      int3 n1_3 = c.getCellTripletNumber(n1);
		      int3 n_Nm1_3 = c.getCellTripletNumber(n_Nm1);
	#endif
		      if(dir != dirE)
		      {
		         E[n1] += E[n_Nm1];


		      }
		      if(dir != 1 || dirE != 1)
		      {
		         E[n_Nm1] =  E[n1];
		      }
		      int n_Nm2 = c.getGlobalBoundaryCellNumber(i,k,dir,N-2);
		      int n0    = c.getGlobalBoundaryCellNumber(i,k,dir,0);
	#ifdef DEBUG_PLASMA
		      int3 n0_3 = c.getCellTripletNumber(n0);
		      int3 n_Nm2_3 = c.getCellTripletNumber(n_Nm2);
	#endif

		      E[n0] += E[n_Nm2];

		      E[n_Nm2] = E[n0];


		      //   E[n0] = E[n_Nm2];
		      //   E[n_Nm1] = E[n1];


		     // }
		  }
	      }
	      return 0;
	}

	void ClearAllParticles(void )
	{
	    for(int n = 0;n < (*AllCells).size();n++)
	    {
	        Cell c = (*AllCells)[n];

		c.ClearParticles();

	    }
	}





	public:

	  Plasma(int nx,int ny,int nz,double lx,double ly,double lz,double ni1,int n_per_cell1,double q_m,double TAU)
	   {
	     Nx = nx;
	     Ny = ny;
	     Nz = nz;

	     Lx = lx;
	     Ly = ly;
	     Lz = lz;

	     ni = ni1;

	     n_per_cell = n_per_cell1;
	     ion_q_m    = q_m;
	     tau        = TAU;
	   }

	   int initControlPointFile()
	   {
		   f_prec_report = fopen("control_points.dat","wt");
		   fclose(f_prec_report);

		   return 0;
	   }




	   int copyCellsWithParticlesToGPU()
	   {
		   Cell c000 = (*AllCells)[0];
		   magf = 1;
		   cudaError_t err;

		   int size = (Nx+2)*(Ny+2)*(Nz+2);

		   cp = (GPUCell **)malloc(size*sizeof(GPUCell *));
		   if((err = cudaGetLastError() ) != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

		   for(int i = 0;i < size;i++)
		   {
		     	GPUCell c,*d_c;
		   	   	// 	z0 = h_CellArray[i];
		   	    d_c = c.allocateCopyCellFromDevice();
		   	 if((err = cudaGetLastError() ) != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

		   	    cp[i] = d_c;
		   }
		   if((err = cudaGetLastError() ) != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
		   return 0;
	   }




	   void Free();

	   virtual void SetInitialConditions(){}

	   virtual void ParticleSort(){}


	   void ListAllParticles(int nt,std::string where)
	   	{
//	   		int cell_sum = 0;
//	   		int part_number = 0;
//	   		double t_hx,t_hy,t_hz,*dbg;
//	   		FILE *f;
//	   		char str[200];
	   		//Cell<Particle> **cp;

#ifndef LIST_ALL_PARTICLES
	   		return;
#endif
	   	}

//
	double TryCheckCurrent(int nt,double *npJx)
	{
//		double *dbg;//,t_hx;//,t_hy,t_hz;



//	  	dbg = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));

	  	// read magnetic field from "nt+1" exlg file - to consider emh2


	    return 1.0;//t_hx;
	}

	double checkNonPeriodicCurrents(int nt)
	{
//		double *dbg,t_hx,t_hy,t_hz;

		printf("CHECKING Non-periodic currents !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");

	  //	dbg = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
	  	TryCheckCurrent(nt,npJx);



		return 1.0;//(t_hx+t_hy+t_hz)/3.0;
	}

	double checkPeriodicCurrents(int nt)
	{
//		double *dbg;//,t_hx;//,t_hy;//,t_hz;

		printf("CHECKING periodic currents !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");

	  	//dbg = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));

	  	TryCheckCurrent(nt,Jx);

		return 1.0;//(t_hx+t_hy+t_hz)/3.0;
	}


int SetCurrentArraysToZero()
{
    cudaError_t err = cudaGetLastError();
    err = cudaGetLastError();
	printf("%s: %d [%d,%d,%d]\n",__FILE__, __LINE__, Nx,Ny,Nz);
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	memset(Jx,0,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
    memset(Jy,0,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	memset(Jz,0,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	cudaMemset(d_Jx,0,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	cudaMemset(d_Jy,0,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	cudaMemset(d_Jz,0,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	return 0;
}

int SetCurrentsInCellsToZero(int nt)
{
	dim3 dimGrid(Nx+2,Ny+2,Nz+2),dimBlockExt(CellExtent,CellExtent,CellExtent);
	char name[100];
	sprintf(name,"before_set_to_zero_%03d.dat",nt);

	write3D_GPUArray(name,d_Jx);

//	GPU_SetAllCurrentsToZero//<<<dimGrid, dimBlockExt,16000>>>(d_CellArray);
	void* args[] = { (void* )&d_CellArray,0};
    cudaError_t cudaStatus = cudaLaunchKernel(
                 (const void*)GPU_SetAllCurrentsToZero, // pointer to kernel func.
                 dimGrid,                       // grid
                 dimBlockExt,                   // block
                 args,                          // arguments
                 16000,
                 0
             );


	return 0;
}

int StepAllCells_fore_diagnostic(int nt)
{
	char name[100];

	memory_monitor("CellOrder_StepAllCells3",nt);

	sprintf(name,"before_step_%03d.dat",nt);
	write3D_GPUArray(name,d_Jx);
	//			printCellCurrents(270,nt,"jx","step");
	ListAllParticles(nt,"bStepAllCells");

	return 0;
}

int StepAllCells(int nt,double mass,double q_mass)
{
	   dim3 dimGrid(Nx+2,Ny+2,Nz+2),dimBlock(512,1,1);
	   cudaDeviceSynchronize();
       puts("begin step");



//	   GPU_StepAllCells//<<<dimGrid, dimBlock,16000>>>(d_CellArray);//,d_Jx);
//	            		     		                 mass,q_mass);

	   void* args[] = { (void* )&d_CellArray,0};
//	   void *d_args;
//	   cudaError_t err = cudaMalloc(&d_args,sizeof(d_CellArray)+sizeof(d_Jx));
//	   cudaError_t err1 = cudaMemcpy(d_args,args,sizeof(d_CellArray)+sizeof(d_Jx),cudaMemcpyHostToDevice);
	   cudaError_t cudaStatus = cudaLaunchKernel(
	                                            (const void*)GPU_StepAllCells, // pointer to kernel func.
	                                            dimGrid,                       // grid
	                                            dimBlock,                      // block
	                                            args,                          // arguments
	                                            16000,
	                                            0
	                                           );
//                mass,q_mass };


       std::cout<<"GPU_StepAllCells returns "<<cudaStatus<<std::endl;
	   dim3 dimBlock1(1,1,1);
	   void *args1[] = { (void* )&d_CellArray,&nt,0};
	   cudaStatus = cudaFuncSetCacheConfig((const void*)GPU_CurrentsAllCells,cudaFuncCachePreferShared);
	   std::cout<<"cudaFuncSetCacheConfig returns "<<cudaStatus<<" "<<cudaGetErrorString(cudaStatus)<<std::endl;
//	   GPU_CurrentsAllCells//<<<dimGrid, dimBlock,16000>>>(d_CellArray);//,0,d_Jx,
	   cudaStatus = cudaLaunchKernel(
	                                            (const void*)GPU_CurrentsAllCells, // pointer to kernel func.
	                                            dimGrid,                       // grid
	                                            dimBlock,                      // block
	                                            args1,                          // arguments
	                                            4000,
	                                            0
	                                           );
	   std::cout<<"GPU_CurrentsAllCells returns "<<cudaStatus<<" "<<cudaGetErrorString(cudaStatus)<<std::endl;
	            		     		                 //mass,q_mass);
	   puts("end step");
	   cudaDeviceSynchronize();

	   puts("end step-12");

	   return 0;
}

int StepAllCells_post_diagnostic(int nt)
{
	  memory_monitor("CellOrder_StepAllCells4",nt);

      ListAllParticles(nt,"aStepAllCells");
      cudaError_t err2 = cudaGetLastError();
	  char err_s[200];
      strcpy(err_s,cudaGetErrorString(err2));

      return (int)err2;
}


int WriteCurrentsFromCellsToArrays(int nt)
{
	char name[100];
	dim3 dimGrid(Nx+2,Ny+2,Nz+2);

	sprintf(name,"before_write_currents_%03d.dat",nt);
	write3D_GPUArray(name,d_Jx);

    dim3 dimExt(CellExtent,CellExtent,CellExtent);
//    GPU_WriteAllCurrents//<<<dimGrid, dimExt,16000>>>(d_CellArray,0,d_Jx,d_Jy,d_Jz,d_Rho);
    int zero = 0;
    void* args[] = { (void *)&d_CellArray,
    		         (void *)&zero,
    		         (void *)&d_Jx,
    		         (void *)&d_Jy,
    		         (void *)&d_Jz,
    		         (void *)&d_Rho,
    		         0};
        cudaError_t cudaStatus = cudaLaunchKernel(
                     (const void*)GPU_WriteAllCurrents, // pointer to kernel func.
                     dimGrid,                       // grid
                     dimExt,                   // block
                     args,                          // arguments
                     16000,
                     0
                 );

    memory_monitor("CellOrder_StepAllCells5",nt);

    sprintf(name,"after_write_currents_%03d.dat",nt);
	write3D_GPUArray(name,d_Jx);

	memory_monitor("CellOrder_StepAllCells6",nt);

	return 0;
}

int MakeParticleList(int nt,int *stage,int *stage1,int **d_stage,int **d_stage1)
{
	dim3 dimGrid(Nx+2,Ny+2,Nz+2),dimGridOne(1,1,1),dimBlock(512,1,1),
	     dimBlockOne(1,1,1),dimBlockGrow(1,1,1),dimBlockExt(CellExtent,CellExtent,CellExtent);
	dim3 dimGridBulk(Nx,Ny,Nz);
	cudaError_t before_MakeDepartureLists,after_MakeDepartureLists;
//      before_ArrangeFlights;//,after_ArrangeFlights;

    before_MakeDepartureLists = cudaGetLastError();
    printf("before_MakeDepartureLists %d %s blockdim %d %d %d\n",before_MakeDepartureLists,
    cudaGetErrorString(before_MakeDepartureLists),dimGrid.x,dimGrid.y,dimGrid.z);

//    int stage[4000],stage1[4000];//,*d_stage,*d_stage1;
    cudaMalloc(d_stage,sizeof(int)*(Nx+2)*(Ny+2)*(Nz+2));

    cudaMalloc(d_stage1,sizeof(int)*(Nx+2)*(Ny+2)*(Nz+2));

//    GPU_MakeDepartureLists//<<<dimGrid, dimBlockOne>>>(d_CellArray,nt,*d_stage);

    void* args[] = {
    		         (void *)&d_CellArray,
    		         (void *)&nt,
    		         (void *)d_stage,
    		         0};
        cudaError_t cudaStatus = cudaLaunchKernel(
                     (const void*)GPU_MakeDepartureLists, // pointer to kernel func.
                     dimGrid,                       // grid
                     dimBlockOne,                   // block
                     args,                          // arguments
                     16000,
                     0
                 );

    after_MakeDepartureLists = cudaGetLastError();
    if(after_MakeDepartureLists != cudaSuccess)
    {
       printf("after_MakeDepartureLists %d %s\n",after_MakeDepartureLists,cudaGetErrorString(after_MakeDepartureLists));
    }

    cudaDeviceSynchronize();

    int err = cudaGetLastError();

    if(err != cudaSuccess)
        {
           printf("MakeParticleList sync error %d %s\n",err,getErrorString(err));
        }
    err = MemoryCopy(stage,*d_stage,sizeof(int)*(Nx+2)*(Ny+2)*(Nz+2),DEVICE_TO_HOST);

    if(err != cudaSuccess)
    {
       printf("MakeParticleList error %d %s\n",err,getErrorString(err));
       exit(0);
    }

    return (int)err;
}

int inter_stage_diagnostic(int *stage,int nt)
{
	   if(stage[0] == TOO_MANY_PARTICLES)
	    {
	       printf("too many particles flying to (%d,%d,%d) from (%d,%d,%d) \n",
		          stage[1],stage[2],stage[3],stage[4],stage[5],stage[6]);
	       exit(0);
	    }

	    ListAllParticles(nt,"aMakeDepartureLists");
#ifdef BALANCING_PRINTS
	    before_ArrangeFlights = cudaGetLastError();
	    printf("before_ArrangeFlights %d %s\n",before_ArrangeFlights,cudaGetErrorString(before_ArrangeFlights));
#endif

	    return 0;
}

int reallyPassParticlesToAnotherCells(int nt,int *stage1,int *d_stage1)
{
    int err;
    dim3 dimGridBulk(Nx,Ny,Nz),dimBlockOne(1,1,1);
	cudaMemset(d_stage1,0,sizeof(int)*(Nx+2)*(Ny+2)*(Nz+2));


//	    GPU_ArrangeFlights//<<<dimGridBulk, dimBlockOne>>>(d_CellArray,nt,d_stage1);
	    void* args[] = {
	    		          (void* )&d_CellArray,
	    		          (void *)&nt,
	    		          (void*)&d_stage1,
	    		          0};

	        cudaError_t cudaStatus = cudaLaunchKernel(
	                     (const void*)GPU_ArrangeFlights, // pointer to kernel func.
	                     dimGridBulk,                       // grid
	                     dimBlockOne,                   // block
	                     args,                          // arguments
	                     16000,
	                     0
	                 );


#ifdef BALANCING_PRINTS
	    CUDA_Errot_t after_ArrangeFlights = cudaGetLastError();


    printf("after_ArrangeFlights %d %s\n",after_ArrangeFlights,cudaGetErrorString(after_ArrangeFlights));
            cudaDeviceSynchronize();
#endif

	err = MemoryCopy(stage1,d_stage1,sizeof(int)*(Nx+2)*(Ny+2)*(Nz+2),DEVICE_TO_HOST);
	if(err != cudaSuccess)
	{
	   puts("copy error");
	   exit(0);
	}
	ListAllParticles(nt,"aArrangeFlights");


	memory_monitor("CellOrder_StepAllCells7",nt);
	return (int)err;

}

int reorder_particles(int nt)
{
    int stage[4000],stage1[4000],*d_stage,*d_stage1,err;

    MakeParticleList(nt,stage,stage1,&d_stage,&d_stage1);

    inter_stage_diagnostic(stage,nt);

    err = reallyPassParticlesToAnotherCells(nt,stage1,d_stage1);

    return (int)err;
}

int Push(int nt,double mass,double q_mass)
{
	StepAllCells_fore_diagnostic(nt);

	StepAllCells(nt,mass,q_mass);
	puts("after StepAllCell");

	return StepAllCells_post_diagnostic(nt);
}

int SetCurrentsToZero(int nt)
{
    cudaError_t err = cudaGetLastError();
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	SetCurrentArraysToZero();

    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

	return SetCurrentsInCellsToZero(nt);
}


	void CellOrder_StepAllCells(int nt,double mass,double q_mass,int first)
	{
    cudaError_t err = cudaGetLastError();
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
		SetCurrentsToZero(nt);

    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

		Push(nt,mass,q_mass);
		puts("Push");
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

        WriteCurrentsFromCellsToArrays(nt);
        puts("writeCut2arr");
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }

        reorder_particles(nt);
    err = cudaGetLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,cudaGetErrorString(err)); }
	}





double checkControlPointParticlesOneSort(int check_point_num,FILE *f,GPUCell **copy_cells,int nt,int sort)
{

    double t = 0.0;
    int size = 1;
#ifdef CPU_DEBUG_RUN
    double q_m,m;

    memory_monitor("checkControlPointParticlesOneSort",nt);

  //  double x,y,z,px,pz,q_m,*buf,tp,m;
    //double dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz;

    Cell c0 = (*AllCells)[0];
    //int pn_min/*,pn_ave,pn_max*/;

    total_particles = readBinaryParticleArraysOneSort(f,&dbg_x,&dbg_y,&dbg_z,
   		                                             &dbg_px,&dbg_py,&dbg_pz,&q_m,&m,nt,sort);
    memory_monitor("checkControlPointParticlesOneSort2",nt);

    size = (*AllCells).size();

   	for(int i = 0;i < size;i++)
   	{
   	 	GPUCell c = *(copy_cells[i]);

#ifdef checkControlPointParticles_PRINT
             printf("cell %d particles %20d \n",i,c.number_of_particles);
#endif

   	 	t += c.checkCellParticles(check_point_num,dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz,q_m,m);
//   	 	if(t < 1.0)
//   	 	{
//   	 	   t += c.checkCellParticles(check_point_num,dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz,q_m,m);
//   	 	}
   	}
   	memory_monitor("checkControlPointParticlesOneSort3",nt);

   	free(dbg_x);
   	free(dbg_y);
   	free(dbg_z);

   	free(dbg_px);
   	free(dbg_py);
   	free(dbg_pz);
   	memory_monitor("checkControlPointParticlesOneSort4",nt);
#endif
	return t/size;
}

double checkControlPointParticles(int check_point_num,FILE *f,char *fname,int nt)
{
	double te = 0.0,ti = 0.0,tb = 0.0;
	struct sysinfo info;
#ifdef CPU_DEBUG_RUN
 //   Cell<Particle> **cp;

	int size = (*AllCells).size();
//  разобраться, как делается девайсрвый массив указателей на ячейки и сделать его копию здесь
//	cp = (GPUCell<Particle> **)malloc(size*sizeof(GPUCell<Particle> *));

//	copyCells(h_CellArray);

	char where[100];
	sprintf(where, "checkpoint%03d",check_point_num);
	copyCells(where,nt);

	//checkParticleNumbers(cp,check_point_num);

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticles %u \n",info.freeram/1024/1024);
#endif
#endif


//	if(check_point_num == 100)
//		{
//				int qq = 0;
//			//	tb  = checkControlPointParticlesOneSort(check_point_num,f,cp);
//		}
	GPUCell c = *(cp[141]);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticlesOneSort cell 141 particles %20d \n",c.number_of_particles);
#endif

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticles0.9 %u \n",info.freeram/1024/1024);
#endif
#endif

	ti  = checkControlPointParticlesOneSort(check_point_num,f,cp,nt,ION);
//	printf("IONS\n");
#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticles1 %u \n",info.freeram/1024/1024);
#endif
#endif

	te  = checkControlPointParticlesOneSort(check_point_num,f,cp,nt,PLASMA_ELECTRON);
//	printf("ELECTRONS\n");

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
#ifdef checkControlPointParticles_PRINTS
    printf("checkControlPointParticles1.5 %u \n",info.freeram/1024/1024);
#endif
#endif

	tb  = checkControlPointParticlesOneSort(check_point_num,f,cp,nt,BEAM_ELECTRON);
//	printf("BEAM\n");

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticles2 %u \n",info.freeram/1024/1024);
#endif
#endif



#endif

//    freeCellCopies(cp);
    memory_monitor("after_free",nt);
	return (te+ti+tb)/3.0;
}

int readControlFile(int nt)
{


#ifndef ATTRIBUTES_CHECK
	return 0;
#else
	FILE *f;
	char fname[100];
	static int first = 1;
	int size;//,jmp1;

	sprintf(fname,"ctrl%05d",nt);

	if((f = fopen(fname,"rb")) == NULL)
		{
		  puts("no ini-file");
		  exit(0);
		}

	fread(&size,sizeof(int),1,f);
	fread(&ami,sizeof(double),1,f);
	fread(&amf,sizeof(double),1,f);
	fread(&amb,sizeof(double),1,f);
	fread(&size,sizeof(int),1,f);

	fread(&size,sizeof(int),1,f);

	if(first == 1)
	{
		first = 0;
        ctrlParticles = (double *)malloc(size);
#ifdef ATTRIBUTES_CHECK
        memset(ctrlParticles,0,size);
        cudaMalloc(&d_ctrlParticles,size);
        cudaMemset(d_ctrlParticles,0,size);
        size_ctrlParticles = size;
#endif
	}
	fread(ctrlParticles,1,size,f);


	jmp = size/sizeof(double)/PARTICLE_ATTRIBUTES/3;


	return 0;
#endif
}



int memory_monitor(std::string legend,int nt)
{
	static int first = 1;
	static FILE *f;

#ifndef FREE_RAM_MONITOR
	return 1;
#endif

	if(first == 1)
	{
		first = 0;
		f = fopen("memory_monitor.log","wt");
	}

	size_t m_free,m_total;
	struct sysinfo info;


	cudaError_t err = cudaMemGetInfo(&m_free,&m_total);

	sysinfo(&info);                                                                //  1   2              3                 4                5
	fprintf(f,"step %10d %50s GPU memory total %10d free %10d free CPU memory %10u \n",nt,legend.c_str(),((int)m_total)/1024/1024,((int)m_free)/1024/1024,((int)info.freeram)/1024/1024);

	return 0;
}

int memory_status_print(int nt)
{
	size_t m_free,m_total;
	struct sysinfo info;


	cudaMemGetInfo(&m_free,&m_total);
	sysinfo(&info);

	#ifdef MEMORY_PRINTS
    printf("before Step  %10d CPU memory free %10u GPU memory total %10d free %10d\n",
 		   nt,info.freeram/1024/1024,m_total/1024/1024,m_free/1024/1024);
#endif

	return 0;
}


int Compute()
{
	   printf("----------------------------------------------------------- \n");
	   size_t m_free,m_total;

	   cudaMemGetInfo(&m_free,&m_total);




	   for(int nt = START_STEP_NUMBER;nt <= TOTAL_STEPS;nt++)
	   {
		   memory_status_print(nt);

	       Step(nt);

	       memory_status_print(nt);
	       printf("step %d ===================\n",nt);
	   }
	   printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");


	   return 0;
}


};







#endif /* GPU_PLASMA_H_ */

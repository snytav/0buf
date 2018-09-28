#include <math.h>
#include <string>
#include <stdio.h>

#include "archAPI.h"

#include <stdlib.h>
#include<string.h>

#include <sys/resource.h>
#include <stdint.h>

#include <sys/sysinfo.h>
#include <sys/time.h>

#include "particle.h"

//#include<cuda.h>

//struct sysinfo {
//       long uptime;             /* Seconds since boot */
//       unsigned long loads[3];  /* 1, 5, and 15 minute load averages */
//       unsigned long totalram;  /* Total usable main memory size */
//       unsigned long freeram;   /* Available memory size */
//       unsigned long sharedram; /* Amount of shared memory */
//       unsigned long bufferram; /* Memory used by buffers */
//       unsigned long totalswap; /* Total swap space size */
//       unsigned long freeswap;  /* swap space still available */
//       unsigned short procs;    /* Number of current processes */
//       unsigned long totalhigh; /* Total high memory size */
//       unsigned long freehigh;  /* Available high memory size */
//       unsigned int mem_unit;   /* Memory unit size in bytes */
//       char _f[20-2*sizeof(long)-sizeof(int)]; /* Padding for libc5 */
//   };

using namespace std;

int setPrintfLimit()
{
	size_t sizeP;

	printf("Particle size %lu %lu CurrentTensor %d short %d\n",sizeof(Particle),sizeof(Particle)/sizeof(double),sizeof(CurrentTensor),sizeof(char));

	cudaDeviceGetLimit(&sizeP,cudaLimitPrintfFifoSize);

	printf("printf default limit %lu \n",sizeP/1024/1024);

	sizeP *= 10000;
	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, sizeP);

	cudaDeviceGetLimit(&sizeP,cudaLimitPrintfFifoSize);

	printf("printf limit set to %lu \n",sizeP/1024/1024);

	return 0;
}

double get_meminfo(void)
{
	FILE *f;
	char str[100];
	int  mem_free;
	double dmem;
   // return 0.0;

	system("free>&free_mem_out.dat");


	if((f = fopen("free_mem_out.dat","rt")) == NULL) return 0.0;

	fgets(str,100,f);
	fgets(str,100,f);

	mem_free = atoi(str + 30);

	dmem = (((double)mem_free)/1024)/1024;

	return dmem;

}

double get_meminfo1(void)
{
	double retval=0;
	char tmp[256]={0x0};
	/* note= add a path to meminfo like /usr/bin/meminfo
	   to match where meminfo lives on your system */
	FILE *shellcommand=popen("meminfo","r");
	while(fgets(tmp,sizeof(tmp),shellcommand)!=NULL)
	{
		if(memcmp(tmp,"Mem:",4)==0)
		{
			int	wordcount=0;
			std::string delimiter=" ";
			char *p=strtok(tmp,delimiter.c_str());
			while(*p)
			{
				wordcount++;
				if(wordcount==3) retval=atof(p);
			}
		}
	}
	pclose(shellcommand);
	return retval;
}


double CheckArraySilent	(double* a, double* dbg_a,int size)
	{
	   // Cell<Particle> c = (*AllCells)[0];
	    double diff = 0.0;

	    for(int n = 0;n < size;n++)
	    {
            diff += pow(a[n] - dbg_a[n],2.0);

//	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
//		    {
//
//		       int3 i = c.getCellTripletNumber(n);
//
//     		}
	    }

	    return pow(diff/(size),0.5);
	}





void get_load_data_file_names(
		string & t_jxfile,
		string & t_jyfile,
		string & t_jzfile,
		string & t_d_jxfile,
		string & t_d_jyfile,
		string & t_d_jzfile,
		string & t_np_jxfile,
		string & t_np_jyfile,
		string & t_np_jzfile,
		string & t_qxfile,
		string & t_qyfile,
		string & t_qzfile,int nt)
{


    char d_exfile[100],d_eyfile[100],d_ezfile[100],d_hxfile[100],d_hyfile[100],d_hzfile[100];
    char d_0exfile[100],d_0eyfile[100],d_0ezfile[100];
    char jxfile[100],jyfile[100],jzfile[100];
    char np_jxfile[100],np_jyfile[100],np_jzfile[100];
    char np_exfile[100],np_eyfile[100],np_ezfile[100];
    char d_jxfile[100],d_jyfile[100],d_jzfile[100];
    char qxfile[100],qyfile[100],qzfile[100];
    char pfile[100],nextpfile[100];
//    char part_name[100];

    sprintf(qxfile,"dnqx%06d.dat",nt);
    sprintf(qyfile,"dnqy%06d.dat",nt);
    sprintf(qzfile,"dnqz%06d.dat",nt);



    sprintf(d_exfile,"dnex%06d.dat",2*nt-1);
    sprintf(d_eyfile,"dney%06d.dat",2*nt-1);
    sprintf(d_ezfile,"dnez%06d.dat",2*nt-1);

    sprintf(d_0exfile,"dnex%06d.dat",2*nt-2);
    sprintf(d_0eyfile,"dney%06d.dat",2*nt-2);
    sprintf(d_0ezfile,"dnez%06d.dat",2*nt-2);

    sprintf(d_hxfile,"dnhx%06d.dat",2*nt-1);
    sprintf(d_hyfile,"dnhy%06d.dat",2*nt-1);
    puts(d_hyfile);
    sprintf(d_hzfile,"dnhz%06d.dat",2*nt-1);

    sprintf(jxfile,"dnjx%06d.dat",2*nt);
    sprintf(jyfile,"dnjy%06d.dat",2*nt);
    sprintf(jzfile,"dnjz%06d.dat",2*nt);

    sprintf(d_jxfile,"npjx%06d.dat",2*nt);
    sprintf(d_jyfile,"npjy%06d.dat",2*nt);
    sprintf(d_jzfile,"npjz%06d.dat",2*nt);

    sprintf(np_jxfile,"npjx%06d.dat",2*nt);
    sprintf(np_jyfile,"npjy%06d.dat",2*nt);
    sprintf(np_jzfile,"npjz%06d.dat",2*nt);

    sprintf(np_exfile,"exlg%03d.dat",2*nt);
    sprintf(np_eyfile,"eylg%03d.dat",2*nt);
    sprintf(np_ezfile,"ezlg%03d.dat",2*nt);

    sprintf(pfile,    "part%06d000.dat",nt);
    sprintf(nextpfile,"part%06d000.dat",nt+2);


    t_jxfile =    jxfile;
    t_jyfile =    jyfile;
    t_jzfile =    jzfile;
    t_d_jxfile =  d_jxfile;
    t_d_jyfile =  d_jyfile;
    t_d_jzfile =  d_jzfile;
    t_np_jxfile = np_jxfile;
    t_np_jyfile = np_jyfile;
    t_np_jzfile = np_jzfile;
    t_qxfile =    qxfile;
    t_qyfile =    qyfile;
    t_qzfile =    qzfile;
}

void cudaMalloc3D(double **X,double **Y,double**Z,int nx,int ny,int nz)
{
	cudaMalloc(X,sizeof(double)*(nx+2)*(ny+2)*(nz+2));
	cudaMalloc(Y,sizeof(double)*(nx+2)*(ny+2)*(nz+2));
	cudaMalloc(Z,sizeof(double)*(nx+2)*(ny+2)*(nz+2));

}



void copyFieldsToGPU(
						double *d_Ex,double *d_Ey,double *d_Ez,
						double *d_Hx,double *d_Hy,double *d_Hz,
						double *d_Jx,double *d_Jy,double *d_Jz,
						double *d_npJx,double *d_npJy,double *d_npJz,
						double *d_Qx,double *d_Qy,double *d_Qz,
						double *Ex,double *Ey,double *Ez,
		        		double *Hx,double *Hy,double *Hz,
		        		double *Jx,double *Jy,double *Jz,
		        		double *npJx,double *npJy,double *npJz,
		                double *Qx,double *Qy,double *Qz,
		                int Nx,int Ny,int Nz
		)
{
	int err;

    err = MemoryCopy(d_Ex,Ex,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
    {
    	printf("1copyFieldsToGPU err %d %s \n",err,getErrorString(err));
    	exit(0);
    }
    err = MemoryCopy(d_Ey,Ey,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
    {
     	printf("2copyFieldsToGPU err %d %s \n",err,getErrorString(err));
    	exit(0);
    }

    err = MemoryCopy(d_Ez,Ez,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("3copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Hx,Hx,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("4copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    err = MemoryCopy(d_Hy,Hy,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("5copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    err = MemoryCopy(d_Hz,Hz,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("6copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Jx,Jx,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("7copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    err = MemoryCopy(d_Jy,Jy,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("8copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Jz,Jz,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("9copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_npJx,npJx,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("10copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_npJy,npJy,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("11copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_npJz,npJz,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("12copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Qx,Qx,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("13copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Qy,Qy,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("14copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Qz,Qz,sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("15copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
}

void InitGPUFields(
		double **d_Ex,double **d_Ey,double **d_Ez,
		double **d_Hx,double **d_Hy,double **d_Hz,
		double **d_Jx,double **d_Jy,double **d_Jz,
		double **d_npJx,double **d_npJy,double **d_npJz,
        double **d_Qx,double **d_Qy,double **d_Qz,
        double *Ex,double *Ey,double *Ez,
		double *Hx,double *Hy,double *Hz,
		double *Jx,double *Jy,double *Jz,
		double *npJx,double *npJy,double *npJz,
		double *Qx,double *Qy,double *Qz,
		int Nx,int Ny,int Nz
        )
{
	cudaMalloc3D(d_Ex,d_Ey,d_Ez,Nx,Ny,Nz);
	cudaMalloc3D(d_Hx,d_Hy,d_Hz,Nx,Ny,Nz);
	cudaMalloc3D(d_Jx,d_Jy,d_Jz,Nx,Ny,Nz);
	cudaMalloc3D(d_npJx,d_npJy,d_npJz,Nx,Ny,Nz);
	cudaMalloc3D(d_Qx,d_Qy,d_Qz,Nx,Ny,Nz);



    copyFieldsToGPU(
    		                        *d_Ex,*d_Ey,*d_Ez,
    								*d_Hx,*d_Hy,*d_Hz,
    								*d_Jx,*d_Jy,*d_Jz,
    								*d_npJx,*d_npJy,*d_npJz,
    								*d_Qx,*d_Qy,*d_Qz,
    								Ex,Ey,Ez,
    				        		Hx,Hy,Hz,
    				        		Jx,Jy,Jz,
    				        		npJx,npJy,npJz,
    				                Qx,Qy,Qz,
    				                Nx,Ny,Nz
    		);
}



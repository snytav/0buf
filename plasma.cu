#include "gpu_plasma.h"
#include <stdlib.h>
#include "mpi_shortcut.h"



int prepare_code_execution()
{
	 GPUPlasma<GPUCell> *plasma;
	  size_t sizeP;

	      printf("oarticle size %d %d \n",sizeof(Particle),sizeof(Particle)/sizeof(double));
	      cudaDeviceGetLimit(&sizeP,cudaLimitPrintfFifoSize);

	      printf("printf default limit %d \n",sizeP/1024/1024);

	      sizeP *= 10000;
	      cudaDeviceSetLimit(cudaLimitPrintfFifoSize, sizeP);

	      cudaDeviceGetLimit(&sizeP,cudaLimitPrintfFifoSize);



	      printf("printf limit set to %d \n",sizeP/1024/1024);


	      int err = cudaSetDevice(0);

	      printf("err %d \n",err);

	   plasma = new GPUPlasma<GPUCell>(100,4,4,1.1424,0.05,0.05,1.0,2000,1.0,0.001);

	   plasma->Initialize();

}

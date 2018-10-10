#include "gpu_plasma.h"
#include <stdlib.h>
#include "mpi_shortcut.h"
//TODO: gpu cell in the global array at copy from there appears to be not initialized

int main(int argc,char*argv[])
{
   Plasma *plasma;

   InitMPI(argc,argv);

   printf("begin Particle size %ld \n", sizeof(Particle));

   plasma = new Plasma(100,4,4,1.1424,0.05,0.05,1.0,2000,1.0,0.001);

   plasma->Initialize();

   plasma->Compute();

   CloseMPI();

   delete plasma;
   

   return 0;
}

#ifndef CELL_H
#define CELL_H


//#define VIRTUAL_FUNCTIONS



#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
//#include <cstdlib>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <vector>
#include "types.h"
#include "particle.h"

#define CellExtent 5
#define PARTICLE_MASS_TOLERANCE 1e-15
#define WRONG_PARTICLE_TYPE -13333
#define PARTICLES_FLYING_ONE_DIRECTION 50
#define TOO_MANY_PARTICLES -513131313


/*double d_sign(double a,double b)
{
    if(b < 0.0) return -a;
    else return a;
}*/

#ifdef GPU_PARTICLE
#ifdef __CUDACC__
   __global__
#endif
#endif
void MoveParticle(void *,int,double);

typedef double* double_pointer;

typedef struct CellDouble {
    double M[CellExtent][CellExtent][CellExtent];
} CellDouble;

typedef struct CellTotalField {

	CellDouble *Ex, *Ey, *Ez, *Hx,*Hy, *Hz;

} CellTotalField;

#define MAX_PPC   5000
const int MAX_particles_per_cell = MAX_PPC;



surface<void, 2> particle_surface;

//__shared__ CellDouble fd[9];


/*
#ifdef GPU_PARTICLE
__device__ void writeParticleToSurface(int n,Particle *p)
{
   	double3 x;
   	int size_p = sizeof(Particle);
   	x = p->GetX();
   	//surf2Dwrite(t,in_surface,ny*8,nx+shift);
   	surf2Dwrite(x.x,particle_surface,1*8,n);
   	surf2Dwrite(x.y,particle_surface,2*8,n);
   	surf2Dwrite(x.z,particle_surface,3*8,n);
}

__device__ void addParticleToSurface(Particle *p,int *number_of_particles)
{
	writeParticleToSurface(*number_of_particles,p);
	number_of_particles++;

}

__device__ void readParticleFromSurfaceDevice(int n,Particle *p)
{
	double3 x;
	int size_p = sizeof(Particle);
	//surf2Dwrite(t,in_surface,ny*8,nx+shift);
	surf2Dread(&(x.x),particle_surface,1*8,n);
	surf2Dread(&(x.y),particle_surface,2*8,n);
	surf2Dread(&(x.z),particle_surface,3*8,n);
	p->SetX(x);

}

__global__ void readParticleFromSurfaceGlobal(int n,Particle *p)
{
	readParticleFromSurfaceDevice(n,p);
}


 void readParticleFromSurfaceHost(int n,Particle *d_p)
{
	readParticleFromSurfaceGlobal<<<1,1>>>(n,d_p);
}

__device__ void removeParticleFromSurfaceDevice(int n,Particle *p,int *number_of_particles)
{
    int i,k;
    char b;

	readParticleFromSurfaceDevice(n,p);

	for(i = n;i < *number_of_particles-1;i++)
	{
		for(k = 0;k < sizeof(Particle);k++)
		{
			surf2Dread(&b,particle_surface,k,i+1);
			surf2Dwrite(&b,particle_surface,k,i);
		}
	}
	(*number_of_particles)--;
}

__global__ void removeParticleFromSurfaceGlobal(int n,Particle *d_p,int *np)
{
	removeParticleFromSurfaceDevice(n,d_p,np);
}

void removeParticleFromSurfaceHost(int n,Particle *d_p,int *np)
{
	removeParticleFromSurfaceGlobal<<<1,1>>>(n,d_p,np);
}
#endif
*/

__host__ __device__ int isNan(double t)
{
    if(t > 0)
    {
		//int i = 0;
    }
    else
    {
	   if(t <= 0)
	   {
		  //int i = 0;
	   }
	   else
	   {
		  return 1;
	   }
    }

    return 0;
}






class Cell
{
public:
  int i,l,k;
  double hx,hy,hz,tau;
  double x0,y0,z0;
  double xm,ym,zm;
  int Nx,Ny,Nz;
//  double y0;
  int jmp;
  double *d_ctrlParticles;

  int flag_wrong_current_cell;
  double *d_wrong_current_particle_attributes,*h_wrong_current_particle_attributes;


  CellDouble *Jx, *Ex, *Hx,*Jy, *Ey, *Hy,*Jz, *Ez, *Hz,*Rho;



#ifdef GPU_PARTICLE
//  cudaArray *cuParticleArray;
  double *doubParticleArray;
  int number_of_particles;
  int busyParticleArray;

  int arrival[3][3][3],departure[3][3][3];
  int departureListLength;
  Particle departureList[3][3][3][PARTICLES_FLYING_ONE_DIRECTION];
  double *doubArrivalArray;

  //get position in arrival matrix of a cell for the particles for cell "c"
#ifdef __CUDACC__
  __host__ __device__
#endif
  int getArrivalMatrixPosition(Cell *c,int *ix,int *iy,int *iz)
  {
	  int ci = c->i,cl = c->l,ck = c->k;

	  getGlobalCellNumberTriplet(&ci,&cl,&ck);

	  int dx = abs(i-ci),
		  dy = abs(l-cl),
		  dz = abs(k-ck);

	  if((dx == 1 || dy == 1 || dz == 1) && (dx < 2 && dy < 2 && dz < 2)) return 1;

	  *ix = 1 + ci - i;
	  *iy = 1 + cl - l;
	  *iz = 1 + ck - k;

	  return 0;
  }

#else
  thrust::host_vector<Particle> all_particles;
#endif

double *getParticles() { return  doubParticleArray; }

#ifdef __CUDACC__
__host__ __device__
#endif
virtual
int ClearCell()
{

	free(Jx);
	free(Ex);
	free(Hx);
	free(Jy);
	free(Ey);
	free(Hy);
	free(Jz);
	free(Ez);
	free(Hz);
	free(Rho);

    free(doubParticleArray);

    return 0;
}

#ifdef __CUDACC__
__host__ __device__
#endif
int AllocParticles()
{
      int size = sizeof(Particle)/sizeof(double);

      doubParticleArray = new double[size*MAX_particles_per_cell];
      number_of_particles = 0;
      busyParticleArray = 0;

//      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc(16, 16, 16, 16, cudaChannelFormatKindUnsigned);
//
//      //CUDA array is column-major. Thus here the FIRST DIMENSION is twice more (actuall it is the SECOND)
//      cudaMallocArray(&cuParticleArray, &channelDesc2, MAX_particles_per_cell, sizeof(Particle), cudaArraySurfaceLoadStore);

      return 0;
}

#ifdef __CUDACC__
__host__ __device__
#endif
double ParticleArrayRead(int n_particle,int attribute)
{
//	printf("n_p %d att %d\n",n_particle,attribute);
	return doubParticleArray[attribute + n_particle*sizeof(Particle)/sizeof(double)];
}

#ifdef __CUDACC__
__host__ __device__
#endif
void ParticleArrayWrite(int n_particle,int attribute,double t)
{
//	printf("n_p %d att %d\n",n_particle,attribute);
	doubParticleArray[attribute + n_particle*sizeof(Particle)/sizeof(double)] = t;
}

double
#ifdef __CUDACC__
 __host__ __device__
 #endif
int2double(int a,int b)
{
    unsigned long long u;
    double *d;

    u = a;
    u = u << 32 | b;
    d = (double*)&u;
    return *d;
}

void
 #ifdef __CUDACC__
 __host__ __device__
 #endif
 double2int(double d,int *a,int *b)
{
    unsigned long long *u,uu;


    u = (unsigned long long *)&d;
    uu = *u;

    *a = (int)(uu >>32);
    *b = (int)(uu & 0xffffffff);

}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
 void writeParticleToSurface(int n,Particle *p)
{

   	ParticleArrayWrite(n,0,p->m);
   	ParticleArrayWrite(n,1,p->x);
   	ParticleArrayWrite(n,2,p->y);
   	ParticleArrayWrite(n,3,p->z);
   	ParticleArrayWrite(n,4,p->pu);
   	ParticleArrayWrite(n,5,p->pv);
   	ParticleArrayWrite(n,6,p->pw);
 	ParticleArrayWrite(n,7,p->q_m);

 	ParticleArrayWrite(n,8,p->x1);
 	ParticleArrayWrite(n,9,p->y1);
 	ParticleArrayWrite(n,10,p->z1);

 	ParticleArrayWrite(n,11,((double)p->sort));
 	ParticleArrayWrite(n,12,((double)p->direction));

}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
 void addParticleToSurface(Particle *p,int *number_of_particles)
{


	writeParticleToSurface(*number_of_particles,p);

//	atomicAdd(number_of_particles, 1);
	(*number_of_particles)++;

	//busyParticleArray = 0;

}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
 Particle readParticleFromSurfaceDevice(int n)
{
    Particle p;

   	p.m = ParticleArrayRead(n,0);
   	p.x = ParticleArrayRead(n,1);
   	p.y = ParticleArrayRead(n,2);
   	p.z = ParticleArrayRead(n,3);
   	p.pu = ParticleArrayRead(n,4);
   	p.pv = ParticleArrayRead(n,5);
   	p.pw = ParticleArrayRead(n,6);
   	p.q_m = ParticleArrayRead(n,7);

  	p.x1 = ParticleArrayRead(n,8);
   	p.y1 = ParticleArrayRead(n,9);
   	p.z1 = ParticleArrayRead(n,10);
   	p.sort = (particle_sorts)ParticleArrayRead(n,11);
   	p.direction = (char)ParticleArrayRead(n,12);

   	return p;
}


public:


 #ifdef __CUDACC__
 __host__ __device__
 #endif
 void removeParticleFromSurfaceDevice(int n,Particle *p,int *number_of_particles)
{
    int i,k;
    double b;

//    do{
//    	busy = atomicCAS(&busyParticleArray,0,1);
//    }while(busy == 1);

	*p = readParticleFromSurfaceDevice(n);

	i = *number_of_particles-1;
	if(this->i == 1 && this->l == 0 && this->k == 0)
	{
#ifdef STRAY_DEBUG_PRINTS
		printf("deleteLAST FN %10d n %d i %d num %d \n",p->fortran_number,n,i,*number_of_particles);
#endif
	}

	if(n < i)
	{

		for(k = 0;k < sizeof(Particle)/sizeof(double);k++)
		{

			b = ParticleArrayRead(i,k);
            ParticleArrayWrite(n,k,b);
		}
	}

	(*number_of_particles)--;
	//atomicAdd(number_of_particles, -1);

	if(this->i == 1 && this->l == 0 && this->k == 0)
	{
#ifdef STRAY_DEBUG_PRINTS
		printf("deleteLAST FN %10d n %d i %d num %d \n",p->fortran_number,n,i,*number_of_particles);
#endif
	}
	//busyParticleArray = 0;

}




 #ifdef __CUDACC__
 __host__ __device__
 #endif

void printCellParticles(char *where,int nt)
{
	Particle p;

	for(int i = 0;i < number_of_particles;i++)
	{
		p = readParticleFromSurfaceDevice(i);
		printf("%s step %d i %5d sort %d FN %10d c ( %d %d %d ) pointInCell %d x %15.5e y %15.5e z %15.5e m %15.5e q_m %15.5e px %15.5e %15.5e %15.5e \n",where,nt,
				i,(int)p.sort,p.fortran_number,this->i,this->l,this->k,isPointInCell(p.GetX()),p.x,p.y,p.z, p.m,p.q_m, p.pu,p.pv,p.pw);
	}

}


#ifdef __CUDACC__
__host__ __device__
#endif
__host__ __device__
  double get_hx(){return hx;}



 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double get_hy(){return hy;}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double get_hz(){return hz;}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  int    get_i(){return i;}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  int    get_l(){return l;}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  int    get_k(){return k;}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  CellDouble& getJx(){return (*Jx);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  CellDouble& getJy(){return (*Jy);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  CellDouble& getJz(){return (*Jz);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  CellDouble& getRho(){return (*Rho);}

#ifdef GPU_PARTICLE
   double  *GetParticles(){return doubParticleArray;}
#else
   thrust::host_vector<Particle>&  GetParticles(){return all_particles;}
#endif

  void writeToArray(double *E);


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  double getCellFraction(double x,double x0,double hx){ double t = (x - x0)/hx;
                                                        return t;}

 #ifdef __CUDACC__
 __host__ __device__
 #endif


  int    getCellNumber(double x,double x0,double hx){   return ((int) (getCellFraction(x,x0,hx) + 1.0) +1);}
//  int    getCellNumberCenterGlobalized(double x,double x0,double hx){double t = (getCellFraction(x,x0,hx) + 1.5);  // C-style numbering for correct global array reference
//                                                           return (int)(t + 1) - 2;}
__host__ __device__

 #ifdef __CUDACC__
 __host__ __device__
 #endif


int    getCellNumberCenter(double x,double x0,double hx)
 {
	 double t = ((getCellFraction(x,x0,hx) + 1.0) + 1.5);         // Fortran-style numbering used for particle shift computation
                                                         return (int)(t);}

 #ifdef __CUDACC__
 __host__ __device__
 #endif

int    getCellNumberCenterY(double x,double x0,double hx){double t = (x-x0)/hx + 1.5;         // Fortran-style numbering used for particle shift computation
                                                         return (int)(t);}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
 int    getCellNumberCenterCurrent(double x,double x0,double hx){double t = (getCellFraction(x,x0,hx) + 1.5);         // Fortran-style numbering used for particle shift computation
                                                                 return (int)(t+1);}

 #ifdef __CUDACC__
 __host__ __device__
 #endif


  int getPointPosition(double x,double x0,double hx){return (int)getCellFraction(x,x0,hx);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif


  int getPointCell(double3 x){                                                                     // for Particle to Cell distribution:
      int i = getPointPosition(x.x,0.0,hx);
      int l = getPointPosition(x.y,0.0,hy);
      int k = getPointPosition(x.z,0.0,hz);
      return getGlobalCellNumber(i,l,k);
  }

 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getNodeX(int i){return (((double)i-0.5)*hx);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getNodeY(int i){return (((double)i-0.5)*hy);}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getNodeZ(int i){return (((double)i-0.5)*hz);}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCoreCell(int code,int i,int l, int k)
  {
         switch(code)
	 {
	   case 0: return (*Ex).M[i][l][k];
	   case 1: return (*Ey).M[i][l][k];
	   case 2: return (*Ez).M[i][l][k];
	   case 3: return (*Hx).M[i][l][k];
	   case 4: return (*Hy).M[i][l][k];
	   case 5: return (*Hz).M[i][l][k];
	   case 6: return (*Jx).M[i][l][k];
	   case 7: return (*Jy).M[i][l][k];
	   case 8: return (*Jz).M[i][l][k];
	   case 9: return (*Rho).M[i][l][k];
	   default: return 0.0;
	 }
  }

 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCellReminder(double x,double x0,double hx)    {return (getCellNumber(x,x0,hx)- getCellFraction(x,x0,hx)) - 1.0;}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCellCenterReminder(double x,double x0,double hx)
  {
    double t  = getCellNumberCenter(x,x0,hx);
    double tf = getCellFraction(x,x0,hx);
    //double q  = t - tf - 1.0;
    return (t - 0.5 - tf - 1.0);
  }


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCellTransitAverage(double hz,int i1, int i2,double x0){return (hz * (((double)i1 + (double)i2) * 0.5 - 2.0) + x0);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCellTransitRatio(double z1,double z,double z2){return (z2 - z) / (z1 - z);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCellTransitProduct(double z1,double z,double z2){return (z2 - z) * (z1 - z);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getRatioBasedX(double x1,double x,double s){      return (x + (x1 - x) * s);}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
  double getCenterRelatedShift(double x,double x1,int i,double hx,double x0){ //double t = 0.50*(x+x1);
                                                                            //  double t1 = hx*((double)i-2.50);
                                                                              return (0.50*(x+x1)-hx*((double)i-2.50) -x0);}







 #ifdef __CUDACC__
 __host__ __device__
 #endif
 void flyDirection(Particle *p,int *dx,int *dy,int *dz)
   {
	     *dx = (((p->x > x0 + hx) && (p->x < x0 + 2*hx)) || ((p->x < hx) && (i == Nx - 1))) ? 2 : ( ((p->x < x0) || ((p->x > xm-hx) && (i == 0))) ? 0 : 1);
 	     *dy = (((p->y > y0 + hy) && (p->y < y0 + 2*hy)) || ((p->y < hy) && (l == Ny - 1))) ? 2 : ( ((p->y < y0) || ((p->y > ym-hy) && (l == 0))) ? 0 : 1);
 	     *dz = (((p->z > z0 + hz) && (p->z < z0 + 2*hz)) || ((p->z < hz) && (k == Nz - 1))) ? 2 : ( ((p->z < z0) || ((p->z > zm-hz) && (k == 0))) ? 0 : 1);

 	     p->direction = (*dx)*9 +(*dy)*3 +(*dz);//(*dx) | ((*dy)<< 2) | ((*dz) << 4);
   }

 #ifdef __CUDACC__
 __host__ __device__
 #endif
   __host__ __device__ void inverseDirection(int *dx,int *dy,int *dz)
   {
	    *dx = (*dx == 2) ? 0 : ( *dx == 0 ? 2 : 1);
	    *dy = (*dy == 2) ? 0 : ( *dy == 0 ? 2 : 1);
	    *dz = (*dz == 2) ? 0 : ( *dz == 0 ? 2 : 1);
   }


 #ifdef __CUDACC__
 __host__ __device__
 #endif
bool isPointInCell(double3 x){return ((x0 <= x.x) && (x.x < x0 + hx) && (y0 <= x.y) && (x.y < y0 + hy) && (z0 <= x.z) && (x.z < z0 + hz));}

  //public:
 #ifdef __CUDACC__
 __host__ __device__
 #endif
    Cell(){}
 #ifdef __CUDACC__
 __host__ __device__
 #endif
   ~Cell(){}
 #ifdef __CUDACC__
 __host__ __device__
 #endif
 Cell(int i1,int l1,int k1,double Lx,double Ly, double Lz,int Nx1, int Ny1, int Nz1,double tau1)
                                                                                       {Cell();
											Nx = Nx1; Ny = Ny1; Nz = Nz1; i = i1; l = l1;k = k1;
                                                                                        hx = Lx/((double)Nx); hy = Ly/((double)Ny); hz = Lz/((double)Nz);
											xm = Lx; ym = Ly; zm = Lz;
                                                                                        x0 = (double)i*hx; y0= (double)l*hy; z0 = (double)k*hz; tau = tau1;
    }

__host__ __device__
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
 #ifdef __CUDACC__
 __host__ __device__
 #endif
double3 GetElectricField(
			    int3 i,int3  i1,
			    double& s1,double& s2,double& s3,double& s4,double& s5,double& s6,
			    double& s11,double& s21,double& s31,double& s41,double& s51,double& s61,
			    Particle *p,CellDouble *Ex1,CellDouble *Ey1,CellDouble *Ez1
)
{
        double3 E;
        int3 cell_num;

    cell_num.x = i.x;
    cell_num.y = i1.y;
    cell_num.z = i1.z;
	E.x =    Interpolate3D(Ex1,&cell_num,s1,s11,s4,s41,s6,s61,p,0);

	cell_num.x = i1.x;
	cell_num.y = i.y;
	cell_num.z = i1.z;
	E.y =    Interpolate3D(Ey1,&cell_num,s2,s21,s3,s31,s6,s61,p,1);

	cell_num.x = i1.x;//i1;
	cell_num.y = i1.y;//l1;
	cell_num.z = i.z; //k;
	E.z =    Interpolate3D(Ez1,&cell_num,s2,s21,s4,s41,s5,s51,p,2);

	return E;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
int3 getCellTripletNumber(int n){int3 i; i.z = n % (Nz+2);
							  i.y = ((n -i.z)/(Nz+2))%(Ny + 2);
							  i.x = (n - i.z - i.y*(Nz+2))/((Ny+2)*(Nz+2));
							  return i;}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
int getGlobalCellNumberTriplet(int *i,int *l,int *k){
	if(*i >= Nx+2) *i -= Nx+2;
	if(*l >= Ny+2) *l -= Ny+2;
	if(*k >= Nz+2) *k -= Nz+2;

	return 0;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
int getGlobalCellNumber        (int i,int l,int k){

	getGlobalCellNumberTriplet(&i,&l,&k);

	return (i*(Ny+2)*(Nz+2) + l*(Nz+2) + k);
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
int getWrapCellNumberTriplet(int *i,int *l,int *k){
	if(*i >= Nx)
	{
		*i -= Nx;
	}
	else
	{
		if(*i < 0)
		{
			*i += Nx;
		}
	}

	if(*l >= Ny)
	{
		*l -= Ny;
	}
	else
	{
		if(*l < 0)
		{
			*l += Ny;
		}
	}

	if(*k >= Nz)
	{
		*k -= Nz;
	}
	else
	{
		if(*k < 0)
		{
			*k += Nz;
		}
	}
	return 0;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
 int getWrapCellNumber        (int i,int l,int k){

	getWrapCellNumberTriplet(&i,&l,&k);

	return (i*(Ny+2)*(Nz+2) + l*(Nz+2) + k);
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
 int getFortranCellNumber       (int i,int l,int k)
{
//	if(i == 102 && l == 2 && k == 0)
//	{
//		int qq = 0;
//		qq =1;
//	}
	return getGlobalCellNumber(i-1,l-1,k-1);}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
void getFortranCellTriplet       (int n,int *i,int *l,int *k)
{
	*i = n/((Ny+2)*(Nz+2));
	*l = (n - (*i)*(Ny+2)*(Nz+2))/Nz;
	*k = n - (*i)*(Ny+2)*(Nz+2) - (*l)*(Nz+2);

	(*i) += 1;
	(*l) += 1;
	(*k) += 1;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif

    int getGlobalBoundaryCellNumber(int i,int k,int dir,int N){ int i1 = (dir==0)*N + (dir == 1)*i + (dir == 2)*i;
                                                                int l1 = (dir==0)*i + (dir == 1)*N + (dir == 2)*k;
							        int k1 = (dir==0)*k + (dir == 1)*k + (dir == 2)*N;
							return getGlobalCellNumber(i1,l1,k1);
                                                      }

 #ifdef __CUDACC__
 __host__ __device__
 #endif
int getPeriodicShift(int dir,int i,int Nx)
{
//	printf("getPeriodicShift dir %d i %d Nx %d result %d \n",dir,i,Nx,((i == -1)*(Nx-1) + (i == Nx)*0 + (i >= 0 && i < Nx)*i));
	return ((i == -1)*(Nx-1) + (i == Nx)*0 + (i >= 0 && i < Nx)*i);
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
int getDirectedPeriodicShift(int dir,int i)
{
	return (dir == 0)*getPeriodicShift(dir,i,Nx)+(dir == 1)*getPeriodicShift(dir,i,Ny) + (dir == 2)*getPeriodicShift(dir,i,Nz);
}




 #ifdef __CUDACC__
 __host__ __device__
 #endif
    int3 getFlyout(Particle p)
{
	int3 n3;

	n3.x = (x0 > p.x)*getDirectedPeriodicShift(0,i-1) + (x0 + hx < p.x)*getDirectedPeriodicShift(0,i+1) + ((x0 < p.x) && (x0 + hx > p.x))*i;
	n3.y = (y0 > p.y)*getDirectedPeriodicShift(1,l-1) + (y0 + hy < p.y)*getDirectedPeriodicShift(1,l+1) + ((y0 < p.y) && (y0 + hy > p.y))*l;
	n3.z = (z0 > p.z)*getDirectedPeriodicShift(2,k-1) + (z0 + hz < p.z)*getDirectedPeriodicShift(2,k+1) + ((z0 < p.z) && (z0 + hz > p.z))*k;

	return n3;
}

//__host__ __device__
//
//    int getPeriodicCellNumber(int i,int k,int boundary_loop_number,int field_dir,int boundary_dir,int N)
//{
//	int i1 = (dir==0)*i + ((dir == 1))*i + (dir == 2)*i;
//
//	int k1 = (dir==0)*k + (dir == 1)*k + (dir == 2)*N;
//
//	return getGlobalBoundaryCellNumber(i1,k1,boundary_dir,N);
//}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
void ClearParticles(){number_of_particles = 0;}

#ifdef __CUDACC__
__host__ __device__
#endif
void Init()
{
//#ifdef GPU_PARTICLE
//    cudaMalloc((void**)&Jx,sizeof(CellDouble));
//    cudaMalloc((void**)&Jy,sizeof(CellDouble));
//    cudaMalloc((void**)&Jz,sizeof(CellDouble));
//    cudaMalloc((void**)&Ex,sizeof(CellDouble));
//    cudaMalloc((void**)&Ey,sizeof(CellDouble));
//    cudaMalloc((void**)&Ez,sizeof(CellDouble));
//    cudaMalloc((void**)&Rho,sizeof(CellDouble));
//    cudaMalloc((void**)&Hx,sizeof(CellDouble));
//    cudaMalloc((void**)&Hy,sizeof(CellDouble));
//    cudaMalloc((void**)&Hz,sizeof(CellDouble));
//#else
    Jx = new CellDouble;
    Jy = new CellDouble;
    Jz = new CellDouble;
    Ex = new CellDouble;
    Ey = new CellDouble;
    Ez = new CellDouble;
    Hx = new CellDouble;
    Hy = new CellDouble;
    Hz = new CellDouble;
    Rho = new CellDouble;

    AllocParticles();
/*#endif*/
}

#ifdef __CUDACC__
 __host__ __device__
#endif
void SetZero()
{

     Init();

     for(int i1 = 0; i1 < CellExtent;i1++)
     {
         for(int l1 = 0; l1 < CellExtent;l1++)
	 {
	     for(int k1 = 0; k1 < CellExtent;k1++)
	     {
	         (*Jx).M[i1][k1][l1]  = 0.0;
	         (*Jy).M[i1][k1][l1]  = 0.0;
	         (*Jz).M[i1][k1][l1]  = 0.0;
	         (*Rho).M[i1][k1][l1] = 0.0;

	         (*Ex).M[i1][k1][l1]  = 0.0;
	         (*Ey).M[i1][k1][l1]  = 0.0;
	         (*Ez).M[i1][k1][l1]  = 0.0;

	         (*Hx).M[i1][k1][l1]  = 0.0;
	         (*Hy).M[i1][k1][l1]  = 0.0;
	         (*Hz).M[i1][k1][l1]  = 0.0;
	     }
	 }
     }
}

#ifdef __CUDACC__
__host__ __device__
#endif
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
void InverseKernel(double3 x,int3 & i,int3 & i1,
			    double& s1,double& s2,double& s3,double& s4,double& s5,double& s6,
			    double& s11,double& s21,double& s31,double& s41,double& s51,double& s61,
			    Particle *p
 			  )
{

//	int2 g;

//	double2 f;
       // s2 = getCellFraction(x.x,0.0,hx);
	i.x =  getCellNumber(x.x,x0,hx);            //(int) (s2 + 1.);  // FORTRAN-StYLE NUMBERING
	i1.x = getCellNumberCenter(x.x,x0,hx);      //(int) (s2 + 1.5);
	s1 = s1_interpolate(x.x);          //i - s2;
	s2 = s2_interpolate(x.x); //getCellCenterReminder(x,0.0,hx);    //i1 - 0.5 - s2;

//	s4   = getCellFraction(x.y,y0,hy);
    i.y  = getCellNumber(x.y,y0,hy);            //(int) (s2 + 1.);
	i1.y = getCellNumberCenter(x.y,y0,hy);      //(int) (s2 + 1.5);

	s3 = s3_interpolate(x.y);//getCellReminder(y,y0,hy);          //i - s2;
	s4 = s4_interpolate(x.y);//   getCellCenterReminder(y,y0,hy);    //i1 - 0.5 - s2;

//	s6   = getCellFraction(x.z,z0,hz);
	i.z  = getCellNumber(x.z,z0,hz);            //(int) (s2 + 1.);
	i1.z = getCellNumberCenter(x.z,z0,hz);      //(int) (s2 + 1.5);
	s5 = s5_interpolate(x.z); //getCellReminder(z,z0,hz);          //i - s2;
	s6 = s6_interpolate(x.z); //getCellCenterReminder(z,z0,hz);    //i1 - 0.5 - s2;

	s11 = 1. - s1;
	s21 = 1. - s2;
	s31 = 1. - s3;
	s41 = 1. - s4;
	s51 = 1. - s5;
	s61 = 1. - s6;


}


 #ifdef __CUDACC__
 __host__ __device__
 #endif
//double Interpolate3D(int n_field,int i,int l,int k,double sx,double sx1,double sy,double sy1,double sz,double sz1)
double Interpolate3D(CellDouble *E,int3 *cell,
		double sx,double sx1,double sy,
		double sy1,double sz,double sz1,
		Particle *p, int n
		)
{
       double t,t1,t2;
       double t_ilk,t_ilk1,t_il1k,t_il1k1,t_i1lk,t_i1lk1,t_i1l1k,t_i1l1k1;
       int i,l,k;

       i = cell->x;
       l = cell->y;
       k = cell->z;



      // printf("jmp %d number %d sort %d  l %d c.y %d \n",jmp,p->fortran_number,p->sort,l,cell->y);
     //  d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,55)] = l;


       // C STYLE NUMBERING, so -1:  i+this->i-1
//       i--;
//       l--;
//       k--;
       if(i < 0 || i > CellExtent) return 0.0;
       if(l < 0 || l > CellExtent) return 0.0;
       if(k < 0 || k > CellExtent) return 0.0;

#ifdef ATTRIBUTES_CHECK
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,110+20*n)] = (sz  * E->M[i][l][k] + sz1 * E->M[i][l][k + 1]);
//       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,122+20*n)] = E->(sz  * E->M[i][l][k] + sz1 * E->M[i][l][k + 1]);
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,111+20*n)] = (sz  * E->M[i][l + 1][k] + sz1 * E->M[i][l + 1][k + 1]);
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,112+20*n)] = (sz  * E->M[i + 1][l][k] + sz1 * E->M[i + 1][l][k + 1]);
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,113+20*n)] = (sz  * E->M[i + 1][l + 1][k]+ sz1 * E->M[i + 1][l + 1][k + 1]);

       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,114+20*n)] = sx;//  * E->M[i + 1][l + 1][k]+ sz1 * E->M[i + 1][l + 1][k + 1]);
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,115+20*n)] = sx1; //(sz  * E->M[i + 1][l + 1][k]+ sz1 * E->M[i + 1][l + 1][k + 1]);
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,116+20*n)] = sy; //+ sz1 * E->M[i + 1][l + 1][k + 1]);
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,117+20*n)] = sy1;

       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,118+20*n)] = i+this->i-1;

       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,119+20*n)] = l+this->l-1;
       //d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,128+20*n)] = k+this->k-1;

       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,120+20*n)] = E->M[i][l][k];
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,121+20*n)] = E->M[i][l + 1][k];
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,122+20*n)] = E->M[i+1][l][k];
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,123+20*n)] = E->M[i+1][l+1][k];

       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,124+20*n)] = E->M[i][l][k+1];
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,125+20*n)] = E->M[i][l+1][k+1];
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,126+20*n)] = E->M[i+1][l][k+1];
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,127+20*n)] = E->M[i+1][l+1][k+1];
#endif
       t_ilk   = E->M[i][l][k];
       t_ilk1  = E->M[i][l][k + 1];
       t_il1k  = E->M[i][l + 1][k];
       t_il1k1 = E->M[i][l + 1][k + 1];

       t_i1lk   = E->M[i+1][l][k];
       t_i1lk1  = E->M[i+1][l][k + 1];
       t_i1l1k  = E->M[i+1][l + 1][k];
       t_i1l1k1 = E->M[i+1][l + 1][k + 1];

       t1 = sx *  (sy *(sz  * t_ilk
               		   		          + sz1 * t_ilk1)
               		   		    +  sy1*(sz  * t_il1k
               		   		          + sz1 * t_il1k1));
#ifdef ATTRIBUTES_CHECK
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,128+20*n)] = t1;
#endif
       t2 = sx1 * (sy *(sz  * t_i1lk
    		   	                  + sz1 * t_i1lk1)
    		   	            +  sy1*(sz  * t_i1l1k
    		   	                  + sz1 * t_i1l1k1));
#ifdef ATTRIBUTES_CHECK
       d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,129+20*n)] = t2;
#endif
//       t =   sx*sy*sz  * E->M[i][l][k]
//           + sx*sy*sz1 * E->M[i][l][k + 1]
//           + sx*sy1*sz  * E->M[i][l + 1][k]
//           + sx*sy1*sz1 * E->M[i][l + 1][k + 1]
//           + sx1*sy*sz* E->M[i + 1][l][k]
//           + sx1*sy*sz1 * E->M[i + 1][l][k + 1]
//           + sx1*sy1*sz  * E->M[i + 1][l + 1][k]
//           + sx1*sy1*sz1 * E->M[i + 1][l + 1][k + 1];
      t  = t2 + t1;
//       t  =    (sx *  (sy *(sz  * fd[n_field].M[i][l][k]
//		          + sz1 * fd[n_field].M[i][l][k + 1])
//		    +  sy1*(sz  * fd[n_field].M[i][l + 1][k]
//		          + sz1 * fd[n_field].M[i][l + 1][k + 1]))
//	      + sx1 * (sy *(sz  * fd[n_field].M[i + 1][l][k]
//	                  + sz1 * fd[n_field].M[i + 1][l][k + 1])
//	            +  sy1*(sz  * fd[n_field].M[i + 1][l + 1][k]
//	                  + sz1 * fd[n_field].M[i + 1][l + 1][k + 1])));
    //   d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,45+10*n)] = t;

       return t;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
double3 GetMagneticField(
			    int3 i,int3 i1,
			    double& s1,double& s2,double& s3,double& s4,double& s5,double& s6,
			    double& s11,double& s21,double& s31,double& s41,double& s51,double& s61,
			    Particle *p,CellDouble *Hx1,CellDouble *Hy1,CellDouble *Hz1
)
{
        double3 H;
        int3 cn;

        cn.x = i1.x;
        cn.y = i.y;
        cn.z = i.z;
    	H.x =    Interpolate3D(Hx1,&cn,s2,s21,s3,s31,s5,s51,p,3);

    	cn.x = i.x;
    	cn.y = i1.y;
    	cn.z = i.z;
    	H.y =    Interpolate3D(Hy1,&cn,s1,s11,s4,s41,s5,s51,p,4);

    	cn.x = i.x;
    	cn.y = i.y;
    	cn.z = i1.z;
    	H.z =    Interpolate3D(Hz1,&cn,s1,s11,s3,s31,s6,s61,p,5);


	return H;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
 double s1_interpolate(double x)
{
	return (int)(x/hx + 1.0) - x/hx;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
double s2_interpolate(double x)
{
	return (int)(x/hx + 1.5) -0.5 - x/hx;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
double s3_interpolate(double y)
{
	return (int)(y/hy + 1.0) - y/hy;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
double s5_interpolate(double z)
{
	return (int)(z/hz + 1.0) - z/hz;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
double s4_interpolate(double y)
{
	return (int)(y/hy + 1.5) -0.5 - y/hy;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
double s6_interpolate(double z)
{
	return (int)(z/hz + 1.5) -0.5 - z/hz;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
void writeWrongCurrentAttribute(int num_along_cell,double a,int num_attr)
{
	d_wrong_current_particle_attributes[num_along_cell*PARTICLE_ATTRIBUTES + num_attr-1] = a;
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
double getWrongCurrentAttribute(int num_along_cell,int num_attr)
{
	return d_wrong_current_particle_attributes[num_along_cell*PARTICLE_ATTRIBUTES + num_attr-1];
}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
Field GetField(Particle *p,CellTotalField cf)
{
    int3 i,i1;
	double s1,s2,s3,s4,s5,s6,s11,s21,s31,s41,s51,s61;
	double3 x = p->GetX();
	Field fd;

	    if(x.x < 0 || x.y < 0 || x.z < 0)
	    {
	    	double3 x;
	    	x.x = 0;
	    	x.y = 0;
	    	x.z = 0;
	    	fd.E = x;
	    	fd.H = x;
	    	return fd;
	    }

        InverseKernel(x,
	  	              i,i1,
		              s1,s2,s3,s4,s5,s6,
	                  s11,s21,s31,s41,s51,s61,p);



        fd.E = GetElectricField(i,i1,
			     s1,s2,s3,s4,s5,s6,
			     s11,s21,s31,s41,s51,s61,p,cf.Ex,cf.Ey,cf.Ez);

	    fd.H = GetMagneticField(i,i1,
			     s1,s2,s3,s4,s5,s6,
			     s11,s21,s31,s41,s51,s61,p,cf.Hx,cf.Hy,cf.Hz);
	    return fd;
}

#ifdef __CUDACC__
 __host__ __device__
#endif
double d_sign(double a, double b)
{
    if(b > 0) return a;
    else return -a;
}



 #ifdef __CUDACC__
 __host__ __device__
 #endif
void CurrentToMesh(double tau,int *cells,DoubleCurrentTensor *dt,Particle *p)
{
      double3 x2;
      double s;
      int3 i2,i1;
      int  m,i,l,k;
      double3 x   = p->GetX();
      double3 x1  = p->GetX1();
      double mass = p->m;
      double q_m  = p->q_m;
//      DoubleCurrentTensor dt;
      CurrentTensor t;

      CurrentTensor *t1 = &(dt->t1);
      CurrentTensor *t2 = &(dt->t2);


      *cells = 1;

      i1.x=getCellNumberCenterCurrent(x.x,x0,hx);
      i1.y=getCellNumberCenterCurrent(x.y,y0,hy);
      i1.z=getCellNumberCenterCurrent(x.z,z0,hz);


      i2.x=getCellNumberCenterCurrent(x1.x,x0,hx);


      i2.y=getCellNumberCenterCurrent(x1.y,y0,hy);


      i2.z=getCellNumberCenterCurrent(x1.z,z0,hz);

      i=abs(i2.x-i1.x);
      l=abs(i2.y-i1.y);
      k=abs(i2.z-i1.z);
      m=4*i+2*l+k;


	switch (m) {
	    case 1:  goto L1;
	    case 2:  goto L2;
	    case 3:  goto L3;
	    case 4:  goto L4;
	    case 5:  goto L5;
	    case 6:  goto L6;
	    case 7:  goto L7;
	}


	pqr(i1,x,x1,mass,tau,t1,0,p);
	goto L18;
L1:
	x2.z = getCellTransitAverage(hz,i1.z,i2.z,z0);                   //hz * ((i1.z + i2.z) * .5 - 1.);
	s    = getCellTransitRatio(x1.z,x.z,x2.z);                    //   (z2 - z__) / (z1 - z__);
	x2.x = getRatioBasedX(x1.x,x.x,s);                            // x + (x1 - x) * s;
	x2.y = getRatioBasedX(x1.y,x.y,s);                            //y + (y1 - y) * s;

	goto L11;
L2:
	x2.y = getCellTransitAverage(hy,i1.y,i2.y,y0);                   //   d_1.h2 * ((l1 + l2) * .5 - 1.) + y0;
	s    = getCellTransitRatio(x1.y,x.y,x2.y);                    // (y2 - y) / (y1 - y);
	x2.x = getRatioBasedX(x1.x,x.x,s);                            //x + (x1 - x) * s;
	x2.z = getRatioBasedX(x1.z,x.z,s);
	goto L11;
L3:
	x2.y = getCellTransitAverage(hy,i1.y,i2.y,y0);                     //d_1.h2 * ((l1 + l2) * .5 - 1.) + y0;
	x2.z = getCellTransitAverage(hz,i1.z,i2.z,z0);                     //d_1.h3 * ((k1 + k2) * .5 - 1.);

	s = (getCellTransitProduct(x1.z,x.z,x2.z) + getCellTransitProduct(x1.y,x.y,x2.y))
			/ (pow(x1.z - x.z,2.0) + pow(x1.y - x.y,2.0));
	x2.x = getRatioBasedX(x1.x,x.x,s);


	goto L11;
L4:
	x2.x = getCellTransitAverage(hx,i1.x,i2.x,x0);
	s    = getCellTransitRatio(x1.x,x.x,x2.x);                       //s = (x2 - x) / (x1 - x);
	x2.y = getRatioBasedX(x1.y,x.y,s);                               //y2 = y + (y1 - y) * s;
	x2.z = getRatioBasedX(x1.z,x.z,s);                               //z2 = z__ + (z1 - z__) * s;
	goto L11;
L5:
	x2.x = getCellTransitAverage(hx,i1.x,i2.x,x0);                      //x2 = d_1.h1 * ((i1 + i2) * .5 - 1.);
	x2.z = getCellTransitAverage(hz,i1.z,i2.z,z0);                      //z2 = d_1.h3 * ((k1 + k2) * .5 - 1.);

	s = (getCellTransitProduct(x1.z,x.z,x2.z) + getCellTransitProduct(x1.x,x.x,x2.x)) / (pow(x1.z - x.z,2.0) + pow(x1.x - x.x,2.0));
	x2.y = getRatioBasedX(x1.y,x.y,s); //	y2 = y + (y1 - y) * s;
	goto L11;
L6:
	x2.x = getCellTransitAverage(hx,i1.x,i2.x,x0);  //x2 = d_1.h1 * ((i1 + i2) * .5 - 1.);
	x2.y = getCellTransitAverage(hy,i1.y,i2.y,y0);  //y2 = d_1.h2 * ((l1 + l2) * .5 - 1.) + y0;

	s = (getCellTransitProduct(x1.y,x.y,x2.y) +  getCellTransitProduct(x1.x,x.x,x2.x)) / (pow(x1.y - x.y,2.0) + pow(x1.x - x.x,2.0));
	x2.z = getRatioBasedX(x1.z,x.z,s); //	z2 = z__ + (z1 - z__) * s;
	goto L11;
L7:
	x2.x = getCellTransitAverage(hx,i1.x,i2.x,x0);  //x2 = d_1.h1 * ((i1 + i2) * .5 - 1.);
	x2.y = getCellTransitAverage(hy,i1.y,i2.y,y0);  //y2 = d_1.h2 * ((l1 + l2) * .5 - 1.) + y0;
	x2.z = getCellTransitAverage(hz,i1.z,i2.z,z0);  //z2 = d_1.h3 * ((k1 + k2) * .5 - 1.);
L11:
	pqr(i1, x, x2,  mass,tau,t1,0,p);
	pqr(i2, x2, x1, mass,tau,t2,1,p);

	*cells = 2;

L18:  p->x = p->x1;
      p->y = p->y1;
      p->z = p->z1;

      Reflect(p);

//      *dt1 = dt;

      t = *t1;
      *t1 = t;

      return;
}

#ifdef __CUDACC__
 __host__ __device__
#endif
void Reflect(Particle *p)
{
        double3 x1 = p->GetX();

        x1.x = (x1.x > xm)*(x1.x - xm) + (x1.x < 0.0)*(xm + x1.x) + (x1.x > 0 && x1.x < xm)*x1.x;
        x1.y = (x1.y > ym)*(x1.y - ym) + (x1.y < 0.0)*(ym + x1.y) + (x1.y > 0 && x1.y < ym)*x1.y;
        x1.z = (x1.z > zm)*(x1.z - zm) + (x1.z < 0.0)*(zm + x1.z) + (x1.z > 0 && x1.z < zm)*x1.z;

     	p->SetX(x1);


}

 #ifdef __CUDACC__
 __host__ __device__
 #endif
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
void  Kernel(CellDouble& Jx,
		   int i11,int i12,int i13,
		   int i21,int i22,int i23,
		   int i31,int i32,int i33,
		   int i41,int i42,int i43,
		   double su,double dy,double dz,double dy1,double dz1,double s1)
{
  // FORTRAN-STYLE NUMBERING ASSIGNMENT !!!!!
      double t1,t2,t3,t4;

      t1 = su*(dy1*dz1+s1);
      t2 = su*(dy1*dz-s1);
      t3 = su*(dy*dz1-s1);
      t4 = su*(dy*dz+s1);

      Jx.M[i11][i12][i13] += t1;//su*(dy1*dz1+s1);
      Jx.M[i21][i22][i23] += t2;//su*(dy1*dz-s1);
      Jx.M[i31][i32][i33] += t3;//su*(dy*dz1-s1);
      Jx.M[i41][i42][i43] += t4;//su*(dy*dz+s1);
}

__host__ __device__
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
void pqr(int3& i,double3& x,double3& x1,double& a1,double tau,CurrentTensor *t1,
		         int num,Particle *p)
{
      double dx,dy,dz,a,dx1,dy1,dz1,su,sv,sw,s1,s2,s3;
     // double xl,yl;//,res,inv_yl;

#ifdef PARTICLE_TRACE
      if(p->fortran_number == 32587 && p->sort == 2)
        		     {
        		    	 printf("pqr 32587 x1 %25.15e \n",x1.x);
        		     }
#endif

      dx=getCenterRelatedShift(x.x,x1.x,i.x,hx,x0); //0.5d0*(x+x1)-h1*(i-1.5d0)
      dy=getCenterRelatedShift(x.y,x1.y,i.y,hy,y0); //0.5d0*(y+y1)-y0-h2*(l-1.5d0)
      dz=getCenterRelatedShift(x.z,x1.z,i.z,hz,z0); //0.5d0*(z+z1)-h3*(k-1.5d0)
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,48+num)] = dx;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,50+num)] = dy;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,52+num)] = dz;
#endif
      a = a1;

      dx1=hx - dx;
      dy1=hy - dy;
      dz1=hz - dz;
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,90+num)] = dx1;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,92+num)] = dy1;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,94+num)] = dz1;
#endif

      su =x1.x - x.x;
      sv =x1.y - x.y;
      sw =x1.z - x.z;
#ifdef PARTICLE_TRACE
      if(p->fortran_number == 32587 && p->sort == 2)
        		     {
        		    	 printf("pqr su 32587 x1 %25.15e %25.15e \n",x1.x,su);
        		     }
#endif

      s1=sv*sw/12.0;
      s2=su*sw/12.0;
      s3=su*sv/12.0;
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,54+num)] = su;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,56+num)] = sv;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,58+num)] = sw;
#endif
//      xl = su*a;
//      yl = tau*hy*hz;
//#ifdef __CUDACC__
   //   res = __ddiv_rn(xl,yl);
//#else
//      res = xl /yl;
//      inv_yl = 1.0/yl;
//#endif
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,106)]     = su*a;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,107)]     = tau*hy*hz;
//      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,108)]     = res;
//      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,109)]     = inv_yl;
#endif
     // d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,110)]     = su*a/(tau*hy*hz);


      su=su*a/(tau*hy*hz);
      sv=sv*a/(tau*hx*hz);
      sw=sw*a/(tau*hx*hy);
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,96+num)]  = su;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,98+num)]  = sv;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,100+num)] = sw;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,102)]     = a;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,103)]     = tau;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,104)]     = hy;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,105)]     = hz;
#endif

      //Kernel(*Jx,i.x,i.y,i.z, i.x,i.y,i.z+1, i.x,i.y+1,i.z, i.x,i.y+1,i.z+1 ,su,dy,dz,dy1,dz1,s1);
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,60+num)] = i.x+this->i -1;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,62+num)] = i.y+this->l -1;
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,64+num)] = i.z+this->k -1;
#endif

      t1->Jx.i11 = i.x;
      t1->Jx.i12 = i.y;
      t1->Jx.i13 = i.z;
      t1->Jx.t[0] = su*(dy1*dz1+s1);

      t1->Jx.i21 = i.x;
      t1->Jx.i22 = i.y;
      t1->Jx.i23 = i.z+1;
      t1->Jx.t[1] = su*(dy1*dz-s1);

      t1->Jx.i31 = i.x;
      t1->Jx.i32 = i.y+1;
      t1->Jx.i33 = i.z;
      t1->Jx.t[2] = su*(dy*dz1-s1);

      t1->Jx.i41 = i.x;
      t1->Jx.i42 = i.y+1;
      t1->Jx.i43 = i.z+1;
      t1->Jx.t[3] = su*(dy*dz+s1);
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,66+num)] = su*(dy1*dz1+s1);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,68+num)] = su*(dy1*dz-s1);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,70+num)] = su*(dy*dz1-s1);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,72+num)] = su*(dy*dz+s1);
#endif
/*      p(i,l,k)=p(i,l,k)+su*(dy1*dz1+s1)
      p(i,l,k+1)=p(i,l,k+1)+su*(dy1*dz-s1)
      p(i,l+1,k)=p(i,l+1,k)+su*(dy*dz1-s1)
      p(i,l+1,k+1)=p(i,l+1,k+1)+su*(dy*dz+s1) */

      //Kernel(*Jy,i.x,i.y,i.z, i.x,i.y,i.z+1, i.x+1,i.y,i.z,    i.x+1,i.y,i.z+1,sv,dx,dz,dx1,dz1,s2);
      t1->Jy.i11 = i.x;
      t1->Jy.i12 = i.y;
      t1->Jy.i13 = i.z;
      t1->Jy.t[0] = sv*(dx1*dz1+s2);

      t1->Jy.i21 = i.x;
      t1->Jy.i22 = i.y;
      t1->Jy.i23 = i.z+1;
      t1->Jy.t[1] = sv*(dx1*dz-s2);

      t1->Jy.i31 = i.x+1;
      t1->Jy.i32 = i.y;
      t1->Jy.i33 = i.z;
      t1->Jy.t[2] = sv*(dx*dz1-s2);

      t1->Jy.i41 = i.x+1;
      t1->Jy.i42 = i.y;
      t1->Jy.i43 = i.z+1;
      t1->Jy.t[3] = sv*(dx*dz+s2);
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,74+num)] = sv*(dx1*dz1+s2);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,76+num)] = sv*(dx1*dz-s2);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,78+num)] = sv*(dx*dz1-s2);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,80+num)] = sv*(dx*dz+s2);
#endif
//      Kernel(Jy,i,l,k,     i,l,k+1,       i+1,l,k,          i+1,l,k+1,    sv,dx,dz,dx,dz1,s2);
/*      q(i,l,k)=q(i,l,k)+sv*(dx1*dz1+s2)
      q(i,l,k+1)=q(i,l,k+1)+sv*(dx1*dz-s2)
      q(i+1,l,k)=q(i+1,l,k)+sv*(dx*dz1-s2)
      q(i+1,l,k+1)=q(i+1,l,k+1)+sv*(dx*dz+s2) */

//      Kernel(*Jz,i.x,i.y,i.z,    i.x,i.y+1,i.z,            i.x+1,i.y,i.z,          i.x+1,i.y+1,i.z,sw,dx,dy,dx1,dy1,s3);
      t1->Jz.i11 = i.x;
      t1->Jz.i12 = i.y;
      t1->Jz.i13 = i.z;
      t1->Jz.t[0] = sw*(dx1*dy1+s3);

      t1->Jz.i21 = i.x;
      t1->Jz.i22 = i.y+1;
      t1->Jz.i23 = i.z;
      t1->Jz.t[1] = sw*(dx1*dy-s3);

      t1->Jz.i31 = i.x+1;
      t1->Jz.i32 = i.y;
      t1->Jz.i33 = i.z;
      t1->Jz.t[2] = sw*(dx*dy1-s3);

      t1->Jz.i41 = i.x+1;
      t1->Jz.i42 = i.y+1;
      t1->Jz.i43 = i.z;
      t1->Jz.t[3] = sw*(dx*dy+s3);
#ifdef ATTRIBUTES_CHECK
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,82+num)] = sw*(dx1*dy1+s3);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,84+num)] = sw*(dx1*dy-s3);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,86+num)] = sw*(dx*dy1-s3);
      d_ctrlParticles[ParticleAttributePosition(jmp,p->fortran_number,p->sort,88+num)] = sw*(dx*dy+s3);
#endif
//      Kernel(Jx,i,l,k,          i,l+1,k,                  i+1,l,k,              i+1,l+1,k,      sw,dx,dy,dx,dy1,s3);
/*      r(i,l,k)=r(i,l,k)+sw*(dx1*dy1+s3)
      r(i,l+1,k)=r(i,l+1,k)+sw*(dx1*dy-s3)
      r(i+1,l,k)=r(i+1,l,k)+sw*(dx*dy1-s3)
      r(i+1,l+1,k)=r(i+1,l+1,k)+sw*(dx*dy+s3) */


}

#ifdef __CUDACC__
 __host__ __device__
 #endif
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
void pqr(int3& i,double3& x,double3& x1,double& a1,double tau)
{
      double dx,dy,dz,a,dx1,dy1,dz1,su,sv,sw,s1,s2,s3;

      dx=getCenterRelatedShift(x.x,x1.x,i.x,hx,x0); //0.5d0*(x+x1)-h1*(i-1.5d0)
      dy=getCenterRelatedShift(x.y,x1.y,i.y,hy,y0); //0.5d0*(y+y1)-y0-h2*(l-1.5d0)
      dz=getCenterRelatedShift(x.z,x1.z,i.z,hz,z0); //0.5d0*(z+z1)-h3*(k-1.5d0)

      a = a1;

      dx1=hx - dx;
      dy1=hy - dy;
      dz1=hz - dz;
      su =x1.x - x.x;
      sv =x1.y - x.y;
      sw =x1.z - x.z;

      s1=sv*sw/12.0;
      s2=su*sw/12.0;
      s3=su*sv/12.0;

      su=su*a/(tau*hy*hz);
      sv=sv*a/(tau*hx*hz);
      sw=sw*a/(tau*hx*hy);

      Kernel(*Jx,i.x,i.y,i.z,i.x,i.y,i.z+1,i.x,i.y+1,i.z,i.x,i.y+1,i.z+1,su,dy,dz,dy1,dz1,s1);

/*      p(i,l,k)=p(i,l,k)+su*(dy1*dz1+s1)
      p(i,l,k+1)=p(i,l,k+1)+su*(dy1*dz-s1)
      p(i,l+1,k)=p(i,l+1,k)+su*(dy*dz1-s1)
      p(i,l+1,k+1)=p(i,l+1,k+1)+su*(dy*dz+s1) */

      Kernel(*Jy,i.x,i.y,i.z, i.x,i.y,i.z+1, i.x+1,i.y,i.z,    i.x+1,i.y,i.z+1,sv,dx,dz,dx1,dz1,s2);
//      Kernel(Jy,i,l,k,     i,l,k+1,       i+1,l,k,          i+1,l,k+1,    sv,dx,dz,dx,dz1,s2);
/*      q(i,l,k)=q(i,l,k)+sv*(dx1*dz1+s2)
      q(i,l,k+1)=q(i,l,k+1)+sv*(dx1*dz-s2)
      q(i+1,l,k)=q(i+1,l,k)+sv*(dx*dz1-s2)
      q(i+1,l,k+1)=q(i+1,l,k+1)+sv*(dx*dz+s2) */

      Kernel(*Jz,i.x,i.y,i.z,    i.x,i.y+1,i.z,            i.x+1,i.y,i.z,          i.x+1,i.y+1,i.z,sw,dx,dy,dx1,dy1,s3);
//      Kernel(Jx,i,l,k,          i,l+1,k,                  i+1,l,k,              i+1,l+1,k,      sw,dx,dy,dx,dy1,s3);
/*      r(i,l,k)=r(i,l,k)+sw*(dx1*dy1+s3)
      r(i,l+1,k)=r(i,l+1,k)+sw*(dx1*dy-s3)
      r(i+1,l,k)=r(i+1,l,k)+sw*(dx*dy1-s3)
      r(i+1,l+1,k)=r(i+1,l+1,k)+sw*(dx*dy+s3) */


}


#ifdef __CUDACC__
 __host__ __device__
 #endif
bool Insert(Particle& p)
{
//	 if(p.fortran_number == 14536 && ((int)p.sort == 2))
//	 {
//		 int qq = 0;
//	 }
     if(isPointInCell(p.GetX()))
     {
         addParticleToSurface(&p,&number_of_particles);
         return true;
     }
     else return false;

     //return false;
}



//Particle *GetParticle(int n)
//{
//    Particle *p;
//
////#ifdef GPU_PARTICLE
////#else
//	p = &(all_particles[n]);
////#endif
//	return p;
//}


void copyCellFromHostToDevice(Cell *d_p,Cell *h_p)
{
     int err;
	 cudaMalloc((void **)&d_p,sizeof(Cell));
     err = MemoryCopy(d_p,h_p,sizeof(Cell),HOST_TO_DEVICE);
     if(err != cudaSuccess)
	 {
		printf("copyCellFromHostToDevice err %d %s \n",err,getErrorString(err));
		exit(0);
	 }
}


int checkParticleType(Particle p,double mass,double q_mass)
{
	return ((fabs(p.m-mass) < PARTICLE_MASS_TOLERANCE) &&
	    	    		 (fabs(p.q_m-q_mass) < PARTICLE_MASS_TOLERANCE));
}

int checkParticleType(int i,double mass,double q_mass)
{
	Particle p;

	p = readParticleFromSurfaceDevice(i);

	return checkParticleType(p,mass,q_mass);
}

int getParticleTypeNumber(double mass,double q_mass)
{
	Particle p;
	int num = 0;

	for(int i = 0; i < number_of_particles;i++)
	{
		p = readParticleFromSurfaceDevice(i);
		num += checkParticleType(p,mass,q_mass);
	}
	return num;
}



#ifdef __CUDACC__
 __host__ __device__
 #endif
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
void MoveSingleParticle(unsigned int i, CellTotalField cf)
{
     Particle p;
     Field fd;

     if(i >= number_of_particles) return;
     p = readParticleFromSurfaceDevice(i);
	 fd = GetField(&p,cf);

	 p.Move(fd,tau);
	 writeParticleToSurface(i,&p);
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
 DoubleCurrentTensor AccumulateCurrentSingleParticle(unsigned int i,int *cells,DoubleCurrentTensor *dt)
 {
	 Particle p;
//	 DoubleCurrentTensor dt;
	 if(i >= number_of_particles)
	 {
		 dt->t1.Jx.t[0] = 0;
		 dt->t1.Jx.t[1] = 0;
		 dt->t1.Jx.t[2] = 0;
		 dt->t1.Jx.t[3] = 0;

		 dt->t1.Jy = dt->t1.Jx;
		 dt->t1.Jz = dt->t1.Jx;
		 dt->t2    = dt->t1;

		 return (*dt);
	 }

	 p = readParticleFromSurfaceDevice(i);
	 CurrentToMesh(tau,cells,dt,&p);

     writeParticleToSurface(i,&p);

//     dt.t1 = *t1;
//     dt.t2 = *t2;

     return (*dt);
}


#ifdef __CUDACC__
 __host__ __device__
 #endif
void SetAllCurrentsToZero(uint3 threadIdx)
{
	int i1,l1,k1;


//    for(int i1 = 0; i1 < CellExtent;i1++)
//    {
//        for(int l1 = 0; l1 < CellExtent;l1++)
//	    {
//	        for(int k1 = 0; k1 < CellExtent;k1++)
//	        {
	  i1 = threadIdx.x;
	  l1 = threadIdx.y;
	  k1 = threadIdx.z;

	        	Jx->M[i1][l1][k1] = 0;
	        	Jy->M[i1][l1][k1] = 0;
	        	Jz->M[i1][l1][k1] = 0;
//	        }
//	    }
//    }
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
 void writeToArray(double *E, CellDouble Ec,int pnum)
{
    // double g;
     int count = 0;

#ifdef GPU_PARTICLE
     count = number_of_particles;
#else
     count = all_particles.size();
#endif

     if(count > 0)
     {
#ifdef DEBUG_PLASMA
       printf("cell number (%d,%d,%d) begins writing current ============================\n",i,l,k);
#endif
     }
     else return;

  //   g = 0;

     for(int i1 = 0; i1 < CellExtent;i1++)
     {
         for(int l1 = 0; l1 < CellExtent;l1++)
	 {
	     for(int k1 = 0; k1 < CellExtent;k1++)
	     {
//	         int n = getGlobalCellNumber(i+i1-1,l+l1-1,k+k1-1);
//	         int n = getGlobalCellNumber(i+i1-2,l+l1-2,k+k1-2);
	         int n = getFortranCellNumber(i+i1-1,l+l1-1,k+k1-1);

		 double t_b = E[n];//,t_a,t;
		 //int i_f,l_f,k_f;

		 //getFortranCellTriplet(n,&i_f,&l_f,&k_f);

//		 t = Ec.M[i1][l1][k1];
		 E[n] += Ec.M[i1][l1][k1];
//		 t_a = E[n];

//#ifdef DEBUG_PLASMA
//		 if((fabs(t) > 1e-15) && (t_b > TOLERANCE) && (n> 0))
//		 {
//		     int g = 0;
//		     g = 1;
//		 }
//
		 if((fabs(Ec.M[i1][l1][k1]) > 1e-15))
		 {
		    printf("n %5d i %3d l %3d k %3d Ebefore %15.5e E %15.5e Ec %15.5e non-zero %5d \n",n,i+i1-1,l+l1-1,k+k1-1,t_b,E[n],Ec.M[i1][l1][k1],count++);
		 }
//#endif
	     }
	 }
     }


}

#ifdef __CUDACC__
 __host__ __device__
 #endif
 double getFortranArrayValue(double *E,int i,int l,int k)
{
	int n = getFortranCellNumber(i,l,k);

	//if(i == 102 && l == 2 && k == 0) printf("n %d i %d l %d k %d\n",n,i,l,k);

	if(n < 0 || n > (Nx +  2)*(Ny + 2)*(Nz + 2)) return 0.0;

	return E[n];
}

// MAPPING: fORTRAN NODE i GOES TO 2nd NODE OF C++ CELL i-1
// Simple : C++ (i+i1) ----->>>> F(i+i1-1)
#ifdef __CUDACC__
 __host__ __device__
 #endif
 void readField(double *E, CellDouble & Ec,uint3 threadIdx)
{
	int i_f,l_f,k_f;
	int i1,l1,k1;

//     for(int i1 = 0; i1 < CellExtent;i1++)
//     {
//         for(int l1 = 0; l1 < CellExtent;l1++)
//	 {
//	     for(int k1 = 0; k1 < CellExtent;k1++)
//	     {
//	         int n = getGlobalCellNumber(i+i1-1,l+l1-1,k+k1-1);
//	         int n = getGlobalCellNumber(i+i1-2,l+l1-2,k+k1-2);
//	    	 if((i+i1-1 == 102) && (l+l1-1 == 2) && (k+k1-1 == 0))
//	    	 {
//	    		 printf("QQQQ i,i1 %d,%d");
//	    	 }
	         i1 = threadIdx.x;
	         l1 = threadIdx.y;
	         k1 = threadIdx.z;
	    	 int n = getFortranCellNumber(i+i1-1,l+l1-1,k+k1-1);
#ifdef DEBUG_PLASMA
//		 if(((i == 0) ) && ((k == 9)) && ((l == 7)))
//		 {
////		    double t = E[n];
////		   t = 1;
//
//		 }
#endif
		    getFortranCellTriplet(n,&i_f,&l_f,&k_f);
		    double t = getFortranArrayValue(E,i+i1-1,l+l1-1,k+k1-1);

		 Ec.M[i1][l1][k1] = t;
//	     }
//	 }
//     }
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
void readField(double *E, CellDouble & Ec)
{
	int i_f,l_f,k_f;

     for(int i1 = 0; i1 < CellExtent;i1++)
     {
         for(int l1 = 0; l1 < CellExtent;l1++)
	 {
	     for(int k1 = 0; k1 < CellExtent;k1++)
	     {
//	         int n = getGlobalCellNumber(i+i1-1,l+l1-1,k+k1-1);
//	         int n = getGlobalCellNumber(i+i1-2,l+l1-2,k+k1-2);
//	    	 if((i+i1-1 == 102) && (l+l1-1 == 2) && (k+k1-1 == 0))
//	    	 {
//	    		 printf("QQQQ i,i1 %d,%d");
//	    	 }

	    	 int n = getFortranCellNumber(i+i1-1,l+l1-1,k+k1-1);
#ifdef DEBUG_PLASMA
//		 if(((i == 0) ) && ((k == 9)) && ((l == 7)))
//		 {
////		    double t = E[n];
////		   t = 1;
//
//		 }
#endif
		    getFortranCellTriplet(n,&i_f,&l_f,&k_f);
		    double t = getFortranArrayValue(E,i+i1-1,l+l1-1,k+k1-1);

		 Ec.M[i1][l1][k1] = t;
	     }
	 }
     }
}


#ifdef __CUDACC__
 __host__ __device__
 #endif
void writeAllToArrays(double *glob_Jx,double *glob_Jy,double *glob_Jz,double *glob_Rho,int pnum)
{
     writeToArray(glob_Jx,*Jx,pnum);
     writeToArray(glob_Jy,*Jy,pnum);
     writeToArray(glob_Jz,*Jz,pnum);
    // writeToArray(glob_Rho,*Rho,pnum);
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
void readFieldsFromArrays(double *glob_Ex,double *glob_Ey,double *glob_Ez,double *glob_Hx,double *glob_Hy,double *glob_Hz,uint3 threadIdx)
{
     readField(glob_Ex,*Ex,threadIdx);
     readField(glob_Ey,*Ey,threadIdx);
     readField(glob_Ez,*Ez,threadIdx);
     readField(glob_Hx,*Hx,threadIdx);
     readField(glob_Hy,*Hy,threadIdx);
     readField(glob_Hz,*Hz,threadIdx);
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
void readFieldsFromArrays(double *glob_Ex,double *glob_Ey,double *glob_Ez,double *glob_Hx,double *glob_Hy,double *glob_Hz)
{
     readField(glob_Ex,*Ex);
     readField(glob_Ey,*Ey);
     readField(glob_Ez,*Ez);
     readField(glob_Hx,*Hx);
     readField(glob_Hy,*Hy);
     readField(glob_Hz,*Hz);
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
thrust::host_vector< Particle >  getFlyList()
{
   thrust::host_vector<Particle> fl;

     int count;

     count = number_of_particles;

     for(int n = 0;n < count;n++)
     {
         Particle p;
         p  = readParticleFromSurfaceDevice(n);
         if(!isPointInCell(p.GetX()))
	     {
	        removeParticleFromSurfaceDevice(n,&p,&number_of_particles);
	        fl.push_back(p);
	     }
     }
     int s = fl.size();
     return fl;
}

#ifdef __CUDACC__
 __host__ __device__
 #endif
int getFortranParticleNumber(int n)
{
    Particle p;
    p = readParticleFromSurfaceDevice(n);
    return p.fortran_number;
}


void memcpy(unsigned char *tgt,unsigned char *src,int size)
{
	int i;
	for(i = 0;i < size;i++) tgt[i] = src[i];
}


#ifdef __CUDACC__
 __host__ __device__
 #endif
Cell & operator=(Cell const & src)
{
	i   = src.i;
	l   = src.l;
	k   = src.k;
	hx  = src.hx;
	hy  = src.hy;
	hz  = src.hz;
	tau = src.tau;
	x0  = src.x0;
	y0  = src.y0;
	z0  = src.z0;
	xm  = src.xm;
	ym  = src.ym;
	zm  = src.zm;

	Nx  = src.Nx;
	Ny  = src.Ny;
	Nz  = src.Nz;

	number_of_particles = src.number_of_particles;

    Jx = src.Jx;
    Jy = src.Jy;
    Jz = src.Jz;
    Ex = src.Ex;
    Ey = src.Ey;
    Ez = src.Ez;
    Hx = src.Hx;
    Hy = src.Hy;
    Hz = src.Hz;
    Rho = src.Rho;

    doubParticleArray = src.doubParticleArray;

    return (*this);
}

double compareParticleLists(Cell *c)
{
	   char s[100];
	   double res;

	   if((number_of_particles == 0) && number_of_particles == c->number_of_particles) return 1.0;

       if(number_of_particles != c->number_of_particles) return 0;

#ifdef PARTICLE_PRINTS
       printf("cell %5d %5d %5d with  particles %5d check particle begins =============================================================\n",i,l,k,number_of_particles);
       for(int i = 0;i < number_of_particles;i++)
       {
    	   Particle p,p1;

    	  p =  readParticleFromSurfaceDevice(i,&p);
    	   c->readParticleFromSurfaceDevice(i,&p1);

    	   printf("particle %5d X %10.3e %10.3e Y %10.3e %10.3e Z %10.3e %10.3e \n \
    			                Ex %10.3e %10.3e Ey %10.3e %10.3e Ez %10.3e %10.3e \n \
             	                Hx %10.3e %10.3e Hy %10.3e %10.3e Hz %10.3e %10.3e \n",
    			   i,
    			   p.x,p1.x,p.y,p1.y,p.z,p1.z,
    			   p.ex,p1.ex,p.ey,p1.ey,p.ez,p1.ez,
    			   p.hx,p1.hx,p.hy,p1.hy,p.hz,p1.hz
    			   );

       }
#endif

  //     sprintf(s,"cell %5d %5d %5d particles %5d",i,l,k,number_of_particles);

       res = compare(doubParticleArray,c->doubParticleArray,
    		   number_of_particles*sizeof(Particle)/sizeof(double),s,ABSOLUTE_TOLERANCE);

    //   printf("%s %.2f check particle ended ==============================================================================================\n",s,res);

       return res;
}

double
#ifdef VIRTUAL_FUNCTIONS
virtual
#endif
compareToCell(Cell & src)
{
	//double t,t1,t_jx,t_jy,t_jz;
	double t_ex,t_ey,t_ez,t_hx,t_hy,t_hz,t_rho,t2,t;
	int size = sizeof(CellDouble)/sizeof(double);

	t = (i   == src.i) +
	    (l   == src.l) +
	    (k   == src.k) +
	    comd(hx,src.hx) +
	    comd(hy,src.hy) +
	    comd(hz,src.hz) +
	    comd(tau,src.tau) +
	    comd(x0,src.x0) +
	    comd(y0,src.y0) +
	    comd(z0,src.z0) +
	    comd(xm,src.xm) +
	    comd(ym,src.ym) +
	    comd(zm,src.zm) +
	    (Nx  == src.Nx) +
	    (Ny  == src.Ny) +
	    (Nz  == src.Nz);
	t /= 16.0;

	double t1 = compareParticleLists(&src);

	if(isNan(t1))
	{
		t1 = compareParticleLists(&src);
	}
	if((t1 < 1.0))
	{
		t1 = compareParticleLists(&src);
	}

//	t_jy =    compare((double*)Jy,(double*)src.Jy,size,"Jy");
//	t_jx =    compare((double*)Jx,(double*)src.Jx,size,"Jx");
//    t_jz =    compare((double*)Jz,(double*)src.Jz,size,"Jz");

    t_ex =    compare((double*)Ex,(double*)src.Ex,size,"Ex",TOLERANCE);
    t_ey =    compare((double*)Ey,(double*)src.Ey,size,"Ey",TOLERANCE);
	t_ez =    compare((double*)Ez,(double*)src.Ez,size,"Ez",TOLERANCE);

    t_hx =    compare((double*)Hx,(double*)src.Hx,size,"Hx",TOLERANCE);
    t_hy =    compare((double*)Hy,(double*)src.Hy,size,"Hy",TOLERANCE);
	t_hz =    compare((double*)Hz,(double*)src.Hz,size,"Hz",TOLERANCE);

	t_rho =    compare((double*)Rho,(double*)src.Rho,size,"Rho",TOLERANCE);

	t2 = t_ex+t_ey+t_ez+t_hx+t_hy+t_hz+t_rho;
	if(isNan(t2) || (t2/7.0 < 1.0))
	{
		t_ex =    compare((double*)Ex,(double*)src.Ex,size,"Ex",TOLERANCE);
		    t_ey =    compare((double*)Ey,(double*)src.Ey,size,"Ey",TOLERANCE);
			t_ez =    compare((double*)Ez,(double*)src.Ez,size,"Ez",TOLERANCE);

		    t_hx =    compare((double*)Hx,(double*)src.Hx,size,"Hx",TOLERANCE);
		    t_hy =    compare((double*)Hy,(double*)src.Hy,size,"Hy",TOLERANCE);
			t_hz =    compare((double*)Hz,(double*)src.Hz,size,"Hz",TOLERANCE);

			t_rho =    compare((double*)Rho,(double*)src.Rho,size,"Rho",TOLERANCE);
	}

   // return (t+t1+t_jx+t_jy+t_jz+t_ex+t_ey+t_ez+t_hx+t_hy+t_hz+t_rho)/12.0;
    return (t+t1+t_ex+t_ey+t_ez+t_hx+t_hy+t_hz+t_rho)/9.0;
}

double checkCellParticles(int check_point_num,double *x,double *y,double *z,
		                  double *px,double *py,double *pz,double q_m,double m)
{
	int i,num = 0,j,num_sort = 0,correct_particle;
	double t,dm,dqm,dx,dy,dz,dpx,dpy,dpz;
	if(number_of_particles < 0 || number_of_particles > MAX_particles_per_cell)
	{
//		int qq = 0;
	}

	for(i = 0; i < number_of_particles;i++)
	{
		Particle p;

		p = readParticleFromSurfaceDevice(i);

		j = p.fortran_number - 1;
		dm  = fabs(p.m   -   m);
		dqm = fabs(p.q_m - q_m);
		if((dm    > ABSOLUTE_TOLERANCE) || (dqm   > ABSOLUTE_TOLERANCE)) continue;


		dx  = fabs(p.x - x[j]);
		dy  = fabs(p.y - y[j]);
		dz  = fabs(p.z - z[j]);
		dpx = fabs(p.pu - px[j]);
		dpy = fabs(p.pv - py[j]);
		dpz = fabs(p.pw - pz[j]);

		num_sort    +=  (dm    < ABSOLUTE_TOLERANCE) &&
						(dqm   < ABSOLUTE_TOLERANCE);

		if((fabs(m) < 1e-3) && (check_point_num >= 270))
		{
			//int qq = 0;
//		   if((dm    < ABSOLUTE_TOLERANCE) &&
//					(dqm   < ABSOLUTE_TOLERANCE) &&
//					(dx    < PARTICLE_TOLERANCE) &&
//					(dy    < PARTICLE_TOLERANCE) &&
//					(dz    < PARTICLE_TOLERANCE) &&
//					(dpx   < PARTICLE_TOLERANCE) &&
//					(dpy   < PARTICLE_TOLERANCE) &&
//					(dpz   < PARTICLE_TOLERANCE))
//		   {
//			   int qq = 0;
//		   }
		}

		correct_particle =  (dm    < ABSOLUTE_TOLERANCE) &&
			            	(dqm   < ABSOLUTE_TOLERANCE) &&
				            (dx    < PARTICLE_TOLERANCE) &&
				            (dy    < PARTICLE_TOLERANCE) &&
				            (dz    < PARTICLE_TOLERANCE) &&
				            (dpx   < PARTICLE_TOLERANCE) &&
				            (dpy   < PARTICLE_TOLERANCE) &&
				            (dpz   < PARTICLE_TOLERANCE);

		if(!correct_particle && check_point_num == 270 && (int)p.sort == 1)
		{
//			double x0 = x[j];
//			int qq1 = 0;
//			qq1 = 1;
		}
		num += correct_particle;
	}

//	if((fabs(m) < 1e-3) && (check_point_num >= 270) && (num > 0))
//			{
//		      // int qq = 0;
//			}



	if(num_sort > 0)
	{
	   t = ((double)num)/num_sort;
	}
	else
	{
		t = 1.0;
	}

	if((t > 0) && (fabs(m) < 1e-3) && (check_point_num >= 270) && (num > 0))
	{
		//int qq = 0;
	}

	return t;
}

void SetControlSystem(int j,double *c)
{
	Particle p;
 //   int i;

	jmp = j;
#ifdef ATTRIBUTES_CHECK
	d_ctrlParticles = c;
#endif

//	 for(i = 0;i < number_of_particles;i++)
//	 {
//          readParticleFromSurfaceDevice(i,&p);
//          p.SetControlSystem(j,c);
//	 }
}

#ifdef __CUDACC__
 __host__ __device__
#endif
void SetControlSystemToParticles()
{
	Particle p;
    int i;

    for(i = 0;i < number_of_particles;i++)
    {
        p = readParticleFromSurfaceDevice(i);
#ifdef ATTRIBUTES_CHECK
   //     p.SetControlSystem(jmp,d_ctrlParticles);
#endif
    }
}

};




//typedef double (Cell::*cell_work_function)(void);


#endif

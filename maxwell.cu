#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//#include "rnd.h"

#include "read_particles.h"

#include "run_control.h"

#include <vector>

#include "particle.h"

#include "maxwell.h"


#include "f2c.h"
#include <stdio.h>
//#include "rnd.h"
// #include "run_control.h"

/* Common Block Declarations */

ParticleArraysGroup initial;
ParticleFloatArraysGroup diagnostics;

struct cag05b_1_ {
    doublereal store1, store2;
};
struct cag05b_2_ {
    doublereal normal, gamma;
};

#define cag05b_1 (*(struct cag05b_1_ *) &cag05b_)
#define cag05b_2 (*(struct cag05b_2_ *) &cag05b_)

struct cag05a_1_ {
    integer ix, iy, iz;
};

#define cag05a_1 (*(struct cag05a_1_ *) &cag05a_)

/* Initialized data */

struct cag05b {
    doublereal e_1[2];
    } cag05b_ = { 1., -1. };

struct cag05a {
    integer e_1[3];
    } cag05a_ = { 1, 255, 25555 };


/* Table of constant values */

//static integer c__9 = 9;
//static integer c__1 = 1;
//static integer c__3 = 3;
//static integer c__5 = 5;
//static integer c__0 = 0;
//static integer c__2 = 2;
//static doublereal c_b172 = 1.;
//static integer c__12 = 12;
//static doublereal c_b417 = 1.1424;
//static doublereal c_b419 = .5712;
//static doublereal c_b429 = 1.5;
//static integer c__20000 = 20000;
//static doublereal c_b581 = 0.;
//static doublereal c_b582 = .11200000000000002;
//static doublereal c_b587 = .14;
//static doublereal c_b589 = .8;
//static doublereal c_b614 = .001;

/* ------------------------------------------------------ */
/* ������� �.�., */
doublereal g05dde_(doublereal *a, doublereal *b,int dbg_print)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal half = .5;
    static doublereal d__[41] = { 0.,.674489750196082,1.150349380376008,
	    1.534120544352546,1.862731867421652,2.153874694061456,
	    2.417559016236505,2.66006746861746,2.885634912426757,
	    3.097269078198785,3.297193345691964,3.487104104114431,
	    3.668329285121323,3.841930685501911,4.008772594168585,
	    4.169569323349106,4.324919040826046,4.475328424654204,
	    4.621231001499247,4.763001034267814,4.900964207963193,
	    5.035405969463927,5.166578119728753,5.294704084854598,
	    5.419983174916868,5.54259405780294,5.662697617459439,
	    5.780439324478935,5.895951216739571,6.009353565530745,
	    6.120756285971941,6.230260137989044,6.33795775455379,
	    6.443934526538564,6.548269367831731,6.651035379893011,
	    6.752300431407015,6.852127665896068,6.95057594791675,
	    7.047700256664409,7.14355203435219 };

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer n;
    static doublereal t, u, v, w, x;
    extern doublereal wrapg05cae_(doublereal *,int);

    u = cag05b_1.store1;
    for (n = 1; n <= 39; ++n) {
	if (u > half) {
	    goto L40;
	}
	u += u;
/* L20: */
    }
    n = 40;
L40:
    t = d__[n - 1];
    u = wrapg05cae_(&x,dbg_print);
L60:
    w = (d__[n] - t) * u;
    v = w * (w * half + t);
L80:
    u = wrapg05cae_(&x,dbg_print);
    if (v <= u) {
	goto L100;
    }
    v = wrapg05cae_(&x,dbg_print);
    if (u > v) {
	goto L80;
    }
    u = (v - u) / (one - u);
    goto L60;
L100:
    u = (u - v) / (one - v);
    if (u > half) {
	goto L120;
    }
    cag05b_1.store1 = u + u;
    ret_val = *a + *b * (w + t);
    return ret_val;
L120:
    cag05b_1.store1 = u + u - one;
    ret_val = *a - *b * (w + t);
    return ret_val;
} /* g05dde_ */

/* ------------------------------------------------------------------ */
doublereal g05cae_(doublereal *x,int dbg_print)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal ai;
    static integer ii;
    static doublereal ax, ay, az;

    cag05a_1.ix = (cag05a_1.ix - cag05a_1.ix / 177 * 177) * 171 - (
	    cag05a_1.ix / 177 << 1);
    cag05a_1.iy = (cag05a_1.iy - cag05a_1.iy / 176 * 176) * 172 - (
	    cag05a_1.iy / 176 << 1);
    cag05a_1.iz = (cag05a_1.iz - cag05a_1.iz / 178 * 178) * 170 - (
	    cag05a_1.iz / 178 << 1);
    if (cag05a_1.ix < 0) {
	cag05a_1.ix += 30269;
    }
    if (cag05a_1.iy < 0) {
	cag05a_1.iy += 30307;
    }
    if (cag05a_1.iz < 0) {
	cag05a_1.iz += 30323;
    }
    ax = (doublereal) cag05a_1.ix;
    ay = (doublereal) cag05a_1.iy;
    az = (doublereal) cag05a_1.iz;
    ai = ax / 30269. + ay / 30307. + az / 30323.;
    ii = (integer) ai;
    ret_val = ai - ii;
    return ret_val;
} /* g05cae_ */

doublereal wrapg05cae_(doublereal *x,int dbg_print)
{
    static int n = 0;
    double t  = g05cae_(x,dbg_print);
// #ifdef DEBUG_PLASMA
    n++;
    if(dbg_print == 1)
    {
     printf("%10d %25.15e \n",n,t);
    }
// #endif

    return t;
}

doublereal wrapg05dde_(doublereal *a, doublereal *b,int dbg_print)
{
    return g05dde_(a,b,dbg_print);
}

double rnd_uniform(int dbg_print)
{
    doublereal x;

    return (double)wrapg05cae_(&x,dbg_print);
}

double rnd_gaussian(double a,double b,int dbg_print)
{
    return (double)wrapg05dde_(&a,&b,dbg_print);
}


int in_range(double z0,double z,double z1)
{
	return ((z > z0) && (z < z1)) || ((fabs(z - z0) < 1e-13) && (fabs(z - z1) < 1e-13));
}

int AllocateBinaryParticleArraysOneSort(
    		  double **dbg_x,
    		  double **dbg_y,
    		  double **dbg_z,
    		  double **dbg_px,
    		  double **dbg_py,
    		  double **dbg_pz,
    		  double **m,
    		  int total_particles
    		  )
      {
	         *dbg_x = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);

	         *dbg_y = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);

	         *dbg_z = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);

	         *dbg_px = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);

	         *dbg_py = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);

	         *dbg_pz = (double *)malloc(sizeof(double)*total_particles);

	         *m      = (double *)malloc(sizeof(double)*total_particles);

	         //debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);

		 	return 0;
      }

int AllocateBinaryParticleArraysOneSortFloat(
    		  float **dbg_x,
    		  float **dbg_y,
    		  float **dbg_z,
    		  float **dbg_px,
    		  float **dbg_py,
    		  float **dbg_pz,
    		  int total_particles
    		  )
      {
	         *dbg_x = (float*)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);

	         *dbg_y = (float*)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);

	         *dbg_z = (float *)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);

	         *dbg_px = (float *)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);

	         *dbg_py = (float *)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);

	         *dbg_pz = (float *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);

		 	return 0;
      }




void AllocateBinaryParticlesArrays(
		  ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons)
{
   AllocateBinaryParticleArraysOneSort(&(ions->dbg_x), &(ions->dbg_y) ,&(ions->dbg_z),
  		                                         &(ions->dbg_px),&(ions->dbg_py),
  		                                         &(ions->dbg_pz),
  		                                         &(ions->m),
  		                                          ions->total);

   AllocateBinaryParticleArraysOneSort(&(electrons->dbg_x),&(electrons->dbg_y), &(electrons->dbg_z),
  		                                             &(electrons->dbg_px),&(electrons->dbg_py),&(electrons->dbg_pz),
  		                                             &(electrons->m),
  		                                             electrons->total);

   AllocateBinaryParticleArraysOneSort(&(beam_electrons->dbg_x),&(beam_electrons->dbg_y),&(beam_electrons->dbg_z),
  		                                             &(beam_electrons->dbg_px),&(beam_electrons->dbg_py),&(beam_electrons->dbg_pz),
  		                                             &(beam_electrons->m),
  		                                             beam_electrons->total);

   //magf = 1;
}

void AllocateBinaryParticlesArraysFloat(
		  ParticleFloatArrays *ions,ParticleFloatArrays *electrons,ParticleFloatArrays *beam_electrons)
{
   AllocateBinaryParticleArraysOneSortFloat(&(ions->dbg_x), &(ions->dbg_y) ,&(ions->dbg_z),
  		                                         &(ions->dbg_px),&(ions->dbg_py),&(ions->dbg_pz),
  		                                          ions->total);


   AllocateBinaryParticleArraysOneSortFloat(&(electrons->dbg_x),&(electrons->dbg_y), &(electrons->dbg_z),
  		                                             &(electrons->dbg_px),&(electrons->dbg_py),&(electrons->dbg_pz),
  		                                             electrons->total);

   AllocateBinaryParticleArraysOneSortFloat(&(beam_electrons->dbg_x),&(beam_electrons->dbg_y),&(beam_electrons->dbg_z),
  		                                             &(beam_electrons->dbg_px),&(beam_electrons->dbg_py),&(beam_electrons->dbg_pz),
  		                                             beam_electrons->total);

   //magf = 1;
}


int writeParamsFile(double tex0,double tey0,double tez0,
                         double Tb,double rimp,
                         double rbd,double ni,
	                      double lx,double ly,double lz,
	                      int lp,int nx,int ny,int nz,
	                      double tau,double B0,
	                      double bx,double by,double bz,
	                      double py,double pz,
	                      int beam_plasma,int start_from_file,
	                      int ts,int ms,int phase
	                      )
	  {
		  FILE *f;

		  if((f = fopen("000_params.dat","wt")) == NULL) return 1;

//		  fprintf(f,"\n");
		  fprintf(f,"%15.5e plasma electron temperature along X \n ",tex0);
		  fprintf(f,"%15.5e plasma electron temperature along Y \n ",tey0);
		  fprintf(f,"%15.5e plasma electron temperature along Z \n ",tez0);
		  fprintf(f,"%15.5e beam impulse\n ",rimp);
		  fprintf(f,"%15.5e beam velocity dispersion \n ",Tb);
		  fprintf(f,"%15.5e beam and plasma density ratio \n ",rbd);
		  fprintf(f,"%15.5e plasma density \n ",ni);
		  fprintf(f,"%15.5e external magnetic field (along X) \n ",B0);
		  fprintf(f,"%15.5e domain size X \n ",lx);
		  fprintf(f,"%15.5e domain size Y \n ",ly);
		  fprintf(f,"%15.5e domain size Z \n ",lz);
		  fprintf(f,"%15.5e plasma size Y \n ",py);
		  fprintf(f,"%15.5e plasma size Z \n ",pz);
		  fprintf(f,"%15.5e beam size X \n ",bx);
		  fprintf(f,"%15.5e beam size Y \n ",by);
		  fprintf(f,"%15.5e beam size Z \n ",bz);
		  fprintf(f,"%15d   average number of particles in cell \n ",lp);
		  fprintf(f,"%15d   number of mesh nodes along X \n ",nx);
		  fprintf(f,"%15d   number of mesh nodes along Y \n ",ny);
		  fprintf(f,"%15d   number of mesh nodes along Z \n ",nz);
		  fprintf(f,"%15.5e timestep \n ",tau);
		  fprintf(f,"%15d   1 if beam-plasma interaction, 0 if beam-beam \n",beam_plasma);
		  fprintf(f,"%15d   moment to start from saved\n ",start_from_file);
		  fprintf(f,"%15d   phase to start from save\n ",phase);
		  fprintf(f,"%15d   total steps \n ",ts);
		  fprintf(f,"%15d   number of steps between diagnostic files \n ",ms);


         fclose(f);
         return 0;
	  }

int InitUniformMaxwellianParticles(int beamf,int jmb,
				   double tex0,double tey0,double tez0,
				   double beam_lx, double beam_ly,double beam_lz,int *jmb_real,
                   double lx,double ly,double lz, int meh,double Tb,double rimp,double rbd,
				   double *xi,double *yi, double *zi,double *ui,double *vi, double *wi,
				   double *xb,double *yb, double *zb,double *ub,double *vb, double *wb,
				   double *xf,double *yf, double *zf,double *uf,double *vf, double *wf
				  )
{
    double x,y,z,vb0,d__1,d__2,d__3,vy,vz,termx,gb0;
    double vf01,vf02,pinv1,pinv2,mfrq = 0.0;
//     double *ux,*uy,*uz;
    double *ux,*uy,*uz;
    double beam_y_max,beam_y_min, beam_sh;
    double beam_z_max,beam_z_min, beam_shz;

    beam_sh = (ly - beam_ly)/2;
    beam_y_max = ly - beam_sh;
    beam_y_min = beam_sh;

    beam_shz   = (lz - beam_lz)/2;
    beam_z_max = lz - beam_shz;
    beam_z_min = beam_shz;

    int j;
    
    ux = (double *)malloc(jmb*sizeof(double));
    uy = (double *)malloc(jmb*sizeof(double));
    uz = (double *)malloc(jmb*sizeof(double));
    
    for (j = 1; j <= jmb; j++) 
    {
	z =   lz * rnd_uniform(0);
	y =   meh * ly + ly * rnd_uniform(0);
	x =   lx * rnd_uniform(0);
	
	xi[j - 1] = x;
	yi[j - 1] = y;
	zi[j - 1] = z;
	ui[j - 1] = 0.0;
	vi[j - 1] = 0.0;
	wi[j - 1] = 0.0;
    }
    

//    parasetrandombeam_();
/*     END RANDOM GENERATOR */
/* ****************** BEAM **************************************** */
    *jmb_real = 0;
	for (j = 1; j <= jmb; j++) 
	{
		double y = yi[j- 1];
		double z = zi[j- 1];
		if((xi[j - 1] < beam_lx) &&
				(y < beam_y_max) && (y > beam_y_min) &&
				in_range(beam_z_min,z,beam_z_max)
		  )
		{
	        xb[*jmb_real] = xi[j - 1];
	        yb[*jmb_real] = yi[j - 1];
	        zb[*jmb_real] = zi[j - 1];
	        vb0       = rnd_gaussian(0.0, Tb*rimp,0);
	        ux[*jmb_real] = vb0 + rimp;
	        uy[*jmb_real] = rnd_gaussian(0.0, Tb*rimp,0);
	        uz[*jmb_real] = rnd_gaussian(0.0, Tb*rimp,0);
#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
            	printf("ion %10d %25.15e %25.15e %25.15e \n",j,xi[j - 1],yi[j - 1],zi[j - 1]);
#endif
	        (*jmb_real)++;
		}
	    
	}
	
	//1st beam particle impulse:    0.20296063288436139
	for (j = 1; j <= *jmb_real; j++)
	{
//	    double uxt,ubt;
	    d__1 = ux[j - 1];
	    d__2 = uy[j - 1];
	    d__3 = uz[j - 1];
	    
	    vb0 = sqrt(1.0 - d__1 * d__1 - d__2 * d__2 - d__3 * d__3);

	    ub[j - 1] = ux[j - 1] / vb0;

	    double t = fabs(ub[j - 1]-0.20296063288436139);

	    vb[j - 1] = uy[j - 1] / vb0;
	    wb[j - 1] = uz[j - 1] / vb0;
#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
	       printf("beam %10d %25.15e  %25.15e    %25.15e %25.15e %25.15e   %25.15e \n",
	    		         j,  xb[j - 1],yb[j - 1],vb0,    ub[j-1],vb[j - 1],wb[j - 1]);
#endif
	}
/*     MAKING THE RANDOM GENERATOR WORK THE SAME */
//    parasetrandomelectrons_();
/*     END RANDOM GENERATOR */

//     vy = rnd_gaussian(0.0,tey0,1);
//     vz = rnd_gaussian(0.0,tez0,1);
//  
//     termx = rnd_gaussian(0.0,tex0,1); 
	j = 1;
    for (j = 1; j <= jmb;j++) 
    {
    	   if((2*j-1) == 24933)
    	   {
//    		   int qq = 0;
    	   }
           xf[2*j-1-1] = xi[j-1];
           yf[2*j-1-1] = yi[j-1];
           zf[2*j-1-1] = zi[j-1];

           xf[2*j-1]   = xi[j-1];
           yf[2*j-1]   = yi[j-1];
           zf[2*j-1]   = zi[j-1];
	   
//         FIRST SETTING TRANVERSE
// razbros v skorostyax
           
           vy=rnd_gaussian(0.0,tey0,0);    
           vz=rnd_gaussian(0.0,tez0,0);        

//          INVERSE CURRENT
           
           termx = rnd_gaussian(0.0,tex0,0);
// 	   printf("termx %15.5e vx %15.5e vy %15.5e \n",termx,vy,vz);
	   
           gb0 = pow(1.0+pow(ub[j-1],2)+pow(vb[j-1],2)+pow(wb[j-1],2),-0.5);
	   
           vb0=ub[j-1]*gb0;
      	   if ((beamf == 1) && ((xi[j - 1] < beam_lx) && (yi[j - 1] < beam_y_max) && (yi[j - 1] > beam_y_min)))
     	   {
        	    vf01=-rbd*vb0+termx;
                vf02=-rbd*vb0-termx;
	       }
       	   else
	       {
      	        vf01=+termx;
      	        vf02=-termx;
	       }
       	   
           pinv1= vf01*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5); 
           pinv2= vf02*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
            
          vf[2*j-2] =  vy*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
// 	  write(37,127) 2*j-1,vy,vf01,vz,                                                                                                                                                                        
//      +    vy/dsqrt((1d0-vf01**2-vy**2-vz**2))
//           printf("%10d vy %15.5e vf01 %15.5e vz %15.5e vf %15.5e \n",2*j-1,vy,vf01,vz,vy*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5));
          vf[2*j-1] = -vy*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
	  
          wf[2*j-2] =  vz*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
          wf[2*j-1] = -vz*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);

	      uf[2*j-2] = pinv1+0.01*sin(mfrq*2.0*M_PI*xf[2*j-2]/lx);
          uf[2*j-1] = pinv2+0.01*sin(mfrq*2.0*M_PI*xf[2*j-1]/lx);    
          
// c my correct end
//           printf("j %10d ub(j) %15.5e uf(j) %15.5e gb0 %15.5e vb0 %15.5e vf0 %15.5e pinv  %15.5e termx %15.5e\n",
// 		  j,     ub[j-1],     uf[j-1],     gb0,       vb0,       vf01,      pinv1,       termx);
	  
#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
	       printf("electron %10d %25.15e %25.15e %25.15e %25.15e \n",2*j-2,yi[j - 1],uf[2*j-2],vf[2*j-2],wf[2*j-2]);
	       printf("electron %10d %25.15e %25.15e %25.15e %25.15e \n",2*j-1,yi[j - 1],uf[2*j-1],vf[2*j-1],wf[2*j-1]);
#endif
	  }

    free(ux);
    free(uy);
    free(uz);

#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
//    exit(0);
#endif

    return 0;
} /* start_ */

int AddBeamParticles(int jmb,
				   double tex0,double tey0,double tez0,
				   double beam_lx, double beam_ly,int *jmb_real,
                   double lx,double ly,double lz, int meh,double Tb,double rimp,double rbd,
				   double *xb,double *yb, double *zb,double *ub,double *vb, double *wb
				  )
{
    double x,y,z;//,vb0,d__1,d__2,d__3,vy,vz,termx,gb0;
//    double vf01,vf02,pinv1,pinv2,mfrq;
//     double *ux,*uy,*uz;
//    double *ux,*uy,*uz;
//    double beam_y_max;
    double beam_y_min, beam_sh;

    beam_sh = (ly - beam_ly)/2;
//    beam_y_max = ly - beam_sh;
    beam_y_min = beam_sh;

    for (int j = 1; j <= jmb; j++)
    {
    	z =   beam_ly * rnd_uniform(0) + beam_y_min;
	    y =   meh * ly + beam_ly * rnd_uniform(0) + beam_y_min;
	    x =   beam_lx * rnd_uniform(0);

	    xb[j - 1] = x;
	    yb[j - 1] = y;
	    zb[j - 1] = z;

	    double vb0,ux,uy,uz;
        vb0       = rnd_gaussian(0.0, Tb*rimp,0);
        ux        = vb0 + rimp;
        uy        = rnd_gaussian(0.0, Tb*rimp,0);
        uz        = rnd_gaussian(0.0, Tb*rimp,0);

        vb0 = sqrt(1.0 - ux*ux - uy * uy - uz * uz);

        ub[j - 1] = ux / vb0;
        vb[j - 1] = uy / vb0;
        wb[j - 1] = uz / vb0;


#ifdef ADD_BEAM_INITIAL_PARTICLE_PRINTS
	       printf("add beam %10d %25.15e  %25.15e    %25.15e %25.15e %25.15e   %25.15e \n",
	    		         j,  xb[j - 1],yb[j - 1],vb0,    ub[j-1],vb[j - 1],wb[j - 1]);
#endif
    }
	return 0;
}

int getMassCharge(ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons,
		double ni,double rbd,int lp)
{
    //int lp = ((double)N)/(Nx*Ny*Nz);
	electrons->m[0]      = -ni/lp/2.0;                 //!!!!!!
	ions->m[0]           =  (ni+rbd)/lp;
	beam_electrons->m[0] =  -rbd/lp;

	electrons->q_m        = -1.0;
	ions->q_m             =  1.0/1836.0;
	beam_electrons->q_m   = -1.0;

	return 0;
}

int AllocateMemoryForArrays(int N,ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons)
{

	ions->total           = N;
    electrons->total      = 2*N;
    beam_electrons->total = N;

//    sorts = 3;

    AllocateBinaryParticlesArrays(ions,electrons,beam_electrons);
//    AllocateBinaryParticlesArraysFloat(&(diagnostics[0]),&(diagnostics[1]),&(diagnostics[2]));
    return 0;
}

int convertParticleArraysToSTLvector(
		  double *dbg_x,
		  double *dbg_y,
		  double *dbg_z,
		  double *dbg_px,
		  double *dbg_py,
		  double *dbg_pz,
		  double q_m,
		  double m,
		  int total_particles,
		  particle_sorts sort,
		  std::vector<Particle> & vp
		  )
{
	  double x,y,z,px,py,pz;

	  for(int i = 0; i < total_particles;i++)
	  {
		  x   = dbg_x[i];
		  y   = dbg_y[i];
		  z   = dbg_z[i];
		  px   = dbg_px[i];
		  py   = dbg_py[i];
		  pz   = dbg_pz[i];


		  Particle p(x,y,z,px,py,pz,m,q_m);

		  p.fortran_number = i+1;
		  p.sort = sort;
		  p.direction = 0;

		  vp.push_back(p);

	  }
	  int size = vp.size();

	  return 0;
}


int getUniformMaxwellianParticles(std::vector<Particle>  & ion_vp,
		                           std::vector<Particle>  & el_vp,
		                           std::vector<Particle>  & beam_vp)
{
	ParticleArrays ions,electrons,beam;

    int total = 1600000,jmb;

    double tex0 = 1e-3;
    double tey0 = 1e-3;
    double tez0 = 1e-3;
//    double tol  = 1e-15;

    double Tb   = 0.14;
    double rimp = 0.2;
    double rbd  = 2.0e-3;
    double ni   = 1.0;
    int    meh  = 0;
    int    lp   = 1000;
    double lx   = 1.1424;
    double ly   = 0.05;
    double lz   = 0.05;

    AllocateMemoryForArrays(total,&ions,&electrons,&beam);

    getMassCharge(&ions,&electrons,&beam,ni,rbd,lp);

	InitUniformMaxwellianParticles(1,total,tex0,tey0,tez0,
					  lx,ly,lz,
					  &jmb,
					  lx,ly,lz,
					  meh,Tb,rimp,rbd,
						 ions.dbg_x,ions.dbg_y,ions.dbg_z,
						 ions.dbg_px,ions.dbg_py,ions.dbg_pz,
						  beam.dbg_x,beam.dbg_y,beam.dbg_z,
						  beam.dbg_px,beam.dbg_py,beam.dbg_pz,
						  electrons.dbg_x,electrons.dbg_y,electrons.dbg_z,
						  electrons.dbg_px,electrons.dbg_py,electrons.dbg_pz);

	 convertParticleArraysToSTLvector(
			 beam.dbg_x,
			 beam.dbg_y,
			 beam.dbg_z,
			 beam.dbg_px,
			 beam.dbg_py,
			 beam.dbg_pz,
			 beam.q_m,
			 *(beam.m),
	 		 beam.total,
	 		 BEAM_ELECTRON,
	 		 beam_vp);

	 convertParticleArraysToSTLvector(
			 ions.dbg_x,
			 ions.dbg_y,
			 ions.dbg_z,
			 ions.dbg_px,
			 ions.dbg_py,
			 ions.dbg_pz,
			 ions.q_m,
			 *(ions.m),
	 		 ions.total,
	 		 ION,
	 		 ion_vp);

	 convertParticleArraysToSTLvector(
	 			 electrons.dbg_x,
	 			 electrons.dbg_y,
	 			 electrons.dbg_z,
	 			 electrons.dbg_px,
	 			 electrons.dbg_py,
	 			 electrons.dbg_pz,
	 			 electrons.q_m,
	 			 *(electrons.m),
	 	 		 electrons.total,
	 	 		 PLASMA_ELECTRON,
	 	 		 el_vp);

	return 0;


}

/*
 * read_particle.cxx
 *
 *  Created on: Jun 9, 2018
 *      Author: snytav
 */

#include <stdio.h>
#include <stdlib.h>

#include <sys/resource.h>
#include <stdint.h>

#include <sys/sysinfo.h>
#include <sys/time.h>

#include <string>


#include "load_data.h"

std::string getMumuFileName(int nt)
	  {
		  char part_name[100];
		  std::string s;

		  sprintf(part_name,"mumu000%08d.dat",nt);

		  s = part_name;

		  return s;
	  }


int readFortranBinaryArray(FILE *f, double* d)
	{
//	    char str[100];
//	    Cell<Particle>  c = (*AllCells)[0];
	    int t,err;//,n;
//	    double t0;


	    //sprintf(fname,"%s_fiel3d.dat",name);
	    fread(&t,sizeof(int),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

	    fread(d,1,t,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

//	    t0 = d[269];
//	    t0 = d[270];
	    fread(&t,sizeof(int),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }
   	    return t;
	}

FILE *readPreliminary3Darrays(std::string fn,int nt,int nx,int ny,int nz)
{
	     double *buf;
	     FILE *f;

	     buf = (double *)malloc(sizeof(double)*(nx+2)*(ny+2)*(nz+2));

	     if((f = fopen(fn.c_str(),"rb")) == NULL) return NULL;
	     struct sysinfo info;

	     sysinfo(&info);
	     printf("before1  %d free %lu \n",nt,info.freeram/1024/1024);

	     //reading electric field Ex
	     readFortranBinaryArray(f,buf);
	     sysinfo(&info);
	     printf("before1  %d free %lu \n",nt,info.freeram/1024/1024);

	     //reading electric field Ey
	     readFortranBinaryArray(f,buf);
	     sysinfo(&info);
	     printf("before1  %d free %lu \n",nt,info.freeram/1024/1024);

	     //reading electric field Ez
	     readFortranBinaryArray(f,buf);
	     sysinfo(&info);
	     printf("before1  %d free %lu \n",nt,info.freeram/1024/1024);

	     //reading magnetic field Hx
	     readFortranBinaryArray(f,buf);
	     sysinfo(&info);
	     printf("before1  %d free %lu \n",nt,info.freeram/1024/1024);

	     //reading magnetic field Hy
	     readFortranBinaryArray(f,buf);

	     //reading magnetic field Hz
	     readFortranBinaryArray(f,buf);

	     //reading array of current, Jx
	     readFortranBinaryArray(f,buf);
	     //reading array of current, Jy
	     readFortranBinaryArray(f,buf);
	     //reading array of current, Jz
	     readFortranBinaryArray(f,buf);

	     //reading magnetic field at halfstep, Qx
	     readFortranBinaryArray(f,buf);

	     //reading magnetic field at halfstep, Qy
	     readFortranBinaryArray(f,buf);

	     //reading magnetic field at halfstep, Qz
	     readFortranBinaryArray(f,buf);
	     int err;
//--------------------------------------------
	     err = ferror(f);

	     if(err != 0) return NULL;

	     return f;
}

void debugPrintParticleCharacteristicArray(double *p_ch,int np,int nt,std::string name,int sort)
{
//	   char fname[200];
//	   FILE *f;

#ifndef PRINT_PARTICLE_INITIALS
	   return;

#else
	   sprintf(fname,"particle_init_%s_%05d_sort%02d.dat",name.c_str(),nt,sort);

	   if((f = fopen(fname,"wt")) == NULL) return;

	   for (int i = 0;i < np;i++)
	   {
		   fprintf(f,"%10d %10d %25.16e \n",i,i+1,p_ch[i]);
	   }
	   fclose(f);
#endif
}

int readBinaryParticleArraysOneSort(
		  FILE *f,
		  double **dbg_x,
		  double **dbg_y,
		  double **dbg_z,
		  double **dbg_px,
		  double **dbg_py,
		  double **dbg_pz,
		  double *qq_m,
		  double *mm,
		  int nt,
		  int sort
		  )
{
	//     char str[1000];
	     double /*x,y,z,px,py,pz,*/q_m,/* *buf,*/tp,m;
	     int t;
//		     Cell<Particle> c0 = (*AllCells)[0];
	     int total_particles;
	     int err;

	     if((err = ferror(f)) != 0)
	     {
	     	 return err ;
	     }
	     //Reading extra number placed by Fortran
	     fread(&t,sizeof(int),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }
	     //Reading number of particles of sort "sort"
	     fread(&tp,sizeof(double),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

	     //Reading charge for sort "sort"
	     total_particles = (int)tp;
	     fread(&q_m,sizeof(double),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }
         //Reading mass for sort "sort"
	     fread(&m,sizeof(double),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

	     // Reading extra number placed by Fortran
	     fread(&t,sizeof(int),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

       *dbg_x = (double *)malloc(sizeof(double)*total_particles);

       *dbg_y = (double *)malloc(sizeof(double)*total_particles);

       *dbg_z = (double *)malloc(sizeof(double)*total_particles);


       *dbg_px = (double *)malloc(sizeof(double)*total_particles);

       *dbg_py = (double *)malloc(sizeof(double)*total_particles);


       *dbg_pz = (double *)malloc(sizeof(double)*total_particles);

        //Reading X coordinates for particles of sort "sort"
	 	readFortranBinaryArray(f,*dbg_x);
	 	//Reading Y coordinates for particles of sort "sort"
	 	readFortranBinaryArray(f,*dbg_y);
	 	//Reading Z coordinates for particles of sort "sort"
	 	readFortranBinaryArray(f,*dbg_z);
	 	//Reading X impulses for particles of sort "sort"
	 	readFortranBinaryArray(f,*dbg_px);
	 	//Reading Y impulses for particles of sort "sort"
	 	readFortranBinaryArray(f,*dbg_py);
	 	//Reading Z impulses for particles of sort "sort"
	 	readFortranBinaryArray(f,*dbg_pz);
	 	debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);
        debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);
        debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);
        debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);
        debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);
        debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);

	 	//printf("particle 79943 %25.15e \n",(*dbg_x)[79943]);

	 	*qq_m = q_m;
	 	*mm   = m;

	 	if((err = ferror(f)) != 0)
      {
 	 	    return err ;
		}

	 	return total_particles;
}

int getParticlesOneSortFromFile(
		                          FILE *f,
                                particle_sorts sort,
                                int nt,
                                std::vector<Particle> & vp,
                                double *q_m,
                                double *m
                                )
{
//	     double x,y,z,px,py,pz;
	     double *dbg_x,*dbg_y,*dbg_z,*dbg_px,*dbg_py,*dbg_pz;
	     int total_particles;


//	     int err;

	     if((ferror(f)) != 0) return 1;

	     total_particles = readBinaryParticleArraysOneSort(f,&dbg_x,&dbg_y,&dbg_z,
	    		                                             &dbg_px,&dbg_py,&dbg_pz,q_m,m,nt,
	    		                                             sort);

	    // real_number_of_particle[(int)sort] = total_particles;

	     if((ferror(f)) != 0) return 1;
	     convertParticleArraysToSTLvector(dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz,*q_m,*m,
	    			  total_particles,sort,vp);
	     return 0;

}
std::vector<Particle> readBinaryParticlesOneSortSTL(FILE *f, particle_sorts sort,int nt)

	      	  {
	      		    double q_m,m;
//	      		    int err;
	      		    std::vector<Particle> vp;
	      		    getParticlesOneSortFromFile(f,sort,nt,vp,&q_m,&m);

//	      		    err = ferror(f);


	      			struct sysinfo info;
	                 sysinfo(&info);
	      			printf("before1  %d free %lu \n",nt,info.freeram/1024/1024);
//	      			err = ferror(f);
                    return vp;
	      //			printPICstatitstics(m,q_m,total_particles);
	      	  }

int readBinaryParticlesAllSorts(FILE *f,int nt,
			                          std::vector<Particle> & ion_vp,
                                      std::vector<Particle> & el_vp,
                                      std::vector<Particle> & beam_vp)
	  {
		  ion_vp = readBinaryParticlesOneSortSTL(f,ION,nt);
		  el_vp = readBinaryParticlesOneSortSTL(f,PLASMA_ELECTRON,nt);
    	  beam_vp = readBinaryParticlesOneSortSTL(f,BEAM_ELECTRON,nt);

    	  return 0;
	  }

//loading particle data and 12 mesh 3D arrays - everything
int LoadParticleData(int nt,    // number of timestep
		               std::vector<Particle> & ion_vp,     //ion particles
		               std::vector<Particle> & el_vp,      //electron particles
		               std::vector<Particle> & beam_vp,
		               int nx, //mesh size X,Y,Z
		               int ny,
		               int nz)
{

	 FILE *f;

	 std::string part_name = getMumuFileName(nt);

	 //reading 3D arrays: electric and magnetic fields and currents
	 if((f = readPreliminary3Darrays(part_name,nt,nx,ny,nz)) == NULL) return 1;

     //reading three sorts of particles
	 readBinaryParticlesAllSorts(f,nt,ion_vp,el_vp,beam_vp);

	 fclose(f);

	 return 0;
}



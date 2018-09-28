/*
 * particle_target.h
 *
 *  Created on: Aug 30, 2016
 *      Author: snytav
 */

#ifndef PARTICLE_TARGET_H_
#define PARTICLE_TARGET_H_

#include "point.h"
#include "particle.h"
#include "read_particles.h"

#define PARTICLES_RUN_ONE_TIMESTEP 10
#define DIMENSIONS 3


template<class Particle,int dims>
class ParticleTarget
{
public:
	  Point<int,DIMENSIONS,1> mesh;
	  double x;
	  Particle *particle_list_one_timestep;
	  int      *particle_number_one_timestep;
	  double    energy_lost;
//  vector<Particle> particle_history;
	  ParticleArrays *initial;

//	  template<class Particle,int dims>
	  int ncell(Point<int,DIMENSIONS,0> *f){return (f->y*mesh.z() + f->z());}

//	  template<class Particle,int dims>
	  int put(Point<int,DIMENSIONS,0> *from,Particle *p)
	  {
		  int nc = this->ncell(from);

          int num = particle_number_one_timestep[nc];
          Particle *p_array = &(particle_list_one_timestep[ncell(from)]);
          p_array[num] = *p;
          particle_number_one_timestep[ncell(from)] += 1;

          return 0;
 	  }

//	  template<class Particle,int dims>
	  int Diagnose(int nt)
	  {
		  char fname[100];
		  char fname_en[100];

		  sprintf(fname,"bound%e_%08d.dat",x,nt);
		  sprintf(fname_en,"bound%e_energy.dat",x);

		  FILE *f    = fopen(fname,"wt");


		  for(int i = 0;i < mesh.y;i++)
		  {

			  if(i == 10)
			  {
//				  int qq = 0;
			  }
			  for(int k = 0;k < mesh.z();k++)
			  {
				  Point<int,DIMENSIONS,0> from;
				  from.y = i;
				  from.setZ(k);

				  int num = particle_number_one_timestep[ncell(&from)];
				  Particle *p_array = &(particle_list_one_timestep[ncell(&from)]),p,p_init;

				  for(int n = 0;n < num;n++)
				  {
     				  p = p_array[n];

     				  double e, e0;
     				  p.getValue(&e,KINETIC_ENERGY);

     				  int fn = p.fortran_number;
     				  p_init.X.x = initial->dbg_x[fn];
     				  p_init.X.y = initial->dbg_y[fn];
     				  p_init.X.setZ(initial->dbg_z[fn]);
     				  p_init.pu = initial->dbg_px[fn];
     				  p_init.pv = initial->dbg_py[fn];
     				  p_init.pw = initial->dbg_pz[fn];

     				  p_init.getValue(&e0,KINETIC_ENERGY);

				      if(f != NULL)
				      {
				    	 fprintf(f,"%10d %10d %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e init  %15.5e %15.5e %15.5e %15.5e \n",nt,p.fortran_number,
				    			 p.X.x,p.X.y,p.X.z(),p.pu,p.pv,p.pw,
				    			 p_init.pu,p_init.pv,p_init.pw,
				    			 e/e0
				    			 );
				      }

				      energy_lost += e;
				  }
			  }
		  }

		  FILE *f_en;
		  if((f_en = fopen(fname_en,"at")) != NULL)
		  {
			  fprintf(f_en,"%10d %25.15e\n",nt,energy_lost);

			  fclose(f_en);
		  }

		  memset(particle_number_one_timestep,0,sizeof(int)*mesh.y*mesh.z());
          fclose(f);

		  return 0;
	  }

//	  template<class Particle,int dims>
	  ParticleTarget(){}

//	  template<class Particle,int dims>
	  int Init(ParticleArrays *pg,double bnd_x,Point<int,DIMENSIONS,1> mesh1)
	  {

		  energy_lost = 0.0;
		  initial  = pg;
		  x = bnd_x;

		  mesh.x = mesh1.x;
		  mesh.y = mesh1.y;
		  mesh.setZ(mesh1.z());

		  particle_list_one_timestep = (Particle *)malloc(sizeof(Particle)*mesh.y*mesh.z()*PARTICLES_RUN_ONE_TIMESTEP);
		  particle_number_one_timestep = (int *)malloc(sizeof(int)*mesh.y*mesh.z());


		  memset(particle_number_one_timestep,0,sizeof(int)*mesh.y*mesh.z());

		  return 0;
	  }
//	  template<class Particle,int dims>
	  ~ParticleTarget(){}
};




#endif /* PARTICLE_TARGET_H_ */

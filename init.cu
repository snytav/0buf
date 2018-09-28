//#include "load_data.h"


int InitializeGPU()
{
	int err = getLastError();
	 if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }
    InitGPUParticles();
    err = getLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }
    InitGPUFields(&d_Ex,&d_Ey,&d_Ez,
    	          &d_Hx,&d_Hy,&d_Hz,
    		      &d_Jx,&d_Jy,&d_Jz,
    		      &d_npJx,&d_npJy,&d_npJz,
                  &d_Qx,&d_Qy,&d_Qz,
                  Ex,Ey,Ez,
				  Hx,Hy,Hz,
				  Jx,Jy,Jz,
				  npJx,npJy,npJz,
				  Qx,Qy,Qz,
				  Nx,Ny,Nz
            );
    err = getLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }

    setPrintfLimit();

    err = cudaSetDevice(0);

    printf("InitializeGPU error %d \n",err);

    return 0;
}

int initMeshArrays()
{
	   initControlPointFile();

	   Alloc();

	   Cell c000;

	   InitCells();
	   c000 = (*AllCells)[0];

	   InitFields();
	   c000 = (*AllCells)[0];
	   InitCurrents();

	   return 0;
}


void LoadTestData(int nt,
		            int part_nt,
		            std::vector<Particle> & ion_vp,
		            std::vector<Particle> & el_vp,
		            std::vector<Particle> & beam_vp)
{

   if(nt > 1)
   	 {
   		 ClearAllParticles();
   	 }


   LoadParticleData(nt,ion_vp,el_vp,beam_vp,Nx,Ny,Nz);

   magf = 1;
}

void AssignArraysToCells()
{
   for(int n = 0;n < (*AllCells).size();n++)
   {

       Cell c = (*AllCells)[n];

	     c.readFieldsFromArrays(Ex,Ey,Ez,Hx,Hy,Hz);
   }
}

int compare(Particle p,Particle p1)
{
	int tx = comd(p.x,p1.x);
	int ty = comd(p.y,p1.y);
	int tz = comd(p.z,p1.z);
	int tpx = comd(p.pu,p1.pu);
	double dpx = fabs(p.pu - p1.pu);
	int tpy = comd(p.pv,p1.pv);
	int tpz = comd(p.pw,p1.pw);

	return (tx && ty && tz && tpx && tpy && tpz);
}

double compareParticleList(std::vector<Particle> v,std::vector<Particle> v1)
{
	double t = 0.0,s = v.size(),s1 = v1.size();

	if(v.size() != v1.size()) return 0.0;

	for (int i = 0;i < v.size();i++)
	{

			t += compare(v[i],v1[i]);
	}

	return (t/v.size());
}



virtual void InitializeCPU()
{
   std::vector<Particle> ion_vp,el_vp,beam_vp;
   std::vector<Particle> ion_vp1,el_vp1,beam_vp1;

   initMeshArrays();

   int flag_from_file = 0;

   if(flag_from_file == 1)
   {
      LoadTestData(START_STEP_NUMBER,START_STEP_NUMBER, ion_vp,el_vp,beam_vp);

   }
   else
   {
	   getUniformMaxwellianParticles(ion_vp1,el_vp1,beam_vp1);

   }


   addAllParticleListsToCells(ion_vp1,el_vp1,beam_vp1);

   AssignArraysToCells();


}

void Initialize()
{
	int err = getLastError();
	InitializeCPU();
	copyCellsWithParticlesToGPU();
	err = getLastError();
	 if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }
	InitializeGPU();

	err = getLastError();
	 if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }
}





void InitGPUParticles()
 //   :InitParticles(fname,vp)
{
	int size;
	GPUCell *d_c,*h_ctrl;
	GPUCell *n;
//	GPUCell<Particle> *h_c;//*h_copy,
//	double t;
	dim3 dimGrid(Nx+2,Ny+2,Nz+2),dimBlockOne(1,1,1);
	int err = getLastError();

	if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }
	err = getLastError();
	           if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
	 readControlFile(START_STEP_NUMBER);
	 err = getLastError();
	            if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}

	 err = getLastError();
     if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }

	size = (*AllCells).size();
	err = getLastError();
	           if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
	 size_t m_free,m_total;

	h_ctrl = new GPUCell;
	n = new GPUCell;

	err = getLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }
    err = getLastError();
               if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}

    h_CellArray = (GPUCell **)malloc(size*sizeof(Cell*));
    err = cudaMalloc(&d_CellArray,size*sizeof(Cell *));

    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }

//    h_controlParticleNumberArray = (int*)malloc(size*sizeof(int));
    err = getLastError();
               if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
    printf("%s : size = %d\n", __FILE__, size);
    for(int i = 0;i < size;i++)
    {
        //printf("GPU cell %d begins******************************************************\n",i);
    	GPUCell c;
    	c = (*AllCells)[i];

    	h_controlParticleNumberArray[i] = c.number_of_particles;
    	/////////////////////////////////////////
    	*n = c;
    	err = getLastError();
    	           if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
#ifdef ATTRIBUTES_CHECK
    	c.SetControlSystem(jmp,d_ctrlParticles);
#endif


        d_c = c.copyCellToDevice();
        err = getLastError();
                   if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
        err = getLastError();
        if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }

        cudaMemGetInfo(&m_free,&m_total);
//        double mtot;
//        mfree = m_free;
        err = getLastError();
                   if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
#ifdef COPY_CELL_PRINTS
        double mfree,mtot;
        mtot  = m_total;
        mfree = m_free;
        printf("cell %10d Device cell array allocated error %d %s memory: free %10.2f total %10.2f\n",i,err,getErrorString(err),
        		                                                mfree/1024/1024/1024,mtot/1024/1024/1024);
        puts("");

	  dbgPrintGPUParticleAttribute(d_c,50,1," CO2DEV " );
	  puts("COPY----------------------------------");
#endif


#ifdef PARTICLE_PRINTS

        if(t < 1.0)
        {
        	t = c.compareToCell(*h_copy);
        }
#endif
        ////////////////////////////////////////.
        err = getLastError();
                   if(err != cudaSuccess) { printf("%s:%d - error %d %s cell %d \n",__FILE__,__LINE__,err,getErrorString(err),i); exit(0);}

        h_CellArray[i] = d_c;
        err = MemoryCopy(h_ctrl,d_c,sizeof(Cell),DEVICE_TO_HOST);

      //  err = getLastError();
           if(err != cudaSuccess)
           {
        	   printf("%s:%d - error %d %s cell %d\n",__FILE__,__LINE__,err,getErrorString(err),i);
        	   exit(0);
           }
#ifdef InitGPUParticles_PRINTS
	    dbgPrintGPUParticleAttribute(d_c,50,1," CPY " );

       cudaPrintfInit();

        testKernel<<<1,1>>>(h_ctrl->d_ctrlParticles,h_ctrl->jmp);
        cudaPrintfDisplay(stdout, true);
        cudaPrintfEnd();

        printf("i %d l %d k n %d %d %e src %e num %d\n",h_ctrl->i,h_ctrl->l,h_ctrl->k,i,
        		c.ParticleArrayRead(0,7),c.number_of_particles
        		);
	printf("GPU cell %d ended ******************************************************\n",i);
#endif
	err = getLastError();
	           if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); exit(0);}
    }

    err = getLastError();
       if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }

    //int err;
    err = MemoryCopy(d_CellArray,h_CellArray,size*sizeof(Cell *),HOST_TO_DEVICE);
    if(err != cudaSuccess)
        {
         	printf("bGPU_WriteControlSystem err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = getLastError();
    if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }


#ifdef ATTRIBUTES_CHECK
    GPU_WriteControlSystem<<<dimGrid, dimBlockOne,16000>>>(d_CellArray);
#endif
	size = 0;

	err = getLastError();
	     if(err != cudaSuccess) { printf("%s:%d - error %d %s\n",__FILE__,__LINE__,err,getErrorString(err)); }

}


virtual void Alloc()
	  {

		  AllCells = new std::vector<GPUCell>;

	     Ex  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Ey  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Ez  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Hx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Hy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Hz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Jx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Jy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Jz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Rho = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];

	     npJx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     npJy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     npJz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];

	     npEx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     npEy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     npEz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];

	     Qx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Qy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     Qz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];

	#ifdef DEBUG_PLASMA

	     dbgEx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgEy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgEz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgEx0  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgEy0  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgEz0  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];

	     dbgHx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgHy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgHz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgJx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgJy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbgJz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];

	     dbg_Qx  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbg_Qy  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	     dbg_Qz  = new double[(Nx + 2)*(Ny + 2)*(Nz + 2)];
	#endif
	  }

	  virtual void InitFields()
	  {
	     for(int i = 0;i < (Nx+2)*(Ny+2)*(Nz+2);i++)
	     {
	         Ex[i] = 0.0;
	         Ey[i] = 0.0;
	         Ez[i] = 0.0;
	         Hx[i] = 0.0;
	         Hy[i] = 0.0;
	         Hz[i] = 0.0;

	         dbgEx[i] = 0.0;
	         dbgEy[i] = 0.0;
	         dbgEz[i] = 0.0;
	         dbgHx[i] = 0.0;
	         dbgHy[i] = 0.0;
	         dbgHz[i] = 0.0;
	     }
	  }

	  virtual void InitCells()
	  {
	     for(int i = 0;i < Nx+2;i++)
	     {
	         for(int l = 0;l < Ny+2;l++)
		 {
		     for(int k = 0;k < Nz+2;k++)
		     {
	                 GPUCell* c = new GPUCell(i,l,k,Lx,Ly,Lz,Nx,Ny,Nz,tau);
	                 c->Init();
			         (*AllCells).push_back(*c);
#ifdef INIT_CELLS_DEBUG_PRINT
	                printf("%5d %5d %5d size %d \n",i,l,k,(*AllCells).size());
#endif
		     }

		 }

	     }
	  }

	  virtual void InitCurrents()
	  {
	     for(int i = 0;i < (Nx+2)*(Ny+2)*(Nz+2);i++)
	     {
	         Jx[i]  = 0.0;
	         Jy[i]  = 0.0;
	         Jz[i]  = 0.0;
	         Rho[i] = 0.0;

	         dbgJx[i]  = 0.0;
	         dbgJy[i]  = 0.0;
	         dbgJz[i]  = 0.0;

	     }
	  }

	  void InitCurrents(string fnjx,string fnjy,string fnjz,
	                    string dbg_fnjx,string dbg_fnjy,string dbg_fnjz,
	                    string np_fnjx,string np_fnjy,string np_fnjz,
			            int dbg)
	  {

	     read3Darray(np_fnjx, npJx);
	     read3Darray(np_fnjy, npJy);
	     read3Darray(np_fnjz, npJz);

	     if(dbg == 0)
	     {
	        read3Darray(fnjx, Jx);
	        read3Darray(fnjy, Jy);
	        read3Darray(fnjz, Jz);
	     }
	#ifdef DEBUG_PLASMA
	     read3Darray(dbg_fnjx, dbgJx);
	     read3Darray(dbg_fnjy, dbgJy);
	     read3Darray(dbg_fnjz, dbgJz);

	#endif
	   }

	  void InitFields(char *fnex,char *fney,char *fnez,
			          char *fnhx,char *fnhy,char *fnhz,
			          char *dbg_fnex,char *dbg_fney,char *dbg_fnez,
			          char *dbg_0fnex,char *dbg_0fney,char *dbg_0fnez,
			          char *np_ex,char *np_ey,char *np_ez,
			          char *dbg_fnhx,char *dbg_fnhy,char *dbg_fnhz)
	  {
	     InitFields();

	      //double t271 = Hy[27];

	     read3Darray(fnex, Ex);
	     read3Darray(fney, Ey);
	     read3Darray(fnez, Ez);
	     read3Darray(fnhx, Hx);
	     read3Darray(fnhy, Hy);
	     read3Darray(fnhz, Hz);

	#ifdef DEBUG_PLASMA
	     read3Darray(dbg_fnex, dbgEx);
	     read3Darray(dbg_fney, dbgEy);
	     read3Darray(dbg_fnez, dbgEz);

	     read3Darray(dbg_0fnex, dbgEx0);
	     read3Darray(dbg_0fney, dbgEy0);
	     read3Darray(dbg_0fnez, dbgEz0);

	     read3Darray(dbg_fnhx, dbgHx);
	     read3Darray(dbg_fnhy, dbgHy);
	     read3Darray(dbg_fnhz, dbgHz);

	     read3DarrayLog(np_ex, npEx,50,8);
	     read3DarrayLog(np_ey, npEy,50,8);
	     read3DarrayLog(np_ez, npEz,50,8);
	#endif

	  //   double t27 = Hy[27];

	  }


	  virtual void InitParticles(thrust::host_vector<Particle> & vp)
	  {
	     InitIonParticles(n_per_cell,ion_q_m,vp);
	  }

	  virtual void InitParticles(char *fname,thrust::host_vector<Particle>& vp)
	  {
	     FILE *f;
	     char str[1000];
	     double x,y,z,px,py,pz,q_m,m;
	     int n = 0;

	     if((f = fopen(fname,"rt")) == NULL) return;

	     while(fgets(str,1000,f) != NULL)
	     {
	          x   = atof(str);
	          y   = atof(str + 25);
	          z   = atof(str + 50);
	          px  = atof(str + 75);
	          py  = atof(str + 100);
	          pz  = atof(str + 125);
	          m   = fabs(atof(str + 150));
	          q_m = atof(str + 175);
	#undef GPU_PARTICLE
		  Particle *p = new Particle(x,y,z,px,py,pz,m,q_m);
//		      if(n == 829)
//		      {
//		    	  int qq = 0;
//		    	  qq = 1;
//		      }
		  p->fortran_number = ++n;
		  vp.push_back(*p);
	#define GPU_PARTICLE

	     }


         dbg_x = (double *)malloc(sizeof(double)*vp.size());
         dbg_y = (double *)malloc(sizeof(double)*vp.size());
         dbg_z = (double *)malloc(sizeof(double)*vp.size());
         dbg_px = (double *)malloc(sizeof(double)*vp.size());
         dbg_py = (double *)malloc(sizeof(double)*vp.size());
         dbg_pz = (double *)malloc(sizeof(double)*vp.size());

         total_particles = vp.size();

	     magf = 1;
	  }



      void printPICstatitstics(double m,double q_m, int total_particles)
      {
    	  int pn_min,pn_ave,pn_max,pn_sum;//,err;

              pn_min = 1000000000;
              pn_max = 0;
              pn_ave = 0;
 		     for(int n = 0;n < (*AllCells).size();n++)
 		     {
 		    	 Cell & c = (*AllCells)[n];

 		    	 pn_ave += c.number_of_particles;
 		    	 if(pn_min > c.number_of_particles) pn_min = c.number_of_particles;
 		    	 if(pn_max < c.number_of_particles) pn_max = c.number_of_particles;

 		     }

 		     pn_sum = pn_ave;
 		     pn_ave /= (*AllCells).size();

 		     printf("SORT m %15.5e q_m %15.5e %10d (sum %10d) particles in %8d cells: MIN %10d MAX %10d average %10d \n",
 		    		 m,            q_m,       total_particles,pn_sum,
 		    		 (*AllCells).size(),
 		    		 pn_min,pn_max,pn_ave);


      }


      int addParticleListToCells(std::vector<Particle>& vp)
      {
    	  Cell c0 = (*AllCells)[0];
    	  int n;

    	  for(int i = 0; i < vp.size();i++)
    	  {
    	      Particle p = vp[i]; // = new Particle(x,y,z,px,py,pz,m,q_m);

    	  	  double3 d;
    	      d.x = p.x;
    	      d.y = p.y;
    	      d.z = p.z;

    	      n = c0.getPointCell(d);

    	      Cell & c = (*AllCells)[n];


    	      if(c.Insert(p) == true)
    	      {
#ifdef PARTICLE_PRINTS1000
    	  		             if((i+1)%1000 == 0 )
    	  		             {
    	  		        	     printf("particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d)\n",
    	  		        	    		 i+1,
    	  				    		x,y,z,c.number_of_particles,c.i,c.l,c.k);
    	  		             }
#endif
    	  			      }
   		      }// END total_particles LOOP

    	  return 0;

      }









	  int addAllParticleListsToCells(std::vector<Particle> & ion_vp,
			                         std::vector<Particle> & el_vp,
			                         std::vector<Particle> & beam_vp)
	  {
			 addParticleListToCells(ion_vp);
			 addParticleListToCells(el_vp);
			 addParticleListToCells(beam_vp);

			 return 0;
	  }



//
	  int readParticles(FILE *f,int nt)
	  {
		 std::vector<Particle> ion_vp,el_vp,beam_vp;

		 readBinaryParticlesAllSorts(f,nt,ion_vp,el_vp,beam_vp);


		 addAllParticleListsToCells(ion_vp,el_vp,beam_vp);

		 return 0;
	  }



	  virtual void InitBinaryParticles(int nt)
	  {
	     FILE *f;
	     std::string part_name  = getMumuFileName(nt);

		 if((f = readPreliminary3Darrays(part_name,nt,Nx,Ny,Nz)) == NULL) return;

		 std::vector<Particle> ion_vp,el_vp,beam_vp;

		 readBinaryParticlesAllSorts(f,nt,ion_vp,el_vp,beam_vp);


		 addAllParticleListsToCells(ion_vp,el_vp,beam_vp);

	     fclose(f);

	     magf = 1;
	  }



	  virtual void InitElectronParticles(){}
	  virtual void InitIonParticles(int n_per_cell1,double q_m,thrust::host_vector<Particle> &vecp)
	  {
	     int total_ions = Nx*Ny*Nz*n_per_cell;
	     Particle *p;
	     //double ami = ni /((double)n_per_cell);
	     double x,y,z;

	     for(int j = 0;j < total_ions;j++)
	     {
		z = 0.0;
		y = 0.0;
		x = 0.0;

		p = new Particle(x,y,z,0.0,0.0,0.0,ni,q_m);

	#ifdef DEBUG_PLASMA
//		printf("particle %d \n",j);
	#endif

		vecp.push_back(*p);
	     }
	  }

	  virtual void InitBeamParticles(int n_per_cell1){}
	  void Distribute(thrust::host_vector<Particle> &vecp)
	  {
	     Cell c0 = (*AllCells)[0],c111;
	     int    n;//,i;
	     int  vec_size = vecp.size();

	     for(int j = 0;j < vecp.size();j++)
	     {
		 Particle p = vecp[j];
		 double3 d;
		 d.x = p.x;
		 d.y = p.y;
		 d.z = p.z;

		 n = c0.getPointCell(d);

		 Cell & c = (*AllCells)[n];


		 if(c.Insert(p) == true)
		 {
	#ifdef PARTICLE_PRINTS1000
         if((vec_size-vecp.size())%1000 == 0 )	printf("particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d)\n",vec_size-vecp.size(),
		    		p.x,p.y,p.z,c.number_of_particles,c.i,c.l,c.k);
         if((vec_size-vecp.size()) == 10000) exit(0);
	#endif
		    vecp.erase(vecp.begin()+j);
		    j--;
		 }
	     }
	     int pn_min = 1000000,pn_max = 0,pn_ave = 0;
	     for(int n = 0;n < (*AllCells).size();n++)
	     {
	    	 Cell & c = (*AllCells)[n];

	    	 pn_ave += c.number_of_particles;
	    	 if(pn_min > c.number_of_particles) pn_min = c.number_of_particles;
	    	 if(pn_max < c.number_of_particles) pn_max = c.number_of_particles;

	     }
	     pn_ave /= (*AllCells).size();

	     printf("%10d particles in %8d cells: MIN %5d MAX %5d average %5d \n",vec_size,(*AllCells).size(),
	    		                                                              pn_min,pn_max,pn_ave);
	  }


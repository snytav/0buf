







__global__
void GPU_SetAllCurrentsToZero(GPUCell  **cells)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	//int i,l,k;
	Cell *c,*c0 = cells[0],nc;
	//double t;
//	__shared__ extern CellDouble fd[9];
	//double *src;//,*dst;

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];

	nc = *c;

	nc.SetAllCurrentsToZero(threadIdx);
}





__global__
void GPU_SetFieldsToCells(GPUCell  **cells,
        double *Ex,double *Ey,double *Ez,
        double *Hx,double *Hy,double *Hz
		)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	//int i,l,k;
	Cell  *c,*c0 = cells[0];
	//double t;

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];


	c->readFieldsFromArrays(Ex,Ey,Ez,Hx,Hy,Hz,threadIdx);
}

__host__ __device__
double CheckArraySize(double* a, double* dbg_a,int size)
	{
	//    Cell<Particle> c = (*AllCells)[0];
	    int wrong = 0;
#ifdef CHECK_ARRAY_SIZE_DEBUG_PRINTS
	    printf("begin array checking1=============================\n");
#endif
	    for(int n = 0;n < size;n++)
	    {
	        //double t  = a[n];
	//	double dt = dbg_a[n];

	        if(fabs(a[n] - dbg_a[n]) > SIZE_TOLERANCE)
		{

		   //int3 i = c.getCellTripletNumber(n);
#ifdef CHECK_ARRAY_SIZE_DEBUG_PRINTS
		   printf("n %5d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
				   n,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
#endif
		}
	    }
#ifdef CHECK_ARRAY_SIZE_DEBUG_PRINTS
	    printf("  end array checking=============================\n");
#endif

	    return (1.0-((double)wrong/(size)));
	}



__global__ void GPU_WriteAllCurrents(GPUCell **cells,int n0,
		      double *jx,double *jy,double *jz,double *rho)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell  *c,*c0 = cells[0];
//	__shared__ extern CellDouble fd[9];

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];

	int i1,l1,k1;
	        	 i1 = threadIdx.x;
	        	 l1 = threadIdx.y;
	        	 k1 = threadIdx.z;
    	         int n = c->getFortranCellNumber(c->i+i1-1,c->l+l1-1,c->k+k1-1);

    	         if (n < 0 ) n = -n;
        		 double t,t_x,t_y;
		         t_x = c->Jx->M[i1][l1][k1];
		         int3 i3 = c->getCellTripletNumber(n);


		         cuda_atomicAdd(&(jx[n]),t_x);
		         t_y= c->Jy->M[i1][l1][k1];
		         cuda_atomicAdd(&(jy[n]),t_y);
		         t = c->Jz->M[i1][l1][k1];
		         cuda_atomicAdd(&(jz[n]),t);

}

__global__ void GPU_WriteControlSystem(Cell **cells)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
//	int i,l,k;
	Cell  *c,*c0 = cells[0],nc;
	//double t;
//	__shared__ extern CellDouble fd[9];
	//double *src; //,*dst;
//	int pqr2;
	//CurrentTensor t1,t2;

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];
//	c = cells[ n ];

	 nc = *c;

	 nc.SetControlSystemToParticles();

}

//TODO : 1. 3 separate kernels :
//            A. form 3x3x3 array with number how many to fly and departure list with start positions in 3x3x3 array
//            B. func to get 3x3x3 indexes from a pair of cell numbers, to and from function
//            C. 2nd kernel to write arrival 3x3x3 matrices
///           D. 3rd kernel to form arrival positions in the particle list
//            E. 4th to write arriving particles




__global__ void GPU_MakeDepartureLists(GPUCell  **cells,int nt,int *d_stage)
{
	    unsigned int nx = blockIdx.x;
		unsigned int ny = blockIdx.y;
		unsigned int nz = blockIdx.z;
		int ix,iy,iz;//,n;

		Particle p;
		Cell  *c,*c0 = cells[0],nc;//,*new_c;
		c = cells[c0->getGlobalCellNumber(nx,ny,nz)];






		c->departureListLength = 0;
		for(ix = 0;ix < 3;ix++)
		{
			for(iy = 0;iy < 3;iy++)
			{
				for(iz = 0;iz < 3;iz++)
				{
					c->departure[ix][iy][iz]      = 0;
				//	c->departureIndex[ix][iy][iz] = 0;
				}
			}
		}
		c->departureListLength  = 0;
		for(int num = 0;num < c->number_of_particles; num++)
			{
			p = c->readParticleFromSurfaceDevice(num);

			c->flyDirection(&p,&ix,&iy,&iz);
			c->writeParticleToSurface(num,&p);
				if(!c->isPointInCell(p.GetX()))   //check Paricle = operator !!!!!!!!!!!!!!!!!!!!!!!!!!!
				{

//					c->flyDirection(&p,&ix,&iy,&iz);
//					printf("fly %d:(%d,%d,%d) %d \n",p.direction,ix,iy,iz,ix*9 +iy*3 +iz);


					c->writeParticleToSurface(num,&p);


//					if(p.direction == 0) printf("Blin-hren'\n");
//					f(p.direction != (ix | (iy << 2) |(iz << 4))) printf("Blin-hren'^2\n");
                    if(c->departureListLength == PARTICLES_FLYING_ONE_DIRECTION)
                    {
                    	d_stage[0] = TOO_MANY_PARTICLES;
                    	d_stage[1] = c->i;
                    	d_stage[2] = c->l;
                    	d_stage[3] = c->k;
                    	d_stage[1] = ix;
                    	d_stage[2] = iy;
                    	d_stage[3] = iz;
                    	return;
                    }
					c->departureListLength++;
					int num1 = c->departure[ix][iy][iz];

					c->departureList[ix][iy][iz][num1] = p;

					c->departure[ix][iy][iz] += 1;
				}
				else
				{
//					printf("home %d:(%d,%d,%d) %d \n",p.direction,ix,iy,iz,ix*9 +iy*3 +iz);
				}
			}
}


__global__ void GPU_MakeDepartureLists_rm(GPUCell  **cells,int nt,int *d_stage)
{
	    unsigned int nx = blockIdx.x;
		unsigned int ny = blockIdx.y;
		unsigned int nz = blockIdx.z;
		int ix,iy,iz;//,n;

		Particle p;
		Cell  *c,*c0 = cells[0],nc;//,*new_c;
		c = cells[c0->getGlobalCellNumber(nx,ny,nz)];



	     	for(int num = 0;num < c->number_of_particles; num++)
			{
	     		p = c->readParticleFromSurfaceDevice(num);

				if(p.direction != 13)   //check Paricle = operator !!!!!!!!!!!!!!!!!!!!!!!!!!!
				{
					c->removeParticleFromSurfaceDevice(num,&p,&(c->number_of_particles));
										num--;
				}
			}
}

__global__ void GPU_RemoveDepartureParticles(GPUCell  **cells,int nt,int *d_stage)
{
	    unsigned int nx = blockIdx.x;
		unsigned int ny = blockIdx.y;
		unsigned int nz = blockIdx.z;
		int ix,iy,iz;//,n;

		Particle p;
		Cell  *c,*c0 = cells[0],nc;//,*new_c;
		c = cells[c0->getGlobalCellNumber(nx,ny,nz)];








		for(int num = 0;num < c->number_of_particles; num++)
			{
			p = c->readParticleFromSurfaceDevice(num);


				if(!c->isPointInCell(p.GetX()))   //check Paricle = operator !!!!!!!!!!!!!!!!!!!!!!!!!!!
				{
					c->removeParticleFromSurfaceDevice(num,&p,&(c->number_of_particles));
//					c->flyDirection(&p,&ix,&iy,&iz);
//
//					//if(p.direction == 0) printf("Blin-hren'\n");
//					//f(p.direction != (ix | (iy << 2) |(iz << 4))) printf("Blin-hren'^2\n");
//                    if(c->departureListLength == PARTICLES_FLYING_ONE_DIRECTION)
//                    {
//                    	d_stage[0] = TOO_MANY_PARTICLES;
//                    	d_stage[1] = c->i;
//                    	d_stage[2] = c->l;
//                    	d_stage[3] = c->k;
//                    	d_stage[1] = ix;
//                    	d_stage[2] = iy;
//                    	d_stage[3] = iz;
//                    	return;
//                    }
//					c->departureListLength++;
//					int num1 = c->departure[ix][iy][iz];
//
//					c->departureList[ix][iy][iz][num1] = p;
//
//					c->departure[ix][iy][iz] += 1;
//
//					num--;
				}
			}
}

__device__ int d_comd(double a,double b)
{
	return (fabs(a - b) < TOLERANCE);
}

__device__ int d_compare(Particle p,Particle p1)
{
	int tx = d_comd(p.x,p1.x);
	int ty = d_comd(p.y,p1.y);
	int tz = d_comd(p.z,p1.z);
	int tpx = d_comd(p.pu,p1.pu);
	double dpx = fabs(p.pu - p1.pu);
	int tpy = d_comd(p.pv,p1.pv);
	int tpz = d_comd(p.pw,p1.pw);

	return (tx && ty && tz && tpx && tpy && tpz);
}


__global__ void GPU_ArrangeFlights(GPUCell  **cells,int nt, int *d_stage)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	int ix,iy,iz,snd_ix,snd_iy,snd_iz,num,n;
	Particle p,p1;

	Cell  *c,*c0 = cells[0],nc,*snd_c;


		c = cells[ n = c0->getGlobalCellNumber(nx,ny,nz)];

		for(ix = 0;ix < 3;ix++)
			for(iy = 0;iy < 3;iy++)
				for(iz = 0;iz < 3;iz++)
				{

					n = c0->getWrapCellNumber(nx+ix-1,ny+iy-1,nz+iz-1);

		            snd_c  = cells[ n ];

		            for(int j = 0;j < snd_c->number_of_particles;j++)
		            {
		            	p1 = snd_c->readParticleFromSurfaceDevice(j);
//					snd_ix = ix;
//					snd_iy = iy;
//					snd_iz = iz;
//
//					c->inverseDirection(&snd_ix,&snd_iy,&snd_iz);
//
//					num = snd_c->departure[snd_ix][snd_iy][snd_iz];
//
//
//
//					    for(int i = 0;i < num;i++)
//					    {
//						    p = snd_c->departureList[snd_ix][snd_iy][snd_iz][i];
//
//
//
							if(p1.direction != 13)
							{
								p1.direction = 13;
								c->Insert(p1);
							}
//
//						}


					}


				}
}



__global__ void GPU_CollectStrayParticles(Cell **cells,int nt
//		                         int n,
//		                         int i,
//		                         double mass,
//		                         double q_mass,
//		                         double *p_control,
//		                         int jmp
		                         )
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;

//	int busy;
	Particle p;
	int n;
//	int i,l,k;
	Cell  *c,*c0 = cells[0],nc,*new_c;
	//int first = 1;

	c = cells[ n = c0->getGlobalCellNumber(nx,ny,nz)];

	for(int i = 0;i < c->number_of_particles; i++)
	{
		p = c->readParticleFromSurfaceDevice(i);
#ifdef STRAY_DEBUG_PRINTS
		//if((p.fortran_number == 2498) && (p.sort == BEAM_ELECTRON))
	//	{
		    printf("STRAY-BASIC step %d cell %3d %d %d sort %d particle %d FORTRAN %5d X: %15.5e < %15.5e < %15.5e Y: %15.5e < %15.5e < %15.5e Z: %15.5e < %15.5e < %15.5e \n",
		    		nt,c->i,c->l,c->k,(int)p.sort,i,p.fortran_number,
		    		c->x0,p.x,c->x0+c->hx,
		    		c->y0,p.y,c->y0+c->hy,
		    		c->z0,p.z,c->z0+c->hz
		    		);
		//}
#endif
		if(!c->isPointInCell(p.GetX()))// || (p.fortran_number == 753) )//|| (p.fortran_number == 10572))
		{
#ifdef STRAY_DEBUG_PRINTS

   			    printf("STRAY-OUT step %3d cell %3d %d %d sort %d particle %d FORTRAN %5d X: %15.5e < %25.17e < %15.5e \n",
   			    		nt,c->i,c->l,c->k,(int)p.sort,i,p.fortran_number,c->x0,p.x,c->x0+c->hx);
#endif
            int new_n = c->getPointCell(p.GetX());
            new_c = cells[new_n];


            if(c->i == 99 && c->l == 0 && c->k == 3)
            {
         //      c->printCellParticles();

//               if(c->i >= c->Nx-1)
//               {
#ifdef STRAY_DEBUG_PRINTS
//       		if((p.fortran_number == 2498) && (p.sort == BEAM_ELECTRON))
//       		{
//       		   printf("c %3d (%d,%d,%d)->(%d,%d,%d) s %d p %d FN %5d X: %15.5e<%23.16e<%15.5e \n",n,c->i,c->l,c->k,new_c->i,new_c->l,new_c->k,(int)p.sort,i,p.fortran_number,c->x0,p.x,c->x0+c->hx);
//   			   printf("c %3d (%d,%d,%d)->(%d,%d,%d) s %d p %d FN %5d Y: %15.5e<%23.16e<%15.5e \n",n,c->i,c->l,c->k,new_c->i,new_c->l,new_c->k,(int)p.sort,i,p.fortran_number,c->y0,p.y,c->y0+c->hy);
//   			   printf("c %3d (%d,%d,%d)->(%d,%d,%d) s %d p %d FN %5d Z: %15.5e<%23.16e<%15.5e \n",n,c->i,c->l,c->k,new_c->i,new_c->l,new_c->k,(int)p.sort,i,p.fortran_number,c->z0,p.z,c->z0+c->hz);
//       		}
#endif
//               }
            }
//            if(first == 1)
//            {

  //          	do{
  //          	    	busy = atomicCAS(&(c->busyParticleArray),0,1);
            	    	// busy = ((c->busyParticleArra == 0 )? 1: c->busyParticleArra)

            	    	// busy = ((c->busyParticleArra == 1 )? 0: c->busyParticleArra)

   //         	  }while(busy == 1);

            while (atomicCAS(&(c->busyParticleArray),0,1)) {}
               c->removeParticleFromSurfaceDevice(i,&p,&(c->number_of_particles));
               //c->busyParticleArray = 0;
              atomicExch(&(c->busyParticleArray),0u);
               i--;
//               first = 0;
//            }
//            if(c->i == 99 && c->l == 0 && c->k == 3)
//            {
//               c->printCellParticles();
//            }
              //do{
               // 	busy = atomicCAS(&(new_c->busyParticleArray),0,1);
              //}while(busy == 1);

               while (atomicCAS(&(new_c->busyParticleArray),0,1)) {}

              new_c->Insert(p);
#ifdef STRAY_DEBUG_PRINTS

   			    printf("STRAY-INSERT step %d %3d %d %d sort %d particle %d FORTRAN %5d X: %15.5e < %25.17e < %15.5e \n",
   			    		nt,
   			    		new_c->i,new_c->l,new_c->k,(int)p.sort,i,p.fortran_number,new_c->x0,p.x,new_c->x0+new_c->hx);
#endif
              //new_c->busyParticleArray = 0;
              atomicExch(&(new_c->busyParticleArray),0u);
       		if((p.fortran_number == 2498) && (p.sort == BEAM_ELECTRON))
      		{

          //    new_c->printCellParticles();
      		}
            if(c->i == 99 && c->l == 0 && c->k == 3)
            {
            //   new_c->printCellParticles();
            }
		}
	}
//	c->printCellParticles("STRAY-FINAL",nt);

}

__device__ double checkCurrentComponentImpact(
		CurrentTensorComponent *t1,CurrentTensorComponent *t2,
		int i,int l, int k,int pqr2
		)
{
	double res = 0;
    if((t1->i11 == i && t1->i12 == l && t1->i13 == k && (fabs(t1->t[0]) > 1e-15))) res = t1->t[0];
    if((t1->i21 == i && t1->i22 == l && t1->i23 == k && (fabs(t1->t[1]) > 1e-15))) res = t1->t[1];
    if((t1->i31 == i && t1->i32 == l && t1->i33 == k && (fabs(t1->t[2]) > 1e-15))) res = t1->t[2];
    if((t1->i41 == i && t1->i42 == l && t1->i43 == k && (fabs(t1->t[3]) > 1e-15))) res = t1->t[3];


    if(pqr2 == 2)
    {
    	if((t2->i11 == i && t2->i12 == l && t2->i13 == k && (fabs(t2->t[0]) > 1e-15))) res = t2->t[0];
    	if((t2->i21 == i && t2->i22 == l && t2->i23 == k && (fabs(t2->t[1]) > 1e-15))) res = t2->t[1];
    	if((t2->i31 == i && t2->i32 == l && t2->i33 == k && (fabs(t2->t[2]) > 1e-15))) res = t2->t[2];
    	if((t2->i41 == i && t2->i42 == l && t2->i43 == k && (fabs(t2->t[3]) > 1e-15))) res = t2->t[3];
    }
    return res;
}


__device__ void writeCurrentComponent(CellDouble *J,
		CurrentTensorComponent *t1,CurrentTensorComponent *t2,int pqr2)
{
    cuda_atomicAdd(&(J->M[t1->i11][t1->i12][t1->i13]),t1->t[0]);
    cuda_atomicAdd(&(J->M[t1->i21][t1->i22][t1->i23]),t1->t[1]);
    cuda_atomicAdd(&(J->M[t1->i31][t1->i32][t1->i33]),t1->t[2]);
    cuda_atomicAdd(&(J->M[t1->i41][t1->i42][t1->i43]),t1->t[3]);

    if(pqr2 == 2)
    {
        cuda_atomicAdd(&(J->M[t2->i11][t2->i12][t2->i13]),t2->t[0]);
        cuda_atomicAdd(&(J->M[t2->i21][t2->i22][t2->i23]),t2->t[1]);
        cuda_atomicAdd(&(J->M[t2->i31][t2->i32][t2->i33]),t2->t[2]);
        cuda_atomicAdd(&(J->M[t2->i41][t2->i42][t2->i43]),t2->t[3]);
    }

}

__device__ void copyCellDouble(CellDouble *dst,CellDouble *src,unsigned int n,uint3 block)
{
	if(n < CellExtent*CellExtent*CellExtent)
	{
		double *d_dst,*d_src;//,t;

		d_dst = (double *)(dst->M);
		d_src = (double *)(src->M);

//		t = d_dst[n];
//
//		if(fabs(d_dst[n] - d_src[n]) > 1e-15)
//		{
//     		printf("block %5d %3d %3d thread %5d CCD t %15.5e dst %15.5e src %15.5e dst %p src %p d_dst[n] %p d_src[n] %p \n",
//     				block.x,block.y,block.z,threadIdx.x,t,d_dst[n],d_src[n],dst,src,&(d_dst[n]),&(d_src[n]));
//		}
		d_dst[n] = d_src[n];
	}
}

__device__ void addCellDouble(CellDouble *dst,CellDouble *src,unsigned int n,uint3 block,int dim)
{
	if(n < CellExtent*CellExtent*CellExtent)
	{
		double *d_dst,*d_src;//,t;

		d_dst = (double *)(dst->M);
		d_src = (double *)(src->M);

//		t = d_dst[n];
//
//		if(fabs(d_dst[n] - d_src[n]) > 1e-15)
//		{
//     		printf("block %5d %3d %3d thread %5d CCD t %15.5e dst %15.5e src %15.5e dst %p src %p d_dst[n] %p d_src[n] %p \n",
//     				block.x,block.y,block.z,threadIdx.x,t,d_dst[n],d_src[n],dst,src,&(d_dst[n]),&(d_src[n]));
//		}
		for (int i = 0;i < dim;i++)
		{
			CellDouble *t = &(src[i]);
			d_src         = (double *)(t->M);

		    d_dst[n] += d_src[n];
		}
	}
}

__device__ void setCellDoubleToZero(CellDouble *dst,unsigned int n)
{
	if(n < CellExtent*CellExtent*CellExtent)
	{
		double *d_dst;//*d_src;//,t;

		d_dst = (double *)(dst->M);
//		d_src = (double *)(src->M);

//		t = d_dst[n];
//
//		if(fabs(d_dst[n] - d_src[n]) > 1e-15)
//		{
//     		printf("block %5d %3d %3d thread %5d CCD t %15.5e dst %15.5e src %15.5e dst %p src %p d_dst[n] %p d_src[n] %p \n",
//     				block.x,block.y,block.z,threadIdx.x,t,d_dst[n],d_src[n],dst,src,&(d_dst[n]),&(d_src[n]));
//		}
		d_dst[n] = 0.0;
	}
}


__global__ void GPU_GetCellNumbers(Cell **cells,
		                         int *numbers)
{
		Cell  *c;//,nc;
		c = cells[blockIdx.x];

		numbers[blockIdx.x] = (*c).number_of_particles;
}

__device__ void assignSharedWithLocal(
		                         CellDouble **c_jx,
		                         CellDouble **c_jy,
		                         CellDouble **c_jz,
		                         CellDouble **c_ex,
		                         CellDouble **c_ey,
		                         CellDouble **c_ez,
		                         CellDouble **c_hx,
		                         CellDouble **c_hy,
		                         CellDouble **c_hz,
		                         CellDouble *fd)
{
	*c_ex = &(fd[0]);
	*c_ey = &(fd[1]);
	*c_ez = &(fd[2]);

	*c_hx = &(fd[3]);
	*c_hy = &(fd[4]);
	*c_hz = &(fd[5]);

	*c_jx = &(fd[6]);
	*c_jy = &(fd[7]);
	*c_jz = &(fd[8]);
}

__device__ void copyFieldsToSharedMemory(
		 CellDouble *c_jx,
		 CellDouble *c_jy,
		 CellDouble *c_jz,
		 CellDouble *c_ex,
		 CellDouble *c_ey,
		 CellDouble *c_ez,
		 CellDouble *c_hx,
		 CellDouble *c_hy,
		 CellDouble *c_hz,
		 Cell *c,
		 int index,
		 dim3 blockId,
		 int blockDimX
		)
{
	//int index  = threadIdx.x;


	while(index < CellExtent*CellExtent*CellExtent)
	{
//		if(index < 125) {

		copyCellDouble(c_ex,c->Ex,index,blockId);
		copyCellDouble(c_ey,c->Ey,index,blockId);
		copyCellDouble(c_ez,c->Ez,index,blockId);

		copyCellDouble(c_hx,c->Hx,index,blockId);
		copyCellDouble(c_hy,c->Hy,index,blockId);
		copyCellDouble(c_hz,c->Hz,index,blockId);

		copyCellDouble(c_jx,c->Jx,index,blockId);
		copyCellDouble(c_jy,c->Jy,index,blockId);
		copyCellDouble(c_jz,c->Jz,index,blockId);
		//}
		index += blockDimX;
	}

	__syncthreads();

}

__device__ void set_cell_double_arrays_to_zero(
		 CellDouble *m_c_jx,
		 CellDouble *m_c_jy,
		 CellDouble *m_c_jz,
		 int size,
		 int index,
		 int blockDimX
		)
{
	//int index  = threadIdx.x;
	for(int i = 0;i < size;i++)
	        {
	setCellDoubleToZero(&(m_c_jx[i]),index);
	setCellDoubleToZero(&(m_c_jy[i]),index);
	setCellDoubleToZero(&(m_c_jz[i]),index);
	        }
//	setCellDoubleToZero(&(m_c_jx[1]),index,blockId);
//	setCellDoubleToZero(&(m_c_jx[0]),index,blockId);
//	setCellDoubleToZero(&(m_c_jx[0]),index,blockId);

	while(index < CellExtent*CellExtent*CellExtent)
	{
//		if(index < 125) {
        for(int i = 0;i < size;i++)
        {
//		    setCellDoubleToZero(&(m_c_jx[i]),index,blockId);
//		    setCellDoubleToZero(&(m_c_jy[i]),index,blockId);
//		    setCellDoubleToZero(&(m_c_jz[i]),index,blockId);
        }
		index += blockDimX;
	}

	__syncthreads();

}

__device__ void set_cell_double_array_to_zero(CurrentTensorComponent *ca,int length)
{
     for(int i = 0; i<= 100;i++)
     {
    	ca[i].t[0] = 0.0;


     }
}

__device__ void MoveParticlesInCell(
									 CellDouble *c_ex,
									 CellDouble *c_ey,
									 CellDouble *c_ez,
									 CellDouble *c_hx,
									 CellDouble *c_hy,
									 CellDouble *c_hz,
									 Cell  *c,
		                             int index,
		                             int blockDimX//,
//		                             double mass,
//		                             double q_mass
		                             )
{
//	CurrentTensor t1,t2;
    int pqr2;
//	Particle p;
    CellTotalField cf;

    while(index < c->number_of_particles)
    {
    	cf.Ex = c->Ex;
    	cf.Ey = c->Ey;
    	cf.Ez = c->Ez;
    	cf.Hx = c->Hx;
    	cf.Hy = c->Hy;
    	cf.Hz = c->Hz;

        c->MoveSingleParticle(index,cf);


        index += blockDimX;
    }



    __syncthreads();
}

__device__ void assign_cell_double(CellDouble *a,CellDouble *b)
{
	int i,l,k;

	for(i = 0;i < CellExtent;i++)
	{
		for(l = 0;l < CellExtent;l++)
		{
			for(k = 0;k < CellExtent;k++)
			{
				a->M[i][l][k] = b->M[i][l][k];
			}
		}
	}
}

__device__ void add_cell_double(CellDouble *a,CellDouble *b)
{
	int i,l,k;

	for(i = 0;i < CellExtent;i++)
	{
		for(l = 0;l < CellExtent;l++)
		{
			for(k = 0;k < CellExtent;k++)
			{
				a->M[i][l][k] += b->M[i][l][k];
			}
		}
	}
}

__device__ void set_cell_double_array_to_zero(CellDouble *arr,int size)
{
	int i,l,k;
	CellDouble *a;

	for(int num = 0;num < size;num++)
	{
		a = &(arr[num]);
		for(i = 0;i < CellExtent;i++)
		{
			for(l = 0;l < CellExtent;l++)
			{
				for(k = 0;k < CellExtent;k++)
				{
					a->M[i][l][k] = 0.0;
				}
			}
		}
	}
}

__device__ void AccumulateCurrentWithParticlesInCell(
									 CellDouble *c_jx,
									 int CellDouble_array_dim,
									 CellDouble *c_jy,
									 CellDouble *c_jz,
									 Cell  *c,
		                             int index,
		                             int blockDimX,
		                             int nt
		                             )
{
	CurrentTensor t1,t2;
	DoubleCurrentTensor dt,dt1;;
    int pqr2;


    while(index < c->number_of_particles)
    {
        c->AccumulateCurrentSingleParticle    (index,&pqr2,&dt);
        writeCurrentComponent(&(c_jx[index%CellDouble_array_dim]),&(dt.t1.Jx),&(dt.t2.Jx),pqr2);
        writeCurrentComponent(&(c_jy[index%CellDouble_array_dim]),&(dt.t1.Jy),&(dt.t2.Jy),pqr2);
        writeCurrentComponent(&(c_jz[index%CellDouble_array_dim]),&(dt.t1.Jz),&(dt.t2.Jz),pqr2);

//        if(index%2 == 0)
//        {
//
//        }
//        else
//        {
//        	if(index%2 == 1)
//        	{
//        		writeCurrentComponent(&(c_jx[1]),&(dt.t1.Jx),&(dt.t2.Jx),pqr2);
//        	}
//        	else
//        	{
//        		writeCurrentComponent(&(c_jx[2]),&(dt.t1.Jx),&(dt.t2.Jx),pqr2);
//        	}
//        }

//        writeCurrentComponent(c_jy,&(dt.t1.Jy),&(dt.t2.Jy),pqr2);
        //writeCurrentComponent(c_jz,&(dt.t1.Jz),&(dt.t2.Jz),pqr2);

        index += blockDimX;
    }
    __syncthreads();


}





__device__ void copyFromSharedMemoryToCell(
		                                     CellDouble *c_jx,
											 CellDouble *c_jy,
											 CellDouble *c_jz,
											 Cell  *c,
				                             int index,
				                             int blockDimX,
				                             dim3 blockId
		)
{
	while(index < CellExtent*CellExtent*CellExtent)
	{
      	copyCellDouble(c->Jx,c_jx,index,blockIdx);
    	copyCellDouble(c->Jy,c_jy,index,blockIdx);
    	copyCellDouble(c->Jz,c_jz,index,blockIdx);

    	index += blockDim.x;
    }
    c->busyParticleArray = 0;
}


__global__ void GPU_StepAllCells(GPUCell  **cells//,
//		                         int i,
//		                         double *global_jx
//		                         double mass,
//		                         double q_mass
		                         )
{
	Cell  *c,*c0 = cells[0];
	__shared__  CellDouble fd[9];
	CellDouble *c_jx,*c_jy,*c_jz,*c_ex,*c_ey,*c_ez,*c_hx,*c_hy,*c_hz;
//	CurrentTensor t1,t2;
//	int pqr2;
	Particle p;

	c = cells[ c0->getGlobalCellNumber(blockIdx.x,blockIdx.y,blockIdx.z)];

	assignSharedWithLocal(&c_jx,&c_jy,&c_jz,&c_ex,&c_ey,&c_ez,&c_hx,&c_hy,&c_hz,fd);




	copyFieldsToSharedMemory(c_jx,c_jy,c_jz,c_ex,c_ey,c_ez,c_hx,c_hy,c_hz,c,
			threadIdx.x,blockIdx,blockDim.x);


	MoveParticlesInCell(c_ex,c_ey,c_ez,c_hx,c_hy,c_hz,
						 c,threadIdx.x,blockDim.x);//,mass,q_mass);
//	MoveAccCurrent(c_ex,c_ey,c_ez,c_hx,c_hy,c_hz,c_jx,c_jy,c_jz,
//							 c,threadIdx.x,blockDim.x,mass,q_mass);


//    WriteCurrents(c_jx,c_jy,c_jz,c_jx,c_jy,c_jz,
//						 c,threadIdx.x,blockDim.x,mass,q_mass);



    copyFromSharedMemoryToCell(c_jx,c_jy,c_jz,c,threadIdx.x,blockDim.x,blockIdx);

}


__global__ void GPU_CurrentsAllCells(GPUCell  **cells,int nt
//		                         int i,
//		                         double *global_jx,
//		                         double mass,
//		                         double q_mass
//		                         CellDouble *
		                         )
{
	Cell  *c,*c0 = cells[0];
	__shared__  CellDouble fd[9];
	CellDouble *c_jx,*c_jy,*c_jz,*c_ex,*c_ey,*c_ez,*c_hx,*c_hy,*c_hz;
//	CurrentTensor t1,t2;
//	int pqr2;
//	Particle p;
	__shared__ CellDouble m_c_jx[CURRENT_SUM_BUFFER_LENGTH];
	__shared__ CellDouble m_c_jy[CURRENT_SUM_BUFFER_LENGTH];
	__shared__ CellDouble m_c_jz[CURRENT_SUM_BUFFER_LENGTH];


	c = cells[ c0->getGlobalCellNumber(blockIdx.x,blockIdx.y,blockIdx.z)];

	assignSharedWithLocal(&c_jx,&c_jy,&c_jz,&c_ex,&c_ey,&c_ez,&c_hx,&c_hy,&c_hz,fd);




	copyFieldsToSharedMemory(c_jx,c_jy,c_jz,c_ex,c_ey,c_ez,c_hx,c_hy,c_hz,c,
			threadIdx.x,blockIdx,blockDim.x);

	set_cell_double_arrays_to_zero(m_c_jx,m_c_jy,m_c_jz,CURRENT_SUM_BUFFER_LENGTH,
			threadIdx.x,blockDim.x);

	AccumulateCurrentWithParticlesInCell(m_c_jx,CURRENT_SUM_BUFFER_LENGTH,m_c_jy,m_c_jz,
		  					 c,threadIdx.x,blockDim.x,nt);

	addCellDouble(c_jx,&(m_c_jx[0]),threadIdx.x,blockIdx,CURRENT_SUM_BUFFER_LENGTH);
	addCellDouble(c_jy,&(m_c_jy[0]),threadIdx.x,blockIdx,CURRENT_SUM_BUFFER_LENGTH);
	addCellDouble(c_jz,&(m_c_jz[0]),threadIdx.x,blockIdx,CURRENT_SUM_BUFFER_LENGTH);
//	addCellDouble(c_jx,&(m_c_jx[1]),threadIdx.x,blockIdx,3);
//	addCellDouble(c_jx,&(m_c_jx[2]),threadIdx.x,blockIdx,3);

    copyFromSharedMemoryToCell(c_jx,c_jy,c_jz,c,threadIdx.x,blockDim.x,blockIdx);

//    if((c->i == 28) && (c->l == 3) && (c->k == 2) && (nt == 1))
//    	{
//    		for(int i = 0;i < CellExtent;i++)
//    		{
//    			for(int l = 0;l < CellExtent;l++)
//    			{
//    				for(int k = 0;k < CellExtent;k++)
//    				{
//    					if(fabs(c_jx->M[i][l][k]) > 1e-15)
//    					{
//    						printf("c_jx %5d %3d %3d %25.15e thread %5d %3d %3d block %5d %3d %3d nt %3d\n",
//    								i,l,k,
//    								c_jx->M[i][l][k],
//    								threadIdx.x,threadIdx.y, threadIdx.z,
//    								blockIdx.x,blockIdx.y,blockIdx.z,nt);
//    					}
//    				}
//    			}
//    		}
//    	}

}


__global__ void GPU_ControlAllCellsCurrents(Cell  **cells,int n,int i,CellDouble *jx,CellDouble *jy,CellDouble *jz)
{
//	unsigned int nx = blockIdx.x;
//	unsigned int ny = blockIdx.y;
//	unsigned int nz = blockIdx.z;
//	int i,l,k;
	Cell  *c,*c0 = cells[0],nc;
	//double t;
//	__shared__ extern CellDouble fd[9];
	//double *src;
	//int pqr2;
//	CurrentTensor t1,t2;

//	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];
	c = cells[ n ];

	nc = *c;

	// double cjx,cjy,cjz;

//	              cjx = CheckArraySize((double *)jx,(double *)(nc.Jx),sizeof(CellDouble)/sizeof(double));
//	              cjy = CheckArraySize((double *)jy,(double *)(nc.Jy),sizeof(CellDouble)/sizeof(double));
//	              cjz = CheckArraySize((double *)jz,(double *)(nc.Jz),sizeof(CellDouble)/sizeof(double));
#ifdef GPU_CONTROL_ALL_CELLS_CURRENTS_PRINT
//	              printf("cell (%d,%d,%d) particle %d currents %.5f %.5f %.5f \n",nc.i,nc.l,nc.k,i,cjx,cjy,cjz);
#endif


}

__host__ __device__
void emh2_Element(
		Cell *c,
		int i,int l,int k,
		double *Q,double *H)
{
	int n  = c->getGlobalCellNumber(i,l,k);

	H[n] += Q[n];
}


__global__
void GPU_emh2(
		 GPUCell  **cells,
				            int i_s,int l_s,int k_s,
							double *Q,double *H
		)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell  *c0 = cells[0];

	emh2_Element(c0,i_s+nx,l_s+ny,k_s+nz,Q,H);
}


__host__ __device__
void emh1_Element(
		Cell *c,
		int3 i,
		double *Q,double *H,double *E1, double *E2,
		double c1,double c2,
		int3 d1,int3 d2)
{

    int n  = c->getGlobalCellNumber(i.x,i.y,i.z);
	int n1 = c->getGlobalCellNumber(i.x+d1.x,i.y+d1.y,i.z+d1.z);
	int n2 = c->getGlobalCellNumber(i.x+d2.x,i.y+d2.y,i.z+d2.z);

	double e1_n1 = E1[n1];
	double e1_n  = E1[n];
	double e2_n2 = E2[n2];
	double e2_n  = E2[n];

	double t  = 0.5*(c1*(e1_n1 - e1_n)- c2*(e2_n2 - e2_n));
    Q[n] = t;
    H[n] += Q[n];
}


__global__
void GPU_emh1(
		 GPUCell  **cells,
							double *Q,double *H,double *E1, double *E2,
							double c1,double c2,
							int3 d1,int3 d2
		)
{

	int3 i3 = make_int3(blockIdx.x,blockIdx.y,blockIdx.z);
	Cell  *c0 = cells[0];

	emh1_Element(c0,i3,Q,H,E1,E2,c1,c2,d1,d2);
}

__host__ __device__
	void emeElement(Cell *c,int3 i,double *E,double *H1, double *H2,
			double *J,double c1,double c2, double tau,
			int3 d1,int3 d2
			)
	{
	   int n  = c->getGlobalCellNumber(i.x,i.y,i.z);
	  int n1 = c->getGlobalCellNumber(i.x+d1.x,i.y+d1.y,i.z+d1.z);
	  int n2 = c->getGlobalCellNumber(i.x+d2.x,i.y+d2.y,i.z+d2.z);

	  E[n] += c1*(H1[n] - H1[n1]) - c2*(H2[n] - H2[n2]) - tau*J[n];
	}

__host__ __device__
void periodicElement(Cell *c,int i,int k,double *E,int dir, int to,int from)
{
    int n   = c->getGlobalBoundaryCellNumber(i,k,dir,to);
	int n1  = c->getGlobalBoundaryCellNumber(i,k,dir,from);
	E[n]    = E[n1];
}

__global__ void GPU_periodic(GPUCell  **cells,
                             int i_s,int k_s,
                             double *E,int dir, int to,int from)
{
	unsigned int nx = blockIdx.x;
	//unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell  *c0 = cells[0];

	periodicElement(c0,nx+i_s,nz+k_s,E, dir,to,from);
}

__host__ __device__
void periodicCurrentElement(Cell *c,int i,int k,double *E,int dir, int dirE,int N)
{
    int n1    = c->getGlobalBoundaryCellNumber(i,k,dir,1);
    int n_Nm1 = c->getGlobalBoundaryCellNumber(i,k,dir,N-1);
    if(dir != dirE)
    {
       E[n1] += E[n_Nm1];
    }
    if(dir != 1 || dirE != 1)
    {
       E[n_Nm1] =  E[n1];
    }

    int n_Nm2 = c->getGlobalBoundaryCellNumber(i,k,dir,N-2);
    int n0    = c->getGlobalBoundaryCellNumber(i,k,dir,0);

#ifdef PERIODIC_CURRENT_PRINTS
    printf("%e %e \n",E[n0],E[n_Nm2]);
#endif
    E[n0] += E[n_Nm2];
    E[n_Nm2] = E[n0];
}


__global__ void GPU_CurrentPeriodic(GPUCell  **cells,double *E,int dirE, int dir,
                             int i_s,int k_s,int N)
{
	unsigned int nx = blockIdx.x;
	//unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell  *c0 = cells[0];


	periodicCurrentElement(c0,nx+i_s,nz+k_s,E, dir,dirE,N);
}


__global__ void GPU_eme(

		            GPUCell  **cells,
		            int3 s,
					double *E,double *H1, double *H2,
					double *J,double c1,double c2, double tau,
					int3 d1,int3 d2
		)
{
	unsigned int nx = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int ny = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int nz = blockIdx.z*blockDim.z + threadIdx.z;
	Cell  *c0 = cells[0];

    s.x += nx;
    s.y += ny;
    s.z += nz;

    emeElement(c0,s,E,H1,H2,J,c1,c2,tau,d1,d2);
}


__global__ void copy_pointers(Cell  **cells,int *d_flags,double_pointer *d_pointers)
{
	Cell  *c = cells[blockIdx.x];

	c->flag_wrong_current_cell = d_flags[blockIdx.x];
	c->d_wrong_current_particle_attributes = d_pointers[blockIdx.x];

}

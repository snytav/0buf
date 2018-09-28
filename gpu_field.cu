
__global__ void GPU_eme(
		            Plasma *gp,
		            GPUCell  **cells,
		            int i_s,int l_s,int k_s,
					double *E,double *H1, double *H2,
					double *J,double c1,double c2, double tau,
					int dx1,int dy1,int dz1,int dx2,int dy2,int dz2
		)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell<Particle>  *c0 = cells[0];




	gp->emeElement(*c0,i_s+nx,l_s+ny,k_s+nz,E,H1,H2,
			    	  		J,c1,c2,tau,
			    	  		dx1,dy1,dz1,dx2,dy2,dz2);
}





void Plasma::copyCells(std::string where,int nt)
{
	static int first = 1;
	size_t m_free,m_total;
	int size = (*AllCells).size();
	struct sysinfo info;
//	unsigned long c1,c2;

    if(first == 1)
    {
    	cp = (GPUCell **)malloc(size*sizeof(GPUCell *));
    }

	unsigned long m1,m2,delta,accum;
	memory_monitor("beforeCopyCells",nt);

	for(int i = 0;i < size;i++)
	{
//		if(i == 141)
//		{
//			int qq = 0;
//		}
		cudaError_t err = cudaMemGetInfo(&m_free,&m_total);
		sysinfo(&info);
		m1 = info.freeram;
	 	GPUCell c,*d_c,*z0;
	 	z0 = h_CellArray[i];
	 	if(first == 1)
	 	{
	       d_c = c.allocateCopyCellFromDevice();
     	   cp[i] = d_c;
	 	}
	 	else
	 	{
	 	   d_c = cp[i];
	 	}
	    c.copyCellFromDevice(z0,d_c,where,nt);
		m2 = info.freeram;

	    delta = m1-m2;
        accum += delta;

	}

	if(first == 1)
	{
		first = 0;
	}

	memory_monitor("afterCopyCells",nt);
}

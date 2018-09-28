

__device__ double cuda_atomicAdd(double *address, double val)
{
    double assumed,old=*address;
    do {
        assumed=old;

        old=

        		__longlong_as_double
        		(
        		atomicCAS(
        				(unsigned long long int*)
        				address,
                    __double_as_longlong(assumed),
                    __double_as_longlong(val+assumed)));
    }while (assumed!=old);

    //printf("NEW ATOMIC ADD\n");

    return old;
}

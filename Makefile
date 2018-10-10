LD=mpicxx
CUDACC=nvcc
CXX=mpicxx
CPP=g++
CUDA=/usr/local/cuda
CUDALIB=$(CUDA)/lib64

LDFLAGS= -lm -L$(CUDALIB) -arch=sm_61
#CUDAFLAGS= --maxrregcount=128  -arch=sm_35 --ptxas-options=-v -I/usr/local/cuda-7.5/include/
CUDAFLAGS= -arch=sm_61 -lineinfo --maxrregcount=128 -g -I$(CUDA)/include/ #--relocatable-device-code=true
CUDALIBS=  -g -L$(CUDALIB) -lcuda -lcudart #-lthrust 
MPIFLAGS=
CFLAGS=

OBJ = main.o mpi_shortcut.o service_functions.o compare.o maxwell.o load_data.o archAPI.o 
#plasma.o
            
main.o: main.cu $(DEPS)
	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS) 
	
kernels.o: kernels.cu $(DEPS)
	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS) 	--relocatable-device-code=true
	
#compare.o: compare.cu $(DEPS)
#	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS)
	
load_data.o: load_data.cu $(DEPS)
	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS)	
	

service_functions.o: service_functions.cu $(DEPS)
	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS)	

archAPI.o: archAPI.cu $(DEPS)
	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS)	
#kernels.o: kernels.cu $(DEPS)
#	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS) 		 	
                    
#plasma.o: plasma.cu $(DEPS)
#	$(CUDACC) -g -c -o $@ $< $(CUDAFLAGS)                     
                    
%.o: %.cxx $(DEPS)
	$(CXX) -g -c -o $@ $< $(MPIFLAGS) 

%.o: %.cpp $(DEPS)
	$(CPP) -g -c -o $@ $< $(CBFLAGS) 
                            
all: $(OBJ)
	$(LD) -g -o $@ $^ $(CFLAGS) $(DBFLAGS) $(CUDALIBS) 

clean:
	rm *.o all    

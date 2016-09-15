CXX = g++
CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) -I/usr/include/cuda -I/usr/local/cuda/include/ -I. -D GPU
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libh5nb6xx.a

CODEOBJS = h5nb6xx_helper.o cuda_util.o

AR = ar ruv
RANLIB = ranlib
RM = rm

cuda_util.o: cuda_util.cu
	nvcc -arch=sm_21 -c cuda_util.cu -I/usr/include/cuda -I/usr/local/cuda/include -I. -D GPU

h5nb6xx_helper.o: h5nb6xx_helper.cc
	g++ -c h5nb6xx_helper.cc $(CXXFLAGS)

all: $(CODELIB) 


clean:
	$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

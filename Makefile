# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

HDF5_LIBS += -lhdf5 -lz -ldl -L/usr/local/cuda/lib -lcuda -lcudart -stdlib=libstdc++ 
HDF5_FLAGS=-I/opt/local/include -I/Developer/NVIDIA/CUDA-6.5/samples/common/inc -I.

OBJS = interface.o 

CODELIB = src/libh5nb6xx.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: h5nb6xx_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ h5nb6xx_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py H5nb6xxInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py H5nb6xxInterface -o $@

h5nb6xx_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@ $(HDF5_FLAGS) $(HDF5_LIBS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(HDF5_FLAGS) $(HDF5_LIBS)

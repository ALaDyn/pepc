CC=xlc
SIONLIB_INC=`sionconfig --cflags --ser -be`
SIONLIB_LIB=`sionconfig --libs --ser -be`

all: pepc_sion2vtk

pepc_sion2vtk: pepc_sion2vtk.c
	$(CC) $(SIONLIB_INC) pepc_sion2vtk.c -o pepc_sion2vtk $(SIONLIB_LIB)
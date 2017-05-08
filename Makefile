CC = mpic++
### config on the stampede cluster
#BDIR = ${TACC_BOOST_MPI_DIR}
### config on the spear cluster
#BDIR = /lustre/home-2/hyang3/workdir/boost
#LIBS = -I${BDIR}/include/ -L${BDIR}/lib/ -lboost_mpi -lboost_serialization
### config on the spear cluster
#OptizelleDIR = /lustre/home-2/hyang3/workdir/Optizelle/install_release
#LIBS += -I$(OptizelleDIR)/include/ -L$(OptizelleDIR)/lib/ -loptizelle -ljson

### config on laptop
LIBS =  -lboost_mpi -lboost_serialization
##OptizelleDIR = /home/huanhuan/Packages/Optizelle-1.1.2/install-release
##LIBS += -I$(OptizelleDIR)/include/ -L$(OptizelleDIR)/lib/ -loptizelle -ljson

ifeq ($(NETCDF),yes)
NCDFDIR = /lustre/home-2/hyang3/share/third-party/netcdf-install
LIBS += -I$(NCDFDIR)/include/ -L$(NCDFDIR)/lib/ -lnetcdf -lnetcdf_c++
NCDF_FLAGS = -DUSE_NETCDF
endif

FLAGS = -std=c++11 -m64 $(NCDF_FLAGS)
DFLAGS = -g -m64 -D_DEBUG
TRISRC=Triangle/
XMLSRC=Pugixml/
BFGSSRC=HLBFGS/

#PLATFORM=_MACOS
PLATFORM=_LINUX

ifeq ($(PLATFORM),_LINUX)
	FLAGS = -std=c++11 -m64 -DLINUX $(NCDF_FLAGS)
	DFLAGS = -g -m64 -D_DEBUG -DLINUX
endif

ifeq ($(PLATFORM),_MACOS)
	FLAGS = -std=c++11 -m64 $(NCDF_FLAGS)
	DFLAGS = -g -m64 -D_DEBUG
endif

TRILIBDEFS= -DTRILIBRARY


all: trilibrary HLBFGS
	${CC} scvt_HLBFGS.cpp  ${TRISRC}triangle.o $(BFGSSRC)HLBFGS.o $(BFGSSRC)LineSearch.o $(BFGSSRC)HLBFGS_BLAS.o $(BFGSSRC)ICFS.o ${LIBS} ${FLAGS} -o scvt_HLBFGS.exe

#QNewton_Optiz: trilibrary pugixml-library
#	${CC} QNewton_scvt_Optizelle.cpp ${TRISRC}triangle.o ${XMLSRC}pugixml.o ${LIBS} ${FLAGS} -o QNewton_scvt_Optizelle.exe

#Lloyd_scvt: trilibrary pugixml-library
#	${CC} Lloyd_scvt.cpp ${TRISRC}triangle.o ${XMLSRC}pugixml.o ${LIBS} ${FLAGS} -o Lloyd_scvt.exe

scvt_Lloyd: trilibrary
	${CC} scvt_Lloyd.cpp ${TRISRC}triangle.o ${LIBS} ${FLAGS} -o scvt_Lloyd.exe

refine: trilibrary
	${CC} refine.cpp ${TRISRC}triangle.o  ${LIBS} ${FLAGS} -o refine.exe

partition_gen: trilibrary
	${CC} partition_gen.cpp ${TRISRC}triangle.o  ${LIBS} ${FLAGS} -o partition_gen.exe

num_estimate:
	${CC} num_estimate.cpp ${LIBS} ${FLAGS} -o num_estimate.exe

#debug: trilibrary-debug pugixml-library-debug
#	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${XMLSRC}pugixml.o ${LIBS} ${DFLAGS} -o ${EXE}

trilibrary:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${FLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

trilibrary-debug:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${DFLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

#pugixml-library:
#	$(CC) $(CSWITCHES) $(FLAGS) -c -o $(XMLSRC)pugixml.o $(XMLSRC)pugixml.cpp

#pugixml-library-debug:
#	$(CC) $(CSWITCHES) $(DFLAGS) -c -o $(XMLSRC)pugixml.o $(XMLSRC)pugixml.cpp

Lite_Sparse_Matrix:
	$(CC) $(CSWITCHES) $(FLAGS) -c -o $(BFGSSRC)Lite_Sparse_Matrix.o $(BFGSSRC)Lite_Sparse_Matrix.cpp
HLBFGS: LineSearch ICFS
	$(CC) $(CSWITCHES) ${LIBS} $(FLAGS) -c -o $(BFGSSRC)HLBFGS.o $(BFGSSRC)HLBFGS.cpp
LineSearch: HLBFGS_BLAS
	$(CC) $(CSWITCHES) ${LIBS} $(FLAGS) -c -o $(BFGSSRC)LineSearch.o $(BFGSSRC)LineSearch.cpp
HLBFGS_BLAS:
	$(CC) $(CSWITCHES) ${LIBS} $(FLAGS) -c -o $(BFGSSRC)HLBFGS_BLAS.o $(BFGSSRC)HLBFGS_BLAS.cpp
ICFS:
	$(CC) $(CSWITCHES) $(FLAGS) -c -o $(BFGSSRC)ICFS.o $(BFGSSRC)ICFS.cpp

clean:
	rm -f *.dat ${EXE} ${TRISRC}triangle.o $(BFGSSRC)HLBFGS.o $(BFGSSRC)LineSearch.o $(BFGSSRC)HLBFGS_BLAS.o $(BFGSSRC)ICFS.o

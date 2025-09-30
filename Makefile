HEADER   = $(HOME)/WORK/PROGRAMMATION/SuperLU/SuperLU_4.3/SRC
SuperLUroot	= $(HOME)/WORK/PROGRAMMATION/SuperLU/SuperLU_4.3
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_4.3.a
LIBS = $(SUPERLULIB) -isysroot $(SDKROOT) -framework Accelerate -lm
#LIBS	= $(SUPERLULIB) -framework Accelerate -lm
#LIBS	= $(SUPERLULIB)
SDKROOT = $(shell xcrun --sdk macosx --show-sdk-path)


CC           = gcc-mp-11
CFLAGS = -O2 -isysroot $(SDKROOT)
FORTRAN	     = gfortran-mp-11
#FFLAGS       = -c -O2 -fcheck=all
FFLAGS = -c -O2 -fallow-argument-mismatch -isysroot $(SDKROOT)


LOADOPTS     =
CDEFS        = -DAdd_



OBJ = main.f Vag3D.o HFV3D.o Geometry.o StructureCreuse.o JacSm-PT.o Meca.o Lois.o NewtonMecaVBulle.o NewtonMecaNitsche.o ProbingRelaxFracVB.o ProbingRelaxFracNitsche.o NitscheUtils.o SolvePT.o CouplagePoroMeca.o JacSmTraceur.o Visu.o  NumerotationInc.o quicksort1.o c_fortran_dgssv.o SolveurSuperLU.o SolveurCPRAMG.o amg1r5.o dslugm.o dgmres.o  dchkw.o ds2y.o dsilus.o dslui.o dsmv.o dhels.o dheqr.o dorth.o dpigmr.o drlcal.o dxlcal.o isdgmr.o qs2i1d.o dslui2.o daxpy.o dcopy.o ddot.o dnrm2.o dscal.o



frac: $(OBJ) $(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(OBJ) $(LIBS) -o $@


.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o frac





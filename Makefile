## Copy f90 error codes to IPOPT
.SUFFIXES: .f90

SF =  -lblas -llapack -lm -lpthread -lm
CF =  -O3 -g  -ffast-math -finit-local-zero


GCC = gcc -std=c99 -I. -O3 -g
GFORTRAN = gfortran -I.  -g

############################################################################
#Include stuff for IPOPT
############################################################################ -fdefault-real-8 $(SF)

# Fotran Compiler options
FFLAGS = -rdynamic -g -ffast-math -funroll-loops -fomit-frame-pointer -mtune=native -march=native -mfpmath=sse -pipe -ffunction-sections -fdata-sections -mavx -ldl

INCL = `PKG_CONFIG_PATH=/Users/lucamazzone/Desktop/Work/CoinIpopt/Ipopt/lib/pkgconfig/:/Users/lucamazzone/Desktop/Work/CoinIpopt/Ipopt/lib/pkgconfig/:/Users/lucamazzone/Desktop/Work/CoinIpopt/Ipopt/lib/pkgconfig/: pkg-config --cflags ipopt` $(ADDINCFLAGS)

#additional Fortran Compiler options for linking
FCLINKFLAGS = #-Wl,--rpath -Wl,../Ipopt-3.12.5/build_alpha/lib

IPOPT_LIBS = `PKG_CONFIG_PATH=../Ipopt-3.12.5/build_alpha/lib64/pkgconfig:../Ipopt-3.12.5/build_alpha/lib/pkgconfig:../Ipopt-3.12.5/build_alpha/share/pkgconfig: pkg-config --libs ipopt` -lstdc++ -lm


############################################################################
all: growth_model.exec
############################################################################

params.o: params.f90
	$(GFORTRAN) -c $< -o $@

KS_econ.o: KS_econ.f90
	$(GFORTRAN) -c $< -o $@

#ipopt_lib.o: add_lib.f90 KS_econ.o
#	$(GFORTRAN) -c $< -o $@

KS_lib.o: KS_lib.f90 KS_econ.o
	$(GFORTRAN) -c $< -o $@

hs071.o: hs071.f90 KS_lib.o  params.o KS_econ.o  #ipopt_lib.o
	$(GFORTRAN) $(INCL) -c $< -o $@


############################################################################
##executable
growth_model.exec: hs071.o KS_lib.o  KS_econ.o params.o    # ipopt_lib.o
	$(GFORTRAN) $^ -o $@ -lm $(FFLAGS) $(IPOPT_LIBS) $(FCLINKFLAGS)
############################################################################  $(FCLINKFLAGS)


clean:
	rm -rf *.o *.mod *.exec *.f90.*  *.c.*
############################################################################
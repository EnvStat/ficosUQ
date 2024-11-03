#
# FICOS water quality core build recipe
#
# Kai Rasmus 2016
# Janne Ropponen 2019
#

#------------------------------
# Environment specific settings
#------------------------------

# (Karel Kaurila, 2024): Changed LIBDIR to relative path. hdf5 libraries need to be installed there first.
LIBDIR=../hdf5/hdf5-1.8.19/hdf5
COMP=$(LIBDIR)/bin/h5fc # hdf5 build with 64-bit default integers

#-------------
# Flags
#-------------

# Debug compiling
FFLAGS=-fdefault-integer-8 -ffpe-trap=invalid,zero,overflow -ggdb -Og -Wall -fbounds-check -static

# Optimized compiling
#FFLAGS=-fdefault-integer-8 -ffpe-trap=invalid,zero,overflow -O3 -Wall -static -march=native -pipe

#----------------
# Targets
#----------------
TARGETficos=wqficos
OBJECTSficos=prec.o wqficos_helper.o wqficos_names.o errors.o hdfio.o wqficos_tc.o rk4.o equations_ficos.o ini.o

#----------------
# Compiling
#----------------
all: $(TARGETficos)
.PHONY: all

$(TARGETficos): $(OBJECTSficos)
	$(COMP) -o $(TARGETficos) $(OBJECTSficos)

# Compile modules first
prec.o: prec.f90
	$(COMP) $(FFLAGS) -c prec.f90

wqficos_names.o: prec.o wqficos_names.f90
	$(COMP) $(FFLAGS) -c wqficos_names.f90

wqficos_helper.o: prec.o wqficos_names.o wqficos_helper.f90
	$(COMP) $(FFLAGS) -c wqficos_helper.f90

ini.o: ini.f90
	$(COMP) $(FFLAGS) -c ini.f90

hdfio.o: prec.o wqficos_helper.o hdfio.f90
	$(COMP) $(FFLAGS) -c hdfio.f90

errors.o: errors.f90
	$(COMP) $(FFLAGS) -c errors.f90

equations_ficos.o: prec.o equations_ficos.f90
	$(COMP) $(FFLAGS) -c equations_ficos.f90

wqficos_tc.o: prec.o wqficos_helper.o wqficos_names.o errors.o hdfio.o wqficos_tc.f90
	$(COMP) $(FFLAGS) -c wqficos_tc.f90

rk4.o: prec.o rk4.f90
	$(COMP) $(FFLAGS) -c rk4.f90

$(TARGET): $(OBJECTS)
	$(COMP) -o $(TARGET) $(OBJECTS)

.PHONY: clean
clean:
	rm -f *.o *.mod $(TARGETficos)

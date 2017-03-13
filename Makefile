# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = mpifort
# Flags for debugging or for maximum performance
FCFLAGS = -g -O0 -Waliasing -Wampersand -Wcharacter-truncation -Winteger-division -Wsurprising -fimplicit-none -fcheck=no-array-temps -fbacktrace -fdefault-real-8
FLIBS   = -llapack -lblas

# List of executables to be built within the package
PROGRAMS = codice

# "make" builds all
all: $(PROGRAMS)

# ======================================================================
# Here comes the most interesting part: the rules for the various files
# ======================================================================
# [...] when an executable is built from many
# sources, ALL the object files have to be specified in the executable
# dependencies, not just those containing SUBROUTINEs, FUNCTIONs or
# MODULEs that are directly CALLed or USEd in the main program
# file. However among the dependencies of the main program object
# file, only those that are directly CALLed or USEd inside it have to
# be specified, the other object files have to be specified as
# dependencies of the relevant object files that require them and will
# be built by a chain rule.

codice:      codice.o  MPI_module_mod.o grid_mod.o essential_mod.o compact_mod.o bandedmatrix_mod.o input_mod.o ic_and_bc_mod.o set_pressure_mod.o solve_pressure_mod.o library_mod.o thomas_mod.o SPIKE_mod.o variables_mod.o diffusive_term_mod.o convective_term_mod.o time_advancement_mod.o output_mod.o
codice.o:              MPI_module_mod.o grid_mod.o essential_mod.o compact_mod.o                    input_mod.o ic_and_bc_mod.o set_pressure_mod.o solve_pressure_mod.o                            SPIKE_mod.o variables_mod.o diffusive_term_mod.o convective_term_mod.o time_advancement_mod.o output_mod.o
grid_mod.o:            MPI_module_mod.o            essential_mod.o
compact_mod.o:         MPI_module_mod.o grid_mod.o 	                        bandedmatrix_mod.o                                                                     library_mod.o
input_mod.o:           MPI_module_mod.o grid_mod.o                 compact_mod.o                                ic_and_bc_mod.o set_pressure_mod.o
ic_and_bc_mod.o:       MPI_module_mod.o grid_mod.o                                                                                                                                                             variables_mod.o convective_term_mod.o
MPI_module_mod.o:                                  essential_mod.o
set_pressure_mod.o:    MPI_module_mod.o grid_mod.o essential_mod.o compact_mod.o bandedmatrix_mod.o                                                                     library_mod.o                          variables_mod.o
solve_pressure_mod.o:  MPI_module_mod.o grid_mod.o                               bandedmatrix_mod.o                             set_pressure_mod.o
SPIKE_mod.o:           MPI_module_mod.o grid_mod.o                 compact_mod.o bandedmatrix_mod.o                                                                                   thomas_mod.o
variables_mod.o:       MPI_module_mod.o            essential_mod.o compact_mod.o
diffusive_term_mod.o:  MPI_module_mod.o variables_mod.o compact_mod.o bandedmatrix_mod.o essential_mod.o SPIKE_mod.o
convective_term_mod.o: MPI_module_mod.o variables_mod.o compact_mod.o bandedmatrix_mod.o essential_mod.o SPIKE_mod.o
time_advancement_mod.o: variables_mod.o bandedmatrix_mod.o solve_pressure_mod.o convective_term_mod.o diffusive_term_mod.o
output_mod.o:	MPI_module_mod.o variables_mod.o

# ======================================================================
# And now the general rules
# ======================================================================
# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(FLIBS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) $(FLIBS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f *~ $(PROGRAMS)

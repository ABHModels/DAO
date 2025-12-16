FC = gfortran

# If not mods in current path, use mkdir -p mods
MODDIR = mods

FLAGS = -O3 -I$(MODDIR) -J$(MODDIR)

# Path to the heasfot
HEADAS = /data/heasoft/heasoft-6.31.1/x86_64-pc-linux-gnu-libc2.17/

# Please find specific version of these lib
LIBS = $(HEADAS)/lib/libcfitsio.so 	\
       $(HEADAS)/lib/libxanlib_6.31.so \
       $(HEADAS)/lib/libape_2.9.so \
	   -lm -lc -fopenmp

NAME = dao
CODE = $(NAME).f90
EXEC = $(NAME).x

MODSRC = global
MOD_F90 = $(shell find $(MODSRC) -name '*.f90')  
MODOBJ = $(MOD_F90:.f90=.o)

SRCDIR = modellib
SRC_F90 = $(shell find $(SRCDIR) -name '*.f90') 
SRC_F   = $(shell find $(SRCDIR) -name '*.f')
SRC_OBJ = $(SRC_F90:.f90=.o) $(SRC_F:.f=.o) 

OBJECTS = $(MODOBJ) $(SRC_OBJ) $(NAME).o

$(EXEC): $(OBJECTS)
	$(FC) $(FLAGS) $(LIBS) -o $@ $^ $(EXT_LIBS)
	@echo "Compilation successful: $(EXEC) has been built."

$(MODSRC)/%.o: $(MODSRC)/%.f90
	$(FC) $(FLAGS) -c $< -o $@

# Subroutine
$(SRCDIR)/%.o: $(SRCDIR)/%.f90 $(MODOBJ)
	$(FC) $(FLAGS) -c $< -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.f $(MODOBJ)
	$(FC) $(FLAGS) -c $< -o $@

# Main 
$(NAME).o: $(CODE) $(MODOBJ)
	$(FC) $(FLAGS) -c $< -o $@

# Clean
clean:
	find $(MODSRC) -name '*.o' -delete
	find $(MODSRC) -name '*.mod' -delete
	find $(SRCDIR) -name '*.o' -delete
	rm $(NAME)*.o $(NAME)*.x


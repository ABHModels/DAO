FC = gfortran
FLAGS = -fallow-argument-mismatch -O3 -fopenmp

HEADAS = /data/heasoft/heasoft-6.31.1/x86_64-pc-linux-gnu-libc2.17/


LIBS = $(HEADAS)/lib/libcfitsio.so 	\
       $(HEADAS)/lib/libxanlib_6.31.so \
       $(HEADAS)/lib/libape_2.9.so \
	   -lm -lc -fopenmp


NAME = drive_SRF
CODE = $(NAME).f
EXEC = $(NAME).x


SRCDIR = ./my-routines
SRC_F90 = $(shell find $(SRCDIR) -name '*.f90')  
SRC_OBJ = $(SRC_F90:.f90=.o)

SRC_F = $(shell find $(SRCDIR) -name '*.f')  
SRC_OBJf = $(SRC_F:.f=.o)

OBJECTS =  $(SRC_OBJ) $(SRC_OBJf) $(NAME).o

$(EXEC): $(OBJECTS)
	$(FC) $(FLAGS) $(LIBS) -o $@ $^ $(EXT_LIBS)
	@echo "Compilation successful: $(EXEC) has been built."


# Model lib
$(SRCDIR)/%.o: $(SRCDIR)/%.f90 
	$(FC) $(FLAGS) -c $< -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.f
	$(FC) $(FLAGS) -c $< -o $@

# Main
$(NAME).o: $(CODE)
	$(FC) $(FLAGS) -c $< -o $@

# 清理
clean:
	find $(SRCDIR) -name '*.o' -delete
	find $(SRCDIR) -name '*.mod' -delete
	rm -f $(NAME)*.o $(NAME)*.x
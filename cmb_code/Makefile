#####################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
#####################################################

# Set compiler
# CC = g++-mp-4.8 

# DEFINE COSMOLOGICAL MODEL
MODEL    = LCDM
RECMODEL = STANDARD

##########################################################################
############################## OPTIONS ###################################
##########################################################################
## i)   LU Options available are -DUSELAPACK 	
##      Link ot lapack with -llapack -lblas -lf2c
##      This can speed the code up, but doesnt change anything else
##
## ii)  To turn on messages in bash use -DVERBOSE			
##
## iii) For output of various quantities use -DOUTPUT_XXX		
##      XXX = BACKGROUND, PERTURBATIONS, RECOMBINATION			
##      -DOUTPUT_BACKGROUND -DOUTPUT_PERTURBATIONS -DOUTPUT_RECOMBINATION
##
## iv)  To turn off output of matter- and/or TT power-spectrum	   	
##      Use -DNO_OUTPUT_MATTERPOW and/or -DNO_OUTPUT_CL			
##
##  v) -DBESSELLIB to use complex_bessel library to calc. sph. j_nu
## 		 Link with -lgfortran -lcomplex_bessel -lm
##########################################################################

# OPTIONS
OPTIONS = -DVERBOSE -DBESSELLIB
C = -I$(HOME)/local/include -O3 $(OPTIONS)
L = -L$(HOME)/local/lib -lm

# COMPLEX BESSELLIB
L += -lgfortran -lcomplex_bessel

# TARGETS
TARGETS := CMFB_$(MODEL)
all: $(TARGETS)

# PATH TO SEARCH FOR FILES
VPATH = source:source/source_general:source/source_model_$(MODEL):source/source_recmodel_$(RECMODEL)

# OBJECT FILES
OBJS = driver_$(MODEL).o Spline.o LUdcmp.o Background.o Background$(MODEL).o Recombination.o Recombination$(RECMODEL).o Perturbations.o Perturbations$(MODEL).o Powerspectrum.o

# DEPENDENCIES
driver_$(MODEL)		         : Parameters.h Background.h Recombination.h Perturbations.h Powerspectrum.h
Background.o               : Parameters.h Background.h
Background$(MODEL).o       : Parameters.h Background.h Background$(MODEL).h
Recombination.o            : Parameters.h Background.h
Recombination$(RECMODEL).o : Parameters.h Background.h Recombination.h Recombination$(RECMODEL).h
Perturbations.o            : Parameters.h Background.h Recombination.h
Perturbations$(MODEL).o    : Parameters.h Background.h Recombination.h Perturbations.h Perturbations$(MODEL).h
Powerspectrum.o            : Parameters.h Background.h Recombination.h Perturbations.h Powerspectrum.h Sphbessel.h

# HEADERS
HEADERS = Parameters.h Spline.h Matrix.h Odeint.h Background.h Recombination.h Perturbations.h Powerspectrum.h Background$(MODEL).h Recombination$(RECMODEL).h Perturbations$(MODEL).h LUdcmp.h

CMFB_$(MODEL): $(OBJS)
  ${CC} -o $@ $^ $C $L $(SFML)

%.o: %.cpp $(HEADERS)
	${CC} -c -o $@ $< $C

clean:
	rm -rf $(TARGETS) *.o

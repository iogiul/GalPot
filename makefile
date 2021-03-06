
#include make.GalPot


all: testGalPot.exe findOrbit.exe findOrbitProperties.exe \
		findManyOrbitProperties.exe findManyOrbitPropertiesfromGalactic.exe \
		findManyOrbitPropertiesfromEquatorial.exe \
		findManyOrbitPropertiesfromEquatorialwErrors.exe \
		findManyOrbitPropertiesfromRAVEMultiGaussian.exe \
		findManyOrbitPropertiesfromGaiaSimple.exe Coord_converter.exe RotCurve.exe \
		findOrbitMultiPot.exe Orbit_parallel.exe Orbit_parallel_worbit_sample.exe Orbit_parallel.exe  \
		GalPot Other

CPP		= g++
LIBPOT		= obj/libPot.a
LIBOTHER	= obj/libOther.a

CFLAGS        = -c -o $@.o -O3 -ffast-math -fopenmp -Isrc/
CFLAGSKEEP    = -c -o $@ -O3 -ffast-math -fopenmp -Isrc/

#CMPFLAGS        = -c -o $@.o -O3 -ffast-math -fopenmp -Isrc/

MFLAGS	= -O3 -ffast-math -fopenmp -Isrc/

LDFLAGS      = -o $@ -Lobj/ -lPot -lOther -lm

# commands to put file into library
AR            = ar r
RL            = ranlib

ARPOT     = $(AR) $(LIBPOT) $@.o; $(RL) $(LIBPOT); touch $@
AROTHER   = $(AR) $(LIBOTHER) $@.o; $(RL) $(LIBOTHER); touch $@
AWAY	  = rm $@.o

AUXIL_H		   = src/Pi.h src/Inline.h src/FreeMemory.h src/Vector.h \
		     src/Matrix.h src/Numerics.templates src/Numerics.h \
		     src/Pspline.h

GalPot_h	   = src/GalPot.h src/Potential.h

obj/Numerics:		src/Numerics.cc $(AUXIL_H)
			$(CPP) $(CFLAGS) src/Numerics.cc;$(ARPOT);$(AWAY)
obj/WDMath:		src/WDMath.cc $(AUXIL_H)
			$(CPP) $(CFLAGS) src/WDMath.cc;$(ARPOT);$(AWAY)

obj/GalPot: 		src/GalPot.cc $(GalPot_h) $(AUXIL_H)
	      		$(CPP) $(CFLAGS) src/GalPot.cc;$(ARPOT);$(AWAY)

obj/KeplerPot: 		src/KeplerPot.cc src/KeplerPot.h src/Potential.h $(AUXIL_H)
					$(CPP) $(CFLAGS) src/KeplerPot.cc;$(ARPOT);$(AWAY)
obj/MiyamotoNagaiPot: 		src/MiyamotoNagaiPot.cc src/MiyamotoNagaiPot.h src/Potential.h $(AUXIL_H)
			$(CPP) $(CFLAGS) src/MiyamotoNagaiPot.cc;$(ARPOT);$(AWAY)

obj/MultiPot: 		src/MultiPot.cc src/MultiPot.h src/Potential.h $(AUXIL_H)
			$(CPP) $(CFLAGS) src/MultiPot.cc;$(ARPOT);$(AWAY)
GalPot: 		obj/Numerics obj/WDMath obj/GalPot obj/KeplerPot \
						obj/MiyamotoNagaiPot obj/MultiPot



obj/OrbitIntegrator: src/OrbitIntegrator.cc src/OrbitIntegrator.h $(GalPot_h) $(AUXIL_H)
		     $(CPP) $(CFLAGS) $< ; $(AROTHER); $(AWAY)

obj/PJMCoords:	src/PJMCoords.cc src/PJMCoords.h src/Pi.h
	$(CPP) $(CFLAGS) $< ; $(AROTHER); $(AWAY)

obj/Random:	src/Random.cc src/Random.h $(AUXIL_H)
	$(CPP) $(CFLAGS) $< ; $(AROTHER); $(AWAY)

Other: obj/OrbitIntegrator obj/PJMCoords obj/Random


%.exe:	%.cc GalPot Other
	$(CPP) $(MFLAGS) -o $@ $<  -Lobj -lPot -lOther -lm

clean:
	rm -f obj/* *.exe

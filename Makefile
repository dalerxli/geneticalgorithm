CC          = icpc
CLINKER     = icpc

#CFLAGS      =   -Wall -O3 -march=pentium
CFLAGS      =   -Wall -O3 -xHost
#CFLAGS      = -i-fast -lm  
#LIBS        = -lm
DEPEND= makedepend

SRC        = GA.c ran_uniform.c system.c boxmuller.c readinput.c initialization.c write.c potential.c vectoroperation.c latticereduction.c rotation.c minimization.c force.c mnbrak.c brent.c overlap.c store.c
OBJS       = GA.o ran_uniform.o system.o boxmuller.o readinput.o initialization.o write.o potential.o vectoroperation.o latticereduction.o rotation.o minimization.o force.o mnbrak.o brent.o overlap.o store.o
EXECS      = GA

default: GA

all: $(EXECS)

GA:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o GA $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

GA.o: system.h ran_uniform.h 
ran_uniform.o: system.h ran_uniform.h
boxmuller.o: system.h ran_uniform.h
system.o: system.h
readinput.o: system.h ran_uniform.h
initialization.o: system.h ran_uniform.h
write.o: system.h
potential.o: system.h
vectoroperation.o: system.h
latticereduction.o: system.h
rotation.o: system.h
minimization.o: system.h
force.o: system.h
overlap.o: system.h
store.o: system.h

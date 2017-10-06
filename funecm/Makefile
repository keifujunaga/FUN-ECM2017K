# Makefile

PROGRAM = funecm
OBJS =  bsgs.o atkin.o point.o double_add.o normal_add.o Montgomery_ecm.o scalar.o Montgomery_main.o
CC = gcc 
LIBRARY = -lgmp -lrt

ifeq ($(ARCHITECTURE),Xeonphi) # for Xeonphi
	CFLAGS = -mmic -fopenmp -O2
	LIBDIR = -L/usr/local/lib
else # for Xeon
	CFLAGS = -openmp -O2
	LIBDIR = -L/home/project8/gmp/lib
endif

CFLAGS = -fopenmp -O2 -lm -O4 -Wall -Wunused-parameter -g 
LIBDIR = -L/usr/local/lib


$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $^ $(LIBRARY) $(LIBDIR)

Montgomery_main.o: Montgomery_main.c
	$(CC) $(CFLAGS) -c $<
point.o: point.c
	$(CC) $(CFLAGS) -c $<
double_add.o: double_add.c
	$(CC) $(CFLAGS) -c $<
normal_add.o: normal_add.c
	$(CC) $(CFLAGS) -c $<
Montgomery_ecm.o: Montgomery_ecm.c
	$(CC) $(CFLAGS) -c $<
scalar.o: scalar.c
	$(CC) $(CFLAGS) -c $<
atkin.o: atkin.c
	$(CC) $(CFLAGS) -c $<
bsgs.o: bsgs.c
	$(CC) $(CFLAGS) -c $<
all: clean $(PROGRAM)

clean:
	-rm -f $(PROGRAM) $(OBJS)


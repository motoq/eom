# CRD Sample code
#   common_c

#---------------------------------------------------------------------------

WRTSRC= echo_crd.c read_crd.c
WRTOBJ= echo_crd.o read_crd.o

echo_crd: $(WRTOBJ)
	gcc -g $(WRTOBJ) -lm -o echo_crd

#---------------------------------------------------------------------------

all:	read_crd.o  echo_crd

read_crd.o: read_crd.c
	gcc -g -c -o read_crd.o read_crd.c

clean:
	rm -f *.o


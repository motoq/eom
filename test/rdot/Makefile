CC = g++
CPPFLAGS = -g -Wall -I. -I../../include \
                    -I/usr/include/eigen3
LFLAGS = -L../../build -L../../../../../lib -leom -lgen -lsofa_c

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

rdot : $(OBJECTS)
	$(CC) $(CFLAGS) -o rdot $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm rdot $(OBJECTS)

CC = g++
CPPFLAGS = -g -Wall -I../../include \
                    -I/usr/include/eigen3
LFLAGS = -L../../build -L../../../../../lib -leom -lgen -lsofa_c
#                    -I../../../../../include/ \
#LFLAGS = -L../../../../../lib -L../../build -lsofa_c -leom

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

gauss_fg : $(OBJECTS)
	$(CC) $(CFLAGS) -o gauss_fg $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm gauss_fg $(OBJECTS)

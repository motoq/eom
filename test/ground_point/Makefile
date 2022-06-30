CC = g++
CPPFLAGS = -g -Wall -I../../include \
                    -I../../../../../include/ \
                    -I/usr/include/eigen3
LFLAGS = -L../../../../../lib -L../../build -leom
#LFLAGS = -L../../../../../lib -L../../build -lsofa_c -leom

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

ground_point : $(OBJECTS)
	$(CC) $(CFLAGS) -o ground_point $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm ground_point $(OBJECTS)

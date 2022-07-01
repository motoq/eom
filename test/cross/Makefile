CC = g++
CPPFLAGS = -g -Wall -I../../include \
                    -I/usr/include/eigen3
LFLAGS = -L../../build -leom

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

cross : $(OBJECTS)
	$(CC) $(CFLAGS) -o cross $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm cross $(OBJECTS)

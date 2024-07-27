CC = g++
CPPFLAGS = -g -Wall -I.                 \
                    -I../../include
LFLAGS = -L../../build -leom

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

sgp4_app : $(OBJECTS)
	$(CC) $(CFLAGS) -o sgp4_app $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm sgp4_app $(OBJECTS)

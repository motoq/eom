CC = g++
CPPFLAGS = -g -Wall -I../../include

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

fact : $(OBJECTS)
	$(CC) $(CFLAGS) -o fact $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm fact $(OBJECTS)

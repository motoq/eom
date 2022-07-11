CC = g++
CPPFLAGS = -g -Wall -I../../include
LFLAGS = -L../../build -leom

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

egm : $(OBJECTS)
	$(CC) $(CFLAGS) -o egm $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm egm $(OBJECTS)

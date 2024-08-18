CC = g++
CPPFLAGS = -g -Wall -I../../include \
                    -I/usr/include/eigen3
LFLAGS = -L../../build -L../../../../../lib -leom -lsofa_c
#                    -I../../../../../include/ \
#LFLAGS = -L../../../../../lib -L../../build -lsofa_c -leom

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

meme : $(OBJECTS)
	$(CC) $(CFLAGS) -o meme $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm meme $(OBJECTS)

CC = g++
CPPFLAGS = -g -Wall -I. -I../../include \
                    -I/usr/include/eigen3
LFLAGS = -L../../build -L../../../../../lib -leom -lgen -lsofa_c

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

rhat_dot : $(OBJECTS)
	$(CC) $(CFLAGS) -o rhat_dot $(OBJECTS) $(LFLAGS)

.PHONY : clean
clean :
	rm rhat_dot $(OBJECTS)

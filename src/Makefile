PROD =		fluid
LIB =		-lm -lGL -lGLU -lglut
LIBDIR =	-L/usr/lib -L/usr/local/lib -L/usr/X11R6/lib -L/usr/local/ -L/opt/X11/lib
INCDIR =	-I/usr/include -I/usr/local/include -I/usr/X11R6/include -I/usr/local/ -I/opt/X11/include
CC =		cc
CFLAGS = 	-Wall $(INCDIR) -O3
LD =		cc
LDFLAGS =	-Wall $(LIBDIR) $(LIB) -O3 # -Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib
OBJ =		main.o \
		fluid.o

$(PROD): $(OBJ)
	$(LD) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) -c $< $(CFLAGS) -o $@
clean:
	rm -f *.o $(PROD)
run:
	./$(PROD)
debug:
	gdb $(PROD)

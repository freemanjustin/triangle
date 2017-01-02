###################################################################
#
#  triangle art project
#
# freeman.justin@gmail.com
#
##################################################################


OBJ=	./src/main.o \
		./src/shader_utils.o \
		./src/tr.o \
		./src/triangle.o \
		./src/save_image.o \
		./external/glew/glew.o

# Compliler flags

# jillong
#INC=   -D_OS_X_ -I./include -I/Users/jfreeman/libpng-1.4.1/include

# lust with modules
INC=   -D_OS_X_ \
	-I./include \
	-I./external/glew \
	-I/apps/libpng/1.6.21/include

CFLAGS=	-O3 -g -Wall

CC=	clang $(CFLAGS) $(INC)

# Libraries

# jillong
#LFLAGS= -framework OpenGL -framework GLUT -L/Users/jfreeman/libpng-1.4.1/lib -lpng `/usr/bin/xml2-config --libs`

# lust with modules
LFLAGS= /apps/zlib/1.2.9/lib/libz.a /apps/libpng/1.6.21/lib/libpng.a -framework GLUT -framework OpenGL

# Executable

EXEC=	./bin/tri

$(EXEC):$(OBJ)
	$(CC) -o $(EXEC) $(OBJ) $(LFLAGS)

clean:
	rm $(OBJ)
	rm $(EXEC)

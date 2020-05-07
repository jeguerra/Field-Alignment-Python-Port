# Makefile for examples, Sai Ravela 2/28/13


TOP = ..
#FA2DC4V2LIB
PROJECT =  FA2DImNoHLIB
SRCDIR =  $(TOP)/examples
COREDIR = $(TOP)/core
CORELIB = $(COREDIR)/codegen/lib
COREPROJ = $(CORELIB)/$(PROJECT)
# get a list of all the C sources files
SRCS = $(wildcard $(COREPROJ)/*.c)
INCLUDES = $(COREPROJ)
CC=gcc
# Wall: show all warning
# -c: generate the object file
# -I: gcc will look in the directory for the include file.
CFLAGS=-c -Wall -I$(INCLUDES)
all: $(COREPROJ)/$(PROJECT).a
	
# AR:create a static library from object files previously compiled
# $@: filename representing the target
# $<: The filename of the first prerequisite
$(COREPROJ)/$(PROJECT).a: $(SRCS:.c=.o) 
	$(AR) rv $@ $(SRCS:.c=.o)
	ranlib $@
.c.o: 
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -rf $(COREPROJ)/*.o $(COREPROJ)/$(PROJECT).a


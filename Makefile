### Makefile for integral equation solver “ppie_solver4b”

# Set compiler and linker flags
INCLUDE = -I. -I/gscratch/riekesheabrown/local/include/
#IQUOTE = -iquote/gscratch/riekesheabrown/local/include/
CC = icc
CFLAGS = -Wall -ggdb
LIBS = -lgsl -lgslcblas -lm
LDFLAGS = -L/gscratch/riekesheabrown/local/lib/
#DEPS = functions_method7.h
LINKOPTS = -Wl,-rpath,/gscratch/riekesheabrown/local/lib/

#all: ppie_solver4b

ppie_solver4: ppie_solver4b.c functions_method7.h
	$(CC) $(INCLUDE) $(LDFLAGS) $(LINKOPTS) $(CFLAGS) ppie_solver4b.c -o ppie_solver4 $(LIBS) 

clean:
	\rm -f ppie_solver4


CC      =	gcc
CPP     = 	g++
CFLAGS  =	-O3 -Wall -g --std=c++11


NVCC    = 	nvcc
CUDIR   = 	/usr/local/cuda
NVFLAGS = 	-O3 -I$(CUDIR)/include -m64 -arch=compute_75 -code=sm_75 -Xptxas -v -rdc=true

TARGET	=	culsm

SRC_DIRS	?=	./src

SRCS 	:=	$(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.cu -or -name *.s)

OBJS	= 	$(SRCS:.c=.o)


INC		= 	-I/src
LFLAGS	= 
LIBS    =	-lm -lcudart

MKFILE	= 	Makefile 


$(TARGET) : $(OBJS)
		$(NVCC) $(NVFLAGS) $(INC) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
		$(NVCC) $(NVFLAGS) $(INC) -c $< -o $@


.PHONY: clean

clean:
		rm -f $(TARGET) $(patsubst %,$(SRC_DIRS)/%.o,$(basename $(SRCS)))

depend: $(SRCS)
		makedepend $(INC) $^

#

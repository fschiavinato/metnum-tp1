CC=gcc
CXX=g++ -std=c++11
RM=rm -f
CPPFLAGS=-g 
LDFLAGS=-g 
LDLIBS=-g

SRCS= main.cpp lu.cpp ppmloader/ppmloader.cpp OperacionesMatriciales.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: tp1 calibracion

tp1: $(OBJS)
	$(CXX) $(LDFLAGS) -o tp1 $(OBJS) $(LDLIBS)

lu: lu.o
	$(CXX) $(LDFLAGS) -o lu $< $(LDLIBS) 

calibracion: calibracion.o ppmloader.o
	$(CXX) $(LDFLAGS) -o calibracion calibracion.o ppmloader.o

calibracion.o: calibracion.cpp
	$(CXX) $(LDFLAGS) -c "calibracion.cpp" "ppmloader/ppmloader.cpp" "macros.h"

ppmloader.o: ppmloader/ppmloader.cpp
	$(CXX) $(LDFLAGS) -c "ppmloader/ppmloader.cpp"

optest: OperacionesMatricialesTest.o ./test/OperacionesMatricialesTest.h OperacionesMatriciales.cpp OperacionesMatriciales.h
	$(CXX) $(LDFLAGS) -o optest OperacionesMatricialesTest.o OperacionesMatriciales.o

OperacionesMatricialesTest.o: ./test/OperacionesMatricialesTest.cpp
	$(CXX) $(LDFLAGS) -c "./test/OperacionesMatricialesTest.cpp" "OperacionesMatriciales.cpp"

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)
	$(RM) *o calibracion
	$(RM) *o tp1
	$(RM) *o optest

distclean: clean
	$(RM) *~ .depend

include .depend

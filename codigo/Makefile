CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-g 
LDFLAGS=-g 
LDLIBS=-g

SRCS=calibracion.cpp main.cpp lu.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

calibracion: calibracion.o
	$(CXX) $(LDFLAGS) -o calibracion $< $(LDLIBS) 

lu: lu.o
	$(CXX) $(LDFLAGS) -o lu $< $(LDLIBS) 

all: $(OBJS)
	$(CXX) $(LDFLAGS) -o tool $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend

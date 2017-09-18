CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-g 
LDFLAGS=-g 
LDLIBS=-g

SRCS=calibracion.cpp main.cpp lu.cpp ppmloader/ppmloader.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: $(OBJS)
	$(CXX) $(LDFLAGS) -o run $(OBJS) $(LDLIBS) 

lu: lu.o
	$(CXX) $(LDFLAGS) -o lu $< $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend

CC=gcc
CXX=g++ -std=c++11
RM=rm -f
CPPFLAGS=-g 
LDFLAGS=-g 
LDLIBS=-g

SRCS= main.cpp lu.cpp ../../codigo/ppmloader/ppmloader.cpp OperacionesMatriciales.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: tp1

tp1: $(OBJS)
	$(CXX) $(LDFLAGS) -o tp1 $(OBJS) $(LDLIBS)

lu: lu.o
	$(CXX) $(LDFLAGS) -o lu $< $(LDLIBS) 

ppmloader.o: ././codigo/ppmloader/ppmloader.cpp
	$(CXX) $(LDFLAGS) -c "../../codigo/ppmloader/ppmloader.cpp"

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)
	$(RM) *o tp1

distclean: clean
	$(RM) *~ .depend

include .depend

#OBJS: main.o OperacionesMatriciales.o
#MAIN: tp1

#OperacionesMatriciales.o: ./codigo/OperacionesMatriciales.h ./codigo/OperacionesMatriciales.cpp ./codigo/test/OperacionesMatricialesTest.h ./codigo/test/OperacionesMatricialesTest.cpp
#main.o: ./codigo/main.cpp

#all: $(OBJS)
#    g++ $(OBJS) -o $@ $^

all: 
	g++ ./codigo/OperacionesMatriciales.h ./codigo/OperacionesMatriciales.cpp ./codigo/test/OperacionesMatricialesTest.h ./codigo/test/OperacionesMatricialesTest.cpp ./codigo/main.cpp -o tp1 -I ./codigo/ -I./codigo/test $^

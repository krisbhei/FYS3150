CC = c++
FLAGS = -c -fPIC -o
FLAGS1 = -shared -W1,-soname,
OBJ = jacobi.o divlib.o problem.o
TARGET = jacobi

all: $(TARGET)
	python project2.py

jacobi: $(OBJ)
	$(CC) $(FLAGS1) libjacobi.so -o libjacobi.so $(OBJ)

jacobi.o:jacobi.cpp
	$(CC) $(FLAGS) jacobi.o jacobi.cpp

divlib.o:divlib.cpp
	$(CC) -c divlib.cpp

problem.o:problem.cpp
	$(CC) -c problem.cpp

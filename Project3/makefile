CC = c++ -std=c++11 #-O3 #-Rpass-missed=loop-vectorize
OBJ = solarsim.o solarsys.o celestialbody.o differential.o vec3.o
TARGET = solarsim
FLAGS = -c -fPIC -o
FLAGS1 = -shared -W1,-soname,

all: $(TARGET)
		python project3.py

solarsim: $(OBJ)
					$(CC) $(FLAGS1) libsolarsim.so -o libsolarsim.so  $(OBJ)

solarsim.o: solarsim.cpp
						$(CC) $(FLAGS) solarsim.o solarsim.cpp

solarsys.o: solarsys.cpp
						$(CC) -c solarsys.cpp
celestialbody.o: celestialbody.cpp
					$(CC) -c celestialbody.cpp
differential.o: differential.cpp
						$(CC) -c differential.cpp
vec3.o: vec3.cpp
	$(CC) -c vec3.cpp

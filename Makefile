CC=g++
CFLAGS=-c

all: hello

hello: main.o Matr.o SobV.o FHT.o
	$(CC) main.o Matr.o SobV.o FHT.o -o Cvetkor -lsfml-graphics -lsfml-window -lsfml-system -lGL -lglut -lGLU

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

Matr.o: Matr.cpp
	$(CC) $(CFLAGS) Matr.cpp

FHT.o: FHT.cpp
	$(CC) $(CFLAGS) FHT.cpp

SobV.o: SobV.cpp
	$(CC) $(CFLAGS) SobV.cpp

clean:
	rm -rf *.o Cvetkor Test*.jpg Преобразован*.jpg

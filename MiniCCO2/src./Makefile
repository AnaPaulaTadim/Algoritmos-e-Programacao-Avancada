# Compilador
CXX = g++

# Flags
CXXFLAGS = -I ../include -std=c++11

# Arquivos fonte
SRC = $(wildcard *.cpp)

# Executavel (./minicco1 Nterm R gamma M)
TARGET = minicco1

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET) arvore.vtk raio.txt segmentos.txt

run:
	./$(TARGET) 100 10.0 3.0 5

view:
	python3 visualizador.py arvore.vtk

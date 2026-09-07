# Compilador
CXX = g++

# Flags
CXXFLAGS = -I ../include -std=c++11

# Arquivos fonte
SRC = $(wildcard *.cpp)

# Executavel (./minicco1 Nterm R gamma modo seed [M])
TARGET = minicco1

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET) arvore.vtk raio.txt segmentos.txt execucao.log

# modo=bin: arvore binaria classica (MiniCCO-1 original)
run:
	./$(TARGET) 100 10.0 3.0 bin 42

# modo=tricco: variante com multifurcacao (ate 3 filhos por no)
run-tricco:
	./$(TARGET) 100 10.0 3.0 tricco 42

view:
	python3 visualizador.py arvore.vtk raio

view-fluxo:
	python3 visualizador.py arvore.vtk fluxo

view-profundidade:
	python3 visualizador.py arvore.vtk profundidade


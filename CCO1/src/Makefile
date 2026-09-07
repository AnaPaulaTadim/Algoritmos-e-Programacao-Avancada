# Compilador
CXX = g++

# Flags
CXXFLAGS = -I ../include

# Bibliotecas
LIBS = -lGL -lGLU -lglut

# Arquivos fonte
SRC = $(wildcard *.cpp)

# Executável
TARGET = programa

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)

run:
	./$(TARGET) 10 10.0

view:
	python3 visualizar.py arvore.vtk
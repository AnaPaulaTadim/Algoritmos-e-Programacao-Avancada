#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include "point.hpp"
#include <string>

using namespace std;

//======================Estrutura do segmento(aresta)==================================
typedef struct{

    Point a; //ponto de inicio do segmento
    Point b; //ponto de fim do segmento

    //Indica se o segmento já foi bifurcado
    bool bifurcado;

    string nome;

} Segment;

#endif

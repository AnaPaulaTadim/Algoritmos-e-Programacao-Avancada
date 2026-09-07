#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include "point.hpp"

//======================Estrutura do segmento(aresta)==================================

struct Segment {

    int id;      //identificador unico do segmento/no
    int paiId;   //id do segmento pai (-1 = raiz da arvore)

    Point a; //ponto proximal (inicio do segmento)
    Point b; //ponto distal   (fim do segmento)

    //==================Atributos fisicos (MiniCCO-1)==================
    double raio;                //raio do vaso (segmento)
    double comprimento;         //comprimento do segmento
    double fluxo;                //vazao que passa pelo segmento
    double resistencia;          //resistencia hidraulica do segmento
    double volume;                //volume intravascular do segmento
    int    qtdTerminaisDistais;  //quantidade de terminais na sub-arvore distal a este segmento

};

#endif

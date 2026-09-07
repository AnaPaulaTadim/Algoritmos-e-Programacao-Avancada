#ifndef NODE_HPP
#define NODE_HPP

#include "point.hpp"

//==========================Estrutura de cada nó=======================================
typedef struct No {

    struct No* esquerda;//ponteiro para o filho da esquerda 

    struct No* direita;//ponteiro para o filho da direita

    struct No* pai;//ponteiro para o pai

    Point ponto;//cooordenada do nó(posição)
    
    bool terminal; //indica se o nó é um terminal

    int id;//identificados do nó

} No;

//==========================Ponteiro para o Nó==================================
typedef No* ptrNo;

#endif
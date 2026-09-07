#ifndef CANDIDATE_HPP
#define CANDIDATE_HPP

#include "point.hpp"
#include "segment.hpp"

//==========================Estrutura para armazenar os candidatos a conexão=======================================

typedef struct Candidato {

    Point novoPonto;

    Segment segmento;

    double custo;

    bool eValido;

} Candidato;

#endif
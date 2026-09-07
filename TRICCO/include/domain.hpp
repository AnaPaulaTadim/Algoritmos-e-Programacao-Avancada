#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "point.hpp"

//==========================Funcoes para manipulacao do dominio circular====================================

//verifica se um ponto esta dentro do dominio circular de centro "centro" e raio R
bool pontoDentroCirculo(Point p, Point centro, double R);
//gera um ponto aleatorio valido (dentro do dominio circular)
Point geraPontoValido(Point centro, double R);
//calcula a distancia minima exigida entre um novo terminal e os terminais ja existentes
double calculaDistanciaMinimaTerminal(double R, int nTerm);

#endif

#ifndef DOMAIN_HPP
#define DOMAIN_HPP


#include "point.hpp"


//==========================Funções para manipulação do domínio circular====================================
bool dentroDominio(Point p);
Point geraPontoValido(double R);
double calculaDistanciaMinimaTerminal(double R,int Nterm);

#endif
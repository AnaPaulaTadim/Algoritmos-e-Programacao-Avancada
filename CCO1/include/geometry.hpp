
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "point.hpp"
#include "segment.hpp"
#include "tree.hpp"

//==========================Funções Matemáticas====================================

//calcula a distancia entre dois pontos
double distancia(Point a, Point b);
//Direção de um ponto com realção a um segmento
double orientacao(Point a, Point b, Point c);
//função para verificar se os segmentos se intersectam
bool intersecaoSegmentos(Segment s1, Segment s2);
//distância de um ponto a um segmento específico
double distanciaPontoSegmento(Point p, Segment s);
//calcula o angulo entre um vetor existende e novo vetor
double calculaAngulo(Point a, Point b, Point novo);
//calcula a distência de um ponto a um segmento
double menorDistanciaSegmentos(Tree& arvore, Point p);
//calcula a distância de um ponto a toda a estrutura da árvore
double calculaDcrit(Point p, Segment seg);
//calcula a distência de um ponto a todos os terminais da árvore 
double menorDistanciaTerminais(Tree& arvore,Point p);

#endif
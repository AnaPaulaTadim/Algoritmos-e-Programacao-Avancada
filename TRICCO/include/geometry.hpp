#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "point.hpp"
#include "tree.hpp"

//==========================Funcoes Matematicas====================================

//---- Parte A (Comprimento, Resistencia e Volume): base para calculaComprimento ----
//calcula a distancia entre dois pontos -> usada em calculaComprimento(seg) = distancia(seg.a, seg.b)
double distancia(Point a, Point b);
//orientacao (produto vetorial) de tres pontos p, q, r
double orientacao(Point p, Point q, Point r);
//verifica se os segmentos p1-p2 e p3-p4 se intersectam (ignorando extremos coincidentes)
bool segmentosSeIntersectam(Point p1, Point p2, Point p3, Point p4);
//distancia de um ponto p ao segmento a-b
double distanciaPontoSegmento(Point p, Point a, Point b);

//---- Parte F (Otimizacao Geometrica da Bifurcacao): X dentro do triangulo A-B-C ----
//verifica se o ponto P esta dentro do triangulo A-B-C (X in triangle ABC)
bool pontoDentroTriangulo(Point A, Point B, Point C, Point P);
//ponto em coordenadas baricentricas: X = alpha*A + beta*B + lambda*C, com alpha+beta+lambda=1
Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda);

//==========================Distancias em relacao a arvore====================================

//menor distancia de um ponto a todos os terminais (pontos distais) ja presentes na arvore
double menorDistanciaTerminais(Tree& arvore, Point p);
//menor distancia de um ponto a qualquer segmento ja existente na arvore
double menorDistanciaSegmentos(Tree& arvore, Point p);
//menor distancia de um ponto a qualquer no (ponto proximal ou distal) da arvore,
//ignorando ate 3 ids especificos (protecao contra trifurcacao durante a insercao)
double menorDistanciaNos3(Tree& arvore, Point p, int idExcluir1, int idExcluir2, int idExcluir3);

#endif

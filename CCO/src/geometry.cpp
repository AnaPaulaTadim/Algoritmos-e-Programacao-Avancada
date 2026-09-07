#include "geometry.hpp"
#include "point.hpp"
#include "segment.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

//calcula a distancia entre dois pontos(Euclidiana)
double distancia(Point a, Point b){

    return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));

}
//========================================================================================
//*=================================================
//Direção de um ponto com realção a um segmento
//Positivo: curvas no sentido horário,
//Negativo: sentido anti horário; 0 colineares
//*=================================================
double orientacao(Point a, Point b, Point c){

    //Direção do ponto C em relação ao segmento AB
    return((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));

}
//=========================================================================================
//Função para verificar se os segmentos se intersectam
bool intersecaoSegmentos(Segment s1, Segment s2){

    //Checando a posição do segundo segmento(C->D) em relação ao priemiro(A->B)
    double orientacao1 = orientacao(s1.a, s1.b, s2.a);
    double orientacao2 = orientacao(s1.a,s1.b, s2.b);

    //Analisando primeiro segmento com relação ao segundo
    double orientacao3 = orientacao(s2.a, s2.b, s1.a); 
    double orientacao4 = orientacao(s2.a, s2.b, s1.b);

    //Verificando se os segmentos se cruzam
    if(orientacao1 * orientacao2 < 0 && orientacao3 * orientacao4 < 0){
        return true; //segmentos se cruzam
    }
    else{
        return false; //segementos não se cruzam
    }
}
//======================================================================================================
//Calculando distância de um ponto a um segmento 
double distanciaPontoSegmento(Point p, Segment s){

    //====================================
    // VETOR DO SEGMENTO
    //====================================

    double dx = s.b.x - s.a.x;
    double dy = s.b.y - s.a.y;

    //====================================
    // VETOR DO PONTO
    //====================================

    double px = p.x - s.a.x;
    double py = p.y - s.a.y;

    //====================================
    // PROJEÇÃO ESCALAR
    //====================================

    double dproj = (px*dx + py*dy) / (dx*dx + dy*dy);

    //====================================
    // PROJEÇÃO DENTRO DO SEGMENTO
    //====================================

    if(dproj >= 0.0 && dproj <= 1.0){

        double projx = s.a.x + dproj * dx;

        double projy = s.a.y + dproj * dy;

        return sqrt((p.x - projx)*(p.x - projx) + (p.y - projy)*(p.y - projy));
    }

    //====================================
    // PROJEÇÃO FORA DO SEGMENTO
    //====================================

    double distA = sqrt((p.x - s.a.x)*(p.x - s.a.x) + (p.y - s.a.y)*(p.y - s.a.y));

    double distB = sqrt((p.x - s.b.x)*(p.x - s.b.x) + (p.y - s.b.y)*(p.y - s.b.y));

    return min(distA, distB);
}

//==========================================================================================
//calcula o angulo entre um vetor existente (A->B) e novo vetor(B->novo)
double calculaAngulo(Point a, Point b, Point novo){

    //u*v = |u|*|v|*cos(theta)
    //|u| = raiz quadrada de (ux^2 + uy^2)->norma
    //|v| = raiz quadrada de (vx^2 + vy^2)->norma

    //Vetor do segmento existente(A->B)
    double ux = b.x - a.x;
    double uy = b.y - a.y;

    //Vetor do novo segmento(B->novo)
    double vx = novo.x - b.x;
    double vy = novo.y - b.y;

    //PRODUTO ESCALAR 
    
    //obtendo U*V + U*V dos dois segmentos para calcular o valor do angulo
    double produtoEscalar = ux * vx + uy * vy;

    //Normalização dos vetores
    double normaU = sqrt(ux * ux + uy * uy);
    double normaV = sqrt(vx * vx + vy * vy);

    //Se não tiver ângulo definido
    if(normaU == 0.0 || normaV == 0.0) {
        return 0.0;
    }

    //Calculando o cosseno
    double cosTheta = produtoEscalar / (normaU * normaV);

    //Evita erros numéricos
    cosTheta = max(-1.0, min(1.0, cosTheta));

    //Retorna o valor do cosseno em radianos
    return acos(cosTheta);

}
//==============================================================================================
//calcula a distância de ponto a toda a estrutura da árvore(segmentos)
double menorDistanciaSegmentos(Tree& arvore, Point p){
    
    double menor = 999999.0;
    
    for(Segment s : arvore.segmentos){
        
        double d = calculaDcrit(p, s); 
        
        if(d < menor){
            menor = d;
        }
    }

    return menor;
}
//======================================================================================
double calculaDcrit(Point p, Segment seg) {
    
    // Extremidade proximal (começo do segmento): A
    double xBj = seg.a.x;
    double yBj = seg.a.y;
    
    // Extremidade distal (fim do segmento): B
    double xj = seg.b.x;
    double yj = seg.b.y;

    //Vetor correto do segmento (Final - Inicial) -> A para B
    double dx = xj - xBj;
    double dy = yj - yBj;
    double lenSq = dx * dx + dy * dy;

    if (lenSq == 0) {
        //retorna a distância do ponto á um segmento(que é um ponto)
        return sqrt((p.x - xj) * (p.x - xj) + (p.y - yj) * (p.y - yj));
    }

    // Vetor do ponto em relação ao início do segmento (P - A)
    double px = p.x - xBj;
    double py = p.y - yBj;

    // Projeção escalar 
    double dproj = (px * dx + py * dy) / lenSq;

    double dcrit = 0.0;

    // Se 0 <= dproj <= 1, o ponto projeta-se perpendicularmente sobre o corpo do segmento
    if (dproj >= 0.0 && dproj <= 1.0) {
        // Distância Ortogonal
        double num = dx * py - dy * px;
        dcrit = abs(num) / sqrt(lenSq);
    } else {
        // Distância até as extremidades mais próximas
        double distToJ = sqrt((p.x - xj) * (p.x - xj) + (p.y - yj) * (p.y - yj));
        double distToBj = sqrt((p.x - xBj) * (p.x - xBj) + (p.y - yBj) * (p.y - yBj));
        dcrit = min(distToJ, distToBj);
    }

    return dcrit;
}
//==================================================================================
double menorDistanciaTerminais(Tree& arvore,Point p){

    double menor = 1e9;

    for(Segment s : arvore.segmentos){

        menor = min(menor,distancia(p,s.b));
    }

    return menor;
}
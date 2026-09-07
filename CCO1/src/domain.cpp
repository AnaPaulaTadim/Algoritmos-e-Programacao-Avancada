#include "domain.hpp"
#include "point.hpp"
#include "random.hpp"


#include <cmath>

//Verifica se o ponto está dentro do domínio circular
bool dentroDominio(Point p) {

    return (p.x * p.x + p.y * p.y) <= 1.0;

}
//Gera um ponto aleatório
double numeroAleatorio(double min, double max){

    double t = (double)rand() / (double)RAND_MAX;

    return min + t*(max-min);
}
//Gera um ponto aleatório válido
Point geraPontoValido(double R){

    Point p;

    do{
        p.x = numeroAleatorio(-R, R);

        p.y = numeroAleatorio(-R, R);

    }while(p.x*p.x + p.y*p.y > R*R);

    return p;
}
//Calcula a distência minima entre um novo terminal e os terminais já existentes
double calculaDistanciaMinimaTerminal(double R,int Nterm){

    double dArea = 0.45 * R * sqrt(M_PI / (double)Nterm);

    double dMin = 0.04 * R;

    double dMax = 0.28 * R;

    if(dArea < dMin)
        return dMin;

    if(dArea > dMax)
        return dMax;

    return dArea;
}




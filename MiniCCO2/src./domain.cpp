#include "domain.hpp"
#include "geometry.hpp"
#include "random.hpp"

#include <cmath>

//Verifica se o ponto esta dentro do dominio circular de centro "centro" e raio R
bool pontoDentroCirculo(Point p, Point centro, double R) {

    return distancia(p, centro) <= R + 1e-9;

}

//Gera um ponto aleatorio valido dentro do dominio circular (rejeicao)
Point geraPontoValido(Point centro, double R) {

    Point p;

    do {
        p.x = centro.x + randomDouble(-R, R);
        p.y = centro.y + randomDouble(-R, R);

    } while (distancia(p, centro) > R);

    return p;
}

//Calcula a distancia minima entre um novo terminal e os terminais ja existentes,
//em funcao da area media "disponivel" por terminal
double calculaDistanciaMinimaTerminal(double R, int nTerm) {

    if (nTerm <= 0) nTerm = 1;

    double dArea = 0.45 * R * sqrt(M_PI / (double) nTerm);

    double dMin = 0.04 * R;

    double dMax = 0.28 * R;

    if (dArea < dMin)
        return dMin;

    if (dArea > dMax)
        return dMax;

    return dArea;
}

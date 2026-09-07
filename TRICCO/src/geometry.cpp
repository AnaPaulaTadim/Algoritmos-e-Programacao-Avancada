#include "geometry.hpp"

#include <cmath>
#include <algorithm>
#include <cfloat>
using namespace std;

//calcula a distancia entre dois pontos
double distancia(Point a, Point b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    return sqrt(dx * dx + dy * dy);
}

//orientacao (produto vetorial) de tres pontos p, q, r
double orientacao(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    return val;
}

//auxiliar: verifica se q esta sobre o segmento p-r (ja sabendo que p, q, r sao colineares)
static bool noSegmento(Point p, Point q, Point r) {
    return (q.x <= max(p.x, r.x) + 1e-9 && q.x >= min(p.x, r.x) - 1e-9 &&
            q.y <= max(p.y, r.y) + 1e-9 && q.y >= min(p.y, r.y) - 1e-9);
}

//verifica se os segmentos p1-p2 e p3-p4 se intersectam, ignorando extremos coincidentes

bool segmentosSeIntersectam(Point p1, Point p2, Point p3, Point p4) {
    const double EPS = 1e-9;

    if (distancia(p1, p3) < EPS || distancia(p1, p4) < EPS ||
        distancia(p2, p3) < EPS || distancia(p2, p4) < EPS) {
        return false;
    }

    double o1 = orientacao(p1, p2, p3);
    double o2 = orientacao(p1, p2, p4);
    double o3 = orientacao(p3, p4, p1);
    double o4 = orientacao(p3, p4, p2);

    bool s1 = (fabs(o1) > EPS), s2 = (fabs(o2) > EPS);
    bool s3 = (fabs(o3) > EPS), s4 = (fabs(o4) > EPS);

    if ((o1 * o2 < 0) && (o3 * o4 < 0)) return true;

    if (!s1 && noSegmento(p1, p3, p2)) return true;
    if (!s2 && noSegmento(p1, p4, p2)) return true;
    if (!s3 && noSegmento(p3, p1, p4)) return true;
    if (!s4 && noSegmento(p3, p2, p4)) return true;

    return false;
}

//distancia de um ponto p ao segmento a-b (projecao ortogonal, com clamp nas extremidades)
double distanciaPontoSegmento(Point p, Point a, Point b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double lenSq = dx * dx + dy * dy;

    if (lenSq == 0.0) return distancia(p, a);

    double dproj = ((p.x - a.x) * dx + (p.y - a.y) * dy) / lenSq;

    if (dproj >= 0.0 && dproj <= 1.0) {
        double projx = a.x + dproj * dx;
        double projy = a.y + dproj * dy;
        return distancia(p, Point{projx, projy});
    }

    return min(distancia(p, a), distancia(p, b));
}

//verifica se o ponto P esta dentro do triangulo A-B-C
bool pontoDentroTriangulo(Point A, Point B, Point C, Point P) {
    double d1 = orientacao(P, A, B);
    double d2 = orientacao(P, B, C);
    double d3 = orientacao(P, C, A);

    bool temNeg = (d1 < -1e-9) || (d2 < -1e-9) || (d3 < -1e-9);
    bool temPos = (d1 > 1e-9) || (d2 > 1e-9) || (d3 > 1e-9);

    return !(temNeg && temPos);
}

//ponto em coordenadas baricentricas (alpha,beta,lambda) dentro do triangulo A-B-C:
//X = alpha*A + beta*B + lambda*C
Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda) {
    Point X;
    X.x = alpha * A.x + beta * B.x + lambda * C.x;
    X.y = alpha * A.y + beta * B.y + lambda * C.y;
    return X;
}

//menor distancia de um ponto a todos os TERMINAIS (pontos distais) ja presentes na arvore
double menorDistanciaTerminais(Tree& arvore, Point p) {
    double menor = DBL_MAX;
    for (const Segment& s : arvore.segmentos) menor = min(menor, distancia(p, s.b));
    return menor;
}

//menor distancia de um ponto a QUALQUER segmento ja existente na arvore
double menorDistanciaSegmentos(Tree& arvore, Point p) {
    double menor = DBL_MAX;
    for (const Segment& s : arvore.segmentos) {
        menor = min(menor, distanciaPontoSegmento(p, s.a, s.b));
    }
    return menor;
}

//menor distancia de um ponto a qualquer NO (ponto proximal ou distal) ja existente na
//arvore, ignorando ate 3 ids especificos - usada para impedir que um ponto de bifurcacao
//se aproxime de um ponto de bifurcação ja existente ou um segmento
double menorDistanciaNos3(Tree& arvore, Point p, int idExcluir1, int idExcluir2, int idExcluir3) {
    double menor = DBL_MAX;
    for (const Segment& s : arvore.segmentos) {
        if (s.id == idExcluir1 || s.id == idExcluir2 || s.id == idExcluir3) continue;
        menor = min(menor, distancia(p, s.a));
        menor = min(menor, distancia(p, s.b));
    }
    return menor;
}

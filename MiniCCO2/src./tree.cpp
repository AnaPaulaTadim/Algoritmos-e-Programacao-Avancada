#include "tree.hpp"
#include "geometry.hpp"
#include "domain.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <iomanip>
using namespace std;

//============================================================================================
//==========================Gerenciamento da arvore==========================================
//============================================================================================
//
// A arvore real deste projeto e o vetor arvore.segmentos : cada Segment guarda seu proprio id e o id do pai
// (paiId). Por isso os utilitarios abaixo trabalham por id, e nao por
// coincidencia de coordenadas.

//Retorna o indice, dentro de arvore.segmentos, do segmento com o id informado
int indicePorId(Tree& arvore, int id) {
    for (int i = 0; i < (int) arvore.segmentos.size(); i++) {
        if (arvore.segmentos[i].id == id) return i;
    }
    return -1;
}

//Monta o mapa id do pai -> indices dos filhos
void construirMapaDeFilhos(Tree& arvore, map<int, vector<int>>& filhosDe) {
    filhosDe.clear();
    for (int i = 0; i < (int) arvore.segmentos.size(); i++) {
        filhosDe[arvore.segmentos[i].paiId].push_back(i);
    }
}

//Conta o numero total de segmentos da arvore
int contaSegmentos(Tree& arvore) { return (int) arvore.segmentos.size(); }

//Conta quantos terminais (folhas, sem filhos) existem na arvore
int contaTerminais(Tree& arvore) {
    map<int, vector<int>> filhosDe;
    construirMapaDeFilhos(arvore, filhosDe);
    int n = 0;
    for (const Segment& s : arvore.segmentos) {
        if (filhosDe.find(s.id) == filhosDe.end() || filhosDe[s.id].empty()) n++;
    }
    return n;
}

//Calcula o comprimento total da arvore
double comprimentoTotalArvore(Tree& arvore) {
    double total = 0.0;
    for (const Segment& s : arvore.segmentos) total += s.comprimento;
    return total;
}

//============================================================================================
//==========================Insercao de bifurcacao============================================
//============================================================================================
//
// Antes:  pai -> old(A..B) -> (filhos de old, se houver)
// Depois: pai -> bif(A..X) -> { old(X..B) com seus filhos originais, term(X..C) }
//
// "bif" assume o lugar de "old" junto ao pai (mesmo paiId de old).
// "old" passa a comecar em X (em vez de A) e aponta para bif como pai.
// "term" e o novo ramo terminal, de X ate C.

void inserirBifurcacao(Tree& arvore, int idxOld, Point X, Point C, int* idBifOut, int* idTermOut) {
    Segment& old = arvore.segmentos[idxOld];

    Segment bif;
    bif.id = arvore.idCounter++;
    bif.paiId = old.paiId;
    bif.a = old.a;
    bif.b = X;
    bif.raio = 0; bif.comprimento = 0; bif.fluxo = 0; bif.resistencia = 0; bif.volume = 0;
    bif.qtdTerminaisDistais = 0;

    Segment term;
    term.id = arvore.idCounter++;
    term.paiId = bif.id;
    term.a = X;
    term.b = C;
    term.raio = 0; term.comprimento = 0; term.fluxo = 0; term.resistencia = 0; term.volume = 0;
    term.qtdTerminaisDistais = 0;

    old.a = X;
    old.paiId = bif.id;

    arvore.segmentos.push_back(bif);
    arvore.segmentos.push_back(term);

    if (idBifOut) *idBifOut = bif.id;
    if (idTermOut) *idTermOut = term.id;
}

//============================================================================================
//==========================Modelo Fisico do CCO (MiniCCO-1)=================================
//============================================================================================

//===================== Parte B - Contagem de Terminais Distais =====================
//Calcula recursivamente (pos-ordem, via pilha explicita) a quantidade de terminais
//distais de cada segmento. Como cada segmento tem no maximo 2 filhos 
// isso e sempre uma soma de 0, 1 ou 2 termos.
//Regra da Parte B: folha -> qtd_term_distal=1 ; senao -> soma dos filhos.
void atualizaQtdTerminaisDistais(Tree& arvore) {
    map<int, vector<int>> filhosDe;
    construirMapaDeFilhos(arvore, filhosDe);

    vector<int> ordem;
    vector<int> pilha;
    for (int i = 0; i < (int) arvore.segmentos.size(); i++) {
        if (arvore.segmentos[i].paiId == -1) pilha.push_back(i);
    }
    vector<bool> visitado(arvore.segmentos.size(), false);
    while (!pilha.empty()) {
        int idx = pilha.back(); pilha.pop_back();
        if (visitado[idx]) continue;
        visitado[idx] = true;
        ordem.push_back(idx);
        auto it = filhosDe.find(arvore.segmentos[idx].id);
        if (it != filhosDe.end()) {
            for (int filhoIdx : it->second) pilha.push_back(filhoIdx);
        }
    }

    for (int i = (int) ordem.size() - 1; i >= 0; i--) {
        int idx = ordem[i];
        auto it = filhosDe.find(arvore.segmentos[idx].id);
        if (it == filhosDe.end() || it->second.empty()) {
            arvore.segmentos[idx].qtdTerminaisDistais = 1;
        } else {
            int soma = 0;
            for (int filhoIdx : it->second) soma += arvore.segmentos[filhoIdx].qtdTerminaisDistais;
            arvore.segmentos[idx].qtdTerminaisDistais = soma;
        }
    }
}

//===================== Parte C - Fluxo em Cada Segmento =====================
//Qj = qtd_term_distal(j) * Qterm  (fluxos terminais iguais)
void atualizaFluxos(Tree& arvore, double Qterm) {
    for (Segment& s : arvore.segmentos) s.fluxo = s.qtdTerminaisDistais * Qterm;
}

//===================== Parte D - Lei de Bifurcacao e Escala dos Raios =====================
//rj = C * Qj^(1/gamma). Com C=1 (valor inicial pedido no enunciado) e o mesmo gamma para
//todo segmento
// Q_pai = Q_esq + Q_dir (conservacao de fluxo, Parte C).
void atualizaRaiosPorFluxo(Tree& arvore, double gamma, double C) {
    for (Segment& s : arvore.segmentos) s.raio = C * pow(s.fluxo, 1.0 / gamma);
}

//===================== Parte A - Comprimento, Resistencia e Volume =====================
//lj = distancia(a,b) = sqrt((x2-x1)^2+(y2-y1)^2)
void atualizaComprimentos(Tree& arvore) {
    for (Segment& s : arvore.segmentos) s.comprimento = distancia(s.a, s.b);
}

//Rj = 8*mu*lj / (pi*rj^4)  (Lei de Poiseuille)
void atualizaResistencias(Tree& arvore, double mu) {
    for (Segment& s : arvore.segmentos) {
        s.resistencia = (s.raio > 0.0) ? (8.0 * mu * s.comprimento) / (M_PI * pow(s.raio, 4)) : HUGE_VAL;
    }
}

//Vj = pi*rj^2*lj  (volume de um segmento cilindrico)
void atualizaVolumes(Tree& arvore) {
    for (Segment& s : arvore.segmentos) s.volume = M_PI * s.raio * s.raio * s.comprimento;
}

//===================== Parte D - funcao "orquestradora" =====================
//1) qtd. terminais distais (B) -> 2) fluxos (C) -> 3) raios (D) -> 4) comprimentos (A)
//-> 5) resistencias (A) -> 6) volumes (A)
void atualizaGeometriaFisica(Tree& arvore, double Qterm, double gamma, double mu) {
    atualizaQtdTerminaisDistais(arvore);   //Parte B
    atualizaFluxos(arvore, Qterm);         //Parte C
    atualizaRaiosPorFluxo(arvore, gamma);  //Parte D
    atualizaComprimentos(arvore);          //Parte A
    atualizaResistencias(arvore, mu);      //Parte A
    atualizaVolumes(arvore);               //Parte A
}

//===================== Parte A - Volume Total (tambem usado na Parte E) =====================
//Vtotal = soma(pi*rj^2*lj) para todos os N segmentos da arvore
double calculaVolumeTotal(Tree& arvore) {
    double total = 0.0;
    for (const Segment& s : arvore.segmentos) total += s.volume;
    return total;
}

//===================== Parte E - Funcao Custo: Volume Intravascular =====================
//J = Vtotal. E esta a funcao minimizada pela otimizacao geometrica das Partes F/G.
double funcaoCustoVolume(Tree& arvore) { return calculaVolumeTotal(arvore); }

//Raio medio de todos os segmentos da arvore
double raioMedio(Tree& arvore) {
    if (arvore.segmentos.empty()) return 0.0;
    double soma = 0.0;
    for (const Segment& s : arvore.segmentos) soma += s.raio;
    return soma / arvore.segmentos.size();
}

//============================================================================================
//==========================Normalizacao por resistencia global (extensao)===================
//============================================================================================
//
// C=1 em atualizaRaiosPorFluxo so fixa as PROPORCOES corretas entre os raios.


static double resistenciaEquivalenteRecursiva(Tree& arvore, int idx, map<int, vector<int>>& filhosDe) {
    Segment& atual = arvore.segmentos[idx];
    auto it = filhosDe.find(atual.id);

    if (it == filhosDe.end() || it->second.empty()) return atual.resistencia;

    double somaInversos = 0.0;
    bool algumValido = false;
    for (int filhoIdx : it->second) {
        double Rfilho = resistenciaEquivalenteRecursiva(arvore, filhoIdx, filhosDe);
        if (Rfilho > 0.0) { somaInversos += 1.0 / Rfilho; algumValido = true; }
    }
    double Rparalelo = algumValido ? (1.0 / somaInversos) : 0.0;
    return atual.resistencia + Rparalelo;
}

//Resistencia hidraulica equivalente de toda a arvore (serie no tronco, paralelo nas bifurcacoes)
double resistenciaEquivalenteArvore(Tree& arvore) {
    if (arvore.segmentos.empty()) return 0.0;
    map<int, vector<int>> filhosDe;
    construirMapaDeFilhos(arvore, filhosDe);
    int idxRaiz = -1;
    for (int i = 0; i < (int) arvore.segmentos.size(); i++) {
        if (arvore.segmentos[i].paiId == -1) { idxRaiz = i; break; }
    }
    if (idxRaiz < 0) return 0.0;
    return resistenciaEquivalenteRecursiva(arvore, idxRaiz, filhosDe);
}

//Multiplica o raio de todos os segmentos pelo mesmo fator k (preserva as proporcoes entre eles)
void escalaRaios(Tree& arvore, double k) {
    for (Segment& s : arvore.segmentos) s.raio *= k;
}

//Reescala os raios da arvore para que a resistencia equivalente bata com deltaP/Qperf,
//e recalcula resistencias e volumes com os raios finais
void normalizaRaiosPorResistenciaGlobal(Tree& arvore, double mu, double deltaP, double Qperf) {
    if (arvore.segmentos.empty() || Qperf <= 0.0) return;

    double Rtarget = deltaP / Qperf;
    double Rraw = resistenciaEquivalenteArvore(arvore);
    if (Rraw <= 0.0 || Rtarget <= 0.0) return;

    double k = pow(Rraw / Rtarget, 1.0 / 4.0);
    escalaRaios(arvore, k);
    atualizaResistencias(arvore, mu);
    atualizaVolumes(arvore);
}

//============================================================================================
//==========================Otimizacao geometrica da bifurcacao (Partes F e G)================
//============================================================================================
//
// Alem do criterio de rejeitar pontos degenerados (X muito perto de A, B ou C), tambem
// rejeita X muito perto de QUALQUER outro no ja existente na arvore. O FALLBACK (quando
// nenhum ponto da grade e valido) e o PONTO MEDIO de A-B (nunca B), o que garante que a
// bifurcacao inserida e sempre um ponto novo e distinto - nunca gera trifurcacao.

//===================== Parte G - Avaliacao de uma Bifurcacao Temporaria =====================
//Para cada candidato X: (1) altera temporariamente a bifurcacao [feito em otimizaBifurcacaoPorGrade,
//via copia da arvore], (2)-(4) atualiza comprimentos/fluxos/raios [atualizaGeometriaFisica],
//(5) calcula o volume total [funcaoCustoVolume, chamado pelo caller], (6) verifica intersecao
//geometrica e demais restricoes - e exatamente o que esta funcao faz:
bool configuracaoValida(Tree& arvoreTemp, int idxBif, int idxOld, int idxTerm,
                         Point domCentro, double domR, double distMinNos) {
    Segment& bif = arvoreTemp.segmentos[idxBif];
    Segment& old = arvoreTemp.segmentos[idxOld];
    Segment& term = arvoreTemp.segmentos[idxTerm];

    Point X = bif.b;
    if (!pontoDentroCirculo(X, domCentro, domR)) return false;

    //protecao extra contra trifurcacao: X nao pode coincidir (nem quase) com nenhum outro
    //no ja existente na arvore, fora do proprio triangulo local (bif, old e term)
    if (menorDistanciaNos3(arvoreTemp, X, bif.id, old.id, term.id) < distMinNos) return false;

    Point A = bif.a, B = old.b, C = term.b;

    for (const Segment& s : arvoreTemp.segmentos) {
        if (s.id == bif.id || s.id == old.id || s.id == term.id) continue;
        if (segmentosSeIntersectam(A, X, s.a, s.b)) return false;
        if (segmentosSeIntersectam(X, B, s.a, s.b)) return false;
        if (segmentosSeIntersectam(X, C, s.a, s.b)) return false;
    }
    return true;
}

//===================== Parte F - Otimizacao Geometrica da Bifurcacao =====================
//Busca em grade sobre coordenadas baricentricas (alpha,beta,lambda), com
//X = alpha*A + beta*B + lambda*C, alpha+beta+lambda=1 (Secao 10.1/10.2). Para cada X:
//guarda em (temp) a bifurcacao "temporaria" [Parte G, passo 1], recalcula a fisica inteira
//[Parte G, passos 2-4 = atualizaGeometriaFisica], mede o volume [Parte G, passo 5 =
//funcaoCustoVolume] e so aceita se configuracaoValida() [Parte G, passo 6] for verdadeiro.
//No final, retorna X* = argmin_{X na vizinhanca local do ponto medio de AB} Vtotal(X) entre
//os candidatos validos.
Point otimizaBifurcacaoPorGrade(Tree& arvore, int idxOld, Point A, Point B, Point C, int M,
                                 double Qterm, double gamma, double mu,
                                 Point domCentro, double domR, double distMinNos,
                                 double* melhorCusto, EstatisticasOtim* stats) {
    double custoMinimo = DBL_MAX;
    //Fallback: ponto medio de A-B (nunca B) - garante que, mesmo se a grade nao encontrar
    //nenhuma posicao valida, a bifurcacao inserida continua sendo um ponto novo e distinto
    Point melhorX = Point{ (A.x + B.x) / 2.0, (A.y + B.y) / 2.0 };
    bool encontrouValido = false;

    double escalaTriangulo = (distancia(A, B) + distancia(B, C) + distancia(A, C)) / 3.0;
    const double EPS_REL = 0.02;
    double epsMin = max(escalaTriangulo * EPS_REL, 1e-9);

    //Janela de busca LOCAL ao redor do ponto medio de A-B (t=0.5 na parametrizacao do
    //segmento AB). t_efetivo = beta/(alpha+beta) mede, para qualquer lambda, em que posicao
    //proporcional entre A (t=0) e B (t=1) o candidato X projeta ao longo de AB. Restringir
    //t_efetivo a [0.5-JANELA, 0.5+JANELA] impede que a otimizacao arraste X para perto de A
    //(ou de B), mantendo o carater LOCAL da otimizacao geometrica descrito na literatura.
    const double JANELA_LOCAL = 0.30;

    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= M - i; j++) {
            double alpha = i / (double) M;
            double beta = j / (double) M;
            double lambda = 1.0 - alpha - beta;

            double somaAB = alpha + beta;
            if (somaAB > 1e-9) {
                double tEfetivo = beta / somaAB;
                if (fabs(tEfetivo - 0.5) > JANELA_LOCAL) continue;
            }

            Point X = pontoBaricentrico(A, B, C, alpha, beta, lambda);

            if (!pontoDentroTriangulo(A, B, C, X)) continue;

            if (distancia(X, A) < epsMin ||
                distancia(X, B) < epsMin ||
                distancia(X, C) < epsMin) continue;

            // nunca altera a arvore real durante a busca em grade
            Tree temp;
            temp.segmentos = arvore.segmentos;
            temp.idCounter = arvore.idCounter;

            int idBif = -1, idTerm = -1;
            inserirBifurcacao(temp, idxOld, X, C, &idBif, &idTerm);

            if (stats) stats->conexoesTestadas++;

            int idxBifTemp = indicePorId(temp, idBif);
            int idxOldTemp = idxOld; //"old" continua no mesmo indice, so foi editado in-place
            int idxTermTemp = indicePorId(temp, idTerm);

            atualizaGeometriaFisica(temp, Qterm, gamma, mu);

            if (!configuracaoValida(temp, idxBifTemp, idxOldTemp, idxTermTemp,
                                     domCentro, domR, distMinNos)) {
                if (stats) stats->conexoesRejeitadas++;
                continue;
            }

            double custo = funcaoCustoVolume(temp);
            if (custo < custoMinimo) {
                custoMinimo = custo;
                melhorX = X;
                encontrouValido = true;
            }
        }
    }

    if (melhorCusto) *melhorCusto = encontrouValido ? custoMinimo : DBL_MAX;
    return melhorX;
}


//============================================================================================
//==========================Insercao de terminais=============================================
//============================================================================================
//
// Geracao de um novo terminal respeitando o CRITERIO DE DISTANCIA (distancia minima aos
// terminais existentes E aos segmentos existentes), aplicado ANTES da escolha do segmento
// candidato - o que impede que o novo ramo seja proposto colado em outra parte da arvore.

Point geraNovoTerminal(Tree& arvore, Point domCentro, double R, Point rootPoint) {
    int tentativas = 0;
    const int MAX_TENTATIVAS = 9000;

    double limiteTerminal = calculaDistanciaMinimaTerminal(R, contaTerminais(arvore) + 1);

    Point melhorPonto = geraPontoValido(domCentro, R);
    double melhorScore = -1.0;
    bool valido = false;

    while (tentativas < MAX_TENTATIVAS) {
        Point p = geraPontoValido(domCentro, R);

        //regra extra para os primeiros nos: afasta da entrada e do centro, evitando
        //um primeiro ramo degenerado
        if (arvore.segmentos.size() <= 1) {
            double dEntrada = distancia(p, rootPoint);
            double dCentro = distancia(p, domCentro);
            if (dEntrada < 0.40 * R || dCentro < 0.35 * R) { tentativas++; continue; }
        }

        double dTerminal = arvore.segmentos.empty() ? DBL_MAX : menorDistanciaTerminais(arvore, p);
        double dSegmento = arvore.segmentos.empty() ? DBL_MAX : menorDistanciaSegmentos(arvore, p);
        double dRaiz = distancia(p, rootPoint);

        if (dTerminal >= limiteTerminal && dSegmento >= 0.01 * R && dRaiz >= 0.01 * R) {

            double score;

            if (arvore.segmentos.empty()) {
                // Para o PRIMEIRO terminal não há nada na árvore ainda para medir
                // distância (dTerminal/dSegmento ficam ambos em DBL_MAX para
                // qualquer candidato válido) 
                // algoritmo aceitar o primeiro ponto que raspa o mínimo de
                // 0.40*R, em vez do mais distante. Por isso, só para este caso,
                // o score maximiza diretamente a distância até a raiz (mesmo
                // critério usado no MiniCCO-0 original), produzindo o tronco
                // longo esperado.
                score = dRaiz;
            } else {
                score = dTerminal + 0.35 * dSegmento;
            }

            if (score > melhorScore) {
                melhorScore = score;
                melhorPonto = p;
                valido = true;
            }
        }
        tentativas++;
    }

    return valido ? melhorPonto : melhorPonto;
}

//Insere um novo terminal na arvore, avaliando TODOS os segmentos como candidatos
bool insereTerminal(Tree& arvore, Point Pnovo, Point rootPoint, int modo, int M,
                     double Qterm, double gamma, double mu, Point domCentro, double domR,
                     double raioFixo, double distMinNos, EstatisticasOtim* statsAcum) {
    if (arvore.segmentos.empty()) {
        //primeira insercao: cria o segmento raiz diretamente (rootPoint -> Pnovo)
        Segment raizSeg;
        raizSeg.id = arvore.idCounter++;
        raizSeg.paiId = -1;
        raizSeg.a = rootPoint;
        raizSeg.b = Pnovo;
        raizSeg.qtdTerminaisDistais = 1;
        arvore.segmentos.push_back(raizSeg);

        if (modo == 0) {
            arvore.segmentos[0].raio = raioFixo;
            atualizaComprimentos(arvore);
            atualizaVolumes(arvore);
        } else {
            atualizaGeometriaFisica(arvore, Qterm, gamma, mu);
        }
        return true;
    }

    int n = (int) arvore.segmentos.size();
    double melhorCustoGlobal = DBL_MAX;
    int melhorIdx = -1;
    Point melhorX = Pnovo;

    for (int idx = 0; idx < n; idx++) {
        Point A = arvore.segmentos[idx].a;
        Point B = arvore.segmentos[idx].b;
        Point C = Pnovo;

        if (modo == 2) {
            double custo;
            EstatisticasOtim statsLocal;
            Point X = otimizaBifurcacaoPorGrade(arvore, idx, A, B, C, M,
                                                 Qterm, gamma, mu,
                                                 domCentro, domR, distMinNos,
                                                 &custo, &statsLocal);
            statsAcum->conexoesTestadas += statsLocal.conexoesTestadas;
            statsAcum->conexoesRejeitadas += statsLocal.conexoesRejeitadas;

            if (custo < melhorCustoGlobal) {
                melhorCustoGlobal = custo;
                melhorIdx = idx;
                melhorX = X;
            }
        } else {
            //modos 0/1 (e a primeira bifurcacao do modo 2): SEM otimizacao
            //geometrica -> bifurca sempre no PONTO MEDIO do segmento candidato
            //(nunca no ponto B, o que geraria trifurcacao)
            Point Xmedio = Point{ (A.x + B.x) / 2.0, (A.y + B.y) / 2.0 };

            Tree temp;
            temp.segmentos = arvore.segmentos;
            temp.idCounter = arvore.idCounter;

            int idBif = -1, idTerm = -1;
            inserirBifurcacao(temp, idx, Xmedio, C, &idBif, &idTerm);

            statsAcum->conexoesTestadas++;

            if (modo == 0) {
                int idxBifT = indicePorId(temp, idBif);
                int idxTermT = indicePorId(temp, idTerm);
                temp.segmentos[idxBifT].raio = raioFixo;
                temp.segmentos[idx].raio = raioFixo; //"old" (mesmo indice)
                temp.segmentos[idxTermT].raio = raioFixo;
                atualizaComprimentos(temp);
                atualizaVolumes(temp);
            } else {
                atualizaGeometriaFisica(temp, Qterm, gamma, mu);
            }

            bool valido = configuracaoValida(temp, indicePorId(temp, idBif), idx,
                                              indicePorId(temp, idTerm),
                                              domCentro, domR, distMinNos);

            double custo = valido ? funcaoCustoVolume(temp) : DBL_MAX;
            if (!valido) statsAcum->conexoesRejeitadas++;

            if (custo < melhorCustoGlobal) {
                melhorCustoGlobal = custo;
                melhorIdx = idx;
                melhorX = Xmedio;
            }
        }
    }

    if (melhorIdx < 0) return false; //nenhuma conexao valida encontrada

    //efetiva a melhor bifurcacao encontrada (permanente)
    int idBif = -1, idTerm = -1;
    inserirBifurcacao(arvore, melhorIdx, melhorX, Pnovo, &idBif, &idTerm);

    if (modo == 0) {
        int idxBif = indicePorId(arvore, idBif);
        int idxTerm = indicePorId(arvore, idTerm);
        arvore.segmentos[idxBif].raio = raioFixo;
        arvore.segmentos[melhorIdx].raio = raioFixo;
        arvore.segmentos[idxTerm].raio = raioFixo;
        atualizaComprimentos(arvore);
        atualizaVolumes(arvore);
    } else {
        atualizaGeometriaFisica(arvore, Qterm, gamma, mu);
    }

    return true;
}

//============================================================================================
//==========================Verificacao estrutural============================================
//============================================================================================
//
// Garante (em tempo de execucao) que nenhum id aparece como pai de mais de 2 segmentos -
// ou seja, que a arvore final e estritamente binaria (somente bifurcacoes, nunca trifurcacoes)

bool verificaSomenteBifurcacoes(Tree& arvore) {
    map<int, int> contagemFilhos;
    for (const Segment& s : arvore.segmentos) contagemFilhos[s.paiId]++;
    for (const auto& par : contagemFilhos) {
        if (par.first == -1) continue; //raiz: "pai" ficticio, nao conta
        if (par.second > 2) return false;
    }
    return true;
}

//============================================================================================
//==========================Exportacao para visualizacao======================================
//============================================================================================

void exportarVTK(Tree& arvore, string nomeArquivo, double R) {

    ofstream arquivo(nomeArquivo);

    if (!arquivo.is_open()) {
        cout << "Erro ao criar arquivo VTK" << endl;
        return;
    }

    //========================================
    // CABECALHO VTK
    //========================================
    arquivo << "# vtk DataFile Version 3.0\n";
    arquivo << "Arvore arterial\n";
    arquivo << "ASCII\n";
    arquivo << "DATASET POLYDATA\n";

    //========================================
    // PONTOS
    //========================================
    int totalPontos = arvore.segmentos.size() * 2;

    arquivo << "POINTS " << totalPontos << " float\n";

    for (Segment s : arvore.segmentos) {
        arquivo << s.a.x << " " << s.a.y << " 0\n";
        arquivo << s.b.x << " " << s.b.y << " 0\n";
    }

    //========================================
    // LINHAS
    //========================================
    int totalLinhas = arvore.segmentos.size();

    arquivo << "LINES " << totalLinhas << " " << totalLinhas * 3 << "\n";

    int indice = 0;
    for (size_t i = 0; i < arvore.segmentos.size(); i++) {
        arquivo << "2 " << indice << " " << indice + 1 << "\n";
        indice += 2;
    }

    //========================================
    // DADOS POR CELULA (RAIO DE CADA SEGMENTO)
    // Necessario para a visualizacao em tubos proporcionais ao raio
    //========================================
    arquivo << "CELL_DATA " << totalLinhas << "\n";
    arquivo << "SCALARS raio float 1\n";
    arquivo << "LOOKUP_TABLE default\n";

    for (Segment s : arvore.segmentos) {
        arquivo << s.raio << "\n";
    }

    arquivo.close();

    cout << "Arquivo VTK exportado com sucesso!" << endl;
}
//------------------------------------------------------------------------------------
//Exporta os segmentos como uma TABELA alinhada (texto simples), com os campos exigidos:
//id, pai, x0, y0, x1, y1, raio, comprimento, fluxo, resistencia, volume
void exportarSegmentosTXT(Tree& arvore, string nomeArquivo) {

    ofstream arquivo(nomeArquivo);

    if (!arquivo.is_open()) {
        cout << "Erro ao criar arquivo de segmentos" << endl;
        return;
    }

    //Larguras de coluna fixas, para as colunas ficarem alinhadas na leitura
    const int W_ID = 6, W_X = 12, W_RAIO = 13, W_LEN = 12, W_FLUXO = 13,
              W_RES = 14, W_VOL = 13;

    arquivo << left
            << setw(W_ID)    << "id"
            << setw(W_ID)    << "pai"
            << setw(W_X)     << "x0"
            << setw(W_X)     << "y0"
            << setw(W_X)     << "x1"
            << setw(W_X)     << "y1"
            << setw(W_RAIO)  << "raio"
            << setw(W_LEN)   << "comprimento"
            << setw(W_FLUXO) << "fluxo"
            << setw(W_RES)   << "resistencia"
            << setw(W_VOL)   << "volume"
            << "\n";

    //Linha separadora, do tamanho total do cabecalho
    int larguraTotal = W_ID*2 + W_X*4 + W_RAIO + W_LEN + W_FLUXO + W_RES + W_VOL;
    arquivo << string(larguraTotal, '-') << "\n";

    arquivo << fixed << setprecision(5);

    for (const Segment& s : arvore.segmentos) {
        arquivo << left
                << setw(W_ID)    << s.id
                << setw(W_ID)    << s.paiId
                << setw(W_X)     << s.a.x
                << setw(W_X)     << s.a.y
                << setw(W_X)     << s.b.x
                << setw(W_X)     << s.b.y
                << setw(W_RAIO)  << s.raio
                << setw(W_LEN)   << s.comprimento
                << setw(W_FLUXO) << scientific << s.fluxo
                << setw(W_RES)   << scientific << s.resistencia
                << setw(W_VOL)   << fixed << s.volume
                << fixed << "\n";
    }

    arquivo.close();
}

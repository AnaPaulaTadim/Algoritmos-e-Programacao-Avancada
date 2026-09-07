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
//==========================TRICCO: Multifurcacao==============================================
//============================================================================================
//
// Anexa o novo terminal C DIRETO no ponto distal do segmento idxPai, sem quebra-lo.
// Antes:  pai(A..B), com filhos existentes {f1, f2, ...}
// Depois: pai(A..B), com filhos {f1, f2, ..., term(B..C)}  -- pai passa a ter mais um filho.
// Isso e o que permite graus de saida > 2 (multifurcacao), sem alterar a estrutura do
// vetor de segmentos nem duplicar nenhuma logica de bifurcacao ja existente.
void inserirComoFilhoDireto(Tree& arvore, int idxPai, Point C, int* idTermOut) {
    Segment& pai = arvore.segmentos[idxPai];

    Segment term;
    term.id = arvore.idCounter++;
    term.paiId = pai.id;
    term.a = pai.b;
    term.b = C;
    term.raio = 0; term.comprimento = 0; term.fluxo = 0; term.resistencia = 0; term.volume = 0;
    term.qtdTerminaisDistais = 0;

    arvore.segmentos.push_back(term);

    if (idTermOut) *idTermOut = term.id;
}

//Validacao equivalente a configuracaoValida (Parte G), adaptada para o caso de multifurcacao:
//nao ha bifurcacao nova (ponto X), so o novo ramo term(pai.b -> C). Reaproveita as mesmas
//funcoes de dominio/distancia/interseccao ja usadas na validacao de bifurcacao.
bool configuracaoValidaMultifurcacao(Tree& arvoreTemp, int idxPai, int idxTerm,
                                      Point domCentro, double domR, double distMinNos) {
    Segment& pai = arvoreTemp.segmentos[idxPai];
    Segment& term = arvoreTemp.segmentos[idxTerm];

    Point A = pai.b;   //ponto de saida, ja existente (nao e um ponto novo)
    Point C = term.b;  //novo terminal

    if (!pontoDentroCirculo(C, domCentro, domR)) return false;

    //C nao pode coincidir (nem quase) com nenhum outro no ja existente na arvore
    if (menorDistanciaNos3(arvoreTemp, C, pai.id, term.id, term.id) < distMinNos) return false;

    //restricao fisica: diametro <= comprimento (2*r <= l) nos segmentos AFETADOS por esta
    //multifurcacao - "pai" ganha um filho a mais (seu raio pode mudar por redistribuicao de
    //fluxo, seu comprimento nao muda) e "term" e o segmento novo criado
    if (2.0 * pai.raio > pai.comprimento) return false;
    if (2.0 * term.raio > term.comprimento) return false;

    for (const Segment& s : arvoreTemp.segmentos) {
        if (s.id == pai.id || s.id == term.id) continue;
        if (segmentosSeIntersectam(A, C, s.a, s.b)) return false;
    }
    return true;
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

    //restricao fisica: diametro <= comprimento (2*r <= l) em todos os segmentos AFETADOS
    //por esta bifurcacao (bif, old e term sao os 3 segmentos criados/remodelados aqui)
    if (2.0 * bif.raio > bif.comprimento) return false;
    if (2.0 * old.raio > old.comprimento) return false;
    if (2.0 * term.raio > term.comprimento) return false;

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
                // o score maximiza diretamente a distância até a raiz 
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
                     double raioFixo, double distMinNos, EstatisticasOtim* statsAcum,
                     int maxFilhos) {
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

    //--------------------------------------------------------------------------------
    // TRICCO: bifurcacao e multifurcacao sao calculadas em PARALELO, cada uma buscando
    // seu proprio melhor candidato. So depois disso e feita a comparacao entre as duas
    // (ver bloco de decisao logo apos os dois loops).
    //--------------------------------------------------------------------------------
    double melhorCustoBifurcacao = DBL_MAX;
    int melhorIdxBifurcacao = -1;
    Point melhorXBifurcacao = Pnovo;

    double melhorCustoMultifurcacao = DBL_MAX;
    int melhorIdxMultifurcacao = -1;

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

            if (custo < melhorCustoBifurcacao) {
                melhorCustoBifurcacao = custo;
                melhorIdxBifurcacao = idx;
                melhorXBifurcacao = X;
            }
        } else {
            //modos 0/1: SEM otimizacao geometrica -> bifurca sempre no PONTO MEDIO
            //do segmento candidato (nunca no ponto B, o que geraria trifurcacao)
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

            if (custo < melhorCustoBifurcacao) {
                melhorCustoBifurcacao = custo;
                melhorIdxBifurcacao = idx;
                melhorXBifurcacao = Xmedio;
            }
        }
    }

    //========================================================================
    // TRICCO: busca, em paralelo, a MELHOR multifurcacao possivel - anexar Pnovo
    // DIRETO a algum segmento que ja e um ponto de ramificacao (>=1 filho) e ainda
    // tem espaco (< maxFilhos). So roda quando maxFilhos>2; com maxFilhos=2
    // (default) este bloco nunca executa e o MiniCCO-1 original fica intacto.
    //========================================================================
    if (maxFilhos > 2) {
        map<int, vector<int>> filhosDe;
        construirMapaDeFilhos(arvore, filhosDe);

        for (int idx = 0; idx < n; idx++) {
            int idSeg = arvore.segmentos[idx].id;
            int numFilhos = filhosDe.count(idSeg) ? (int) filhosDe[idSeg].size() : 0;

            //So aceita anexar direto em nos que JA SAO pontos de ramificacao (>=1 filho
            //existente) e ainda tem espaco (< maxFilhos). Anexar num segmento com ZERO
            //filhos (uma folha de verdade) criaria um no "de passagem" - grau 2, sem
            //ramificacao real, fisiologicamente sem sentido (um vaso so seguindo para um
            //unico descendente). Esse caso ja e coberto pela bifurcacao normal (loop acima).
            if (numFilhos == 0 || numFilhos >= maxFilhos) continue;

            Tree temp;
            temp.segmentos = arvore.segmentos;
            temp.idCounter = arvore.idCounter;

            int idTerm = -1;
            inserirComoFilhoDireto(temp, idx, Pnovo, &idTerm);

            statsAcum->conexoesTestadas++;

            if (modo == 0) {
                int idxTermT = indicePorId(temp, idTerm);
                temp.segmentos[idxTermT].raio = raioFixo;
                atualizaComprimentos(temp);
                atualizaVolumes(temp);
            } else {
                atualizaGeometriaFisica(temp, Qterm, gamma, mu);
            }

            bool valido = configuracaoValidaMultifurcacao(temp, idx, indicePorId(temp, idTerm),
                                                            domCentro, domR, distMinNos);

            double custo = valido ? funcaoCustoVolume(temp) : DBL_MAX;
            if (!valido) statsAcum->conexoesRejeitadas++;

            if (custo < melhorCustoMultifurcacao) {
                melhorCustoMultifurcacao = custo;
                melhorIdxMultifurcacao = idx;
            }
        }
    }

    if (melhorIdxBifurcacao < 0 && melhorIdxMultifurcacao < 0) return false; //nenhuma conexao valida

    //========================================================================
    // DECISAO: so aceita a multifurcacao se ela passar em DUAS regras, comparada
    // contra a melhor bifurcacao (calculada acima, em paralelo):
    //
    // 1) TOLERANCIA DE CUSTO: a multifurcacao pode custar, no maximo, 5% de volume
    //    a mais que a melhor bifurcacao (ligar direto a um no existente nunca e
    //    geometricamente otimizado, entao quase sempre custa um pouco mais).
    //
    // 2) ALCANCE LOCAL: o novo ramo da multifurcacao (pai.b -> Pnovo) nao pode ser
    //    mais que 1.6x mais longo que o ramo terminal que a MELHOR bifurcacao
    //    criaria (X_bifurcacao -> Pnovo) - evita vasos compridos/finos atravessando
    //    o dominio so pra "economizar" uma bifurcacao, mantendo a trifurcacao local.
    //
    // Se nenhuma multifurcacao valida existir, ou se ela nao passar nas duas regras,
    // cai de volta pra melhor bifurcacao
    //========================================================================
    const double TOLERANCIA_CUSTO = 1.05;   //ate 5% de volume a mais
    const double FATOR_ALCANCE_LOCAL = 1.6; //ate 1.6x o comprimento do ramo terminal da bifurcacao

    bool aceitarMultifurcacao = false;

    if (melhorIdxMultifurcacao >= 0 && melhorIdxBifurcacao >= 0) {
        bool dentroTolerancia = melhorCustoMultifurcacao <= TOLERANCIA_CUSTO * melhorCustoBifurcacao;

        double comprimentoRamoMultifurcacao = distancia(arvore.segmentos[melhorIdxMultifurcacao].b, Pnovo);
        double comprimentoRamoBifurcacao = distancia(melhorXBifurcacao, Pnovo);
        bool dentroAlcance = comprimentoRamoMultifurcacao <= FATOR_ALCANCE_LOCAL * comprimentoRamoBifurcacao;

        aceitarMultifurcacao = dentroTolerancia && dentroAlcance;
    } else if (melhorIdxMultifurcacao >= 0 && melhorIdxBifurcacao < 0) {
        //nenhuma bifurcacao valida foi encontrada, mas a multifurcacao foi - aceita direto
        aceitarMultifurcacao = true;
    }

    //bifurcacao OU multifurcacao (TRICCO)
    if (aceitarMultifurcacao) {
        int idTerm = -1;
        inserirComoFilhoDireto(arvore, melhorIdxMultifurcacao, Pnovo, &idTerm);

        if (modo == 0) {
            int idxTerm = indicePorId(arvore, idTerm);
            arvore.segmentos[idxTerm].raio = raioFixo;
            atualizaComprimentos(arvore);
            atualizaVolumes(arvore);
        } else {
            atualizaGeometriaFisica(arvore, Qterm, gamma, mu);
        }

        return true;
    }

    if (melhorIdxBifurcacao < 0) return false; //so aconteceria se so a multifurcacao tivesse sido testada e rejeitada

    int idBif = -1, idTerm = -1;
    inserirBifurcacao(arvore, melhorIdxBifurcacao, melhorXBifurcacao, Pnovo, &idBif, &idTerm);

    if (modo == 0) {
        int idxBif = indicePorId(arvore, idBif);
        int idxTerm = indicePorId(arvore, idTerm);
        arvore.segmentos[idxBif].raio = raioFixo;
        arvore.segmentos[melhorIdxBifurcacao].raio = raioFixo;
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

bool verificaSomenteBifurcacoes(Tree& arvore, int maxFilhos) {
    map<int, int> contagemFilhos;
    for (const Segment& s : arvore.segmentos) contagemFilhos[s.paiId]++;
    for (const auto& par : contagemFilhos) {
        if (par.first == -1) continue; //raiz: "pai" ficticio, nao conta
        if (par.second > maxFilhos) return false;
    }
    return true;
}

//============================================================================================
//==========================TRICCO: Metricas de grafo/estatisticas topologicas================
//============================================================================================
//
// Todas reaproveitam construirMapaDeFilhos (ja existente, usado em varias outras funcoes
// acima) - nenhuma logica de percurso da arvore e duplicada aqui, so o pos-processamento.

//Auxiliar: localiza o indice do segmento raiz (paiId == -1)
static int indiceRaiz(Tree& arvore) {
    for (int i = 0; i < (int) arvore.segmentos.size(); i++) {
        if (arvore.segmentos[i].paiId == -1) return i;
    }
    return -1;
}

//Profundidade (em saltos a partir da raiz) de CADA segmento, indexada por id. Reaproveitada
//por profundidadeMaxima, distribuicaoProfundidades e por exportarVTK (campo CELL_DATA).
map<int, int> profundidadePorSegmento(Tree& arvore) {
    map<int, int> prof;
    int idxRaiz = indiceRaiz(arvore);
    if (idxRaiz < 0) return prof;

    map<int, vector<int>> filhosDe;
    construirMapaDeFilhos(arvore, filhosDe);

    int idRaiz = arvore.segmentos[idxRaiz].id;
    vector<pair<int,int>> fila;
    fila.push_back({idRaiz, 1});
    size_t cursor = 0;

    while (cursor < fila.size()) {
        int idAtual = fila[cursor].first;
        int p = fila[cursor].second;
        cursor++;
        prof[idAtual] = p;

        auto it = filhosDe.find(idAtual);
        if (it != filhosDe.end()) {
            for (int filhoIdx : it->second) fila.push_back({arvore.segmentos[filhoIdx].id, p + 1});
        }
    }
    return prof;
}

//Profundidade maxima: BFS a partir da raiz contando saltos (arestas) ate o terminal mais distante
int profundidadeMaxima(Tree& arvore) {
    map<int, int> prof = profundidadePorSegmento(arvore);
    int maxProf = 0;
    for (const auto& par : prof) if (par.second > maxProf) maxProf = par.second;
    return maxProf;
}

//Numero de nos com exatamente 2 filhos (bifurcacoes classicas)
int numeroBifurcacoes(Tree& arvore) {
    map<int, int> contagem;
    for (const Segment& s : arvore.segmentos) contagem[s.paiId]++;
    int n = 0;
    for (const auto& par : contagem) {
        if (par.first == -1) continue;
        if (par.second == 2) n++;
    }
    return n;
}

//Numero de nos com 3 ou mais filhos
int numeroMultifurcacoes(Tree& arvore) {
    map<int, int> contagem;
    for (const Segment& s : arvore.segmentos) contagem[s.paiId]++;
    int n = 0;
    for (const auto& par : contagem) {
        if (par.first == -1) continue;
        if (par.second >= 3) n++;
    }
    return n;
}

//Grau medio do grafo da arvore: V = Nseg+1 pontos (raiz + ponto distal de cada segmento),
//E = Nseg arestas (cada segmento e uma aresta) -> grau medio = 2E/V (definicao padrao de
//teoria de grafos: soma dos graus = 2*arestas)
double grauMedio(Tree& arvore) {
    int nSeg = (int) arvore.segmentos.size();
    if (nSeg == 0) return 0.0;
    int nNos = nSeg + 1;
    return (2.0 * nSeg) / (double) nNos;
}

//Distribuicao de profundidades: para cada profundidade (saltos a partir da raiz), quantos
//NOS (pontos distais de segmento) existem naquela profundidade
map<int, int> distribuicaoProfundidades(Tree& arvore) {
    map<int, int> distribuicao;
    map<int, int> prof = profundidadePorSegmento(arvore);
    for (const auto& par : prof) distribuicao[par.second]++;
    return distribuicao;
}

//Comprimento (geometrico) da raiz ate cada terminal: soma de s.comprimento subindo por
//paiId ate chegar na raiz. Assume atualizaComprimentos ja foi chamado (sempre e, dentro
//de atualizaGeometriaFisica).
vector<double> caminhosRaizTerminal(Tree& arvore) {
    vector<double> caminhos;

    map<int, vector<int>> filhosDe;
    construirMapaDeFilhos(arvore, filhosDe);
    map<int, int> idParaIdx;
    for (int i = 0; i < (int) arvore.segmentos.size(); i++) idParaIdx[arvore.segmentos[i].id] = i;

    for (int i = 0; i < (int) arvore.segmentos.size(); i++) {
        int idSeg = arvore.segmentos[i].id;
        bool ehFolha = (filhosDe.find(idSeg) == filhosDe.end() || filhosDe[idSeg].empty());
        if (!ehFolha) continue;

        double soma = 0.0;
        int idAtual = idSeg;
        while (idAtual != -1) {
            int idx = idParaIdx[idAtual];
            soma += arvore.segmentos[idx].comprimento;
            idAtual = arvore.segmentos[idx].paiId;
        }
        caminhos.push_back(soma);
    }
    return caminhos;
}

//Ordem de Strahler (generalizada para multifurcacao): folha=1; no interno com filhos de
//ordens o1..ok -> ordem = max(oi), e +1 SE o valor maximo aparecer em 2 ou mais filhos
//(regra classica de Strahler, valida tanto para bifurcacao quanto multifurcacao).
map<int, int> ordemStrahler(Tree& arvore) {
    map<int, int> ordem;

    map<int, vector<int>> filhosDe;
    construirMapaDeFilhos(arvore, filhosDe);

    //ordem de visita pos-ordem via pilha (raiz -> folhas -> processa de tras pra frente)
    int idxRaiz = indiceRaiz(arvore);
    if (idxRaiz < 0) return ordem;

    vector<int> pilha;
    vector<int> ordemVisita;
    pilha.push_back(arvore.segmentos[idxRaiz].id);
    while (!pilha.empty()) {
        int idAtual = pilha.back(); pilha.pop_back();
        ordemVisita.push_back(idAtual);
        auto it = filhosDe.find(idAtual);
        if (it != filhosDe.end()) {
            for (int filhoIdx : it->second) pilha.push_back(arvore.segmentos[filhoIdx].id);
        }
    }

    for (int i = (int) ordemVisita.size() - 1; i >= 0; i--) {
        int idAtual = ordemVisita[i];
        auto it = filhosDe.find(idAtual);
        if (it == filhosDe.end() || it->second.empty()) {
            ordem[idAtual] = 1; //folha
            continue;
        }
        int maxOrdem = 0;
        int quantosNoMax = 0;
        for (int filhoIdx : it->second) {
            int idFilho = arvore.segmentos[filhoIdx].id;
            int ordemFilho = ordem[idFilho];
            if (ordemFilho > maxOrdem) { maxOrdem = ordemFilho; quantosNoMax = 1; }
            else if (ordemFilho == maxOrdem) { quantosNoMax++; }
        }
        ordem[idAtual] = (quantosNoMax >= 2) ? (maxOrdem + 1) : maxOrdem;
    }
    return ordem;
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

    
    // DADOS POR CELULA (RAIO, FLUXO E PROFUNDIDADE DE CADA SEGMENTO)
    
    arquivo << "CELL_DATA " << totalLinhas << "\n";

    arquivo << "SCALARS raio float 1\n";
    arquivo << "LOOKUP_TABLE default\n";
    for (Segment s : arvore.segmentos) {
        arquivo << s.raio << "\n";
    }

    arquivo << "SCALARS fluxo float 1\n";
    arquivo << "LOOKUP_TABLE default\n";
    for (Segment s : arvore.segmentos) {
        arquivo << s.fluxo << "\n";
    }

    map<int, int> profundidade = profundidadePorSegmento(arvore);
    arquivo << "SCALARS profundidade float 1\n";
    arquivo << "LOOKUP_TABLE default\n";
    for (Segment s : arvore.segmentos) {
        auto it = profundidade.find(s.id);
        arquivo << (it != profundidade.end() ? it->second : 0) << "\n";
    }

    arquivo.close();
}
//------------------------------------------------------------------------------------
//Exporta os segmentos como uma TABELA alinhada (texto simples) e
// acrescenta tabelas adicionais bem separadas com as metricas
//globais, distribuicao de profundidades, caminhos raiz-terminal e ordem de Strahler.
void exportarSegmentosTXT(Tree& arvore, string nomeArquivo, int Nterm, double tempoExecucao,
                           long long conexoesTestadas, long long conexoesRejeitadas,
                           string modoStr, unsigned int seed) {

    ofstream arquivo(nomeArquivo);

    if (!arquivo.is_open()) {
        cout << "Erro ao criar arquivo de segmentos" << endl;
        return;
    }

    //========================================================================
    // TABELA 1: SEGMENTOS
    //========================================================================
    arquivo << "===================================================================\n";
    arquivo << " TABELA 1 - SEGMENTOS DA ARVORE\n";
    arquivo << "===================================================================\n";

    //Larguras de coluna fixas, para as colunas ficarem alinhadas na leitura.
    //z0/z1 sao sempre 0.0 (dominio 2D), mas ficam explicitas no arquivo porque o
    //enunciado do TRICCO pede o formato 3D-compativel: id,pai,x0,y0,z0,x1,y1,z1,...
    const int W_ID = 6, W_X = 12, W_Z = 8, W_RAIO = 13, W_LEN = 12, W_FLUXO = 13,
              W_RES = 14, W_VOL = 13;

    arquivo << left
            << setw(W_ID)    << "id"
            << setw(W_ID)    << "pai"
            << setw(W_X)     << "x0"
            << setw(W_X)     << "y0"
            << setw(W_Z)     << "z0"
            << setw(W_X)     << "x1"
            << setw(W_X)     << "y1"
            << setw(W_Z)     << "z1"
            << setw(W_RAIO)  << "raio"
            << setw(W_LEN)   << "comprimento"
            << setw(W_FLUXO) << "fluxo"
            << setw(W_RES)   << "resistencia"
            << setw(W_VOL)   << "volume"
            << "\n";

    int larguraTotal = W_ID*2 + W_X*4 + W_Z*2 + W_RAIO + W_LEN + W_FLUXO + W_RES + W_VOL;
    arquivo << string(larguraTotal, '-') << "\n";

    arquivo << fixed << setprecision(5);

    for (const Segment& s : arvore.segmentos) {
        arquivo << left
                << setw(W_ID)    << s.id
                << setw(W_ID)    << s.paiId
                << setw(W_X)     << s.a.x
                << setw(W_X)     << s.a.y
                << setw(W_Z)     << 0.0
                << setw(W_X)     << s.b.x
                << setw(W_X)     << s.b.y
                << setw(W_Z)     << 0.0
                << setw(W_RAIO)  << s.raio
                << setw(W_LEN)   << s.comprimento
                << setw(W_FLUXO) << scientific << s.fluxo
                << setw(W_RES)   << scientific << s.resistencia
                << setw(W_VOL)   << fixed << s.volume
                << fixed << "\n";
    }

    double raioDaRaiz = 0.0;
    for (const Segment& s : arvore.segmentos) {
        if (s.paiId == -1) { raioDaRaiz = s.raio; break; }
    }
    int nSeg = contaSegmentos(arvore);
    int nNos = nSeg + 1;
    int nTerm = contaTerminais(arvore);
    double compTotal = comprimentoTotalArvore(arvore);
    double volTotal = funcaoCustoVolume(arvore);
    double raioMed = raioMedio(arvore);
    int profMax = profundidadeMaxima(arvore);
    int nBifurc = numeroBifurcacoes(arvore);
    int nMultif = numeroMultifurcacoes(arvore);
    double grauMed = grauMedio(arvore);

    //========================================================================
    // TABELA 2: METRICAS GLOBAIS 
    //========================================================================
    arquivo << "\n===================================================================\n";
    arquivo << " TABELA 2 - METRICAS GLOBAIS\n";
    arquivo << "===================================================================\n";
    arquivo << left << setw(28) << "metrica" << "valor" << "\n";
    arquivo << string(45, '-') << "\n";
    arquivo << left << setw(28) << "Nterm"                 << Nterm << "\n";
    arquivo << left << setw(28) << "Nseg"                  << nSeg << "\n";
    arquivo << left << setw(28) << "numero_nos"            << nNos << "\n";
    arquivo << left << setw(28) << "numero_terminais"      << nTerm << "\n";
    arquivo << left << setw(28) << "comprimento_total"     << compTotal << "\n";
    arquivo << left << setw(28) << "volume_total"          << scientific << volTotal << fixed << "\n";
    arquivo << left << setw(28) << "raio_raiz"              << scientific << raioDaRaiz << fixed << "\n";
    arquivo << left << setw(28) << "raio_medio"             << scientific << raioMed << fixed << "\n";
    arquivo << left << setw(28) << "profundidade_maxima"    << profMax << "\n";
    arquivo << left << setw(28) << "numero_bifurcacoes"     << nBifurc << "\n";
    arquivo << left << setw(28) << "numero_multifurcacoes"  << nMultif << "\n";
    arquivo << left << setw(28) << "grau_medio"             << grauMed << "\n";
    arquivo << left << setw(28) << "conexoes_testadas"      << conexoesTestadas << "\n";
    arquivo << left << setw(28) << "conexoes_rejeitadas"    << conexoesRejeitadas << "\n";
    arquivo << left << setw(28) << "tempo_execucao_s"       << tempoExecucao << "\n";
    arquivo << left << setw(28) << "modo"                   << modoStr << "\n";
    arquivo << left << setw(28) << "seed"                   << seed << "\n";

    //========================================================================
    // TABELA 3: DISTRIBUICAO DE PROFUNDIDADES
    //========================================================================
    arquivo << "\n===================================================================\n";
    arquivo << " TABELA 3 - DISTRIBUICAO DE PROFUNDIDADES (saltos a partir da raiz)\n";
    arquivo << "===================================================================\n";
    arquivo << left << setw(16) << "profundidade" << "quantidade_de_nos" << "\n";
    arquivo << string(35, '-') << "\n";
    map<int, int> distProf = distribuicaoProfundidades(arvore);
    for (const auto& par : distProf) {
        arquivo << left << setw(16) << par.first << par.second << "\n";
    }

    //========================================================================
    // TABELA 4: CAMINHOS RAIZ-TERMINAL (resumo estatistico)
    //========================================================================
    arquivo << "\n===================================================================\n";
    arquivo << " TABELA 4 - CAMINHOS RAIZ-TERMINAL (comprimento geometrico, resumo)\n";
    arquivo << "===================================================================\n";
    vector<double> caminhos = caminhosRaizTerminal(arvore);
    if (!caminhos.empty()) {
        double soma = 0.0, minC = caminhos[0], maxC = caminhos[0];
        for (double c : caminhos) { soma += c; if (c < minC) minC = c; if (c > maxC) maxC = c; }
        double media = soma / caminhos.size();
        double somaQuad = 0.0;
        for (double c : caminhos) somaQuad += (c - media) * (c - media);
        double desvio = (caminhos.size() > 1) ? sqrt(somaQuad / (caminhos.size() - 1)) : 0.0;

        arquivo << left << setw(28) << "quantidade_de_caminhos" << caminhos.size() << "\n";
        arquivo << left << setw(28) << "comprimento_minimo"     << minC << "\n";
        arquivo << left << setw(28) << "comprimento_maximo"     << maxC << "\n";
        arquivo << left << setw(28) << "comprimento_medio"      << media << "\n";
        arquivo << left << setw(28) << "desvio_padrao"          << desvio << "\n";
    } else {
        arquivo << "(arvore vazia)\n";
    }

    //========================================================================
    // TABELA 5: ORDEM DE STRAHLER
    //========================================================================
    arquivo << "\n===================================================================\n";
    arquivo << " TABELA 5 - ORDEM DE STRAHLER\n";
    arquivo << "===================================================================\n";
    map<int, int> strahler = ordemStrahler(arvore);
    int ordemDaRaiz = 0;
    for (const Segment& s : arvore.segmentos) {
        if (s.paiId == -1 && strahler.count(s.id)) { ordemDaRaiz = strahler[s.id]; break; }
    }
    map<int, int> contagemPorOrdem;
    for (const auto& par : strahler) contagemPorOrdem[par.second]++;

    arquivo << left << setw(28) << "ordem_strahler_da_raiz" << ordemDaRaiz << "\n";
    arquivo << "\n" << left << setw(16) << "ordem" << "quantidade_de_segmentos" << "\n";
    arquivo << string(40, '-') << "\n";
    for (const auto& par : contagemPorOrdem) {
        arquivo << left << setw(16) << par.first << par.second << "\n";
    }

    arquivo.close();
}

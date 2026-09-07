#ifndef TREE_HPP
#define TREE_HPP

#include "point.hpp"
#include "segment.hpp"

#include <string>
#include <vector>
#include <map>
using namespace std;

//=======================Estrutura da arvore======================================
//
// A arvore e representada como um vetor de Segmentos (arvore.segmentos), cada
// um identificado por id e apontando para o pai por paiId (-1 = raiz). 
struct Tree {

    vector<Segment> segmentos;

    int idCounter = 0;

};

//==========================Estatisticas de tentativas de conexao=====================================
struct EstatisticasOtim {

    int conexoesTestadas = 0;
    int conexoesRejeitadas = 0;

};

//==========================Gerenciamento da arvore=====================================

//Retorna o indice, dentro de arvore.segmentos, do segmento com o id informado (-1 se nao existir)
int indicePorId(Tree& arvore, int id);
//Monta o mapa id do pai -> indices dos filhos (reconstruido a cada chamada)
void construirMapaDeFilhos(Tree& arvore, map<int, vector<int>>& filhosDe);
//Conta o numero total de segmentos da arvore
int contaSegmentos(Tree& arvore);
//Conta quantos terminais (folhas) existem na arvore
int contaTerminais(Tree& arvore);
//Calcula o comprimento total da arvore
double comprimentoTotalArvore(Tree& arvore);

//Insere uma bifurcacao no segmento idxOld, no ponto X, ligando ao novo terminal C.
//pai -> old(A..B) vira pai -> bif(A..X) -> { old(X..B), term(X..C) }
void inserirBifurcacao(Tree& arvore, int idxOld, Point X, Point C, int* idBifOut, int* idTermOut);

//==========================TRICCO: Multifurcacao (variante do MiniCCO-1)==========================
//
// O MiniCCO-1 so permitia BIFURCACAO: um novo terminal sempre quebra um segmento existente
// em dois (inserirBifurcacao). A variante TRICCO generaliza isso permitindo tambem anexar o
// novo terminal DIRETAMENTE a um ponto (segmento) ja existente, sem quebra-lo - criando um
// no com 3+ filhos (multifurcacao). Nada na estrutura de dados precisou mudar: como os
// filhos de um segmento ja sao identificados via paiId (nao ha limite embutido no vetor), a
// unica coisa nova e ESTA funcao de insercao alternativa + o criterio para escolher entre
// bifurcar ou multifurcar (ver insereTerminal).

//Anexa um novo terminal C diretamente como filho extra do segmento idxPai (sem quebra-lo).
//pai(..) permanece igual; term(pai.b .. C) e adicionado como mais um filho de pai.
void inserirComoFilhoDireto(Tree& arvore, int idxPai, Point C, int* idTermOut);
//Verifica se uma multifurcacao temporaria (pai ganhando o filho term) e valida: dominio,
//distancia minima a outros nos e ausencia de interseccao com o resto da arvore.
bool configuracaoValidaMultifurcacao(Tree& arvoreTemp, int idxPai, int idxTerm,
                                      Point domCentro, double domR, double distMinNos);

//==========================Modelo Fisico do CCO (MiniCCO-1)=====================================

//---- Parte B (Contagem de Terminais Distais) ----
//Calcula recursivamente, em pos-ordem, a quantidade de terminais distais de cada segmento:
//folha -> qtd_term_distal = 1; caso contrario -> soma dos filhos
void atualizaQtdTerminaisDistais(Tree& arvore);
//---- Parte C (Fluxo em Cada Segmento) ----
//Atualiza o fluxo de cada segmento: Qj = qtd_term_distal(j) * Qterm
void atualizaFluxos(Tree& arvore, double Qterm);
//---- Parte D (Lei de Bifurcacao e Escala dos Raios) ----
//Atualiza os raios de cada segmento a partir do fluxo: rj = C * Qj^(1/gamma)
//(satisfaz automaticamente r_pai^gamma = r_esq^gamma + r_dir^gamma, pois Q_pai = Q_esq + Q_dir)
void atualizaRaiosPorFluxo(Tree& arvore, double gamma, double C = 1.0);
//---- Parte A (Comprimento, Resistencia e Volume) ----
//Atualiza o comprimento de cada segmento: lj = distancia(a,b)
void atualizaComprimentos(Tree& arvore);
//Atualiza a resistencia hidraulica de cada segmento (Lei de Poiseuille): Rj = 8*mu*lj/(pi*rj^4)
void atualizaResistencias(Tree& arvore, double mu);
//Atualiza o volume intravascular de cada segmento: Vj = pi*rj^2*lj
void atualizaVolumes(Tree& arvore);
//---- Parte D (funcao "orquestradora" pedida no enunciado) ----
//Executa a atualizacao fisica completa, EXATAMENTE na ordem pedida na Parte D:
//1) qtd. terminais distais -> 2) fluxos -> 3) raios -> 4) comprimentos -> 5) resistencias -> 6) volumes
void atualizaGeometriaFisica(Tree& arvore, double Qterm, double gamma, double mu);
//---- Parte A (Volume total, usado tambem na Parte E) ----
//Soma o volume de todos os segmentos da arvore: Vtotal = soma(pi*rj^2*lj)
double calculaVolumeTotal(Tree& arvore);
//---- Parte E (Funcao Custo: Volume Intravascular) ----
//Funcao custo do CCO, a ser minimizada pela otimizacao geometrica (Parte F/G): J = Vtotal
double funcaoCustoVolume(Tree& arvore);
//Raio medio de todos os segmentos da arvore (usado no relatorio final)
double raioMedio(Tree& arvore);

//==========================Normalizacao por resistencia global====================================

//Resistencia hidraulica equivalente de toda a arvore (serie no tronco, paralelo nas bifurcacoes)
double resistenciaEquivalenteArvore(Tree& arvore);
//Multiplica o raio de todos os segmentos pelo mesmo fator k (preserva as proporcoes entre eles)
void escalaRaios(Tree& arvore, double k);
//Reescala os raios da arvore para que a resistencia equivalente bata com deltaP/Qperf
void normalizaRaiosPorResistenciaGlobal(Tree& arvore, double mu, double deltaP, double Qperf);

//==========================Otimizacao geometrica da bifurcacao (Partes F e G)=========================

//---- Parte G (Avaliacao de uma Bifurcacao Temporaria) ----
//Verifica se uma configuracao de bifurcacao temporaria (bif/old/term) e valida: dentro do
//dominio, sem coincidir com outros nos da arvore e sem interseccao geometrica com outros
//segmentos. Chamada para CADA posicao candidata X testada na busca em grade (Parte F).
bool configuracaoValida(Tree& arvoreTemp, int idxBif, int idxOld, int idxTerm,
                         Point domCentro, double domR, double distMinNos);
//---- Parte F (Otimizacao Geometrica da Bifurcacao) ----
//Busca exaustiva em grade (coordenadas baricentricas alpha,beta,lambda; ver Parte F, Secao
//10.2), dentro do triangulo A-B-C (A,B = segmento antigo escolhido; C = novo terminal), pela
//posicao de bifurcacao X que resulta no MENOR volume total (Parte E) entre as configuracoes
//validas (Parte G). Equivale a X* = argmin_{X em ABC} Vtotal(X).
Point otimizaBifurcacaoPorGrade(Tree& arvore, int idxOld, Point A, Point B, Point C, int M,
                                 double Qterm, double gamma, double mu,
                                 Point domCentro, double domR, double distMinNos,
                                 double* melhorCusto, EstatisticasOtim* stats);

//==========================Insercao de terminais=====================================

//Gera um novo ponto terminal respeitando o criterio de distancia minima aos
//terminais e segmentos ja existentes
Point geraNovoTerminal(Tree& arvore, Point domCentro, double R, Point rootPoint);
//---- Parte H (Comparacao com o MiniCCO-0) ----
//Insere um novo terminal na arvore, avaliando TODOS os segmentos existentes como candidatos
//. O parametro "modo" seleciona qual das
//tres configuracoes da Parte H esta ativa, para permitir comparar as metricas entre elas:
//  modo 0 = raio FIXO em todos os segmentos            (arvore "sem escala dos raios")
//  modo 1 = raio pela Lei de Bifurcacao (Parte D), MAS SEM otimizar a posicao da bifurcacao
//           (sempre no ponto medio de A-B)              (arvore "com escala dos raios")
//  modo 2 = Lei de Bifurcacao (Parte D) + otimizacao geometrica completa (Partes F/G)
//                                                        (arvore "com escala + otimizacao")
//O main.cpp atual sempre roda em modo=2 (configuracao final pedida no enunciado); os modos
//0/1 continuam implementados e podem ser usados para reproduzir a comparacao da Parte H.
//
//maxFilhos: quantos filhos um segmento pode ter no maximo antes de ser obrigado a bifurcar
//(TRICCO). maxFilhos=2 preserva o comportamento binario original do MiniCCO-1 (default,
//compativel com chamadas antigas); maxFilhos=3 (ou mais) habilita multifurcacao: a cada novo
//terminal, alem de testar bifurcar cada segmento (como antes), tambem testa anexa-lo direto
//a cada segmento com menos de maxFilhos filhos, e escolhe a opcao de menor custo entre as duas.
bool insereTerminal(Tree& arvore, Point Pnovo, Point rootPoint, int modo, int M,
                     double Qterm, double gamma, double mu, Point domCentro, double domR,
                     double raioFixo, double distMinNos, EstatisticasOtim* statsAcum,
                     int maxFilhos = 2);

//==========================Verificacao estrutural=====================================

//Garante que nenhum segmento e pai de mais de maxFilhos segmentos.
//maxFilhos=2 (default) preserva o comportamento original (somente bifurcacoes).
bool verificaSomenteBifurcacoes(Tree& arvore, int maxFilhos = 2);

//==========================TRICCO: Metricas de grafo/estatisticas topologicas=====================
//
//Todas reaproveitam construirMapaDeFilhos (ja existente) - nenhum calculo duplicado.

//Profundidade maxima da arvore (numero de arestas da raiz ate o terminal mais distante em saltos)
int profundidadeMaxima(Tree& arvore);
//Profundidade de CADA segmento (id -> profundidade em saltos a partir da raiz). Reaproveitada
//por profundidadeMaxima, distribuicaoProfundidades e exportarVTK (colorir por profundidade).
map<int, int> profundidadePorSegmento(Tree& arvore);
//Numero de nos com exatamente 2 filhos (bifurcacoes classicas)
int numeroBifurcacoes(Tree& arvore);
//Numero de nos com 3 ou mais filhos (multifurcacoes - o que a variante TRICCO introduz)
int numeroMultifurcacoes(Tree& arvore);
//Grau medio do grafo da arvore: 2*Nseg/(Nseg+1) (V=Nseg+1 pontos, E=Nseg arestas)
double grauMedio(Tree& arvore);
//Distribuicao de profundidades: retorna, para cada profundidade (em saltos a partir da
//raiz), quantos NOS (pontos distais de segmento) existem naquela profundidade
map<int, int> distribuicaoProfundidades(Tree& arvore);
//Comprimento (geometrico, soma dos segmentos no caminho) da raiz ate cada terminal (folha)
vector<double> caminhosRaizTerminal(Tree& arvore);
//Ordem de Strahler de cada segmento (map id -> ordem). Generaliza para multifurcacao:
//folha=1; no interno com filhos de ordens o1..ok -> ordem = max(oi), +1 SE o maximo se repetir
//em 2 ou mais filhos.
map<int, int> ordemStrahler(Tree& arvore);

//==========================Exportacao para visualizacao=====================================

//gera um arquivo VTK para a visualizacao da arvore (raio de cada segmento como escalar por celula)
void exportarVTK(Tree& arvore, string nomeArquivo, double R);
//Exporta os segmentos como uma TABELA alinhada (texto simples), com os campos exigidos, e
//ALEM DISSO acrescenta, no mesmo arquivo, tabelas adicionais bem separadas com as metricas
//globais, distribuicao de profundidades, caminhos raiz-terminal e ordem de Strahler -
//consolidando num unico relatorio legivel (em vez de espalhar em varios arquivos).
void exportarSegmentosTXT(Tree& arvore, string nomeArquivo, int Nterm, double tempoExecucao,
                           long long conexoesTestadas, long long conexoesRejeitadas,
                           string modoStr, unsigned int seed);

#endif

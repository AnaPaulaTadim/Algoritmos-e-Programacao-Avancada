#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cfloat>

#include "tree.hpp"
#include "domain.hpp"
#include "geometry.hpp"

using namespace std;

//viscosidade sanguinea padrao (Pa*s), usada apenas para calcular resistencia/normalizacao
#define MU_DEFAULT 3.6e-3

//========================================================
// PARAMETROS FISICOS FIXOS  o enunciado
// pede so ./minicco1 Nterm R gamma M.
//========================================================
#define QPERF_DEFAULT 8.33e-6   //fluxo total de perfusao (m^3/s)
#define PPERF_DEFAULT 1.33e4    //pressao de perfusao (Pa)
#define PTERM_DEFAULT 7.98e3    //pressao terminal (Pa)

int main(int argc, char* argv[]) {

    //========================================================
    // VERIFICACAO DOS ARGUMENTOS DA LINHA DE COMANDO
    //========================================================
    if (argc != 5) {
        cout << "Uso: ./minicco1 Nterm R gamma M" << endl;
        cout << "  Nterm = numero de terminais" << endl;
        cout << "  R     = raio do dominio circular" << endl;
        cout << "  gamma = expoente da lei de bifurcacao (r = C * Q^(1/gamma))" << endl;
        cout << "  M     = resolucao da busca em grade no triangulo (Parte J)" << endl;
        return 1;
    }

    //========================================================
    // LEITURA DOS ARGUMENTOS
    //========================================================
    int Nterm = atoi(argv[1]);
    double R = atof(argv[2]);
    double gamma = atof(argv[3]);
    int M = atoi(argv[4]);

    //Parametros fisicos fixos (nao vem mais da linha de comando)
    double Qperf = QPERF_DEFAULT;
    double deltaP = PPERF_DEFAULT - PTERM_DEFAULT;
    int modo = 2; 

    if (Nterm <= 0 || R <= 0.0 || gamma <= 0.0 || M <= 0) {
        cout << "Parametros invalidos." << endl;
        return 1;
    }

    double mu = MU_DEFAULT;
    double Qterm = Qperf / (double) Nterm;
    double raioFixo = 0.02 * R; //nao usado (modo 0 nao e mais acessivel via CLI, mantido so internamente)

    //========================================================
    // INICIALIZA O GERADOR ALEATORIO
    //========================================================
    srand((unsigned) time(NULL));

    //========================================================
    // CRIA A ESTRUTURA DA ARVORE
    //========================================================
    Tree arvore;
    EstatisticasOtim statsAcum;

    Point domCentro = {0.0, 0.0};
    Point rootPoint = {0.0, R};

    //distancia minima entre nos (protecao anti-trifurcacao na busca em grade) 
   
    double distMinNos = 0.15 * calculaDistanciaMinimaTerminal(R, Nterm);

    clock_t inicio = clock();

    //========================================================
    // GERACAO DA ARVORE (CRITERIO GULOSO DO ARTIGO)
    //========================================================
    int terminaisInseridos = 0;
    int tentativasTotais = 0;

    //Sem limite global de tentativas: segue tentando pontos novos até conseguir
    //inserir os Nterm terminais,ele nunca "desiste" só continua testando até achar um ponto válido segundo
    //os critérios.
    while (terminaisInseridos < Nterm) {

        Point Pnovo = geraNovoTerminal(arvore, domCentro, R, rootPoint);

        bool ok = insereTerminal(arvore, Pnovo, rootPoint, modo, M,
                                  Qterm, gamma, mu, domCentro, R,
                                  raioFixo, distMinNos, &statsAcum);
        tentativasTotais++;

        if (ok) {
            terminaisInseridos++;
            cout << "--------------------------------------------------" << endl;
            cout << " Terminal " << terminaisInseridos << "/" << Nterm
                 << " inserido. Coordenadas: (" << Pnovo.x << ", " << Pnovo.y << ")" << endl;
        }
    }

    //========================================================
    // NORMALIZACAO DOS RAIOS PELA RESISTENCIA HIDRAULICA GLOBAL
    //========================================================
    if (modo != 0 && !arvore.segmentos.empty()) {
        normalizaRaiosPorResistenciaGlobal(arvore, mu, deltaP, Qperf);
    }

    clock_t fim = clock();
    double tempoExecucao = ((double) (fim - inicio)) / CLOCKS_PER_SEC;

    //========================================================
    // RELATORIO FINAL
    //========================================================
    cout << endl;
    cout << "==================================================" << endl;
    cout << "               RELATORIO DA ARVORE                " << endl;
    cout << "==================================================" << endl;
    //Localiza o segmento que toca a raiz (paiId == -1) para reportar o raio
    //correto - o indice 0 do vetor NAO e garantia de ser esse segmento, pois
    //a ordem interna muda a cada insercao 
    double raioDaRaiz = 0.0;
    for (const Segment& s : arvore.segmentos) {
        if (s.paiId == -1) { raioDaRaiz = s.raio; break; }
    }

    cout << "Numero total de nos:       " << contaSegmentos(arvore) + 1 << endl;
    cout << "Numero total de segmentos: " << contaSegmentos(arvore) << endl;
    cout << "Numero de terminais:       " << contaTerminais(arvore) << endl;
    cout << "Comprimento total:         " << comprimentoTotalArvore(arvore) << endl;
    cout << "Volume total (J):          " << funcaoCustoVolume(arvore) << endl;
    cout << "Raio da raiz:              " << raioDaRaiz << endl;
    cout << "Raio medio dos segmentos:  " << raioMedio(arvore) << endl;
    cout << "Conexoes testadas:         " << statsAcum.conexoesTestadas << endl;
    cout << "Conexoes rejeitadas:       " << statsAcum.conexoesRejeitadas << endl;
    cout << "Tempo de execucao:         " << tempoExecucao << " s" << endl;
    cout << "==================================================" << endl;

    if (!verificaSomenteBifurcacoes(arvore)) {
        cout << endl << "[ALERTA] Inconsistencia estrutural detectada: algum no possui "
             << "mais de 2 filhos (trifurcacao)." << endl;
    }

    //========================================================
    // EXPORTACAO PARA VISUALIZACAO 
    //========================================================
    exportarVTK(arvore, "arvore.vtk", R);
    exportarSegmentosTXT(arvore, "segmentos.txt");

    ofstream raioFile("raio.txt");
    raioFile << R;
    raioFile.close();

    return 0;
}

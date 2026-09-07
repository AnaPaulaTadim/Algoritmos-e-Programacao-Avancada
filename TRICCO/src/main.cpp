#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <string>
#include <algorithm>

#include "tree.hpp"
#include "domain.hpp"
#include "geometry.hpp"

using namespace std;

//viscosidade sanguinea padrao (Pa*s), usada apenas para calcular resistencia/normalizacao
#define MU_DEFAULT 3.6e-3

//========================================================
// PARAMETROS FISICOS FIXOS (valores sugeridos no enunciado do MiniCCO-1/TRICCO).
// Continuam fixos internamente - o enunciado do TRICCO define a entrada como
// "Nterm R gamma modo seed", sem espaco para eles na linha de comando.
//========================================================
#define QPERF_DEFAULT 8.33e-6   //fluxo total de perfusao (m^3/s)
#define PPERF_DEFAULT 1.33e4    //pressao de perfusao (Pa)
#define PTERM_DEFAULT 7.98e3    //pressao terminal (Pa)

//Resolucao da grade da otimizacao geometrica (Parte F/G do MiniCCO-1), usada quando o
//parametro opcional M nao e informado na linha de comando
#define M_DEFAULT 6

static void imprimirUso(const char* prog) {
    cout << "Uso: " << prog << " Nterm R gamma modo seed [M]" << endl;
    cout << "  Nterm = numero de terminais" << endl;
    cout << "  R     = raio do dominio circular" << endl;
    cout << "  gamma = expoente da lei de bifurcacao (r = C * Q^(1/gamma))" << endl;
    cout << "  modo  = variante da arvore:" << endl;
    cout << "            bin    -> binaria classica (MiniCCO-1 original, maximo 2 filhos por no)" << endl;
    cout << "            tricco -> multifurcacao (TRICCO, ate 3 filhos por no)" << endl;
    cout << "  seed  = semente do gerador aleatorio (inteiro)" << endl;
    cout << "  M     = (opcional) resolucao da busca em grade da bifurcacao. Default = "
         << M_DEFAULT << endl;
    cout << "Exemplo: " << prog << " 100 0.05 3.0 tricco 123" << endl;
}

int main(int argc, char* argv[]) {

    //========================================================
    // VERIFICACAO E LEITURA DOS ARGUMENTOS DA LINHA DE COMANDO
    //========================================================
    // Formato pedido no enunciado do TRICCO: ./trabalho Nterm R gamma modo seed
    // O parametro M (resolucao da busca geometrica), que so existia no MiniCCO-1, foi
    // incorporado como um 6o argumento OPCIONAL, para nao quebrar a especificacao de 5
    // argumentos e ao mesmo tempo preservar a possibilidade de ajusta-lo.
    if (argc != 6 && argc != 7) {
        imprimirUso(argv[0]);
        return 1;
    }

    int Nterm = atoi(argv[1]);
    double R = atof(argv[2]);
    double gamma = atof(argv[3]);
    string modoStr = argv[4];
    unsigned int seed = (unsigned int) strtoul(argv[5], nullptr, 10);
    int M = (argc == 7) ? atoi(argv[6]) : M_DEFAULT;

    if (Nterm <= 0 || R <= 0.0 || gamma <= 0.0 || M <= 0) {
        cout << "Parametros invalidos." << endl;
        imprimirUso(argv[0]);
        return 1;
    }

    //--------------------------------------------------------------------------
    // "modo" (string, pedido pelo enunciado do TRICCO) seleciona a VARIANTE da
    // arvore: binaria classica (maxFilhos=2, comportamento 100% original do
    // MiniCCO-1) ou multifurcacao TRICCO (maxFilhos=3). Isto e ortogonal ao
    // parametro interno "modoGeometrico" ja existente no MiniCCO-1 (0=raio
    // fixo, 1=Murray sem otimizacao, 2=Murray+otimizacao geometrica) - esse
    // continua fixo em 2 (a configuracao mais completa), como ja era antes.
    //--------------------------------------------------------------------------
    string modoNormalizado = modoStr;
    transform(modoNormalizado.begin(), modoNormalizado.end(), modoNormalizado.begin(), ::tolower);

    int maxFilhos = 2;
    string nomeVariante = "bin (binaria classica, MiniCCO-1 original)";
    if (modoNormalizado == "tricco" || modoNormalizado == "tri" || modoNormalizado == "multi") {
        maxFilhos = 3;
        nomeVariante = "tricco (multifurcacao, ate 3 filhos por no)";
    } else if (modoNormalizado == "bin" || modoNormalizado == "cco" || modoNormalizado == "bifurc") {
        maxFilhos = 2;
        nomeVariante = "bin (binaria classica, MiniCCO-1 original)";
    } else {
        cout << "Modo desconhecido: '" << modoStr << "'. Use 'bin' ou 'tricco'." << endl;
        imprimirUso(argv[0]);
        return 1;
    }

    const int modoGeometrico = 2; //Murray + otimizacao geometrica completa (Partes D/F/G)

    double Qperf = QPERF_DEFAULT;
    double deltaP = PPERF_DEFAULT - PTERM_DEFAULT;
    double mu = MU_DEFAULT;
    double Qterm = Qperf / (double) Nterm;
    double raioFixo = 0.02 * R; //nao usado no modoGeometrico atual (so existe para o modo 0, mantido por compatibilidade)

    //========================================================
    // ARQUIVO DE LOG (substitui os prints de depuracao no terminal)
    //========================================================
    ofstream log("execucao.log");
    log << "===================================================" << "\n";
    log << " EXECUCAO CCO-TRICCO" << "\n";
    log << "===================================================" << "\n";
    log << "Parametros utilizados:" << "\n";
    log << "  Nterm         = " << Nterm << "\n";
    log << "  R             = " << R << "\n";
    log << "  gamma         = " << gamma << "\n";
    log << "  modo          = " << modoStr << "  (" << nomeVariante << ")" << "\n";
    log << "  seed          = " << seed << "\n";
    log << "  M (grade)     = " << M << "\n";
    log << "  Qperf         = " << Qperf << "\n";
    log << "  pperf         = " << PPERF_DEFAULT << "\n";
    log << "  pterm         = " << PTERM_DEFAULT << "\n";
    log << "  deltaP        = " << deltaP << "\n";
    log << "  mu            = " << mu << "\n";
    log << "  Qterm         = " << Qterm << "\n";
    log << "---------------------------------------------------" << "\n";
    log << "Inicio da execucao." << "\n";
    log << "---------------------------------------------------" << "\n";

    //========================================================
    // INICIALIZA O GERADOR ALEATORIO COM A SEMENTE INFORMADA
    //========================================================
    // (antes usava srand(time(NULL)), o que impedia reproduzir experimentos;
    // o enunciado do TRICCO pede a semente como parametro justamente para isso -
    // os "Experimentos Minimos" do trabalho final pedem rodar cada configuracao
    // com 3 sementes diferentes e comparar media/desvio padrao.)
    srand(seed);

    //========================================================
    // CRIA A ESTRUTURA DA ARVORE
    //========================================================
    Tree arvore;
    EstatisticasOtim statsAcum;

    Point domCentro = {0.0, 0.0};
    Point rootPoint = {0.0, R};

    double distMinNos = 0.15 * calculaDistanciaMinimaTerminal(R, Nterm);

    clock_t inicio = clock();

    //========================================================
    // GERACAO DA ARVORE (CRITERIO GULOSO, IGUAL AO MiniCCO-1 - NAO ALTERADO)
    //========================================================
    int terminaisInseridos = 0;
    int tentativasTotais = 0;

    while (terminaisInseridos < Nterm) {

        Point Pnovo = geraNovoTerminal(arvore, domCentro, R, rootPoint);

        long long testadasAntes = statsAcum.conexoesTestadas;
        long long rejeitadasAntes = statsAcum.conexoesRejeitadas;

        bool ok = insereTerminal(arvore, Pnovo, rootPoint, modoGeometrico, M,
                                  Qterm, gamma, mu, domCentro, R,
                                  raioFixo, distMinNos, &statsAcum, maxFilhos);
        tentativasTotais++;

        if (ok) {
            terminaisInseridos++;

            //Segmento recem-inserido = ultimo do vetor (bifurcacao adiciona 2,
            //multifurcacao adiciona 1 - em ambos os casos o terminal novo e o ultimo)
            const Segment& term = arvore.segmentos.back();

            //Localiza o segmento da raiz para registrar fluxo/raio globais apos esta insercao
            double fluxoRaiz = 0.0, raioRaizAtual = 0.0;
            for (const Segment& s : arvore.segmentos) {
                if (s.paiId == -1) { fluxoRaiz = s.fluxo; raioRaizAtual = s.raio; break; }
            }

            log << "Terminal " << terminaisInseridos << "/" << Nterm
                << " inserido em (" << Pnovo.x << ", " << Pnovo.y << ")"
                << " | segmento id=" << term.id << " pai_id=" << term.paiId
                << " | fluxo(raiz)=" << fluxoRaiz
                << " | raio(raiz)=" << raioRaizAtual
                << " | volume_total=" << funcaoCustoVolume(arvore)
                << " | conexoes testadas nesta insercao=" << (statsAcum.conexoesTestadas - testadasAntes)
                << " | rejeitadas nesta insercao=" << (statsAcum.conexoesRejeitadas - rejeitadasAntes)
                << "\n";
        }
    }

    //========================================================
    // NORMALIZACAO DOS RAIOS PELA RESISTENCIA HIDRAULICA GLOBAL (ja existia)
    //========================================================
    if (modoGeometrico != 0 && !arvore.segmentos.empty()) {
        normalizaRaiosPorResistenciaGlobal(arvore, mu, deltaP, Qperf);
        log << "---------------------------------------------------" << "\n";
        log << "Normalizacao por resistencia hidraulica global aplicada." << "\n";
        log << "Volume total apos normalizacao: " << funcaoCustoVolume(arvore) << "\n";
    }

    clock_t fim = clock();
    double tempoExecucao = ((double) (fim - inicio)) / CLOCKS_PER_SEC;

    //========================================================
    // METRICAS CALCULADAS (reaproveitando funcoes existentes - nada duplicado)
    //========================================================
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

    log << "---------------------------------------------------" << "\n";
    log << "Fim da execucao." << "\n";
    log << "Conexoes aceitas (terminais inseridos): " << terminaisInseridos << "\n";
    log << "Conexoes testadas (total): " << statsAcum.conexoesTestadas << "\n";
    log << "Conexoes rejeitadas (total): " << statsAcum.conexoesRejeitadas << "\n";
    log << "Tempo total: " << tempoExecucao << " s" << "\n";
    log << "===================================================" << "\n";
    log.close();

    //========================================================
    // RELATORIO FINAL (UNICA SAIDA NO TERMINAL)
    //========================================================
    cout << "==================================================" << endl;
    cout << "               RELATORIO DA ARVORE                 " << endl;
    cout << "==================================================" << endl;
    cout << "Variante:                  " << nomeVariante << endl;
    cout << "Modo:                      " << modoStr << endl;
    cout << "Seed:                      " << seed << endl;
    cout << "Numero total de nos:       " << nNos << endl;
    cout << "Numero total de segmentos: " << nSeg << endl;
    cout << "Numero de terminais:       " << nTerm << endl;
    cout << "Comprimento total:         " << compTotal << endl;
    cout << "Volume total (J):          " << volTotal << endl;
    cout << "Raio da raiz:              " << raioDaRaiz << endl;
    cout << "Raio medio dos segmentos:  " << raioMed << endl;
    cout << "Profundidade maxima:       " << profMax << endl;
    cout << "Numero de bifurcacoes:     " << nBifurc << endl;
    cout << "Numero de multifurcacoes:  " << nMultif << endl;
    cout << "Grau medio do grafo:       " << grauMed << endl;
    cout << "Conexoes testadas:         " << statsAcum.conexoesTestadas << endl;
    cout << "Conexoes rejeitadas:       " << statsAcum.conexoesRejeitadas << endl;
    cout << "Tempo de execucao:         " << tempoExecucao << " s" << endl;
    cout << "==================================================" << endl;
    cout << "Arquivos gerados: arvore.vtk, segmentos.txt, execucao.log" << endl;

    if (!verificaSomenteBifurcacoes(arvore, maxFilhos)) {
        cout << endl << "[ALERTA] Inconsistencia estrutural detectada: algum no possui "
             << "mais filhos do que o permitido (" << maxFilhos << ")." << endl;
    }

    //========================================================
    // EXPORTACAO PARA VISUALIZACAO (ja existia, nao alterado)
    //========================================================
    // O arquivo metricas.csv separado foi removido: a mesma informacao (Tabela 2)
    // ja fica dentro de segmentos.txt, e o relatorio do terminal acima ja imprime
    // tudo (agora incluindo modo/seed) - manter os dois era redundante.
    exportarVTK(arvore, "arvore.vtk", R);
    exportarSegmentosTXT(arvore, "segmentos.txt", Nterm, tempoExecucao,
                         statsAcum.conexoesTestadas, statsAcum.conexoesRejeitadas,
                         modoStr, seed);

    ofstream raioFile("raio.txt");
    raioFile << R;
    raioFile.close();

    return 0;
}

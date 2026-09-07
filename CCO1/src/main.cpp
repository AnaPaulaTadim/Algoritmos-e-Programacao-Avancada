#include <iostream> 
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#include "tree.hpp"
#include "domain.hpp"
#include "segment.hpp"
#include "geometry.hpp"

using namespace std;

int main(int argc, char* argv[]){

    //========================================================
    // VERIFICAÇÃO DOS ARGUMENTOS DA LINHA DE COMANDO
    //========================================================
    if(argc != 3) {
        cout << "Uso: ./programa Nterm R" << endl;
        return 1;
    }

    //========================================================
    // LEITURA DOS ARGUMENTOS
    //========================================================
    int Nterm = atoi(argv[1]);
    double R = atof(argv[2]);

    //========================================================
    // INICIALIZA O GERADOR ALEATÓRIO
    //========================================================
    srand((unsigned)time(NULL));

    //========================================================
    // CRIA A ESTRUTURA DA ÁRVORE
    //========================================================
    Tree arvore;
    arvore.totalNos = 0;
    arvore.totalFolhas = 0;
    
    //========================================================
    // CONTADORES DE CONTROLE
    //========================================================
    int conexoesRejeitadas = 0;
    int tentativasInvalidas = 0;
    //========================================================
    // CRIA O NÓ RAIZ (Fixo no topo do Domínio: 0, R)
    //========================================================
    ptrNo raiz = new No;
    raiz->ponto.x = 0.0;
    raiz->ponto.y = R;
    raiz->pai = nullptr;
    raiz->esquerda = nullptr;
    raiz->direita = nullptr;
    raiz->id = 0;
    raiz->terminal = false;

    arvore.raiz = raiz;
    arvore.totalNos = 1;

    //========================================================
    // GERAÇÃO DO PRIMEIRO TERMINAL VIA CRITÉRIO GULOSO
    //========================================================
    Point primeiroPonto;
    bool primeiroValido = false;
    int tentativasPrimeiro = 0;
    int maxTentativasPrimeiro = 5000;
    double melhorScorePrimeiro = -1.0;

    while(tentativasPrimeiro < maxTentativasPrimeiro) {
        
        Point candidato = geraPontoValido(R);

        // REGRAS EXTRA DA ESPECIFICAÇÃO PARA O PRIMEIRO TERMINAL:
        double dEntrada = distancia(candidato, raiz->ponto);
        double dCentro  = sqrt(candidato.x * candidato.x + candidato.y * candidato.y);

        // Deve ficar longe da entrada (0.40*R) e longe do centro (0.35*R)
        if(dEntrada >= 0.40 * R && dCentro >= 0.35 * R) {
            
            
            double score = dEntrada; 

            if(score > melhorScorePrimeiro) {
                melhorScorePrimeiro = score;
                primeiroPonto = candidato;
                primeiroValido = true;
            }
        }
        tentativasPrimeiro++;
    }

    //garante um ponto na extremidade oposta
    if(!primeiroValido) {
        primeiroPonto.x = 0.0;
        primeiroPonto.y = -R * 0.99;
    }

    //========================================================
    // CONECTA O PRIMEIRO TERMINAL À RAIZ
    //========================================================
    ptrNo primeiro = new No;
    primeiro->ponto = primeiroPonto;
    primeiro->pai = raiz;
    primeiro->esquerda = nullptr;
    primeiro->direita  = nullptr;
    primeiro->terminal = true;
    primeiro->id = 1;

    raiz->esquerda = primeiro;
    arvore.totalNos = 2;
    arvore.totalFolhas = 1;

    Segment inicial;
    inicial.a = raiz->ponto;       // (0, R)
    inicial.b = primeiroPonto;     // Ponto sorteado no domínio
    inicial.bifurcado = false;
    inicial.nome = "Raiz->A";

    arvore.segmentos.push_back(inicial);
    
    cout << "--------------------------------------------------" << endl;
    cout << "[SUCESSO] Primeiro Terminal Gerado por Criterio!" << endl;
    cout << "  -> Coordenadas: (" << primeiroPonto.x << ", " << primeiroPonto.y << ")" << endl;
    cout << "  -> Tentativas gastas para achar: " << tentativasPrimeiro << endl;
    



    //========================================================
    // VERIFICAÇÃO EXCLUSIVA PARA APENAS 1 TERMINAL
    //========================================================
    if(Nterm <= 1) {
        cout << "\n[INFO] Execucao finalizada com 1 terminal." << endl;
        cout << "Numero total de nos: " << arvore.totalNos << endl;
        cout << "Numero de folhas:    " << contaFolhas(arvore) << endl;
        
        //Salva os arquivos antes de fechar o programa!
        exportarVTK(arvore, "arvore.vtk", R);
        ofstream raioFile("raio.txt");
        raioFile << R;
        raioFile.close();
        return 0;
    }

    //========================================================
    // GERAÇÃO DOS PRÓXIMOS RAMOS (PROCESSO GULOSO DO ARTIGO)
    //========================================================
    char proximaLetra = 'B';

    for(int i = 2; i <= Nterm; i++){
        
        Point novoPonto;
        bool valido = false;

        int tentativas = 0;
        int maxTentativas = 9000; // Aumentado para garantir melhor busca espacial

        double limiteTerminal = calculaDistanciaMinimaTerminal(R, contaFolhas(arvore) + 1);

        Point melhorPonto;
        double melhorScore = -1.0;

        //Loop de seleção por score guloso para encontar um ponto bem espalhado
        while(tentativas < maxTentativas){

            novoPonto = geraPontoValido(R);

            // Regra Extra de afastamento para os primeiros nós periféricos
            if(arvore.segmentos.size() == 1){
                double dEntrada = distancia(novoPonto, raiz->ponto);
                double dCentro = sqrt(novoPonto.x * novoPonto.x + novoPonto.y * novoPonto.y);
                if(dEntrada < 0.40 * R || dCentro < 0.35 * R){
                    tentativas++;
                    continue;
                }
            }

            double dTerminal = menorDistanciaTerminais(arvore, novoPonto);
            double dSegmento = menorDistanciaSegmentos(arvore, novoPonto);

            // Filtro de restrição mínima espacial
            if(dTerminal >= limiteTerminal && dSegmento >= 0.01 * R){
                
                // Fórmula de pontuação estrita do artigo CCO
                double score = dTerminal + 0.35 * dSegmento;

                // Armazena de forma gulosa o candidato com o MAIOR score (mais espalhado)
                if(score > melhorScore){
                    melhorScore = score;
                    melhorPonto = novoPonto;
                    valido = true;
                }
            }

            tentativas++;
        }

        // Se após as tentativas um ponto ótimo foi selecionado
        if(valido){
            novoPonto = melhorPonto;
        } else {
            //assume o último ponto gerado
            novoPonto = geraPontoValido(R);
        }

        cout << "--------------------------------------------------" << endl;
        cout << "[SUCESSO] Iteracao " << i << " -> Novo Ponto Aceito: Ponto " << proximaLetra << "!" << endl;
        cout << "  -> Coordenadas Finais: (" << novoPonto.x << ", " << novoPonto.y << ")" << endl;
        cout << "  -> Tentativas gastas para achar: " << tentativas << endl;

        //====================================================
        // GERA LISTA DE SEGMENTOS CANDIDATOS
        //====================================================
        vector<Candidato> candidatos = geraCandidatos(arvore, novoPonto);

        if(candidatos.empty()) {
            cout << "  [AVISO] Lista de candidatos vazia para este ponto!" << endl;
            conexoesRejeitadas++;
            continue;
        }

        //====================================================
        // ESCOLHE O MELHOR CANDIDATO (MENOR CUSTO DE CONEXÃO)
        //====================================================
        Candidato melhor = candidatos[0];
        for(size_t j = 1; j < candidatos.size(); j++) {
            if(candidatos[j].custo < melhor.custo) {
                melhor = candidatos[j];
            }
        }

        cout << "  -> Escolhido o Melhor Segmento Alvo para Conexao: " << melhor.segmento.nome << endl;
        cout << "     Coordenadas: (" << melhor.segmento.a.x << ", " << melhor.segmento.a.y << ") -> ("
             << melhor.segmento.b.x << ", " << melhor.segmento.b.y << ")" << endl;
        cout << "     Distancia/Custo de conexao ate ele: " << melhor.custo << endl;
            
        conexoesRejeitadas += (candidatos.size() - 1);

        //====================================================
        // INSERE O NOVO TERMINAL NO SEGMENTO ESCOLHIDO
        //====================================================
        string nomeOriginal = melhor.segmento.nome;
        bool inseriu = insereNo(arvore, melhor.segmento, novoPonto);

        if(!inseriu) {
            cout << "  [ERRO] Falha critica ao inserir no nó!" << endl;
            conexoesRejeitadas++;
        } else {
            if (arvore.segmentos.size() >= 3) {
                size_t sz = arvore.segmentos.size();
                arvore.segmentos[sz - 3].nome = nomeOriginal + "_1";
                arvore.segmentos[sz - 2].nome = nomeOriginal + "_2";
                arvore.segmentos[sz - 1].nome = string("Ramo->") + proximaLetra;
            }
            cout << "  [SUCESSO] Arvore atualizada. Total de folhas: " << contaFolhas(arvore) << endl;
        }

        proximaLetra++;
    }

    //========================================================
    // RELATÓRIO FINAL E EXPORTAÇÃO VTK
    //========================================================
    cout << endl;
    cout << "==================================================" << endl;
    cout << "               RELATORIO DA ARVORE                " << endl;
    cout << "==================================================" << endl;
    cout << "Numero total de nos: " << arvore.totalNos << endl;
    cout << "Numero de folhas:    " << contaFolhas(arvore) << endl;
    cout << "Comprimento total:   " << comprimentoTotal(arvore) << endl;
    cout << "Conexoes rejeitadas: " << conexoesRejeitadas << endl;
    cout << "==================================================" << endl;

    exportarVTK(arvore, "arvore.vtk", R);
    
    ofstream raioFile("raio.txt");
    raioFile << R;
    raioFile.close();

    return 0;
}
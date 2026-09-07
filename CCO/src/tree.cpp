#include "tree.hpp"
#include "candidate.hpp"
#include "point.hpp"
#include "node.hpp"
#include "geometry.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


//Insere um novo terminal em um segmento 
bool insereNo(Tree& arvore, Segment& segmentoAlvo, Point novoPonto) {
    
    // 1. Localizar o segmento alvo no vetor para modificar ou remover
    for (auto it = arvore.segmentos.begin(); it != arvore.segmentos.end(); ++it) {
        if (it->a.x == segmentoAlvo.a.x && it->a.y == segmentoAlvo.a.y &&
            it->b.x == segmentoAlvo.b.x && it->b.y == segmentoAlvo.b.y) {
            
            // 2. Calcular o ponto médio exato do segmento alvo (Onde ocorrerá a bifurcação)
            Point pontoMedio;
            pontoMedio.x = (it->a.x + it->b.x) / 2.0;
            pontoMedio.y = (it->a.y + it->b.y) / 2.0;

            // Guardamos os pontos originais para não perdê-los na remoção
            Point pontoA = it->a;
            Point pontoB = it->b;
            string nomeAntigo = it->nome;

            // 3. Remove o segmento antigo que foi partido
            arvore.segmentos.erase(it);

            // 4. Cria os 3 novos segmentos com as coordenadas exatas
            Segment s1; // Parte 1 do segmento antigo (do início até o ponto médio)
            s1.a = pontoA;
            s1.b = pontoMedio;
            s1.bifurcado = true;
            s1.nome = nomeAntigo + "_1";

            Segment s2; // Parte 2 do segmento antigo (do ponto médio até o fim original)
            s2.a = pontoMedio;
            s2.b = pontoB;
            s2.bifurcado = true;
            s2.nome = nomeAntigo + "_2";

            Segment s3; // O novo ramo (do ponto médio até o novo ponto sorteado)
            s3.a = pontoMedio;
            s3.b = novoPonto; // garante que usa o novoPonto com Y negativo correto
            s3.bifurcado = false;
            s3.nome = "Ramo_Novo";

            // 5. Adiciona os novos segmentos à árvore
            arvore.segmentos.push_back(s1);
            arvore.segmentos.push_back(s2);
            arvore.segmentos.push_back(s3);

            // Atualiza a contagem de nós da estrutura da árvore
            arvore.totalNos += 2; // Adicionou o ponto médio e o novo terminal
            return true;
        }
    }
    return false;
}
//====================================================================================
//Percorre todos os nós da árvore
void percorrerArvore(ptrNo raiz){

    //Verifica se a árvore esta vazia(caso base) 
    if(raiz == nullptr){
        return;
    }

    //*==========================================================
    //*Percorrendo á arvore em pré-ordem(raiz, esquerda, direita)
    //===========================================================
    
    //Percorrendo o nó atual(raiz)
    cout << "ID: " << raiz->id << endl;
    //Verifica todos os nós filhos da esquerda
    percorrerArvore(raiz->esquerda);
    //Verifica todos os nós da direita
    percorrerArvore(raiz->direita);
 

}
//===================================================================================
//Conta quantas folhas existem na árvore
int contaFolhas(Tree& arvore){

    int folhas = 0;

    for(Segment s1 : arvore.segmentos){

        bool ehFolha = true;

        //Verifica se o ponto final deste segmento
        //é usado como início de outro segmento
        for(Segment s2 : arvore.segmentos){

            if(s1.b.x == s2.a.x &&
               s1.b.y == s2.a.y){

                ehFolha = false;

                break;
            }
        }

        if(ehFolha){
            folhas++;
        }
    }

    return folhas;
}
//====================================================================================
//Calcula comprimento total da árvore
double comprimentoTotal(Tree& arvore){

    double total = 0.0;

    for(Segment s : arvore.segmentos){

        total += distancia(s.a, s.b);

    }

    return total;
}
//==========================================================================================
//Lista de segmentos candidatos para o novo terminal
vector<Candidato> geraCandidatos(Tree& arvore, Point novoPonto){

    vector<Candidato> candidatos;

    cout << "\n====================================" << endl;
    cout << "CANDIDATOS PARA O PONTO: ("
         << novoPonto.x << ", "
         << novoPonto.y << ")" << endl;
    cout << "====================================" << endl;

    int idCandidato = 1;

    for(Segment& s : arvore.segmentos){

        Candidato c;
        c.novoPonto = novoPonto;
        c.segmento = s;
        c.eValido = true;

        c.custo = distanciaPontoSegmento(novoPonto,s);

        //====================================================
        // NOVO RAMO QUE SERIA CRIADO
        //====================================================

        Point pontoBif;

        pontoBif.x = (s.a.x + s.b.x)/2.0;
        pontoBif.y = (s.a.y + s.b.y)/2.0;

        Segment novoRamo;
        novoRamo.a = pontoBif;
        novoRamo.b = novoPonto;

        bool intersecta = false;

        //====================================================
        // VERIFICA INTERSEÇÃO COM TODOS OS SEGMENTOS EXISTENTES
        //====================================================

        for(Segment& existente : arvore.segmentos){

            // ignora o próprio segmento que será bifurcado
            if(&existente == &s)
                continue;

            if(intersecaoSegmentos(novoRamo, existente)){

                intersecta = true;
                break;
            }
        }

        if(intersecta){

            cout << "[REJEITADO] "
                 << s.nome
                 << " -> cruza outro segmento"
                 << endl;

            continue;
        }

        cout << "[" << idCandidato << "] "
             << s.nome
             << " | custo = "
             << c.custo
             << endl;

        candidatos.push_back(c);

        idCandidato++;
    }

    cout << "Total de candidatos validos: "
         << candidatos.size()
         << endl;

    cout << "===================================="
         << endl;

    return candidatos;
}
//============================================================================================
void exportarVTK(Tree& arvore, string nomeArquivo, double R){

    ofstream arquivo(nomeArquivo);

    if(!arquivo.is_open()){

        cout << "Erro ao criar arquivo VTK" << endl;

        return;
    }

    //========================================
    // CABEÇALHO VTK
    //========================================

    arquivo << "# vtk DataFile Version 3.0\n";

    arquivo << "Arvore arterial\n";

    arquivo << "ASCII\n";

    arquivo << "DATASET POLYDATA\n";


    //========================================
    // PONTOS
    //========================================

    int totalPontos = arvore.segmentos.size() * 2;

    arquivo << "POINTS "
            << totalPontos
            << " float\n";

    for(Segment s : arvore.segmentos){

        arquivo << s.a.x
                << " "
                << s.a.y
                << " 0\n";

        arquivo << s.b.x
                << " "
                << s.b.y
                << " 0\n";
    }

    //========================================
    // LINHAS
    //========================================

    int totalLinhas = arvore.segmentos.size();

    arquivo << "LINES "
            << totalLinhas
            << " "
            << totalLinhas * 3
            << "\n";

    int indice = 0;

    for(size_t i = 0;
        i < arvore.segmentos.size();
        i++){

        arquivo << "2 "
                << indice
                << " "
                << indice + 1
                << "\n";

        indice += 2;
    }

    arquivo.close();

    cout << "Arquivo VTK exportado com sucesso!"
         << endl;
}

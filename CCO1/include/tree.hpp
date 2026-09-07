#ifndef TREE_HPP
#define TREE_HPP


#include "node.hpp"
#include "candidate.hpp"
#include "point.hpp"
#include <string>


#include <vector>
using namespace std;

//=======================Estrutura da árvore======================================
struct Tree {

    ptrNo raiz;

    int totalNos;
    int totalFolhas;
    
    //Lista de segmentos para inserção 
    vector<Segment> segmentos;

};

//==========================Gerenciamento da árvore=====================================

//Insere um novo nó, por segmento
bool insereNo(Tree& arvore, Segment& segmento, Point novoTerminal); 
//Percorre todos os nós da árvore
void percorrerArvore(ptrNo raiz); 
//Conta quantas folhas existem na árvore
int contaFolhas(Tree& arvore);
//Calcula o comprimento total da árvore
//double comprimentoTotal(ptrNo raiz);
double comprimentoTotal(Tree& arvore);
//Lista de segmentos candidatos para o novo terminal
vector<Candidato> geraCandidatos(Tree& arvore, Point novoPonto);
//gera um arquivo VTK para a visualização da árvore
void exportarVTK(Tree& arvore, string nomeArquivo, double R);

#endif
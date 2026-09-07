# 🧬 Modelagem e Simulação de Árvores Vasculares (CCO & TRICCO)

Este repositório contém a implementação do algoritmo **Constrained Constructive
Optimization (CCO)** para construção de árvores vasculares artificiais e sua
extensão **TRICCO (TRIfurcating CCO)**, desenvolvida para permitir
multifurcações de até três filhos por nó.

O projeto combina **modelagem geométrica**, **otimização construtiva** e
**modelagem hemodinâmica**, permitindo gerar, analisar e visualizar árvores
vasculares sintéticas em um domínio tridimensional.

---

## Objetivos

O projeto possui quatro objetivos principais:

- **Construir árvores vasculares sintéticas** utilizando o método CCO;
- **Implementar uma variante com multifurcações**, denominada TRICCO;
- **Minimizar o volume intravascular total** respeitando restrições geométricas
  e fisiológicas;
- **Comparar quantitativamente** as árvores binárias e as árvores geradas com
  suporte a trifurcações.

A implementação também permite analisar características topológicas,
geométricas e computacionais das árvores geradas.

---

#  Fundamentação

## Constrained Constructive Optimization (CCO)

O CCO é um método construtivo para geração de redes vasculares que adiciona
novos terminais progressivamente a uma árvore existente.

A cada novo terminal, o algoritmo:

1. gera um novo ponto terminal dentro do domínio;
2. verifica os segmentos existentes que podem receber uma nova conexão;
3. testa diferentes posições possíveis para a bifurcação;
4. recalcula as propriedades físicas da árvore;
5. verifica as restrições geométricas;
6. calcula o volume total da árvore;
7. seleciona a configuração viável de menor custo.

Dessa forma, a árvore é construída de maneira incremental, utilizando uma
estratégia de otimização gulosa.

---

#  Modelo físico

Cada vaso é representado como um segmento cilíndrico definido por:

- posição inicial;
- posição final;
- comprimento;
- raio;
- fluxo;
- resistência hidráulica;
- volume;
- identificador do segmento;
- identificador do segmento pai.

### Escoamento de Poiseuille

A resistência hidráulica de cada segmento é calculada utilizando a relação de
Poiseuille:

```text
R = 8 μ l / (π r⁴)

#include "random.hpp"

#include<cstdlib> //Para rand() e srand()
#include<ctime> 

//Gera número aleatório entre min e max
double randomDouble(double min, double max){
    
    //Gera numero aleatorio entre 0 e 1
    double r = (double) rand() / RAND_MAX;

    //Calcula o meio do intervalo e retorna o valor aleatorio 
    return min + r * (max - min);

}
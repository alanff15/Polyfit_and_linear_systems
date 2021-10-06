/* 
 * File:   sislinmmq.h
 * Author: Alan
 *
 * Created on 05 de Maio de 2017, 18:29
 */

#ifndef SISLINMMQ_H
#define SISLINMMQ_H

#include <stdlib.h>

//valor mínimo do determinante para considerar válido resultado do SisLin
#define DET_VALOR_MIN 1e-5

//número máximo de coeficientes no polinômio do Mmq
//defina como 0 para alocação dinâmica de memória (apenas se tiver sistema operacional)
//ou o valor calculado para o pior caso de uso (m*o + 1) para alocação estática (baremetal)
#define PARAMS_MAX_MEM_ALOC 0

#ifdef __cplusplus
extern "C" {
#endif    
    /*
    * Resolve o sistema linear A.x=b por escalonamento
    * x  - vetor de tamanho 'n' com o resultado
    * Ab - matriz de tamanho 'n' x 'n+1' que contem a matriz 'A' e o vetor 'b' na ultima coluna
    * n  - é a ordem da matriz A e número de incógnitas
    * retorna 1 se o resultado for ruim (determinante menor que o valor mínimo)
    * retorna 0 caso contrário
    */
    int SisLin(double *x, double *Ab, int n);

    /*
    * Calcula os coeficientes do polinômio de ajuste de curva pelo método dos mínimos quadrados
    * a  - vetor de tamanho 'm*o + 1' com o resultado dos coeficientes calculados
    * yx - matriz de tamanho 'n' x 'm+1' com os dados, cada linha contem o valor de saída 'y' seguido dos valores de entrada 'x's
    * m  - quantidade de incógnitas (x's)
    * n  - quantidade de linhas de dados
    * o  - ordem do polinômio de saída
    * retorna 1 se houver erro
    * retorna 0 caso contrátrio
    */
    int Mmq(double *a, double *yx, int m, int n, int o);

    /*
    * Calcula o erro quadrático do polinômio calculado comparando com os dados
    * yx - matriz de tamanho 'n' x 'm+1' com os dados, cada linha contem o valor de saída 'y' seguido dos valores de entrada 'x's
    * m  - quantidade de incógnitas (x's)
    * n  - quantidade de linhas de dados
    * o  - ordem do polinômio de saída
    * a  - vetor de tamanho 'm*o + 1' com os coeficientes calculados
    */
    double Erro(double *yx, int m, int n, int o, double *a);

#ifdef __cplusplus
}
#endif

#endif /* SISLINMMQ_H */
/* 
 * File:   linsysmmq.h
 * Author: Alan
 *
 * Created on 05 de Maio de 2017, 18:29
 */

#ifndef SISLINMMQ_H
#define SISLINMMQ_H

#include <stdlib.h>

/*
 * EN: Min value of the determinant to consider valid the linear system result
 * 
 * PT: Valor mínimo do determinante para considerar válido o resultado do linsys
 */
#define DET_VAL_MIN 1e-5

/*
 * EN: Coefficients in polyfit polynomial. Set to 0 for dynamic memory
 * allocation or at least m*o + 1 for static allocation
 * 
 * PT: Coeficientes no polinômio do polyfit. Defina 0 para alocação dinâmica de
 * memória ou, no mínimo, m*o + 1 para alocação estática
 */
#define PARAMS_MAX_MEM_ALOC 0

#ifdef __cplusplus
extern "C" {
#endif    
    /*
     * EN:
     * Solve the linear system A.x=b
     * return - 1 if result is valid or 0 if invalid (determinant less than minimum value)
     * x - vector of size 'n' with the result
     * Ab - matrix of size 'n' x 'n+1' that contains the matrix 'A' and the vector 'b' in the last column
     * n - is the order of matrix A and number of unknowns
     * 
     * PT:
     * Resolve o sistema linear A.x=b por escalonamento
     * retorno - 1 se o resultado for válido ou 0 se inválido (determinante menor que o valor mínimo)
     * x  - vetor de tamanho 'n' com o resultado
     * Ab - matriz de tamanho 'n' x 'n+1' que contem a matriz 'A' e o vetor 'b' na ultima coluna
     * n  - é a ordem da matriz A e número de incógnitas
     */
    int linsys(double *x, double *Ab, int n);

    /*
     * EN:
     * Calculates the coefficients of the fit polynomial by the least squares method
     * return - 1 if result is valid or 0 if invalid (insufficient data)
     * a - vector of size 'm*o + 1' with the result of the calculated coefficients
     * yx - array of size 'n' x 'm+1' with the data, each row contains the output value 'y' followed by the input values 'x's
     * m - number of unknowns (x's)
     * n - number of data rows
     * o - order of the output polynomial
     * 
     * PT:
     * Calcula os coeficientes do polinômio de ajuste pelo método dos mínimos quadrados
     * retorno - 1 se o resultado for válido ou 0 se inválido (dados insuficientes)
     * a  - vetor de tamanho 'm*o + 1' com o resultado dos coeficientes calculados
     * yx - matriz de tamanho 'n' x 'm+1' com os dados, cada linha contem o valor de saída 'y' seguido dos valores de entrada 'x's
     * m  - quantidade de incógnitas (x's)
     * n  - quantidade de linhas de dados
     * o  - ordem do polinômio de saída
     */
    int polyfit(double *a, double *yx, int m, int n, int o);

    /*
     * EN:
     * Calculates the squared error of the polynomial compared to the data
     * yx - array of size 'n' x 'm+1' with the data, each row contains the output value 'y' followed by the input values 'x's
     * m - number of unknowns (x's)
     * n - number of data rows
     * o - order of the output polynomial
     * a - vector of size 'm*o + 1' with calculated coefficients
     * 
     * PT:
     * Calcula o erro quadrático do polinômio comparado aos dados
     * yx - matriz de tamanho 'n' x 'm+1' com os dados, cada linha contem o valor de saída 'y' seguido dos valores de entrada 'x's
     * m  - quantidade de incógnitas (x's)
     * n  - quantidade de linhas de dados
     * o  - ordem do polinômio de saída
     * a  - vetor de tamanho 'm*o + 1' com os coeficientes calculados
     */
    double polyfit_err(double *yx, int m, int n, int o, double *a);

#ifdef __cplusplus
}
#endif

#endif /* SISLINMMQ_H */
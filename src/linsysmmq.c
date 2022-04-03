#include "linsysmmq.h"

#define MOD(EXP) ( (EXP)<0?-(EXP):(EXP) )

double Pow(double a, int b) {
    double resp = 1;
    for (; b > 0; b--) {
        resp *= a;
    }
    return resp;
}

int linsys(double *x, double *Ab, int n) {
    if (n < 1) {
        return 1;
    }
    int ret = 1, Abw = n + 1;
    double aux;
    /*escalonar equações*/
    for (int l, k, i, j = 0; j < n - 1; j++) //varrer matriz da esquerda pra direita
    {
        //ordenar linhas pela coluna j, crescente de cima pra baixo,
        //da diagonal princiapal até a ultima linha
        for (i = j; i < n; i++) //seleção direta
        {
            //apontar pra linha atual
            k = i;
            //comparar em módulo com valores das demais linhas, selecionar o menor deles
            for (l = i + 1; l < n; l++) {
                if (MOD(Ab[l * Abw + j]) < MOD(Ab[k * Abw + j])) {
                    k = l;
                }
            }
            //se o menor valor estiver enrte as demais linhas, trocar linhas
            if (k != i) {
                for (l = 0; l <= n; l++) {
                    aux = Ab[k * Abw + l];
                    Ab[k * Abw + l] = Ab[i * Abw + l];
                    Ab[i * Abw + l] = aux;
                }
            }
        }
        //rolar linhas para cima (shift), da ultima linha até a diagonal principal
        //até que o elemento da diagonal principal não seja 0
        for (k = 0; Ab[j * Abw + j] == 0 && k < n; k++) {
            //trocar linha 'i' com linha 'i+1'
            for (i = j; i < n - 1; i++) {
                for (l = 0; l <= n; l++) {
                    aux = Ab[i * Abw + l];
                    Ab[i * Abw + l] = Ab[(i + 1) * Abw + l];
                    Ab[(i + 1) * Abw + l] = aux;
                }
            }
        }
        //zerar elementos
        for (i = n - 1; i > j; i--) {
            if (Ab[i * Abw + j] != 0) {
                //elemento a ser zerado apontado por (i,j)
                //pivo
                aux = -Ab[i * Abw + j] / Ab[j * Abw + j];
                Ab[i * Abw + j] = 0;
                for (l = j + 1; l <= n; l++) {
                    Ab[i * Abw + l] += Ab[j * Abw + l] * aux;
                }
            }
        }
    }
    /*validar sistema pelo determinante*/
    aux = Ab[0 * Abw + 0];
    for (int i = 1; i < n; i++) {
        aux *= Ab[i * Abw + i];
    }
    if (MOD(aux) < DET_VAL_MIN) {
        ret = 0;
    }
    /*resolver sistema triangular superior*/
    for (int j, i = n - 1; i >= 0; i--) {
        x[i] = Ab[i * Abw + n];
        for (j = n - 1; j > i; j--) {
            x[i] -= x[j] * Ab[i * Abw + j];
        }
        x[i] /= Ab[i * Abw + i];
    }
    return ret;
}

int polyfit(double *a, double *yx, int m, int n, int o) {
    /*
    Resolver o sistema linear A.a=b
     */
    int ret, l, c, i, k, kl, j, jl;
    //se dados forem nulos retorna erro
    if (n < 1) {
        return 1;
    }
    int tam = m * o + 1; //ordem da matriz 'A' e tamanho dos vetores 'a' e 'b'
    //se dados não forem suficientes retornar erro
    if (n < tam) {
        return 1;
    }

#if (PARAMS_MAX_MEM_ALOC == 0)
    //alocar memória do sistema linear
    double *Ab = (double*) malloc(tam * (tam + 1) * sizeof (double));
#else
    double Ab[(PARAMS_MAX_MEM_ALOC * (PARAMS_MAX_MEM_ALOC + 1))];
#endif

    int Abw = tam + 1;
    int yxw = m + 1;

    /*preencher matriz A*/
    //elemento 0,0
    Ab[0 * Abw + 0] = n;
    //demais elementos
    l = 1;
    for (jl = 1; jl <= m; jl++) {
        for (kl = 1; kl <= o; kl++) {
            //coluna 0
            c = 0;
            //x_jl^kl
            Ab[l * Abw + c] = 0;
            for (i = 0; i < n; i++) {
                Ab[l * Abw + c] += Pow(yx[i * yxw + jl], kl);
            }
            //copiar para linha 0
            Ab[c * Abw + l] = Ab[l * Abw + c];
            //demais colunas
            for (j = 1; j <= m; j++) {
                for (k = 1; k <= o; k++) {
                    //próxima coluna
                    c++;
                    //x_jl^kl * x_j^k
                    Ab[l * Abw + c] = 0;
                    for (i = 0; i < n; i++) {
                        Ab[l * Abw + c] += Pow(yx[i * yxw + jl], kl) * Pow(yx[i * yxw + j], k);
                    }
                }
            }
            //próxima linha
            l++;
        }
    }
    /*Preencher vetor b*/
    //elemento 0
    Ab[0 * Abw + tam] = 0;
    for (i = 0; i < n; i++) {
        Ab[0 * Abw + tam] += yx[i * yxw + 0];
    }
    //demais elementos
    l = 1;
    for (j = 1; j <= m; j++) {
        for (k = 1; k <= o; k++) {
            Ab[l * Abw + tam] = 0;
            for (i = 0; i < n; i++) {
                Ab[l * Abw + tam] += yx[i * yxw + 0] * Pow(yx[i * yxw + j], k);
            }
            l++;
        }
    }

    /*resolver sistema*/
    ret = linsys(a, Ab, tam);

#if (PARAMS_MAX_MEM_ALOC == 0)
    //liberar memória do sistema linear
    free(Ab);
#endif
    return ret;
}

double polyfit_err(double *yx, int m, int n, int o, double *a) {
    double aux, err = 0;
    int k, j, i, c;
    int yxw = m + 1;
    //calcular valor do polinômio para cada linha de dados
    for (i = 0; i < n; i++) {
        c = 0;
        aux = a[c++];
        for (j = 1; j <= m; j++) {
            for (k = 1; k <= o; k++) {
                aux += a[c++] * Pow(yx[i * yxw + j], k);
            }
        }
        //somar quadrado do erro da linha
        err += Pow(yx[i * yxw + 0] - aux, 2);
    }
    return err;
}

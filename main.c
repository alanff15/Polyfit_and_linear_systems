#include <stdio.h>
#include <stdlib.h>

#include "sislinmmq.h"

int main(int argc, char** argv) {
    int n = 5; //quantidade de dados
    int m = 1; //número de incógnitas
    int o = 2; //ordem do polinômio

    double a[3]; //coeficientes [m*o + 1]
    double yx[5][2]; //dados [n][m+1]

    //dados Y e x
    yx[0][0] = 2.0; yx[0][1] = 2.0;
    yx[1][0] = 2.1; yx[1][1] = 2.1;
    yx[2][0] = 2.2; yx[2][1] = 2.2;
    yx[3][0] = 2.3; yx[3][1] = 2.3;
    yx[4][0] = 2.4; yx[4][1] = 2.4;

    printf("dados:\n%10c%10c\n", 'y', 'x');
    for (int j, i = 0; i < n; i++) {
        for (j = 0; j < m + 1; j++) {
            printf("%10g", yx[i][j]);
        }
        printf("\n");
    }

    //calcular
    Mmq(a, (double*) yx, m, n, o);

    printf("\ncoeficientes:\n");
    for (int i = 0; i < (m * o + 1); i++) printf("a[%d] = %g\n", i, a[i]);
    
    printf("\nerro:\n%g\n", Erro((double*)yx, m, n, o, a));

    getchar();
    return (EXIT_SUCCESS);
}
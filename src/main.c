#include <stdio.h>
#include <stdlib.h>

#include "linsysmmq.h"

void ex_linsys();
void ex_polyfit();

int main(int argc, char** argv) {
    ex_linsys();
    ex_polyfit();
    getchar();
    return (EXIT_SUCCESS);
}

void ex_linsys() {
    printf("-------------------------------------------------\n");
    printf("sislin:\n");

    //A*x = b
    int n = 4;
    double x[4];
    double Ab_eqs[4][5] = //EQs
    {
        1.2, 3.5, -2.0, 2.4, -4.2,
        2.2, -2.5, 1.0, 2.4, 4.2,
        -1.0, -2.0, 1.0, -0.4, -0.6,
        -1.5, -3.0, 2.5, 2.0, -0.6
    };

    //print
    printf(" Eqs:\n");
    for (int j, i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (Ab_eqs[i][j] > 0) {
                printf("  +");
            } else {
                printf("  ");
            }
            printf("%3.1f * x%d", Ab_eqs[i][j], i);
        }
        printf(" = %4.2f", Ab_eqs[i][n]);
        printf("\n");
    }

    //calc
    linsys(x, (double*) Ab_eqs, n);

    //result
    printf(" result:\n");
    for (int i = 0; i < n; i++) {
        printf("  x%d = %4.1f\n", i, x[i]);
    }
    printf("\n");
}

void ex_polyfit() {
    printf("-------------------------------------------------\n");
    printf("polyfit:\n");

    int n = 5; //data size / quantidade de dados
    int m = 1; //unknowns / número de incógnitas
    int o = 2; //order / ordem do polinômio

    double a[3]; //coefs [m*o + 1]
    double yx[5][2] = //data [n][m+1]
    {
        1.40, 1.00,
        1.50, 1.20,
        1.40, 1.40,
        1.10, 1.60,
        0.59, 1.80
    };

    //print
    printf("%6c%6c\n", 'y', 'x');
    for (int j, i = 0; i < n; i++) {
        for (j = 0; j < m + 1; j++) {
            printf("%6.2f", yx[i][j]);
        }
        printf("\n");
    }

    //calc
    polyfit(a, (double*) yx, m, n, o);

    //result
    printf(" coefs:\n");
    for (int i = 0; i < (m * o + 1); i++) printf("  a%d = %5.2f\n", i, a[i]);

    printf(" err:\n  %g\n", polyfit_err((double*) yx, m, n, o, a));
}

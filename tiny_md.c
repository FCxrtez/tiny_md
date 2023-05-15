#define _XOPEN_SOURCE 500 // M_PI
#include "core.h"
#include "parameters.h"
#include "wtime.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

float rxyz[3 * N], vxyz[3 * N], fxyz[3 * N];

int main()
{
    int m = cbrtf(N / 4.0);
    assert(m * m * m * 4 == N);

    // FILE *file_xyz, *file_thermo;
    // file_xyz = fopen("trajectory.xyz", "w");
    // file_thermo = fopen("thermo.log", "w");
    float Ekin, Epot, Temp, Pres; // variables macroscopicas
    float Rho, cell_V, cell_L, tail, Etail, Ptail;

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);
    printf("# densidad, volumen, energía potencial media, presión media\n");
    // fprintf(file_thermo, "# t Temp Pres Epot Etot\n");

    srand(SEED);
    float t = 0.0, sf;
    float Rhob;
    Rho = RHOI;
    init_pos(rxyz, Rho);
    float start = wtime();
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1 * (float)m;
        cell_V = (float)N / Rho;
        cell_L = cbrtf(cell_V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
        Etail = tail * (float)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrtf(Rhob / Rho);
        for (int k = 0; k < 3 * N; k++) { // reescaleo posiciones a nueva densidad
            rxyz[k] *= sf;
        }
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < TEQ; i++) { // loop de equilibracion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrtf(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) { // reescaleo de velocidades
                vxyz[k] *= sf;
            }
        }

        int mes = 0;
        float epotm = 0.0, presm = 0.0;
        for (i = TEQ; i < TRUN; i++) { // loop de medicion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrtf(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) { // reescaleo de velocidades
                vxyz[k] *= sf;
            }

            if (i % TMES == 0) {
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;

                // fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot, Epot + Ekin);
                // fprintf(file_xyz, "%d\n\n", N);
                // for (int k = 0; k < 3 * N; k += 3) {
                //     fprintf(file_xyz, "Ar %e %e %e\n", rxyz[k + 0], rxyz[k + 1], rxyz[k + 2]);
                // }
            }

            t += DT;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (float)mes, presm / (float)mes);
    }

    float elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Natom 40
#define L 10
#define Nsim 30001
#define T 10
#define PI 3.145926

void init_xy(double x[Natom], double y[Natom]){
    
    for(size_t i = 0; i < Natom; i++){
        x[i] = L*(double)rand()/RAND_MAX;
        y[i] = L*(double)rand()/RAND_MAX;
    }
}

double dist2(double x1, double x2, double y1, double y2){
    return (x1-x2)*(x1-x2) + (y1 - y2)*(y1-y2);
}

double config_NRG(double x[Natom], double y[Natom]) {
    double E = 0.0, dist;
    
    for (size_t i = 0; i < Natom; i++) {
        for (size_t j = i + 1; j < Natom; j++) { 
            dist = dist2(x[i], x[j], y[i], y[j]);
            dist = dist * dist;

            E += 4.0 / (dist * dist) - 4.0 / dist;
        }
    }

    return E;
}

void MC_move(double x[Natom], double y[Natom], int idx) {
    double d = 5.0 * (double)rand()/RAND_MAX;
    double arg = 2 * PI * (double)rand()/RAND_MAX;

    double xtmp, ytmp;
    double xcurr = x[idx], ycurr = y[idx];

    xtmp = x[idx] + d * cos(arg);
    ytmp = y[idx] + d * sin(arg);

    // no pbc: reject move if out of bounds
    if (xtmp < 0 || xtmp > L || ytmp < 0 || ytmp > L) {
        return;
    }

    double Ei = config_NRG(x, y);

    x[idx] = xtmp;
    y[idx] = ytmp;

    double Ef = config_NRG(x, y);

    double boltzmann_fac = exp(-(Ef - Ei) / 1.0/ T);
    double n = (double)rand() / RAND_MAX;
    
    //Metropolis step
    if (n < boltzmann_fac) {
        return;
    } else {
        x[idx] = xcurr;
        y[idx] = ycurr;
    }
}


int main(void){
    int seed[3] = {42,0xC0FFEE,12};
    srand(seed[2]);
    double x[Natom] = {0};
    double y[Natom] = {0};
    init_xy(x,y);


    printf("init energy of config: %.2e\n", config_NRG(x, y));
    //printf("seed: %d\n", 0xC0FFEE);
    //FILE *positions = fopen("positionsT10seed12.csv", "w");
    
    FILE *energies = fopen("energiesT10seed12.csv", "w");
    fprintf(energies, "STEP,E\n");

    for(int i = 0; i < Nsim; i++){
        int idx = Natom*(double)rand()/RAND_MAX;
        MC_move(x, y,idx);

        if(i%1000 == 0){
            fprintf(energies, "%d,%.4f\n", i, config_NRG(x,y));
            
        }
    }
    //fclose(positions);
    fclose(energies);

    return EXIT_SUCCESS;
}

//GRAVEYARD:

/*

void init_xy(double x[Natom], double y[Natom]) {
    for (size_t i = 0; i < Natom; i++) {
        int valid;
        do {
            valid = 1;  // Assume the new point is valid
            x[i] = L * (double)rand() / (double)RAND_MAX;
            y[i] = L * (double)rand() / (double)RAND_MAX;

            // Check against all previous atoms
            for (size_t j = 0; j < i; j++) {
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dist = sqrt(dx * dx + dy * dy);

                if (dist < TOL) {
                    valid = 0;  // Too close, retry
                    break;
                }
            }
        } while (!valid);  // Repeat until a valid position is found
    }
}



    while(abs(config_NRG(x,y)) > 1e20){
        int idx = Natom*rand()/RAND_MAX;
        MC_move(x,y,idx);
    }

    for(int i = 0; i < Natom; i++){
        printf("%lf ", x[i]);
    }
    printf("\n");
    for(int i = 0; i < Natom; i++){
        printf("%lf ", y[i]);
    }






*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Natom 40
#define Nsim 30001
#define T 0.5e1
#define p 1e1
#define PI 3.145926
double L = 10.0;
double Lf = 10.0;

void init_xy(double x[Natom], double y[Natom], double z[Natom]){
    
    for(size_t i = 0; i < Natom; i++){
        x[i] = L*(double)rand()/RAND_MAX;
        y[i] = L*(double)rand()/RAND_MAX;
        z[i] = L*(double)rand()/RAND_MAX;
        
    }
}

double dist2(double x1, double x2, double y1, double y2, double z1, double z2){
    return (x1-x2)*(x1-x2) + (y1 - y2)*(y1-y2) + (z1 - z2)*(z1 - z2);
}

double config_NRG(double x[Natom], double y[Natom], double z[Natom]) {
    double E = 0.0, dist;
    
    for (size_t i = 0; i < Natom; i++) {
        for (size_t j = i + 1; j < Natom; j++) { 
            dist = dist2(x[i], x[j], y[i], y[j], z[i], z[j]);
            dist = dist * dist;

            E += 4.0 / (dist * dist) - 4.0 / dist;
        }
    }

    return E;
}


void MC_move(double x[Natom], double y[Natom], double z[Natom], int idx) {
    double d = 5.0 * (double)rand()/RAND_MAX;
    double theta = 2 * PI * ((double)rand()/RAND_MAX);
    double phi = PI * ((double)rand()/RAND_MAX);

    double xtmp, ytmp, ztmp;
    double xcurr = x[idx], ycurr = y[idx], zcurr = z[idx];

    xtmp = x[idx] + d * sin(phi) * cos(theta);
    ytmp = y[idx] + d * sin(phi) * sin(theta);
    ztmp = z[idx] + d * cos(phi);

    // no pbc: reject move if out of bounds
    if (xtmp < 0 || xtmp > L || ytmp < 0 || ytmp > L || ztmp < 0 || ztmp > L) {
        return;
    }

    double Ei = config_NRG(x, y, z);

    x[idx] = xtmp;
    y[idx] = ytmp;
    z[idx] = ztmp;

    double Ef = config_NRG(x, y, z);

    double deltaV = (L*L*L) - (Lf*Lf*Lf);
    //THINK
    double boltzmann_arg = -(1.0/T)*((Ef - Ei) - p*deltaV) + Natom*log(((L*L*L)/(Lf*Lf*Lf)));
    double boltzmann_fac = exp(boltzmann_arg);
    double p_acc = fmin(1, boltzmann_fac);
    double n = (double)rand() / RAND_MAX;
    
    //Metropolis step
    if (n < p_acc) {
        return;
    } else {
        x[idx] = xcurr;
        y[idx] = ycurr;
        z[idx] = zcurr;
    }
}


int main(void){
    //int seed[3] = {42,0xC0FFEE,12};
    srand(time(NULL));
    double x[Natom] = {0};
    double y[Natom] = {0};
    double z[Natom] = {0};
    init_xy(x,y,z);


    printf("init energy of config: %.2e\n", config_NRG(x, y, z));
    //printf("seed: %d\n", 0xC0FFEE);
    //FILE *positions = fopen("positionsT10seed12.csv", "w");
    //fprintf(positions, "STEP,X,Y,Z\n");
    FILE *volumes = fopen("volumesT05e1p1e1.csv", "w");
    fprintf(volumes, "STEP,V\n");
    FILE *energies = fopen("energiesT05e1p1e1.csv", "w");
    fprintf(energies, "STEP,E\n");

    for(int i = 0; i < Nsim; i++){
        int idx = Natom*(double)rand()/RAND_MAX;
        double choice = rand() / (double)RAND_MAX;
        if(choice <= 0.95){
            MC_move(x, y, z, idx);
        } else {
            Lf = L;
            if(choice <= 0.5){
                L *=1.0000001;
            } else {
                L *= (1-0.0000001);
            }
            double factor = pow((L*L*L) / (Lf*Lf*Lf), 1.0/3.0);
            for(size_t i = 0; i < Natom; i++){
                x[i] *= factor;
                y[i] *= factor;
                z[i] *= factor;
            }
        }
        

        if(i%1000 == 0){
            fprintf(energies, "%d,%.4f\n", i, config_NRG(x,y,z));
            fprintf(volumes, "%d,%.4f\n", i, L*L*L);
            /*for(size_t i = 0; i < Natom; i++){
                fprintf(positions,"%zu,%.4f,%.4f,%.4f\n",i,x[i],y[i],z[i]);
            }*/
            printf("L is now: %.4f\n", L);
        }
    }
    //fclose(positions);
    
    fclose(energies);
    fclose(volumes);

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
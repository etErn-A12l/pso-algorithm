#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 10 // Number of particles
#define D 30 // Number of dimensions
#define MAX_ITER 1000

// Structure to hold a particle's position and velocity
typedef struct {
    double position[D];
    double velocity[D];
    double bestPosition[D];
    double bestError;
} Particle;

// Function to initialize a particle's position and velocity
void initializeParticle(Particle *particle) {
    for (int i = 0; i < D; i++) {
        particle->position[i] = (double)rand() / RAND_MAX;
        particle->velocity[i] = 0.0;
        particle->bestPosition[i] = particle->position[i];
    }
    particle->bestError = INFINITY;
}

// Function to evaluate a particle's fitness
double evaluateFitness(double *position) {
    double sum = 0.0;
    for (int i = 0; i < D; i++) {
        sum += pow(position[i], 2.0);
    }
    return sum;
}

// PSO function
void pso() {
    // Initialization
    Particle particles[N];
    Particle globalBest;
    globalBest.bestError = INFINITY;
    for (int i = 0; i < N; i++) {
        initializeParticle(&particles[i]);
        double error = evaluateFitness(particles[i].position);
        if (error < particles[i].bestError) {
            particles[i].bestError = error;
            for (int j = 0; j < D; j++) {
                particles[i].bestPosition[j] = particles[i].position[j];
            }
        }
        if (error < globalBest.bestError) {
            globalBest.bestError = error;
            for (int j = 0; j < D; j++) {
                globalBest.bestPosition[j] = particles[i].position[j];
            }
        }
    }

    // Main loop
    for (int iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < N; i++) {
            // Update particle's velocity
            for (int j = 0; j < D; j++) {
                double r1 = (double)rand() / RAND_MAX;
                double r2 = (double)rand() / RAND_MAX;
                particles[i].velocity[j] = 0.5 * particles[i].velocity[j] + 2.0 * r1 * (particles[i].bestPosition[j] - particles[i].position[j]) + 2.0 * r2 * (globalBest.bestPosition[j] - particles[i].position[j]);
            }

            // Update particle's position
            for (int j = 0; j < D; j++) {
                particles[i].position[j] += particles[i].velocity[j];
                if (particles[i].position[j] < 0.0) {
                    particles[i].position[j] = 0.0;
                } else if (particles[i].position[j] > 1.0) {
                    particles[i].position[j] = 1.0;
                }
            }

            // Update particle's best position and global best
            double error = evaluateFitness(particles[i].position);
            if (error < particles[i].bestError) {
                particles[i].bestError = error;
                for (int j = 0; j < D; j++) {
                    particles[i].bestPosition[j] = particles[i].position[j];
                }
            }
            if (error < globalBest.bestError) {
                globalBest.bestError = error;
                for (int j = 0; j < D; j++) {
                    globalBest.bestPosition[j] = particles[i].position[j];
                }
            }
        }
    }

    // Print the best solution found
    printf("Best solution found: ");
    for (int i = 0; i < D; i++) {
        printf("%f ", globalBest.bestPosition[i]);
    }
    printf("\n");
}

int main() {
    pso();
    return 0;
}

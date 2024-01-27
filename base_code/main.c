#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define SWARM 10 // Number of particles
#define D 1      // Number of dimensions
#define MAX_ITER 100
#define W 0.5 // Inertia Weight
#define C1 2  // Acceleration Factor
#define C2 2  // Acceleration Factor

// Structure to hold a particle's position and velocity
typedef struct
{
    double position[D];
    double velocity[D];
    double bestPosition[D];
    double bestFit;
} Particle;


//======== Function Declarations ===========
/*--------------------------------------*/

void initializeParticle(Particle *);
void rankPrticles(Particle *, Particle *);
double evaluateFitness(double *);
void updateVelocity(Particle *, Particle *);
void updatePosition(Particle *);
void pso();

/*---------------MAIN FUNCTION------------------*/

int main(int argc, char const *argv[])
{
    pso();
    return 0;
}

//======== Function Declarations ===========
/*--------------------------------------*/

// Function to initialize a particle's position and velocity
void initializeParticle(Particle *particle)
{
    for (int i = 0; i < D; i++)
    {
        particle->position[i] = ((double)rand() / (double)(RAND_MAX)) * (5.12 + 5.12) - 5.12; // Between -5 and 5
        particle->velocity[i] = ((double)rand() / (double)(RAND_MAX)) * (1 - 0) + 0;          // Between 0 and 1
        particle->bestPosition[i] = particle->position[i];
    }
    particle->bestFit = INFINITY;
}

// Function to Rank Particles and Best Solutions
void rankPrticles(Particle *particle, Particle *global)
{
    double fit = evaluateFitness(particle->position);
    if (fit < particle->bestFit)
    {
        particle->bestFit = fit;
        for (int j = 0; j < D; j++)
        {
            particle->bestPosition[j] = particle->position[j];
        }
    }
    if (fit < global->bestFit)
    {
        global->bestFit = fit;
        for (int j = 0; j < D; j++)
        {
            global->bestPosition[j] = particle->position[j];
        }
    }
}

// Function to evaluate a particle's fitness
double evaluateFitness(double *position)
{
    /*
        Formulae for 3 Dimension:
        F(x) = 10 * (x1 - 1) ^ 2 + 20 * (x2 - 2) ^ 2 + 30 * (x3 - 3) ^ 2
    */
    double sum = 0.0;
    // for (int i = 0; i < D; i++)
    // {
    // sum += ((i + 1) * 10) * pow((position[i] - (i + 1)), 2.0);
    // }
    // sum += pow(position[0], 2.0) + pow(position[1], 2.0) - cos(18 * position[0]) - cos(18 * position[1]);
    int n = 4;
    for (int i = 0; i < n; i++) // n = 4
    {
        sum += pow(position[0], 2.0) - 10 * cos(2 * 3.17 * position[0]);
    }
    sum = sum + (10 * n);
    return sum;
}

void updateVelocity(Particle *particle, Particle *globalBest)
{
    for (int j = 0; j < D; j++)
    {
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        particle->velocity[j] = W * particle->velocity[j] + C1 * r1 * (particle->bestPosition[j] - particle->position[j]) + C2 * r2 * (globalBest->bestPosition[j] - particle->position[j]);
    }
}

void updatePosition(Particle *particle)
{
    for (int j = 0; j < D; j++)
    {
        particle->position[j] += particle->velocity[j]; // Xi ^ (t+1) = (Xi ^ t) + (Vi ^ t)

        if (particle->position[j] < -5.0) // Checking Minimum Value
            particle->position[j] = 0.0;

        else if (particle->position[j] > 5.0) // Checking Maximum Value
            particle->position[j] = 1.0;
    }
}

// PSO function
void pso()
{
    // Initialization
    Particle particles[SWARM];
    Particle globalBest;
    globalBest.bestFit = INFINITY;

    // Rank the Particles and Find Local & Global Best
    for (int i = 0; i < SWARM; i++)
    {
        initializeParticle(&particles[i]);
        rankPrticles(&particles[i], &globalBest);
    }

// Main loop
#pragma omp parallel for schedule(static) // Using OpenMP for parallel computation
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        for (int i = 0; i < SWARM; i++)
        {
            // Update particle's position
            updatePosition(&particles[i]);

            // Update particle's velocity
            updateVelocity(&particles[i], &globalBest);

            // Rank Particles
            rankPrticles(&particles[i], &globalBest);
        }
        printf("Itr = %d, Best = %f\n", iter, globalBest.bestFit);
    }

    // Print the best solution found
    printf("\nBest solution found with fitness: %lf\n\n\t", globalBest.bestFit);
    for (int i = 0; i < D; i++)
    {
        printf("%.3f ", globalBest.bestPosition[i]);
    }
    printf("\n");
}

/*

1. Initialize position and velocity
2. Calculate fitnesses and global_best
3. Update velocity and position
4. Calculate fitness and find global_best
5. Increase iteration

*/
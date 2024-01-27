#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define SWARM 10 // Number of particles
#define D 10     // Number of dimensions
#define MAX_ITER 100
#define W 0.5 // Inertia Weight
#define C1 2  // Acceleration Factor
#define C2 2  // Acceleration Factor

// File Paths
const char *data_file = "/home/xron/Drives/Docs_And_Media/Documents/Everything/Code-More/Application Development/Programming Languages/C/Others/HK mam TSP Prob 2/Datasets/ftv55.txt";
const char *sol_file = "/home/xron/Drives/Docs_And_Media/Documents/Everything/Code-More/Application Development/Programming Languages/C/Others/HK mam TSP Prob 2/Solutions/ftv55_sol.txt";

// Structure to hold a particle's position and velocity
typedef struct
{
    int position[D];
    double velocity[D];
    int bestPosition[D];
    long bestFit;
} Particle;

// Global variables
int **matrix;

//======== Function Declarations ===========
/*--------------------------------------*/

void initializeParticle(Particle *);
void rankPrticles(Particle *, Particle *);
double evaluateFitness(int *);
void updateVelocity(Particle *, Particle *);
void updatePosition(Particle *);
int **read_matrix(const char *);
long cal_fitness(int []);
void pso();

/*---------------MAIN FUNCTION------------------*/

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    pso();
    return 0;
}

//======== Function Declarations ===========
/*--------------------------------------*/

int **read_matrix(const char *data_file)
{
    FILE *file = fopen(data_file, "r");
    char line[10000];

    if (file == NULL)
    {
        printf("Error opening file\n");
        return NULL;
    }

    // printf("\n  %d x %d Cost Matrix\n\n", D, D);

    // Allocate memory for the 2D array
    int **matrix = (int **)malloc(D * sizeof(int *));
    for (int i = 0; i < D; i++)
        matrix[i] = (int *)malloc(D * sizeof(int));

    // Reset file pointer to beginning
    fseek(file, 0L, SEEK_SET);

    // Read the matrix values
    for (int i = 0; i < D; i++)
    {
        for (int j = 0; j < D; j++)
            fscanf(file, "%d", &matrix[i][j]);
    }

    fclose(file);
    return matrix;
}

// Function to initialize a particle's position and velocity
void initializeParticle(Particle *particle)
{
    bool pos[D + 1] = {false}; // Initially make all false
    int count = 0;
    int upper = D - 1, lower = 0;

    // Randomly Generate Positions
    while (count != D)
    {
        int num = (rand() % (upper - lower + 1)) + lower;
        if (!pos[num])
        {
            particle->position[count] = num;
            pos[num] = true;
            count++;
        }
    }

    for (int i = 0; i < D; i++)
    {
        particle->velocity[i] = ((double)rand() / (double)(RAND_MAX)) * (1 - 0) + 0; // Between 0 and 1
        particle->bestPosition[i] = particle->position[i];
    }

    particle->bestFit = INFINITY;
}

// Function to Rank Particles and Best Solutions
void rankPrticles(Particle *particle, Particle *global)
{
    double fit = cal_fitness(particle->position);
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

// Function to evaluate a particle's fitness
long cal_fitness(int pos[])
{
    long fit = 0;
    for (int i = 0; i < D - 1; i++)
    {
        assert(pos[i] != pos[i + 1]);
        if (pos[i] != pos[i + 1])
            fit += matrix[pos[i]][pos[i + 1]];
    }
    fit += matrix[pos[D - 1]][pos[0]];
    return fit;
}

// PSO function
void pso()
{
    // Read Matrix from file
    matrix = read_matrix(data_file);

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
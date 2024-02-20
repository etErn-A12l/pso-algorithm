#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <assert.h>
#include <string.h>

#define SWARM 10 // Number of particles
#define D 17     // Number of dimensions
#define MAX_ITER 10000
#define W 0.5    // Inertia Weight
#define C1 2     // Acceleration Factor
#define C2 2     // Acceleration Factor
#define TEMP 100 // Initial Temperature

// File Paths
const char *data_file = "/home/eternal/Documents/CodeSpace/bat-algorithm/Datasets/br17.txt";
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
int **matrix; // Cost Matrix

//======== Function Declarations ===========
/*--------------------------------------*/

void initializeParticle(Particle *);
void rankPrticles(Particle *, Particle *);
void updateVelocity(Particle *, Particle *);
void updatePosition(Particle *);
int **read_matrix(const char *);
long cal_fitness(int[]);
void pso();

void reverse(int *, int, int);
void apply_3_opt(int *, int, int, int);
void Simul_Annel_3_Opt(int *, double, double);
void InverseMutation(int *);
void swap(int *, int, int);

/*---------------MAIN FUNCTION------------------*/

int main(int argc, char const *argv[])
{
    srand(1);
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

    particle->bestFit = __INT64_MAX__;
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
        float vel = W * particle->velocity[j] + C1 * r1 * (particle->bestPosition[j] - particle->position[j]) + C2 * r2 * (globalBest->bestPosition[j] - particle->position[j]);
        int i_vel = (int)vel;
        particle->velocity[j] = vel - i_vel;
        // printf("\nnew velocity = %.4f", particle->velocity[j]);
    }
}

void updatePosition(Particle *Particle)
{
    float max_vel = 0;
    for (int k = 0; k < D; k++)
    {
        max_vel = (Particle->velocity[k] > max_vel) ? Particle->velocity[k] : max_vel;
    }

    Simul_Annel_3_Opt(Particle->position, TEMP, max_vel);

    long new_fit = cal_fitness(Particle->position);

    if (Particle->bestFit < new_fit)
    {
        InverseMutation(Particle->position);
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

// ===================================

void reverse(int *tour, int i, int j)
{
    if (i < j)
    {
        while (i < j)
        {
            swap(tour, i, j);
            i++;
            j--;
        }
    }
    else
    {
        int prev_low = i;
        while (prev_low != j)
        {
            swap(tour, i, j);
            i = (i + 1) % D;
            if (j - 1 == -1)
                j = D - 1;
            else
                --j;
        }
    }
}

void swap(int *pos, int i, int j)
{
    int temp = pos[i];
    pos[i] = pos[j];
    pos[j] = temp;
}

void InverseMutation(int *pos)
{
    int lower = rand() % D, upper = rand() % D;

    while (lower == upper)
        upper = rand() % D;

    reverse(pos, lower, upper);
}

void apply_3_opt(int *path, int i, int j, int k)
{
    // Apply the 3-opt move to the solution
    int tmp_path[D];
    memcpy(tmp_path, path, D * sizeof(int));

    reverse(tmp_path, i + 1, j);
    reverse(tmp_path, j + 1, k);

    memcpy(path, tmp_path, D * sizeof(int));
}

void Simul_Annel_3_Opt(int *path, double temp, double cooling_rate)
{

    int tour[D];
    memcpy(tour, path, D * sizeof(int));
    while (temp > 1)
    {
        // Generate a new solution by applying a 3-opt move
        int i = rand() % (D - 1);
        int j = rand() % (D - 1);
        int k = rand() % (D - 1);

        while (i == j || j == k || i == k)
        {
            j = rand() % (D - 1);
            k = rand() % (D - 1);
        }

        int small = (i <= j && i <= k) ? i : ((j <= i && j <= k) ? j : k);
        int big = (i >= j && i >= k) ? i : ((j >= i && j >= k) ? j : k);
        int middle = (i != small && i != big) ? i : ((j != small && j != big) ? j : k);

        // printf("\nEEntering 3-opt: i = %d, j = %d, k = %d", i, j, k);

        apply_3_opt(tour, small, middle, big);

        long current_cost = cal_fitness(path);
        long new_cost = cal_fitness(tour);
        // printf("\nAAAAAA");
        if (new_cost < current_cost)
            memcpy(path, tour, D * sizeof(int));

        else
        {
            double p = exp((current_cost - new_cost) / temp);
            if ((double)rand() / RAND_MAX < p)
                memcpy(path, tour, D * sizeof(int));

            // else
            //     memcpy(tour, path, D * sizeof(int));
        }

        temp *= cooling_rate;
    }
}

// ===================================

// PSO function
void pso()
{
    // Read Matrix from file
    matrix = read_matrix(data_file);

    // Initialization
    Particle particles[SWARM];
    Particle globalBest;
    globalBest.bestFit = __INT64_MAX__;

    // Rank the Particles and Find Local & Global Best
    for (int i = 0; i < SWARM; i++)
    {
        initializeParticle(&particles[i]);
        rankPrticles(&particles[i], &globalBest);
    }

    // Main loop
    for (int iter = 0; iter < MAX_ITER; iter++)
    // int iter = 0;
    // do
    {
#pragma omp parallel for // Using OpenMP for parallel computation
        for (int i = 0; i < SWARM; i++)
        {
            // Update particle's velocity
            updateVelocity(&particles[i], &globalBest);

            // Update particle's position
            updatePosition(&particles[i]);

            // Rank Particles
            rankPrticles(&particles[i], &globalBest);
        }
        printf("Itr = %d, Best = %ld\n", iter, globalBest.bestFit);
        // iter++;
#pragma omp barrier
    }
    // }while (globalBest.bestFit != 39);

    // Print the best solution found
    printf("\nBest solution found with fitness: %ld\n\n\t", globalBest.bestFit);
    for (int i = 0; i < D; i++)
    {
        printf("%d ", globalBest.bestPosition[i]);
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

// gcc main.c -fopenmp -lm -o code
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <omp.h>
#include <assert.h>
#include <string.h>

#define SWARM 10 // Number of ants
#define D 170    // Number of cities
#define MAX_ITER 10000
#define ALPHA 1 // Importance of pheromone trail
#define BETA 5  // Importance of distance heuristic
#define RHO 0.1 // Evaporation rate
#define Q 100   // Amount of pheromone deposited

// File Paths
const char *data_file = "datasets/ftv170.txt";
const char *sol_file = "solutions/ftv170-aco.txt";

// Structure to hold a city's coordinates
typedef struct
{
    int id;
    double x, y;
} City;

// Structure to hold an ant's tour and fitness
typedef struct
{
    int tour[D];
    double fitness;
} Ant;

// Global variables
City **cities; // Array of cities
double pheromone[D][D];

//======== Function Declarations ===========
/*--------------------------------------*/

void initializeCity(City *);
void initializePheromone();
void initializeAnts(Ant *);
void updatePheromone(Ant *);
void printBestTour(Ant *);
double calDistance(City, City);
double calFitness(Ant *);
void aco();

/*---------------MAIN FUNCTION------------------*/
int main(int argc, char const *argv[])
{
    srand(1);
    aco();
    return 0;
}

// Function to initialize a city's coordinates
void initializeCity(City *city)
{
    city->id = -1;
    city->x = 0;
    city->y = 0;
}

// Function to initialize the pheromone trails
void initializePheromone()
{
    int i, j;
    for (i = 0; i < D; i++)
    {
        for (j = 0; j < D; j++)
        {
            pheromone[i][j] = 0.01; // Small positive value
        }
    }
}

// Function to initialize an ant's tour
void initializeAnts(Ant *ants)
{
    int i, j;
    for (i = 0; i < SWARM; i++)
    {
        ants[i].fitness = 0;
        for (j = 0; j < D; j++)
        {
            ants[i].tour[j] = j;
        }
        ants[i].tour[D - 1] = -1; // Mark the end of the tour
    }
}

// Function to update the pheromone trails
void updatePheromone(Ant *bestAnt)
{
    int i;
    for (i = 0; i < D * D; i++)
    {
        pheromone[bestAnt->tour[i % D]][bestAnt->tour[(i + 1) % D]] += Q / bestAnt->fitness;
    }
}

// Function to print the best tour found
void printBestTour(Ant *bestAnt)
{
    int i;
    printf("\nBest solution found with fitness: %f\n\n\t", bestAnt->fitness);
    for (i = 0; i < D; i++)
    {
        printf("%d ", bestAnt->tour[i]);
    }
    printf("\n");
}

// Function to calculate the distance between two cities
double calDistance(City city1, City city2)
{
    return sqrt(pow(city1.x - city2.x, 2) + pow(city1.y - city2.y, 2));
}

// Function to calculate an ant's fitness
double calFitness(Ant *ant)
{
    double fitness = 0;
    int i;
    for (i = 0; ant->tour[i] != -1; i++)
    {
        if (i != 0)
        {
            fitness += pheromone[ant->tour[i - 1]][ant->tour[i]] + calDistance(cities[ant->tour[i - 1]], cities[ant->tour[i]]);
        }
    }
    return 1 / fitness;
}

// Function to implement the ACO algorithm
void aco()
{
    // Read cities from file
    int i;
    FILE *file = fopen(data_file, "r");
    if (file == NULL)
    {
        printf("Error opening file\n");
        return;
    }
    cities = (City **)malloc(D * sizeof(City *));
    for (i = 0; i < D; i++)
    {
        cities[i] = (City *)malloc(sizeof(City));
        fscanf(file, "%d %lf %lf\n", &(cities[i]->id), &(cities[i]->x), &(cities[i]->y));
    }
    fclose(file);

    // Initialization
    Ant ants[SWARM];
    Ant bestAnt;
    initializePheromone();
    initializeAnts(ants);
    bestAnt.fitness = 0;

    // Main loop
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
#pragma omp parallel for private(i)
        for (i = 0; i < SWARM; i++)
        {
            double prob[D];
            double sum = 0;
            int j;
            for (j = 0; j < D; j++)
            {
                if (ants[i].tour[j] == -1)
                    break;
                prob[j] = pow(pheromone[ants[i].tour[j]][ants[i].tour[(j + 1) % D]], ALPHA) * pow(1 / calDistance(cities[ants[i].tour[j]], cities[ants[i].tour[(j + 1) % D]]), BETA);
                sum += prob[j];
            }
            for (j = 0; j < D; j++)
            {
                if (ants[i].tour[j] == -1)
                    break;
                prob[j] /= sum;
            }
            double r = (double)rand() / RAND_MAX;
            int prev = -1;
            for (j = 0; j < D; j++)
            {
                if (ants[i].tour[j] == -1)
                    break;
                r -= prob[j];
                if (r <= 0)
                {
                    prev = j;
                    break;
                }
            }
            int next = -1;
            for (j = 0; j < D; j++)
            {
                if (ants[i].tour[j] == -1)
                    break;
                if (j == prev)
                    continue;
                if (next == -1)
                {
                    next = j;
                    break;
                }
            }
            ants[i].tour[prev + 1] = ants[i].tour[next];
            ants[i].tour[next] = -1;
            ants[i].fitness = calFitness(&ants[i]);
        }

        // Find the best ant
        for (i = 0; i < SWARM; i++)
        {
            if (ants[i].fitness > bestAnt.fitness)
            {
                bestAnt = ants[i];
            }
        }

        // Update pheromone trails
        updatePheromone(&bestAnt);

        // Evaporate pheromone trails
        for (i = 0; i < D; i++)
        {
            for (int j = 0; j < D; j++)
            {
                pheromone[i][j] *= (1 - RHO);
            }
        }

// Print after each iteration
#pragma omp single
        {
            printf("Itr = %d, Best = %f\n", iter, bestAnt.fitness);
        }
    }

    // Print the best solution found
    printBestTour(&bestAnt);

    // Free memory
    for (i = 0; i < D; i++)
    {
        free(cities[i]);
    }
    free(cities);
}
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <omp.h>
#include <assert.h>
#include <string.h>

#define SWARM 10 // Number of particles
#define D 34    // Number of dimensions
#define MAX_ITER 10000
// #define W 0.5              // Inertia Weight
#define C1 2               // Acceleration Factor
#define C2 2               // Acceleration Factor
#define TEMP 100           // Initial Temperature
#define BEST_SOLUTION 1286 // Best Solution

#define W_MAX 0.9
#define W_MIN 0.4
int W;

#define NUM_THREADS 6 // Maximum Threads

// File Paths
const char *data_file = "/content/TSP/ftv33.txt";
const char *sol_file = "/content/TSP/solution.txt";

// Structure to hold a particle's position and velocity
typedef struct
{
    int position[D];
    double velocity[D];
    int bestPosition[D];
    long bestFit;
} Particle;

// Define a structure to pass arguments to the thread function
struct ThreadArgs
{
    Particle *particles;
    Particle *globalBest;
    int startIndex;
    int endIndex;
};

// Global variables
int **matrix; // Cost Matrix

// Ant Colony Variables
double tou[D][D];
float alpha = 2.0, row = 0.0001, ita = 8.0;

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
void Simul_Annel_3_Opt(int *, double, double);
void InverseMutation(int *);
void swap(int *, int, int);

void *updateParticles(void *arg);

void initialize_tou();
void evaporation();
void pheromore_update(Particle *particle);
void construct_path(int *PATH);

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
    int i, j;
    for (i = 0; i < D; i++)
        matrix[i] = (int *)malloc(D * sizeof(int));

    // Reset file pointer to beginning
    fseek(file, 0L, SEEK_SET);

    // Read the matrix values
    for (i = 0; i < D; i++)
    {
        for (j = 0; j < D; j++)
            fscanf(file, "%d", &matrix[i][j]);
    }

    fclose(file);

    // Ant Colony

    // initialize_tou();

    int h, l;
    for (h = 0; h < D; h++)
        for (l = 0; l < D; l++)
        {
            // printf(" %d", matrix[h][l]);
            tou[h][l] = 1.0 / matrix[h][l];
        }

    return matrix;
}

// Function to initialize a particle's position and velocity
void initializeParticle(Particle *particle)
{
    bool pos[D + 1] = {false}; // Initially make all false
    int count = 0, i;
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

    for (i = 0; i < D; i++)
    {
        particle->velocity[i] = ((double)rand() / (double)(RAND_MAX)) * (1 - 0) + 0; // Between 0 and 1
        particle->bestPosition[i] = particle->position[i];
    }

    particle->bestFit = __INT64_MAX__;
}

// Function to Rank Particles and Best Solutions
void rankPrticles(Particle *particle, Particle *global)
{
    long fit = cal_fitness(particle->position);
    int j;
    if (fit < particle->bestFit)
    {
        particle->bestFit = fit;
        for (j = 0; j < D; j++)
        {
            particle->bestPosition[j] = particle->position[j];
        }
    }
    if (fit < global->bestFit)
    {
        FILE *solf = fopen(sol_file, "w");

        global->bestFit = fit;

        // I/O Section
        printf("Best = %ld\n", fit);

        // Write Best Found Solution into Disk
        fprintf(solf, "\n\n\t\t\t\t----* %d x %d Matrix *----\n\n\n", D, D);
        fprintf(solf, "\t  Best Fitness : %6ld | Optimum : %4d \n\n\nPATH:\n\n", fit, BEST_SOLUTION);

        for (j = 0; j < D; j++)
        {
            global->bestPosition[j] = particle->position[j];
            fprintf(solf, "%d  ", particle->position[j]);
        }

        fprintf(solf, "\n\n\n================================\n\n\n");

        // for (j = 0; j < SWARM; j++)
        //     fprintf(solf, "\n  BAT %d : \t%ld  ", j + 1, x[j].fit);

        fclose(solf);
    }
}

void updateVelocity(Particle *particle, Particle *globalBest)
{
    int j;
    for (j = 0; j < D; j++)
    {
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        float vel = W * particle->velocity[j] + C1 * r1 * (particle->bestPosition[j] - particle->position[j]) + C2 * r2 * (globalBest->bestPosition[j] - particle->position[j]);
        int i_vel = (int)vel;
        particle->velocity[j] = vel - i_vel;
        // printf("\nnew velocity = %.4f", particle->velocity[j]);
    }
}

// Optimizing function

void *updateParticles(void *arg)
{
    struct ThreadArgs *args = (struct ThreadArgs *)arg;
    int i;
    for (i = args->startIndex; i < args->endIndex; i++)
    {
        // Update particle's velocity
        updateVelocity(&(args->particles[i]), args->globalBest);

        // Update particle's position
        updatePosition(&(args->particles[i]));

        // Rank Particles
        rankPrticles(&(args->particles[i]), args->globalBest);
    }
    pthread_exit(NULL);
}

void updatePosition(Particle *Particle)
{
    float max_vel = 0;
    int k;
    for (k = 0; k < D; k++)
    {
        max_vel = (Particle->velocity[k] > max_vel) ? Particle->velocity[k] : max_vel;
    }

    // for (int i = 0; i < D; i++)
    // {
    //     for (int j = i + 1; j < D; j++)
    //     {
    //         if (Particle->position[i] == Particle->position[j])
    //         {
    //             printf("\n\n *** ERRROR ***\n");
    //         }
    //     }
    // }

    // printf("\n\n");
    // for (int i = 0; i < D; i++)
    // {
    //     printf("%d\t", Particle->position[i]);
    // }
    // printf("\n\n");

    Simul_Annel_3_Opt(Particle->position, TEMP, max_vel);

    long new_fit = cal_fitness(Particle->position);

    if (Particle->bestFit < new_fit)
    {
       
        InverseMutation(Particle->position);
    }

    new_fit = cal_fitness(Particle->position);

    if (Particle->bestFit < new_fit)
    {
         construct_path(Particle->position);
    }
}

// Function to evaluate a particle's fitness
long cal_fitness(int pos[])
{
    long fit = 0;
    int i;
    for (i = 0; i < D - 1; i++)
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

        // ==== APPLY 3 OPT ALGORITHM ====
        // apply_3_opt(tour, small, middle, big);
        reverse(tour, small + 1, middle);
        reverse(tour, middle + 1, big);
        // ==== APPLY 3 OPT ALGORITHM ====

        long current_cost = cal_fitness(path);
        long new_cost = cal_fitness(tour);
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
    int i;
    // Initialization
    Particle particles[SWARM];
    Particle globalBest;
    globalBest.bestFit = __INT64_MAX__;

    // Rank the Particles and Find Local & Global Best
    for (i = 0; i < SWARM; i++)
    {
        initializeParticle(&particles[i]);
        rankPrticles(&particles[i], &globalBest);
    }

    // Main loop
    for (long long iter = 0; iter < MAX_ITER; iter++)

    /**/
    // long iter = 0;
    // do
    {
        #pragma omp parallel for private(i)
        for (i = 0; i < SWARM; i++)
        {
            // Calculate W
            W = W_MAX - (W_MAX - W_MIN) * ( iter / MAX_ITER );
            // Update particle's velocity
            updateVelocity(&particles[i], &globalBest);

            // Update particle's position
            updatePosition(&particles[i]);

            // Rank Particles
            rankPrticles(&particles[i], &globalBest);
        }

        // Ant Colony
        evaporation();

        for (i = 0; i < SWARM; i++)
            pheromore_update(&particles[i]);

    }
    // } while (globalBest.bestFit != BEST_SOLUTION);

    // pthread_t threads[NUM_THREADS];
    // struct ThreadArgs threadArgs[NUM_THREADS];

    // do
    // {
    //     // Create threads
    //     for (i = 0; i < NUM_THREADS; i++)
    //     {
    //         threadArgs[i].particles = particles;
    //         threadArgs[i].globalBest = &globalBest;
    //         threadArgs[i].startIndex = i * (SWARM / NUM_THREADS);
    //         threadArgs[i].endIndex = (i + 1) * (SWARM / NUM_THREADS);
    //         pthread_create(&threads[i], NULL, updateParticles, (void *)&threadArgs[i]);
    //     }

    //     // Join threads
    //     for (i = 0; i < NUM_THREADS; i++)
    //     {
    //         pthread_join(threads[i], NULL);
    //     }

    //     // Print after each iteration
    //     printf("Itr = %d, Best = %ld\n", iter, globalBest.bestFit);
    //     iter++;
    // } while (globalBest.bestFit != BEST_SOLUTION);

    // Print the best solution found
    printf("\nBest solution found with fitness: %ld\n\n\t", globalBest.bestFit);
    for (i = 0; i < D; i++)
    {
        printf("%d ", globalBest.bestPosition[i]);
    }
    printf("\n");
}

// ===================================

// void initialize_tou()
// {
//     int i, j;
//     for (i = 0; i < D; i++)
//         for (j = 0; j < D; j++) {
//             printf(" %d", matrix[i][j]);
//             tou[i][j] = 1.0 / matrix[i][j];
//         }
// }

void construct_path(int *PATH)
{
    int i, j, l, n_list[D], count = 0;
    float p[D], sum = 0, val;
start:
    for (i = 0; i < D - 1; i++)
        n_list[i] = i + 1;
    PATH[0] = 0;
    for (i = 1; i < D; i++)
    {
        sum = 0;
        for (j = 0; j < D - i; j++)
            sum += pow(tou[PATH[i - 1]][n_list[j]], alpha);
        for (j = 0; j < D - i; j++)
            p[j] = pow(tou[PATH[i - 1]][n_list[j]], alpha) / sum;
        sum = 0;
        for (j = 0; j < D - i; j++)
            sum += p[j];
        for (j = 1; j < D - i; j++)
            p[j] = p[j] + p[j - 1];
        val = (rand() % 1000) / 1000.;
        for (j = 0; j < D - i; j++)
        {
            if (val < p[j])
            {
                PATH[i] = n_list[j];
                for (l = j; l < D - i; l++)
                    n_list[l] = n_list[l + 1];
                goto end;
            }
        }
    end:;
    }
}

void evaporation()
{
    int i, j;
    for (i = 0; i < D; i++)
        for (j = 0; j < D; j++)
            tou[i][j] = (1 - row) * tou[i][j];
}

void pheromore_update(Particle *particle)
{
    int i;
    long val = cal_fitness(particle->position);
    for (i = 0; i < D - 1; i++)
        tou[particle->position[i]][particle->position[i + 1]] += 1. / val;
    tou[particle->position[D - 1]][particle->position[0]] += 1. / val;
}

/*

1. Initialize position and velocity
2. Calculate fitnesses and global_best
3. Update velocity and position
4. Calculate fitness and find global_best
5. Increase iteration

*/

// gcc main.c -fopenmp -lm -o code
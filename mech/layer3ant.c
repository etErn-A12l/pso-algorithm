#include <assert.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SWARM 50
#define D 71 // Number of dimensions

#define BEST_SOLUTION 1286

int count = 0;

typedef struct
{
  int path[D];
  long fit;
} BOM;

BOM best;

#define NUM_THREADS 6 // Maximum Threads

// File Paths
const char *data_file = "/home/eternal/Documents/pso-algorithm/datasets/ftv70.txt";
const char *sol_file = "/home/eternal/Documents/pso-algorithm/mech/solution/custom70.txt";

// Global variables
int **matrix; // Cost Matrix

// Ant Colony Variables
double tou[D][D];
float alpha = 2.0, row = 0.0001, ita = 8.0;

int **read_matrix(const char *);
long cal_fitness(int[]);
void evaporation();
void pheromore_update(int *);
void construct_path(int *);
void initializePath(int *);
void rankPrticles(BOM *, BOM *);
void runAlgorithm(BOM[], BOM *, int);

int compareBOM(const void *, const void *);
void sortBOMArray(BOM *, size_t);

int main(int argc, char const *argv[])
{
  matrix = read_matrix(data_file);

  int i, j;

  BOM S1[SWARM];
  BOM S2[SWARM];

  best.fit = __INT64_MAX__; // Set Best Fitness to Max

  printf("\n======================= Made By 0000 ========================\n\n\n");

  // Initialize Paths
  #pragma omp parallel for private(i)
  for (i = 0; i < SWARM; i++)
  {
    initializePath(S1[i].path);
    S1[i].fit = __INT64_MAX__;
    rankPrticles(&S1[i], &best);
  }

  runAlgorithm(S1, &best, 50000);

  // printf("======================================================\n");

  #pragma omp parallel for private(i)
  for (i = 0; i < SWARM; i++)
  {
    initializePath(S2[i].path);
    S2[i].fit = __INT64_MAX__;
    rankPrticles(&S2[i], &best);
  }

  // for (i = 0; i < SWARM; i++)
  //   printf("S2[%d].fit = %ld\n", i, S2[i].fit);

  runAlgorithm(S2, &best, 50000);

  sortBOMArray(S1, SWARM);
  sortBOMArray(S2, SWARM);

  BOM FINAL[SWARM];

  printf("\n\n===================== COMPUTING FINAL =======================\n\n");

  for (i = 0; i < SWARM; i += 2)
  {
    FINAL[i] = S1[i];
    FINAL[i + 1] = S2[i];
  }

  runAlgorithm(FINAL, &best, 100000);

  sortBOMArray(FINAL, SWARM);

  printf("\n======================================================\n");
  printf("======================================================\n\n\n");
  for (i = 0; i < SWARM; i++)
    printf("FINAL [%d] | -> %ld <- |\n", i, FINAL[i].fit);

  return 0;
}

void runAlgorithm(BOM bom[], BOM *best, int MAX_ITER)
{
  int i, j, iter = 0;

  for (iter = 0; iter < MAX_ITER; iter++)
  {
    // Update paths
    #pragma omp parallel for private(i)
    for (i = 0; i < SWARM; i++)
    {
      // Update path
      construct_path(bom[i].path);
      // Rank Paths
      rankPrticles(&bom[i], best);
    }

    // Update Pheromore
    #pragma omp parallel for private(i)
    for (i = 0; i < SWARM; i++)
      pheromore_update(bom[i].path);
  }
}

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

void pheromore_update(int *path)
{
  int i;
  long val = cal_fitness(path);
  for (i = 0; i < D - 1; i++)
    tou[path[i]][path[i + 1]] += 1 / val;
  tou[path[D - 1]][path[0]] += 1 / val;
}

void initializePath(int *path)
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
      path[count] = num;
      pos[num] = true;
      count++;
    }
  }
}

int **read_matrix(const char *data_file)
{
  FILE *file = fopen(data_file, "r");
  char line[10000];

  if (file == NULL)
  {
    printf("Error opening file\n");
    return NULL;
  }

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

  int h, l;
  for (h = 0; h < D; h++)
    for (l = 0; l < D; l++)
    {
      // printf(" %d", matrix[h][l]);
      tou[h][l] = 1.0 / matrix[h][l];
    }

  return matrix;
}

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

void rankPrticles(BOM *bom, BOM *best)
{
  long fit = cal_fitness(bom->path);
  int j;
  if (fit < bom->fit)
    bom->fit = fit;

  if (fit < best->fit)
  {
    FILE *solf = fopen(sol_file, "w");

    best->fit = fit;

    // I/O Section
    printf("\r| FITNESS -> %ld <- |", fit);
    fflush(stdout); // Flush the output buffer

    // Write Best Found Solution into Disk
    fprintf(solf, "\n\n\t\t\t\t----* %d x %d Matrix *----\n\n\n", D, D);
    fprintf(solf, "\t  Best Fitness : %6ld | Optimum : %4d \n\n\nPATH:\n\n",
            fit, BEST_SOLUTION);

    for (j = 0; j < D; j++)
    {
      best->path[j] = bom->path[j];
      fprintf(solf, "%d  ", bom->path[j]);
    }

    fprintf(solf, "\n\n\n================================\n\n\n");

    fclose(solf);
  }
}

// Utility Functions

int compareBOM(const void *a, const void *b)
{
  BOM *bomA = (BOM *)a;
  BOM *bomB = (BOM *)b;
  if (bomA->fit < bomB->fit)
    return -1;
  if (bomA->fit > bomB->fit)
    return 1;
  return 0;
}

// Function to sort an array of BOM structs
void sortBOMArray(BOM *array, size_t size)
{
  qsort(array, size, sizeof(BOM), compareBOM);
}
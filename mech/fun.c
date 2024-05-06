#include <stdio.h>
#include <stdbool.h> // Include this header for bool type
#include <stdlib.h>

int cost[10][10];

void initialize_cost_matrix()
{
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            if (i == j)

                cost[i][j] = 9999;

            else
                cost[i][j] = (rand() % 20) + 1;
        }
    }
}

/**
 * Checks if a given number exists in an array.
 *
 * @param arr The array to search in.
 * @param size The size of the array.
 * @param num The number to search for.
 * @return true if the number is found in the array, false otherwise.
 */
bool doesNumberExist(int arr[], int size, int num)
{
    for (int i = 0; i < size; i++)
    {
        if (arr[i] == num)
        {
            return true; // Number found
        }
    }
    return false; // Number not found
}

void reArrange(int *arr, int size, int n)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (arr[i] == n)
        {
            while (i > 0)
            {
                arr[i] = arr[i - 1];
                i--;
            }
            arr[0] = n;
            return;
        }
    }
}

void mechs_crossover(int *p1, int *p2, int *child, int size)
{
    int i, j, u;

    // Choose a random city v
    int v = rand() % size;
    child[0] = v;
    int childIndex = 1;

    // Initialize child array with -1
    for (i = 1; i < size; i++)
        child[i] = -1;

    int *parent1 = p1;
    int *parent2 = p2;

    reArrange(parent1, size, v);
    reArrange(parent2, size, v);

    // printf("\nRearranged parent1:\n");
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%d ", parent1[i]);
    // }

    // printf("\nRearranged parent2:\n");
    // for (int i = 0; i < 10; i++)
    // {
    //     printf("%d ", parent2[i]);
    // }

    i = 1;
    j = 1;

    while (i < size && j < size)
    {
        if (doesNumberExist(child, size, parent2[j]) && doesNumberExist(child, size, parent1[i]))
        {
            i++;
            j++;
        }
        else if (doesNumberExist(child, size, parent1[i]))
        {
            child[childIndex] = parent2[j];
            j++;
            childIndex++;
        }
        else if (doesNumberExist(child, size, parent2[j]))
        {
            child[childIndex] = parent1[i];
            i++;
            childIndex++;
        }
        else
        {
            u = child[childIndex - 1];
            if (cost[u, parent1[i]] < cost[u, parent2[j]])
            {
                child[childIndex] = parent1[i];
                i++;
            }
            else
            {
                child[childIndex] = parent2[j];
                j++;
            }
            childIndex++;
        }
    }
}

void my_nech(int *p1, int *p2, int *child, int size)
{
    int i, j, u;
    // Choose a random city v
    int v = rand() % size;
    child[0] = v;
    int childIndex = 1;

    // Initialize child array with -1
    for (i = 1; i < size; i++)
        child[i] = -1;

    int *parent1 = p1;
    int *parent2 = p2;

    reArrange(parent1, size, v);
    reArrange(parent2, size, v);

    i = 1;
    j = 1;

    while (childIndex < size)
    {
        if (i < size && j < size)
        {
            if (doesNumberExist(child, size, parent2[j]) && doesNumberExist(child, size, parent1[i]))
            {
                i++;
                j++;
            }
            else if (doesNumberExist(child, size, parent1[i]))
                child[childIndex++] = parent2[j++];
            else if (doesNumberExist(child, size, parent2[j]))
                child[childIndex++] = parent1[i++];
            else
            {
                u = child[childIndex - 1];
                if (cost[u, parent1[i]] < cost[u, parent2[j]])
                    child[childIndex++] = parent1[i++];
                else
                    child[childIndex++] = parent2[j++];
            }
        }
        else if (i < size && j == size)
        {
            while (i < size)
            {
                if (doesNumberExist(child, size, parent1[i]))
                    i++;
                else
                    child[childIndex++] = parent1[i++];
            }
        }
        else if (i == size && j < size)
        {
            while (j < size)
            {
                if (doesNumberExist(child, size, parent2[j]))
                    j++;
                else
                    child[childIndex++] = parent2[j++];
            }
        }
        else
            break;
    }
}

long cal_fitness(int pos[])
{
    int i;
    long fit = 0;
    for (i = 0; i < 10 - 1; i++)
    {
        if (pos[i] != pos[i + 1])
            fit += cost[pos[i]][pos[i + 1]];
    }
    fit += cost[pos[9]][pos[0]];
    return fit;
}

int main()
{
    int parent1[] = {0, 2, 3, 4, 1, 6, 5, 9, 7, 8};
    int parent2[] = {2, 3, 6, 5, 0, 1, 4, 9, 8, 7};
    int child[10];

    initialize_cost_matrix();

    my_nech(parent1, parent2, child, 10);

    printf("\nParent1 Fitness: %ld\n", cal_fitness(parent1));
    printf("\nParent2 Fitness: %ld\n", cal_fitness(parent2));

    printf("\nChild:\n");

    for (int i = 0; i < 10; i++)
    {
        printf("%d ", child[i]);
    }

    printf("\n");

    printf("\nChild Fitness: %ld\n", cal_fitness(child));

    return 0;
}
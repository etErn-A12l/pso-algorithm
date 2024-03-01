/* Program for Simple ACO Crisp cost and Constraints*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define node 70
#define nk 20
#define nit 4000

void initialize_tou();
void construct_path(int i);
void evaporation();
void pheromore_update();
float cost(int k);
int cheak_constraint(int k);
float randval(float, float);
void createfile();
float rnd();

float TOU[node][node], COST[node][node], safety[node][node];

int PATH[nk][node], optimum[node], opti_cost = 50000;
float alpha = 2.0, row = 0.0001, ita = 8.0;
FILE *fptr, *fptr1;
int main()
{
	int seed, i, k, j, count = 0;
	float val;
	float a, b;
	printf("\n Please enter seed:");
	scanf("%d", &seed);
	srand(seed);
	createfile();
	initialize_tou();
	for (j = 0; j < nit; j++)
	{
		for (i = 0; i < nk; i++)
			construct_path(i);
		evaporation();
		pheromore_update();
	}
	fptr = fopen("result", "w");
	for (k = 0; k < nk; k++)
	{
		val = 0.0;
		for (i = 0; i < node - 1; i++)
			val += COST[PATH[k][i]][PATH[k][i + 1]];
		val += COST[PATH[k][node - 1]][PATH[k][0]];

		for (i = 0; i < node; i++)
		{
			printf(" %d ", PATH[k][i]);
			fprintf(fptr, " %d ", PATH[k][i]);
		}
		printf(" %f ", val);
		fprintf(fptr, " %f ", val);
		fprintf(fptr, "\n");
		printf("\n");
	}
	fprintf(fptr, " optimum path is given below\n ");
	for (i = 0; i < node; i++)
		fprintf(fptr, " %d ", optimum[i]);
	fprintf(fptr, " %d ", opti_cost);

	fclose(fptr);
	// exit(0);
	return 0;
}

void initialize_tou()
{
	int i, j;
	for (i = 0; i < node; i++)
		for (j = 0; j < node; j++)
			TOU[i][j] = 1.0 / COST[i][j];
}
void construct_path(int k)
{
	int i, j, l, n_list[node], count = 0;
	float p[node], sum = 0, val;
start:
	for (i = 0; i < node - 1; i++)
		n_list[i] = i + 1;
	PATH[k][0] = 0;
	for (i = 1; i < node; i++)
	{
		sum = 0;
		for (j = 0; j < node - i; j++)
			sum += pow(TOU[PATH[k][i - 1]][n_list[j]], alpha);
		for (j = 0; j < node - i; j++)
			p[j] = pow(TOU[PATH[k][i - 1]][n_list[j]], alpha) / sum;
		sum = 0;
		for (j = 0; j < node - i; j++)
			sum += p[j];
		for (j = 1; j < node - i; j++)
			p[j] = p[j] + p[j - 1];
		val = rnd();
		for (j = 0; j < node - i; j++)
		{
			if (val < p[j])
			{
				PATH[k][i] = n_list[j];
				for (l = j; l < node - i; l++)
					n_list[l] = n_list[l + 1];
				goto end;
			}
		}
	end:;
	}
	if (cheak_constraint(k) == 0)
	{
		if (count++ < 200)
			goto start;
		else
		{
			printf("\n Exit from path Construct");
			exit(0);
		}
	}
}
int cheak_constraint(int k)
{
	/*	int x[node],i;
		float val=0.0;
	  for(i=0;i<node;i++)
		 x[i]=PATH[k][i];
	 for(i=0;i<node-1;i++)
		val+=safety[x[i]][x[i+1]];
	 val+=safety[x[node-1]][x[0]];
	 if(val>ita)
		return 1;
	  else
		return 0;*/
	return 1;
}
float rnd()
{
	float val;
	val = (rand() % 1000) / 1000.;
	return val;
}
void evaporation()
{
	int i, j;
	for (i = 0; i < node; i++)
		for (j = 0; j < node; j++)
			TOU[i][j] = (1 - row) * TOU[i][j];
}
void pheromore_update()
{
	int i, k;
	float val;
	for (k = 0; k < nk; k++)
	{
		val = cost(k);
		for (i = 0; i < node - 1; i++)
			TOU[PATH[k][i]][PATH[k][i + 1]] += 1. / val;
		TOU[PATH[k][node - 1]][PATH[k][0]] += 1. / val;
	}
}
float randval(float a, float b)
{
	return a + ((rand() % 1000) / 1000.) * (b - a);
}
float cost(int k)
{
	int i;
	float val = 0;
	for (i = 0; i < node - 1; i++)
		val += COST[PATH[k][i]][PATH[k][i + 1]];
	val += COST[PATH[k][node - 1]][PATH[k][0]];
	if (opti_cost > val)
	{
		for (i = 0; i < node - 1; i++)
			optimum[i] = PATH[k][i];
		opti_cost = val;
	}
	return val;
}

void createfile()
{
	int a = 0, j;
	float ta[70], t;
	int i;
	int count = 0;

	fptr = fopen("e:\\ft701.txt", "r");
	fptr1 = fopen("e:\\ft701read.txt", "w");
	while (fscanf(fptr, "%d ", &i) != EOF)
	{
		printf("%d\t", i);
		fprintf(fptr1, "%d\t", i);
		t = (float)i;
		ta[count] = t;
		count++;

		if (count == 70)
		{
			fprintf(fptr1, "\n");
			count = 0;
			for (j = 0; j < 70; j++)
				COST[a][j] = ta[j];
			a++;
		}
	}

	fclose(fptr);
	fclose(fptr1);
}

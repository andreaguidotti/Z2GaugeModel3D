#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<time.h>

#include "../include/random.h"
#include "../include/geometry.h"

void init_lattice(int ***lattice, long int volume, int DIM)
{
    int link;
    *lattice = (int **)malloc((unsigned long int)volume * sizeof(int *));
    if (!*lattice)
    {
        fprintf(stderr, "malloc failed for lattice");
        exit(EXIT_FAILURE);
    }
    for (long int i = 0; i < volume; i++)
    {
        (*lattice)[i] = (int *)malloc((unsigned int)DIM * sizeof(int));
        for (int ii = 0; ii < DIM; ii++)
        {
            link = 2 * (int)(2 * myrand()) - 1;
            (*lattice)[i][ii] = link;
        }
    }
}

int main()
{
    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;

    myrand_init(seed1, seed2);

    int **Lattice;
    long int volume = 10;
    int DIM = 3;

    init_lattice(&Lattice, volume, DIM);
    for (long int  i = 0; i < volume; i++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            printf("%d\n", Lattice[i][ii]);
        }
        
    }
 
}
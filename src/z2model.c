#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/random.h"
#include "../include/geometry.h"

#define DEBUG 1

int DIM = 3;

#if DEBUG
int ****debugLattice = NULL;
int cartcoord[3];
#endif

void init_lattice(int ***lattice, int size, int dim)
{
#if DEBUG
    debugLattice = malloc((unsigned int)size * sizeof(int ***));
    for (int x = 0; x < size; x++)
    {
        debugLattice[x] = malloc((unsigned int)size * sizeof(int **));
        for (int y = 0; y < size; y++)
        {
            debugLattice[x][y] = malloc((unsigned int)size * sizeof(int *));
            for (int z = 0; z < size; z++)
            {
                debugLattice[x][y][z] = malloc((unsigned int)DIM * sizeof(int));
            }
        }
    }
#endif
    int link;
    long int volume = 1;
    for (int i = 0; i < dim; i++)
    {
        volume *= size;
    }
    *lattice = (int **)malloc((unsigned long int)volume * sizeof(int *));
    if (!*lattice)
    {
        fprintf(stderr, "malloc failed for lattice");
        exit(EXIT_FAILURE);
    }
    for (long int i = 0; i < volume; i++)
    {
        (*lattice)[i] = (int *)malloc((unsigned int)dim * sizeof(int));
        if (!(*lattice)[i])
        {
            fprintf(stderr, "malloc failed for lattice[%li]", i);
            exit(EXIT_FAILURE);
        }
        for (int ii = 0; ii < dim; ii++)
        {
            link = 2 * (int)(2 * myrand()) - 1;
        #if DEBUG
            lex_to_cart(i, cartcoord, dim, size);
            debugLattice[cartcoord[0]][cartcoord[1]][cartcoord[2]][ii] = link;
        #endif
            (*lattice)[i][ii] = link;
        }
    }
}
void init_neighbours(long int **nnp, long int **nnm, int size, int dim)
{
    long int volume = 1;
    for (int i = 0; i < dim; i++)
    {
        volume *= size;
    }
    (*nnp) = (long int *)malloc((unsigned long int)(dim * volume) * sizeof(long int));
    if (!(*nnp))
    {
        fprintf(stderr, "failed to initialize neighbours at (%s, %d)", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    (*nnm) = (long int *)malloc((unsigned long int)(dim * volume) * sizeof(long int));
    if (!(*nnm))
    {
        fprintf(stderr, "failed to initialize neighbours at (%s, %d)", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    init_geo(*nnp, *nnm, size, dim);
}
#if DEBUG
void check_lattice(int **Lattice, int ****debugLattice, int size, int dim)
{
    long int volume = 1;
    int linkLattice;
    int linkDebug;
    for (int i = 0; i < dim; i++)
    {
        volume *= size;
    }
    for (long int lex = 0; lex < volume; lex++)
    {
        for (int dir = 0; dir < dim; dir++)
        {
            lex_to_cart(lex, cartcoord, dim, size);

            linkLattice = Lattice[lex][dir];
            linkDebug = debugLattice[cartcoord[0]][cartcoord[1]][cartcoord[2]][dir];

            if (linkLattice != linkDebug)
            {
                printf("mismatch between Lattice and debugLattice");
                exit(EXIT_FAILURE);
            }
        }
    }
    printf("Check Lattice passed\n");
}

void check_neighbours(long int *nnp, long int *nnm, int size, int dim)
{
    int cartCheck[3];

    long int forward;
    long int backward;

    long int volume = 1;

    for (int i = 0; i < dim; i++)
    {
        volume *= size;
    }
    for (long int lex = 0; lex < volume; lex++)
    {
        for (int dir = 0; dir < dim; dir++)
        {
            lex_to_cart(lex, cartCheck, dim, size);

            cartCheck[dir] = (cartCheck[dir] + 1) % size;
            cart_to_lex(&forward, cartCheck, dim, size);

            if (nnp[dirgeo(lex, dir, volume)] != forward)
            {
                printf("Mismatch nnp \n");
                exit(EXIT_FAILURE);
            }
            cartCheck[dir] = (cartCheck[dir] - 2 + size) % size; 
            cart_to_lex(&backward, cartCheck, dim, size);

            if (nnm[dirgeo(lex, dir, volume)] != backward)
            {
                printf("Mismatch nnm \n");
                exit(EXIT_FAILURE);
            }
        }
    }
    printf("Check neighbours passed\n");
}
#endif
int main()
{
    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;

    myrand_init(seed1, seed2);

    int **Lattice;
    long int *nnp, *nnm;

#if DEBUG

    int size = 5;

    init_lattice(&Lattice,size,DIM);
    init_neighbours(&nnp,&nnm,size,DIM);

    check_lattice(Lattice,debugLattice,size,DIM);
    check_neighbours(nnp,nnm,size,DIM);

#endif
}
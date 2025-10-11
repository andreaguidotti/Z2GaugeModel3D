#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/random.h"
#include "../include/geometry.h"

int dim = 3;

/* Initialize the lattice configuration.

   Allocates memory for a 'dim'-dimensional lattice of linear size 'size'.
   Each site contains 'dim' links initialized randomly to ±1.

   The resulting lattice is stored in a dynamically allocated array:
     lattice[lex][dir], where 'lex' is the lexicographic site index and
     'dir' labels the direction.
*/
void init_lattice(int ***lattice, int size)
{
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
            (*lattice)[i][ii] = link;
        }
    }
}
/* Initialize nearest-neighbour lookup tables.

   Allocates and fills the arrays 'nnp' and 'nnm', which store the forward
   and backward nearest-neighbour indices for each lattice site.

   These tables are used to navigate the lattice efficiently in lexicographic
   representation.
*/
void init_neighbours(long int **nnp, long int **nnm, int size)
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

/* Initialize the exponential lookup table.

   Precomputes exp(beta * deltaS) for all possible energy variations 'deltaS'
   in the Metropolis update, and stores them in 'expTable' for faster access.

   The table indices are shifted by 'deltaSmax' to allow negative deltas.
*/
void init_expTable(double ** expTable, double beta)
{
    int deltaSmax = 4 * (dim - 1);
    int range = 2 * deltaSmax + 1;

    *expTable = (double*)malloc((unsigned int)range*sizeof(double));
    for (int deltaS = - deltaSmax; deltaS <= deltaSmax; deltaS += 4)
    {
        (*expTable)[deltaS + deltaSmax] = exp(beta*deltaS);
    }
}

/* Compute the sum of staples around a selected link.

   For the link at site 'lex' in direction 'dir' this function returns the sum
   of all staples of that link, i.e each staple contributes a product of three 
   neighbouring links forming the sides of the plaquette.
*/
int computeStaples(int **restrict Lattice, long int const *restrict nnp,
                   long int const *restrict nnm, long int lex, int dir,
                   long int volume)
{
    int linkProd, sumStaples = 0;
    long int lex_minus_orth, lex_plus_dir, lex_plus_orth;

    /*                      ^ dir
                            |
                       lex_plus_dir
                   +--------+--------+
                   |        |        |
                   |        |        |
                   |        |        |    orth
                ---+--------+--------+--->
         lex_minus_orth    lex      lex_plus_orth
                            |
     */

    lex_plus_dir = nnp[dirgeo(lex, dir, volume)];

    for (int orth = 0; orth < dim; orth++)
    {
        if (orth == dir)
            continue;

        lex_plus_orth = nnp[dirgeo(lex, orth, volume)];

        // ---------- forward staple ----------
        linkProd = Lattice[lex][orth];
        linkProd *= Lattice[lex_plus_dir][orth];
        linkProd *= Lattice[lex_plus_orth][dir];

        sumStaples += linkProd;

        // ---------- backward staple ----------
        lex_minus_orth = nnm[dirgeo(lex, orth, volume)];

        linkProd = Lattice[lex_minus_orth][dir];
        linkProd *= Lattice[lex_minus_orth][orth];
        linkProd *= Lattice[nnp[dirgeo(lex_minus_orth, dir, volume)]][orth];

        sumStaples += linkProd;
    }
    return sumStaples;
}

/* Perform a single Metropolis update for a link variable.

   Attempts to flip the link Lattice[lex][dir] with probability min(1, exp(ΔS)),
   where ΔS is the local change in action computed from the surrounding staples.

   Returns 1 if the link is flipped, 0 otherwise.
*/
int metropolis(int **restrict Lattice, long int const *restrict nnp,
               long int const *restrict nnm, double * expTable,
               int lex, int dir,long int volume, double beta)
{
    int deltaS, deltaSmax, sumStaples;

    deltaSmax = 2*(dim-1);
    sumStaples = computeStaples(Lattice, nnp, nnm, lex, dir, volume);

    deltaS = -2 * Lattice[lex][dir] * sumStaples;
    if (deltaS > 0)
    {
        Lattice[lex][dir] *= -1;
        return 1;
    }
    else if (myrand() < expTable[deltaS + deltaSmax])
    {
        Lattice[lex][dir] *= -1;
        return 1;
    }
    return 0;
}
int main()
{
    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;

    myrand_init(seed1, seed2);

    int **Lattice;
    long int *nnp, *nnm;
}
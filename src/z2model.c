#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/random.h"
#include "../include/geometry.h"
#include "../include/debug.h"

#define DEBUG 0

#define dim 3
#define epsilon 0.05

#define STRING_LENGTH 50
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/* Initialize the lattice configuration.

   Allocates memory for a 'dim'-dimensional lattice of linear size 'size'.
   Each site contains 'dim' links initialized randomly to ±1.

   The resulting lattice is stored in a dynamically allocated array:
     lattice[lex][dir], where 'lex' is the lexicographic site index and
     'dir' labels the direction.
*/
void init_lattice(int ***lattice, long int volume)
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
void init_expTable(double **expTable, double beta)
{
    int deltaSmax = 4 * (dim - 1);
    int range = 2 * deltaSmax + 1;

    *expTable = (double *)malloc((unsigned int)range * sizeof(double));
    for (int deltaS = -deltaSmax; deltaS <= deltaSmax; deltaS += 4)
    {
        (*expTable)[deltaS + deltaSmax] = exp(beta * deltaS);
    }
}

/* Compute the sum of staples around a selected link.

   For the link at site 'lex' in direction 'dir' this function returns the sum
   of all staples of that link where each staple consists of a product of three
   neighbouring links forming the sides of the plaquette.
*/
int computeStaple(int **restrict Lattice,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm,
                  long int lex, int dir, long int volume)
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
        lex_minus_orth = nnm[dirgeo(lex, orth, volume)];

        // ---------- forward staple ----------
        linkProd = Lattice[lex][orth];
        linkProd *= Lattice[lex_plus_dir][orth];
        linkProd *= Lattice[lex_plus_orth][dir];

        sumStaples += linkProd;

        // ---------- backward staple ----------
        linkProd = Lattice[lex_minus_orth][dir];
        linkProd *= Lattice[lex_minus_orth][orth];
        linkProd *= Lattice[nnp[dirgeo(lex_minus_orth, dir, volume)]][orth];

        sumStaples += linkProd;
    }
    return sumStaples;
}

/* Perform a full Metropolis sweep over all link variables of the lattice.

   Each link Lattice[lex][dir] is flipped according to the Metropolis algorithm:
   - If ΔS > 0, the flip is always accepted.
   - If ΔS ≤ 0, the flip is accepted with probability exp(ΔS), looked up from expTable.

   The function iterates over all sites and all directions, computes ΔS from
   the surrounding staples, updates the lattice accordingly and returns the
    accepted flips rate.
*/
double sweepMetro(int **restrict Lattice,
                  double *restrict expTable,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm,
                  long int volume)
{
    long int updateCounter = 0, trialCounter = 0;
    int deltaS, deltaSmax, sumStaples, maxStaple;

    // max possible value returned by computeStaple (all link = 1)
    maxStaple = 2 * (dim - 1);
    // max Action variation, used to match expTable index
    deltaSmax = 2 * maxStaple;

    // ordered update of all lattice sites
    for (long int lex = 0; lex < volume; lex++)
    {
        // update in all possible direction for each site
        for (int dir = 0; dir < dim; dir++)
        {
            if (myrand() < epsilon)
            {
                continue;
            }
            else
            {   trialCounter+=1;
                // compute action variation associated with proposed trial update
                sumStaples = computeStaple(Lattice, nnp, nnm, lex, dir, volume);
                deltaS = -2 * Lattice[lex][dir] * sumStaples;
                // to avoid out of range index in expTable
                if (deltaS < -deltaSmax || deltaS > deltaSmax)
                {
                    fprintf(stderr, "Warning: out of range index(lex=%ld, dir=%d)\n",
                            lex, dir);
                    continue;
                }
                if (deltaS > 0)
                {
                    Lattice[lex][dir] *= -1;
                    updateCounter += 1;
                }
                else if (myrand() < expTable[deltaS + deltaSmax])
                {
                    Lattice[lex][dir] *= -1;
                    updateCounter += 1;
                }
            }
        }
    }
    return (double)updateCounter / (double)(trialCounter);
}

/* Compute the Wilson loop W(Wt, Ws).

   The Wilson loop is defined as the product of link variables along a closed
   rectangular path in the (0, dir) plane, with temporal extent Wt and spatial extent Ws.

   Although a single loop would suffice in principle, the function averages over
   all lattice sites and spatial directions to exploit lattice symmetries and
   improve statistics.
*/
double WilsonLoop(int **Lattice,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm,
                  long int volume, int Wt, int Ws)
{
    double res = 0.0;
    long int check_lex;
    // loop over all lattice sites
    for (long int lex = 0; lex < volume; lex++)
    {
        check_lex = lex;
        // loop over all spatial directions
        for (int dir = 1; dir < dim; dir++)
        {
            /*
             ^ 0-temp direction
             | r1      r2
             +--------+
             |        |
             |        |
             |        | dir-spatial direction
          ---+--------+----->
             r0        r3
            */

            // initialize wilson loop variable
            int loop = 1;

            // starting point r0
            for (int i = 0; i < Wt; i++)
            {
                loop *= Lattice[lex][0];
                lex = nnp[dirgeo(lex, 0, volume)];
            }
            // starting point r1
            for (int i = 0; i < Ws; i++)
            {
                loop *= Lattice[lex][dir];
                lex = nnp[dirgeo(lex, dir, volume)];
            }
            // starting point r2
            for (int i = 0; i < Wt; i++)
            {
                lex = nnm[dirgeo(lex, 0, volume)];
                loop *= Lattice[lex][0];
            }
            // starting point r3
            for (int i = 0; i < Ws; i++)
            {
                lex = nnm[dirgeo(lex, dir, volume)];
                loop *= Lattice[lex][dir];
            }
            // check if loop has correctly been closed
            if (lex != check_lex)
            {
                fprintf(stderr, "Warning: loop not closed (lex=%ld, dir=%d)\n",
                        check_lex, dir);
                continue;
            }
            // at r0 again, loop closed
            res += (double)loop;
        }
    }
    // compute average
    res /= (double)volume;
    res /= (double)(dim - 1);

    return res;
}

int main(int argc, char **argv)
{
    FILE *fp;
    int **Lattice;
    long int *nnp, *nnm;
    double *expTable;
    char datafile[STRING_LENGTH];

    int size;
    long int sample, therm;
    double beta, loop, accRate = 0.;

    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;
    const int measevery = 200;

    if (argc != 6)
    {
        fprintf(stdout, "How to use this program:\n");
        fprintf(stdout, "  %s size beta sample therm datafile\n\n", argv[0]);
        fprintf(stdout, "  size = temporal and spatial size of the lattice\n");
        fprintf(stdout, "  (space-time dimension defined by macro dim)\n");
        fprintf(stdout, "  beta = coupling\n");
        fprintf(stdout, "  sample = number of drawn to be extracted\n");
        fprintf(stdout, "  therm = number of thermalization sweeps\n");
        fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
        fprintf(stdout, "Compiled for:\n");
        fprintf(stdout, "  dimensionality = %d\n\n", dim);
        fprintf(stdout, "Output:\n");
        fprintf(stdout, "  Wilson loops (Wt, Ws) using the following order\n");
        fprintf(stdout, "  for(Ws=1; Ws<=size/4; Ws++){\n");
        fprintf(stdout, "      for(Wt=1; Wt<=MIN(size/4,8); Wt++){\n");
        fprintf(stdout, "         Ws, Wt Wilson loop }}\n");

        return EXIT_SUCCESS;
    }
    else
    {
        // read and check input values
        size = atoi(argv[1]);
        beta = atof(argv[2]);
        sample = atol(argv[3]);
        therm = atol(argv[4]);

        if (strlen(argv[5]) >= STRING_LENGTH)
        {
            fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
        else
        {
            strcpy(datafile, argv[5]);
        }
        if (sample <= 0)
        {
            fprintf(stderr, "'sample' must be positive\n");
            return EXIT_FAILURE;
        }
        if (therm <= 0)
        {
            fprintf(stderr, "'therm' must be positive\n");
            return EXIT_FAILURE;
        }
        if (size <= 0)
        {
            fprintf(stderr, "'size' must be positive\n");
            return EXIT_FAILURE;
        }
    }
    // to get n°sample measures
    sample *= measevery;

    // initialize random number generator
    myrand_init(seed1, seed2);

    // compute lattice space time volume
    long int volume = 1;
    for (int dir = 0; dir < dim; dir++)
    {
        volume *= size;
    }

    // initialize data structures
    init_expTable(&expTable, beta);
    init_lattice(&Lattice, volume);
    init_neighbours(&nnp, &nnm, size);

#if DEBUG
    // debug test of data structure functions
    debugTests(Lattice, nnp, nnm, dim, volume, beta, expTable, size);
#endif

    // open file
    fp = fopen(datafile, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    // thermalization
    for (long int iter = 0; iter < therm; iter++)
    {
        sweepMetro(Lattice, expTable, nnp, nnm, volume);
    }
    for (long int iter = 0; iter < sample; iter++)
    {
        accRate += sweepMetro(Lattice, expTable, nnp, nnm, volume);
        if (iter % measevery == 0)
        {
            for (int Ws = 1; Ws <= size / 4; Ws++)
            {
                // fprintf(fp, "%d ",Ws);
                for (int Wt = 1; Wt <= MIN(size / 4, 8); Wt++)
                {
                    loop = WilsonLoop(Lattice, nnp, nnm, volume, Wt, Ws);
                    fprintf(fp, "%.12f ", loop);
                }
                // fprintf(fp, "\n");
            }
            fprintf(fp, "\n");
        }
    }
    accRate /= (double)sample;

    printf("\n");
    printf("Acceptance rate: %f\n", accRate);

    // close datafile
    fclose(fp);

    // free memory
    for (long int i = 0; i < volume; i++)
    {
        free(Lattice[i]);
    }
    free(Lattice);
    free(nnp);
    free(nnm);
    free(expTable);

    return EXIT_SUCCESS;
}
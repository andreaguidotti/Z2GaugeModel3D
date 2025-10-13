#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/geometry.h"
#include "../include/random.h"
#include "../include/debug.h"

extern int computeStaples(int **Lattice, long int const *nnp,
                          long int const *nnm, long int lex, int dir,
                          long int volume);

/* Check staple correctness by flipping links.

   Flips each link to identify which staples should change,
   then verifies computeStaples() returns the correct values.
*/
void check_staple(int **restrict Lattice, long int const *restrict nnp,
                  long int const *restrict nnm, int dim, long int volume)
{
    printf("Running staple consistency check:\n");

    long int lex_plus_orth, lex_plus_dir, lex_minus_orth, errorCounter;

    int **changed = (int **)malloc((unsigned int)volume * sizeof(int *));
    int **expected = (int **)malloc((unsigned int)volume * sizeof(int *));

    int **beforeStaple = (int **)malloc((unsigned int)volume * sizeof(int *));
    int **afterStaple = (int **)malloc((unsigned int)volume * sizeof(int *));

    for (long int site = 0; site < volume; site++)
    {
        expected[site] = (int *)calloc((unsigned int)dim, sizeof(int));
        beforeStaple[site] = (int *)calloc((unsigned int)dim, sizeof(int));
        afterStaple[site] = (int *)calloc((unsigned int)dim, sizeof(int));
        changed[site] = (int *)calloc((unsigned int)dim, sizeof(int));
    }
    // Compute staples before any flips
    for (long int site = 0; site < volume; site++)
        for (int d = 0; d < dim; d++)
            beforeStaple[site][d] = computeStaples(Lattice, nnp, nnm, site, d, volume);

    // Loop over all sites and directions to flip each link
    errorCounter = 0;
    for (long int lex = 0; lex < volume; lex++)
    {
        for (int dir = 0; dir < dim; dir++)
        {
            // Flip the current link
            Lattice[lex][dir] *= -1;

            // Recompute staples after the flip
            for (long int site = 0; site < volume; site++)
                for (int d = 0; d < dim; d++)
                    afterStaple[site][d] = computeStaples(Lattice, nnp, nnm, site, d, volume);

            // Determine which staples are expected to change
            lex_plus_dir = nnp[dirgeo(lex, dir, volume)];
            for (int d = 0; d < dim; d++)
            {
                if (d == dir)
                    continue;

                expected[lex][d] = 1;
                expected[lex_plus_dir][d] = 1;
                lex_plus_orth = nnp[dirgeo(lex, d, volume)];
                expected[lex_plus_orth][dir] = 1;
                lex_minus_orth = nnm[dirgeo(lex, d, volume)];
                expected[lex_minus_orth][d] = 1;
                expected[lex_minus_orth][dir] = 1;
                expected[nnp[dirgeo(lex_minus_orth, dir, volume)]][d] = 1;
            }
            // Compare actual changes with expected ones
            for (long int site = 0; site < volume; site++)
            {
                for (int d = 0; d < dim; d++)
                {
                    if (beforeStaple[site][d] != afterStaple[site][d])
                        changed[site][d] = 1;
                    if (expected[site][d] != changed[site][d])
                    {
                        errorCounter += 1;
                        printf("Mismatch at site %ld, dir %d (expected %d, got %d)\n",
                               site, d, expected[site][d], changed[site][d]);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            // Reset expected and changed arrays to zero for next flip
            for (long int i = 0; i < volume; i++)
                for (int d = 0; d < dim; d++)
                {
                    expected[i][d] = 0;
                    changed[i][d] = 0;
                }
            // Flip the link back to original state
            Lattice[lex][dir] *= -1;
        }
    }
    if (errorCounter == 0)
    {
        printf("ComputeStaples test passed");
    }

    // Free all dynamically allocated memory
    for (long int i = 0; i < volume; i++)
    {
        free(expected[i]);
        free(beforeStaple[i]);
        free(afterStaple[i]);
        free(changed[i]);
    }
    free(expected);
    free(beforeStaple);
    free(afterStaple);
    free(changed);
}

/* Debug function that checks precomputed expTable values .

   Recomputes exp(beta * ΔS) for all possible ΔS values and compares with expTable.
   Prints an error and exits if a mismatch is found, otherwise confirms test passed.
*/
void check_expTable(double const *expTable, double beta, int dim)
{
    printf("Running expTable consistency check:\n");

    int deltaSmax = 4 * (dim - 1);
    int sampleMax = 2 * (dim - 1);

    for (int sampleSum = -sampleMax; sampleSum <= sampleMax; sampleSum += 2)
    {
        int deltaS = -2 * sampleSum;
        double check = exp(beta * deltaS);

        if (fabs(check - expTable[deltaS + deltaSmax]) > 1e-12)
        {
            printf("expTable test failed at deltaS = %d\n", deltaS);
            exit(EXIT_FAILURE);
        }
    }
    printf("expTable test passed\n");
}

/* Check nearest-neighbor tables.

   Verifies that forward (nnp) and backward (nnm) nearest-neighbor arrays
   correctly map lattice sites in all directions by converting between
   lexicographic and Cartesian coordinates.
*/
void check_neighbours(long int *nnp, long int *nnm, int size, int dim)
{
    printf("Running nnp-nnm consistency check:\n");
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

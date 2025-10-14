#ifndef DEBUG_H
#define DEBUG_H

void check_staple(int **restrict Lattice, 
                  long int const * const restrict nnp,
                  long int const * const restrict nnm, 
                  int dim, long int volume);

void check_expTable(double const *expTable, double beta, int dim);

void check_neighbours(long int const * restrict nnp, 
                      long int const * restrict nnm, 
                      int size, int dim);

void debugTests(int **restrict Lattice, 
                long int const * const restrict nnp,
                long int const * const restrict nnm, 
                int dim, long int volume, double beta, 
                double *expTable, int size);

#endif
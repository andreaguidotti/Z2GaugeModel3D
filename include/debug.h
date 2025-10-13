#ifndef DEBUG_H
#define DEBUG_H

void check_staple(int **restrict Lattice, long int const *restrict nnp,
                  long int const *restrict nnm, int dim, long int volume);

void check_expTable(double const *expTable, double beta, int dim);

void check_neighbours(long int *nnp, long int *nnm, int size, int dim);


#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include"../include/geometry.h"

/* from cartesian to lexicographic coordinates
 e.g (1,2,3) with lattice Size = L is mapped in 1+2*L+3*L^2*/

void cart_to_lex(long int *lex, int *cartcoord, int dim, int L)
{
    long int aux, res;

    aux = 1;
    res = 0;

    for (int i = 0; i < dim; i++)
    {
        res += cartcoord[i] * aux;
        aux *= L;
    }

    *lex = res;
}

/* from lexicographic to cartesian coordinates
e.g n with  lattice Size = L is mapped in: n/L^2 = cartcoord[2],...*/

void lex_to_cart(long int lex,  int *cartcoord, int dim, int L)
{
    long int aux;
    aux = 1;

    for (int i = 0; i < dim-1; i++)
    {
        aux *= L;
    }

    for (int j = dim - 1; j >= 0; j--)
    {
        cartcoord[j] = (int)(lex / aux);
        lex -= cartcoord[j] * aux;
        aux /= L;
    }
}
/*init nearest neighbor sites for each element and save them in nnm and nnp*/
void init_geo(long int *restrict nnp,
           long int *restrict nnm,
           int L,
           int dim)
{
    long int volume;
    int auxp, auxm, value;
    long int rp, rm;

    int *cartcoord = ( int *)malloc((unsigned int)dim * sizeof(int));

    if (cartcoord == NULL)
    {
        fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    volume = 1;
    for (int d = 0; d < dim; d++)
    {
        volume *= L;
    }
  //for each element i, in lex-coord, go to cart-coordinates and maps to lex its neighb
    for (long int i = 0; i < volume; i++)
    {
        lex_to_cart(i, cartcoord, dim, L);
        for (int ii = 0; ii < dim; ii++)
        {
            value = cartcoord[ii];

            auxp = value + 1;
            if (auxp >= L)
            {
                auxp -= L;
            }
            cartcoord[ii] = auxp;
            cart_to_lex(&rp, cartcoord, dim, L);
            nnp[dirgeo(i, ii, volume)] = rp;

            auxm = value - 1;
            if (auxm < 0)
            {
                auxm += L;
            }
            cartcoord[ii] = auxm;
            cart_to_lex(&rm, cartcoord, dim, L);
            nnm[dirgeo(i, ii, volume)] = rm;

            cartcoord[ii]=value;
        }
    }
    free(cartcoord);
}

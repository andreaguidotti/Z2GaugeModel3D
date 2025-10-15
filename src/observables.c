#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define STRING_LENGTH 50
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

// Struct for storing (Wilson loop) data.
typedef struct data
{
    double totSum; // Total sum of array elements
    double aux;    // Auxiliary variable for jackknife calculations
    double *arr;   // Array of data

    double kahanCorrection; // correction to be used in sumKahan

    int Wt; // Temporal size Wilson loop
    int Ws; // Spazial size Wilson loop

} data;

// Struct for storing observable results and jackknife statistics.
typedef struct obs
{
    double avg;         // Simple average of observable
    double std;         // Standard deviation
    double jackavg;     // Jackknife mean
    double jackavgSqrd; // Jackknife mean squared

    double kahanCorrection[2]; //  corrections to be used in sumKahan

} obs;

// Reset all fields in a data structure
void resetData(data *W)
{
    W->totSum = 0.0;
    W->aux = 0.0;
    W->kahanCorrection = 0.0;
    W->Wt = 0;
    W->Ws = 0;
}
// Reset all fields in a obs structure
void resetObs(obs *Potential)
{
    Potential->avg = 0.0;
    Potential->std = 0.0;
    Potential->jackavg = 0.0;
    Potential->jackavgSqrd = 0.0;
    Potential->kahanCorrection[0] = 0.0;
    Potential->kahanCorrection[1] = 0.0;
}

/* Kahan summation algorithm.

   Performs numerically stable summation by tracking a small
   correction term to compensate for floating-point errors.
*/
void sumKahan(double addend, double *sumtot, double *correction)
{
    double correctedAddend, tempSum;

    correctedAddend = addend + *correction;
    tempSum = *sumtot + correctedAddend;
    *correction = (*sumtot - tempSum) + correctedAddend;
    *sumtot = tempSum;
}

/* Extract data from a file for a specific column.

   Reads 'sampleEff' lines from file 'fp', extracts the value in
   column 'target_col' and stores it in the array 'W->arr'.

   Uses Kahan summation to accumulate the total sum.
   Updates the Wilson loop dimensions Wt and Ws.
*/
void extractData(FILE *fp, data *W, long int sampleEff, long pos,
                 int target_col, int Wt_counter, int Ws_counter)
{
    char buffer[1024];
    double value;

    int read;

    for (int row = 0; row < sampleEff; row++)
    {
        if (fgets(buffer, sizeof(buffer), fp) == NULL)
        {
            fprintf(stderr, "Error: unexpected end of file"
                            "or read error at row %d\n", row);
            free(W->arr);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
        char *ptr = buffer;
        for (int col = 0; col <= target_col; col++)
        {
            sscanf(ptr, "%lf%n", &value, &read);
            ptr += read;
        }
        W->arr[row] = value;
        // Accumulate sum with Kahan correction for numerical accuracy
        sumKahan(value, &(W->totSum), &(W->kahanCorrection));

        W->Wt = Wt_counter;
        W->Ws = Ws_counter;
    }
    // Add residual Kahan correction after the loop
    W->totSum += W->kahanCorrection;

    if (fseek(fp, pos, SEEK_SET) != 0)
    {
        fprintf(stderr,
                "Error: unable to reset file pointer to saved position, "
                "extraction of data in column %d failed\n",
                target_col);

        fclose(fp);
        free(W->arr);
        exit(EXIT_FAILURE);
    }
}

/* Compute Wilson loop potential.

   Given the average Wilson loop value and temporal size Wt,
   returns the string potential.
*/
double StringU(double Wilson, int Wt)
{
    double res;
    res = -log(Wilson);
    res /= (double)Wt;

    return res;
}

/* Jackknife leave-one-out.

   Computes the mean of the data excluding the i-th block.
   'blockdim' is the block size, 'nblocks' is the total number of blocks.
*/
static inline void jackLeaveOneOut(data *data, long int i, long int blockdim, long int nblocks)
{
    long int idx;

    data->aux = data->totSum;
    for (long int ii = 0; ii < blockdim; ii++)
    {
        idx = i * blockdim + ii;
        data->aux -= data->arr[idx];
    }
    data->aux /= (double)((nblocks - 1) * blockdim);
}

// Accumulate mean and mean squared using Kahan summation
static inline void accumulate(obs *o, double value)
{
    sumKahan(value, &o->jackavg, &o->kahanCorrection[0]);
    sumKahan(value * value, &o->jackavgSqrd, &o->kahanCorrection[1]);
}

/* Perform jackknife analysis.

   Computes the jackknife mean and standard deviation for the
   Wilson loop potential over 'nblocks' blocks of size 'blockdim'.
*/
void jackknife(obs *Potential, data *W,
               long int nblocks, long int blockdim)
{
    double valuePotential;

    for (long int i = 0; i < nblocks; i++)
    {
        jackLeaveOneOut(W, i, blockdim, nblocks);

        if (W->aux <= 0)
        {
            fprintf(stderr,
                    "Error:negative log argument, "
                    "larger statistics needed for "
                    "Wt-%d Ws-%d Wilson loop",
                    W->Wt, W->Ws);

            free(W->arr);
            exit(EXIT_FAILURE);
        }

        // mean and mean squared of jacksamples

        valuePotential = StringU(W->aux, W->Wt);

        accumulate(Potential, valuePotential);
    }

    // Add the Kahan residual from the last addition

    Potential->jackavg += Potential->kahanCorrection[0];
    Potential->jackavgSqrd += Potential->kahanCorrection[1];

    // Normalization

    Potential->jackavg /= nblocks;
    Potential->jackavgSqrd /= nblocks;

    // Standard deviation

    Potential->std = sqrt((nblocks - 1) * (Potential->jackavgSqrd - pow(Potential->jackavg, 2)));
}
int main(int argc, char **argv)
{
    if (argc != 8)
    {
        printf("How to use this program: \n");
        printf("%s inputfile, outputfile, sample, blocksize, beta, size, therm\n\n", argv[0]);
        printf("  inputfile  : path to input data file\n");
        printf("  outputfile : path to output results file\n");
        printf("  sample     : number of measurements\n");
        printf("  blocksize  : size of each jackknife block (>=2)\n");
        printf("  beta       : inverse temperature parameter\n");
        printf("  size       : lattice size (>=4)\n");
        printf("  therm      : thermalization steps (>=0)\n");
        printf("\n");

        return EXIT_SUCCESS;
    }

    // initialize variables and data structures

    obs Potential = {0};
    data W = {0};

    long int sample, sampleEff, blockdim, nblocks, therm;
    int size;
    double beta;

    char infile[STRING_LENGTH];
    char outfile[STRING_LENGTH];

    sample = atol(argv[3]);
    blockdim = atol(argv[4]);
    beta = atof(argv[5]);
    size = atoi(argv[6]);
    therm = atol(argv[7]);

    // Check input

    if (blockdim < 2)
    {
        fprintf(stderr, "blockdim has to be at least 2");
        return EXIT_FAILURE;
    }
    if (size < 1)
    {
        fprintf(stderr, "size has to be at least 1");
        return EXIT_FAILURE;
    }
    if (therm < 0)
    {
        fprintf(stderr, "thermalization has to be at least 0");
        return EXIT_FAILURE;
    }
    nblocks = (sample - therm) / blockdim;
    sampleEff = nblocks * blockdim;

    // Open  input and output file

    strcpy(infile, argv[1]);
    strcpy(outfile, argv[2]);

    FILE *fp = fopen(infile, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file in %s,%d", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    FILE *out = fopen(outfile, "a");
    if (out == NULL)
    {
        fprintf(stderr, "Error defining array %s,%d", __FILE__, __LINE__);
        fclose(fp);
        return EXIT_FAILURE;
    }

    // Allocate dynamic arrays

    W.arr = (double *)malloc((unsigned)sampleEff * sizeof(double));
    if (W.arr == NULL)
    {
        fprintf(stderr, "Error defining array %s,%d", __FILE__, __LINE__);
        fclose(fp);

        return EXIT_FAILURE;
    }

    // Discard thermalization

    char buffer[1024];
    for (long int i = 0; i < therm; i++)
    {
        fgets(buffer, sizeof(buffer), fp);
    }

    // Save pointer position in FILE after thermalization

    long pos = ftell(fp);

    // Determine lattice loop parameters

    int Wt_max = MIN(size / 4, 8), Wt_counter = 1, Ws_counter = 1;
    int ncolumns = MIN(size / 4, 8) * size / 4;

    // Loop over each Wilson loop column

    for (int target_col = 0; target_col < ncolumns; target_col++)
    {
        resetData(&W);
        resetObs(&Potential);

        extractData(fp, &W, sampleEff, pos, target_col, Wt_counter, Ws_counter);

        // Compute average potential for this column

        if (W.totSum / (double)sampleEff <= 0)
        {
            fprintf(stderr,
                    "Error: negative log argument, "
                    "larger statistics needed for "
                    "Wt-%d Ws-%d avg Wilson loop value",
                    W.Wt, W.Ws);

            free(W.arr);
            exit(EXIT_FAILURE);
        }

        Potential.avg = StringU(W.totSum / (double)sampleEff, W.Wt);

        // Perform jackknife analysis

        jackknife(&Potential, &W, nblocks, blockdim);

        // Export results

        fprintf(out, "%lf %d %d %lf %lf \n", beta, W.Ws, W.Wt, Potential.avg, Potential.std);

        // Update lattice counters

        Wt_counter += 1;
        if (Wt_counter > Wt_max)
        {
            Ws_counter += 1;
            Wt_counter = 1;
        }
    }

    // Close files and free allocated memory

    fclose(fp);
    fclose(out);
    free(W.arr);

    return EXIT_SUCCESS;
}
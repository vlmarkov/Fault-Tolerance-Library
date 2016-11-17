/*
 * heat-2d.c: MPI implementation of Laplace equation solver by Jacobi iteration method.
 * 
 * 2D Laplace equation: 
 *   \Delta u = 0
 *   \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0
 * 
 * Domain: x in [0, 1],  y in [0, 1]
 * Boundary conditions: 
 *   u(x, 0) = sin(pi * x)
 *   u(x, 1) = sin(pi * x) * exp(-pi)
 *   u(0, y) = u(1, y) = 0
 * Initial value for interior points is 0
 * Analytical solution: 
 *   u(x, y) = sin(pi * x) * exp(-pi * y)
 * 
 * Parallel implementation: 
 * 2D domain decomposition of grid [0..rows - 1][0..cols -1]
 * Each process is assigned a subgrid [rows / nprocs][cols / nprocs]
 *
 * Input parameters: rows, cols, EPS
 *
 * Usage: mpiexec -np <p> ./heat-2d <rows> <cols>
 *
 * (C) Mikhail Kurnosov, 2015
 */
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include <mpi.h>

#include <string.h>
#include "../checkpoint_lib/checkpoint_lib.h"

#define EPS 0.001
#define PI 3.14159265358979323846

#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))
#define IND(i, j) ((i) * (nx + 2) + (j))


/*************************************************************************/
/* Global variables                                                      */
/*************************************************************************/

int options    = CHECKPOINT_MODE;

int niters     = 0;

int ny, nx;

double ttotal  = 0.0;
double thalo   = 0.0;
double treduce = 0.0;

double *local_grid    = NULL;
double *local_newgrid = NULL;


void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);
    if (p == NULL) {
        fprintf(stderr, "No enough memory\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    return p;    
}

int get_block_size(int n, int rank, int nprocs)
{
    int s = n / nprocs;
    if (n % nprocs > rank)
        s++;
    return s;
}

int get_sum_of_prev_blocks(int n, int rank, int nprocs)
{
    int rem = n % nprocs;
    return n / nprocs * rank + ((rank >= rem) ? rem : rank);
}

/*************************************************************************/
/* Additional functions helps save and get data                          */
/*************************************************************************/
inline static void checkpoint_save(int phase)
{
    if (phase > 0) {
        CHECKPOINT_TIMER_STOP();
    }

    MPI_File local_snapshot;
        
    CHECKPOINT_FILE_OPEN(&local_snapshot, phase);
        
    CHECKPOINT_SAVE(local_snapshot, local_grid, ((ny + 2) * (nx + 2)), MPI_DOUBLE);
    CHECKPOINT_SAVE(local_snapshot, &ttotal, 1, MPI_DOUBLE);
    CHECKPOINT_SAVE(local_snapshot, &thalo, 1, MPI_DOUBLE);
    CHECKPOINT_SAVE(local_snapshot, &treduce, 1, MPI_DOUBLE);
    CHECKPOINT_SAVE(local_snapshot, &niters, 1, MPI_INT);
    
    CHECKPOINT_FILE_CLOSE(&local_snapshot);
}


inline static int checkpoint_get(double *local_grid,
                                 int size,
                                 double *ttotal, 
                                 double *thalo, 
                                 double *treduce, 
                                 int *niters)
{
    char last_chechkpoint[] = { "0_0_0.000000" };
    int phase = CHECKPOINT_GET(last_chechkpoint);

    FILE * file = fopen(last_chechkpoint, "rb");
    if (file) {
        // copy the file into the buffer:
        fread(local_grid, sizeof(double), size, file);
        fread(ttotal, sizeof(double), 1, file);
        fread(thalo, sizeof(double), 1, file);
        fread(treduce, sizeof(double), 1, file);
        fread(niters, sizeof(int), 1, file);
        
        fclose(file);
    }

    return phase;
}

void time_handler(int sig)
{
    checkpoint_save(0);
}

int main(int argc, char *argv[]) 
{
    /*************************************************************************/
    /* Initialize checkpoint library                                         */
    /*************************************************************************/
    CHECKPOINT_LIB_INIT(2, 5); // 2 checkpoints, 5 seconds for timer
    CHECKPOINT_ASSIGN(&&start);
    CHECKPOINT_ASSIGN(&&end);

    CHECKPOINT_SIGNAL_HANDLER(time_handler);


    /*************************************************************************/
    /* Local  variables                                                      */
    /*************************************************************************/
    char procname[MPI_MAX_PROCESSOR_NAME];

    int commsize, rank, rankx, ranky, rows, cols, px, py, namelen;

    int coords[2];
    int dims[2]     = { 0, 0 };
    int periodic[2] = { 0, 0 };

    int left, right, top, bottom;


    /*************************************************************************/
    /* Initialize                                                            */
    /*************************************************************************/

    /* Initialize args, rank, commsize */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Get processor name */
    MPI_Get_processor_name(procname, &namelen);

    /* Create 2D grid of processes: commsize = px * py */
    MPI_Comm cartcomm;
    MPI_Dims_create(commsize, 2, dims);
    
    px = dims[0];
    py = dims[1];

    /* Make start time point */
    ttotal = MPI_Wtime();

    if (px < 2 || py < 2) {
        fprintf(stderr, "Invalid number of processes %d: px %d and py %d must be greater than 1\n", commsize, px, py);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &cartcomm);
    MPI_Cart_coords(cartcomm, rank, 2, coords);

    rankx = coords[0];
    ranky = coords[1];


    /*************************************************************************/
    /* Broadcast command line arguments                                      */
    /*************************************************************************/
    if (rank == 0) {
        rows = (argc > 1) ? atoi(argv[1]) : py * 100;
        cols = (argc > 2) ? atoi(argv[2]) : px * 100;

        if (rows < py) {
            fprintf(stderr, "Number of rows %d less then number of py processes %d\n", rows, py);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        if (cols < px) {
            fprintf(stderr, "Number of cols %d less then number of px processes %d\n", cols, px);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        int args[2] = { rows, cols };
        MPI_Bcast(&args, NELEMS(args), MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        int args[2];
        MPI_Bcast(&args, NELEMS(args), MPI_INT, 0, MPI_COMM_WORLD);
        rows = args[0];
        cols = args[1];
    }


    /*************************************************************************/
    /* Allocate memory for local 2D subgrids with halo cells                 */
    /* [0..ny + 1][0..nx + 1]                                                */
    /*************************************************************************/
    ny = get_block_size(rows, ranky, py);
    nx = get_block_size(cols, rankx, px);

    local_grid    = xcalloc((ny + 2) * (nx + 2), sizeof(*local_grid));
    local_newgrid = xcalloc((ny + 2) * (nx + 2), sizeof(*local_newgrid));

    // Fill boundary points: 
    //   - left and right borders are zero filled
    //   - top border: u(x, 0) = sin(pi * x)
    //   - bottom border: u(x, 1) = sin(pi * x) * exp(-pi)    
    double dx = 1.0 / (cols - 1.0); 
    int sj    = get_sum_of_prev_blocks(cols, rankx, px);

    if (ranky == 0) {
        // Initialize top border: u(x, 0) = sin(pi * x)
        for (int j = 1; j <= nx; j++) {
            // Translate col index to x coord in [0, 1]
            double x = dx * (sj + j - 1);
            int ind  = IND(0, j);
            local_newgrid[ind] = local_grid[ind] = sin(PI * x);
        }
    }

    if (ranky == (py - 1)) {
        // Initialize bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
        for (int j = 1; j <= nx; j++) {
            // Translate col index to x coord in [0, 1]
            double x = dx * (sj + j - 1);
            int ind  = IND(ny + 1, j);
            local_newgrid[ind] = local_grid[ind] = sin(PI * x) * exp(-PI);
        }
    }

    // Neighbours
    MPI_Cart_shift(cartcomm, 0, 1, &left, &right);
    MPI_Cart_shift(cartcomm, 1, 1, &top, &bottom);

    // Left and right borders type
    MPI_Datatype col;
    MPI_Type_vector(ny, 1, nx + 2, MPI_DOUBLE, &col);
    MPI_Type_commit(&col);

    // Top and bottom borders type
    MPI_Datatype row;        
    MPI_Type_contiguous(nx, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    MPI_Request reqs[8];


    if (options == CHECKPOINT_MODE) {
        checkpoint_save(0);
    } else if (options == RECOVERY_MODE) {
        int phase = checkpoint_get(local_grid, ((ny + 2) * (nx + 2)), &ttotal, &thalo, &treduce, &niters);
        CHECKPOINT_GOTO(phase);
    }

    CHECKPOINT_SET(start);
    CHECKPOINT_TIMER_INIT();

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    for (;;) {
        niters++;
        
        // Update interior points
        for (int i = 1; i <= ny; i++) {
            for (int j = 1; j <= nx; j++) {
                local_newgrid[IND(i, j)] = 
                    (local_grid[IND(i - 1, j)] + local_grid[IND(i + 1, j)] +
                     local_grid[IND(i, j - 1)] + local_grid[IND(i, j + 1)]) * 0.25;
            }
        }

        // Check termination condition
        double maxdiff = 0;
        for (int i = 1; i <= ny; i++) {
            for (int j = 1; j <= nx; j++) {
                int ind = IND(i, j);
                maxdiff = fmax(maxdiff, fabs(local_grid[ind] - local_newgrid[ind]));
            }
        }

        // Swap grids (after termination local_grid will contain result)
        double *p     = local_grid;
        local_grid    = local_newgrid;
        local_newgrid = p;    

        treduce -= MPI_Wtime();
        MPI_Allreduce(MPI_IN_PLACE, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        treduce += MPI_Wtime();               
        
        if (maxdiff < EPS) {
            break;
        }
        
        // Halo exchange: T = 4 * (a + b * (rows / py)) + 4 * (a + b * (cols / px))
        thalo -= MPI_Wtime();

        MPI_Irecv(&local_grid[IND(0, 1)], 1, row, top, 0, cartcomm, &reqs[0]);         // top
        MPI_Irecv(&local_grid[IND(ny + 1, 1)], 1, row, bottom, 0, cartcomm, &reqs[1]); // bottom
        MPI_Irecv(&local_grid[IND(1, 0)], 1, col, left, 0, cartcomm, &reqs[2]);        // left
        MPI_Irecv(&local_grid[IND(1, nx + 1)], 1, col, right, 0, cartcomm, &reqs[3]);  // right
        MPI_Isend(&local_grid[IND(1, 1)], 1, row, top, 0, cartcomm, &reqs[4]);         // top
        MPI_Isend(&local_grid[IND(ny, 1)], 1, row, bottom, 0, cartcomm, &reqs[5]);     // bottom
        MPI_Isend(&local_grid[IND(1, 1)], 1, col, left, 0, cartcomm, &reqs[6]);        // left
        MPI_Isend(&local_grid[IND(1, nx)], 1, col, right, 0, cartcomm, &reqs[7]);      // right
        
        MPI_Waitall(8, reqs, MPI_STATUS_IGNORE);
        
        thalo += MPI_Wtime();
    }

    checkpoint_save(1);
    CHECKPOINT_SET(end);


    MPI_Type_free(&row);
    MPI_Type_free(&col);
       
    free(local_newgrid);
    free(local_grid);

    ttotal += MPI_Wtime();

    if (rank == 0) {
        printf("# Heat 2D (mpi): grid: rows %d, cols %d, procs %d (px %d, py %d)\n", rows, cols, commsize, px, py);
    }
    
    printf("# P %4d (%2d, %2d) on %s: grid ny %d nx %d, total %.6f,"
        " mpi %.6f (%.2f) = allred %.6f (%.2f) + halo %.6f (%.2f)\n", 
        rank, rankx, ranky, procname, ny, nx, ttotal, treduce + thalo, 
        (treduce + thalo) / ttotal, treduce, treduce / (treduce + thalo), 
        thalo, thalo / (treduce + thalo)); 
        
    double prof[3] = { ttotal, treduce, thalo };    

    if (rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, prof, NELEMS(prof), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        printf("# procs %d : grid %d %d : niters %d : total time %.6f :" 
            " mpi time %.6f : allred %.6f : halo %.6f\n", 
            commsize, rows, cols, niters, prof[0], prof[1] + prof[2], prof[1], prof[2]);
    } else {
        MPI_Reduce(prof, NULL, NELEMS(prof), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

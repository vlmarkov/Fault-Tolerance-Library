/***************************************************************************************/
/* heat-2d.c: MPI implementation of Laplace equation solver by Jacobi iteration method.*/
/*                                                                                     */ 
/* 2D Laplace equation:                                                                */
/*   \Delta u = 0                                                                      */
/*   \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0         */
/*                                                                                     */
/* Domain: x in [0, 1],  y in [0, 1]                                                   */
/* Boundary conditions:                                                                */
/*   u(x, 0) = sin(pi * x)                                                             */
/*   u(x, 1) = sin(pi * x) * exp(-pi)                                                  */
/*   u(0, y) = u(1, y) = 0                                                             */
/* Initial value for interior points is 0                                              */
/* Analytical solution:                                                                */
/*   u(x, y) = sin(pi * x) * exp(-pi * y)                                              */
/*                                                                                     */
/* Parallel implementation:                                                            */
/* 2D domain decomposition of grid [0..rows - 1][0..cols -1]                           */
/* Each process is assigned a subgrid [rows / nprocs][cols / nprocs]                   */
/*                                                                                     */
/* Input parameters: rows, cols, EPS                                                   */
/*                                                                                     */
/* Usage: mpiexec -np <p> ./heat-2d <rows> <cols>                                      */
/*                                                                                     */
/* (C) Mikhail Kurnosov, 2015                                                          */
/***************************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include <mpi.h>

#include "utils.h"
#include "grid-task.h"

#define EPS       0.001
#define PI        3.14159265358979323846
#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))
#define IND(i, j) ((i) * (nx + 2) + (j))


int main(int argc, char *argv[])
{
    int rank;
    int commsize;

    MPI_Init(&argc, &argv);

    double ttotal = -MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /*
     * Create 2D grid of processes: commsize = px * py
     */
    MPI_Comm cartcomm;
    int dims[2]     = {0, 0};
    int periodic[2] = {0, 0};

    MPI_Dims_create(commsize, 2, dims);
    int px = dims[0];
    int py = dims[1];

    if (px < 2 || py < 2)
    {
        fprintf(stderr, "Invalid number of processes %d: px %d and py %d"
                        "must be greater than 1\n", commsize, px, py);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &cartcomm);
    int coords[2];
    MPI_Cart_coords(cartcomm, rank, 2, coords);
    int rankx = coords[0];
    int ranky = coords[1];

    int rows, cols;

    /* 
     * Broadcast command line arguments
     */
    if (rank == 0)
    {
        rows = (argc > 1) ? atoi(argv[1]) : py * 100;
        cols = (argc > 2) ? atoi(argv[2]) : px * 100;

        if (rows < py)
        {
            fprintf(stderr, "Number of rows %d less then number of py processes %d\n", rows, py);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (cols < px)
        {
            fprintf(stderr, "Number of cols %d less then number of px processes %d\n", cols, px);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        int args[2] = {rows, cols};
        MPI_Bcast(&args, NELEMS(args), MPI_INT, 0, MPI_COMM_WORLD);
    }
    else
    {
        int args[2];
        MPI_Bcast(&args, NELEMS(args), MPI_INT, 0, MPI_COMM_WORLD);
        rows = args[0];
        cols = args[1];
    }

    /* 
     * Allocate memory for local 2D subgrids with halo cells [0..ny + 1][0..nx + 1]
     */
    int ny = get_block_size(rows, ranky, py);
    int nx = get_block_size(cols, rankx, px);

    grid_task_t *grid_task = grid_task_allocate(cols, rows, nx, ny, commsize); // TODO

    grid_task_init(grid_task);

    double *local_grid    = grid_task_local_grid_get(grid_task, rank);
    double *local_newgrid = grid_task_local_newgrid_get(grid_task, rank);

    int    r_ranks[commsize];          // Redundancy ranks
    double *r_local_grid[commsize];    // Redundancy local_grids
    double *r_local_newgrid[commsize]; // Redundancy local_newgrids

    memset(r_ranks, 0, sizeof(int) * commsize);
    memset(r_local_grid, 0, sizeof(double *) * commsize);
    memset(r_local_newgrid, 0, sizeof(double *) * commsize);

    // Get redundancy ranks
    int r_size = grid_task_redundancy_ranks_get(grid_task, rank, r_ranks);

    // Get redundancy local_grid and local_newgrid
    for (int i = 0; i < r_size; i++)
    {
        r_local_grid[i]    = grid_task_redundancy_local_grid_get(grid_task, r_ranks[i]);
        r_local_newgrid[i] = grid_task_redundancy_local_newgrid_get(grid_task, r_ranks[i]);
    }

    int local_grid_size = (ny + 2) * (nx + 2);

/*
    for (int i = 0; i < r_size; i++)
    {
        if (r_local_grid[i] && r_local_newgrid[i])
        {
            printf("I'am rank %04d redundancy local_grid or local_newgrid - good\n", rank);
        }
        else
        {
            fprintf(stderr, "I'am rank %04d redundancy local_grid or local_newgrid - failed\n", rank);
        }
    }
*/

    /*
     * Fill boundary points: 
     *   - left and right borders are zero filled
     *   - top border: u(x, 0) = sin(pi * x)
     *   - bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
     */
    double dx = 1.0 / (cols - 1.0); 
    int    sj = get_sum_of_prev_blocks(cols, rankx, px);

    if (ranky == 0)
    {
        // Initialize top border: u(x, 0) = sin(pi * x)
        for (int j = 1; j <= nx; j++)
        {
            // Translate col index to x coord in [0, 1]
            double x = dx * (sj + j - 1);
            int ind = IND(0, j);
            local_newgrid[ind] = local_grid[ind] = sin(PI * x);
        }
    }

    if (ranky == py - 1)
    {
        // Initialize bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
        for (int j = 1; j <= nx; j++)
        {
            // Translate col index to x coord in [0, 1]
            double x = dx * (sj + j - 1);
            int ind = IND(ny + 1, j);
            local_newgrid[ind] = local_grid[ind] = sin(PI * x) * exp(-PI);
        }
    }

    /*
     * Neighbours
     */
    int left, right, top, bottom;
    grid_task_neighbors_get(grid_task, rank, &top, &bottom, &left, &right);

    /*
     * Left and right borders type
     */
    MPI_Datatype col;
    MPI_Type_vector(ny, 1, nx + 2, MPI_DOUBLE, &col);
    MPI_Type_commit(&col);

    /*
     * Top and bottom borders type
     */
    MPI_Datatype row;
    MPI_Type_contiguous(nx, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    MPI_Request reqs[8];
    double thalo   = 0;
    double treduce = 0;

    int niters = 0;

    while (1)
    {
        niters++;

        // Update interior points
        for (int i = 1; i <= ny; i++)
        {
            for (int j = 1; j <= nx; j++)
            {
                local_newgrid[IND(i, j)] = 
                    (local_grid[IND(i - 1, j)] + local_grid[IND(i + 1, j)] +
                     local_grid[IND(i, j - 1)] + local_grid[IND(i, j + 1)]) * 0.25;
            }
        }

        // Check termination condition
        double maxdiff = 0;
        for (int i = 1; i <= ny; i++)
        {
            for (int j = 1; j <= nx; j++)
            {
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

        if (maxdiff < EPS)
        {
            break;
        }

        // Halo exchange: T = 4 * (a + b * (rows / py)) + 4 * (a + b * (cols / px))
        thalo -= MPI_Wtime();

        MPI_Irecv(&local_grid[IND(0, 1)],      1, row, top,    0, cartcomm, &reqs[0]); // top
        MPI_Irecv(&local_grid[IND(ny + 1, 1)], 1, row, bottom, 0, cartcomm, &reqs[1]); // bottom
        MPI_Irecv(&local_grid[IND(1, 0)],      1, col, left,   0, cartcomm, &reqs[2]); // left
        MPI_Irecv(&local_grid[IND(1, nx + 1)], 1, col, right,  0, cartcomm, &reqs[3]); // right

        MPI_Isend(&local_grid[IND(1, 1)],      1, row, top,    0, cartcomm, &reqs[4]); // top
        MPI_Isend(&local_grid[IND(ny, 1)],     1, row, bottom, 0, cartcomm, &reqs[5]); // bottom
        MPI_Isend(&local_grid[IND(1, 1)],      1, col, left,   0, cartcomm, &reqs[6]); // left
        MPI_Isend(&local_grid[IND(1, nx)],     1, col, right,  0, cartcomm, &reqs[7]); // right

        MPI_Waitall(8, reqs, MPI_STATUS_IGNORE);

        thalo += MPI_Wtime();

        // Exchnage redundancy local_grids
        int r = 0;
        MPI_Request r_reqs[2 *r_size];
        for (int i = 0; i < r_size; i++)
        {
            double *r_grid = r_local_grid[i];
            if (rank == 0)
            {
                //printf("I'am rank %04d redundancy irecv %04d\n", rank, r_ranks[i]);
            }
            MPI_Irecv(r_grid, local_grid_size, MPI_DOUBLE, r_ranks[i], 0, cartcomm, &r_reqs[r++]);
        }

        for (int i = 0; i < r_size; i++)
        {
            if (rank == 0)
            {
                //printf("I'am rank %04d redundancy isend %04d\n", rank, r_ranks[i]);
            }
            MPI_Isend(local_grid, local_grid_size, MPI_DOUBLE, r_ranks[i], 0, cartcomm, &r_reqs[r++]);
        }
        MPI_Waitall(2 *r_size, r_reqs, MPI_STATUS_IGNORE);
    }

    MPI_Type_free(&row);
    MPI_Type_free(&col);

    free(local_newgrid); // todo
    free(local_grid); // todo

    ttotal += MPI_Wtime();

    if (rank == 0)
    {
        printf("# Heat 2D (mpi): grid: rows %d, cols %d, procs %d (px %d, py %d)\n", 
                rows, cols, commsize, px, py);
    }

    int namelen;
    char procname[MPI_MAX_PROCESSOR_NAME];

    MPI_Get_processor_name(procname, &namelen);

    printf("# P %4d (%2d, %2d) on %s: grid ny %d nx %d, total %.6f,"
           " mpi %.6f (%.2f) = allred %.6f (%.2f) + halo %.6f (%.2f)\n", 
           rank, rankx, ranky, procname, ny, nx, ttotal, treduce + thalo, 
           (treduce + thalo) / ttotal, treduce, treduce / (treduce + thalo), 
           thalo, thalo / (treduce + thalo)); 

    double prof[3] = {ttotal, treduce, thalo};
    if (rank == 0)
    {
        MPI_Reduce(MPI_IN_PLACE, prof, NELEMS(prof), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        printf("# procs %d : grid %d %d : niters %d : total time %.6f :" 
               " mpi time %.6f : allred %.6f : halo %.6f\n", 
               commsize, rows, cols, niters, prof[0], prof[1] + prof[2], prof[1], prof[2]);
    }
    else
    {
        MPI_Reduce(prof, NULL, NELEMS(prof), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

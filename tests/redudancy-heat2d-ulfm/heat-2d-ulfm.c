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
#include <mpi-ext.h>

#include "utils.h"
#include "grid-task.h"

#define EPS       0.001
#define PI        3.14159265358979323846
#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))
#define IND(i, j) ((i) * (nx + 2) + (j))

typedef struct
{
    double      *buf;
    int          count;
    MPI_Datatype row;
    MPI_Datatype col;
    int          top;
    int          bottom;
    int          left;
    int          right;
    MPI_Comm     comm;
    MPI_Request *reqs;
    int nx;
    int ny;
} halo_cookie_t;

typedef struct
{
    double      *send;
    double     **receive;
    int         *dest;
    int          ranks;
    int          count;
    MPI_Comm     comm;
    MPI_Request *reqs;
} redundancy_cookie_t;

void update_interior_points(double *newgrid, double *grid, const int ny, const int nx)
{
    for (int i = 1; i <= ny; i++)
    {
        for (int j = 1; j <= nx; j++)
        {
            newgrid[IND(i, j)] = (grid[IND(i - 1, j)] + grid[IND(i + 1, j)] +
                                  grid[IND(i, j - 1)] + grid[IND(i, j + 1)]) * 0.25;
        }
    }
}

double check_termination_condition(double *newgrid, double *grid, const int ny, const int nx)
{
    double maxdiff = 0;
    
    for (int i = 1; i <= ny; i++)
    {
        for (int j = 1; j <= nx; j++)
        {
            int ind = IND(i, j);
            maxdiff = fmax(maxdiff, fabs(grid[ind] - newgrid[ind]));
        }
    }

    return maxdiff;
}

void halo_exchange(halo_cookie_t *cookie)
{
    int tag = 0;

    double *buf = cookie->buf;
    int count   = cookie->count;

    MPI_Datatype row = cookie->row;
    MPI_Datatype col = cookie->col;
    
    // Neighbors
    int top = cookie->top;
    int bottom = cookie->bottom;
    int left   = cookie->left;
    int right  = cookie->right;
    
    // Communicator
    MPI_Comm comm = cookie->comm;

    MPI_Request *reqs = cookie->reqs;

    int nx = cookie->nx;
    int ny = cookie->ny;

    MPI_Irecv(&buf[IND(0, 1)],      count, row, top,    tag, comm, &reqs[0]); // top
    MPI_Irecv(&buf[IND(ny + 1, 1)], count, row, bottom, tag, comm, &reqs[1]); // bottom
    MPI_Irecv(&buf[IND(1, 0)],      count, col, left,   tag, comm, &reqs[2]); // left
    MPI_Irecv(&buf[IND(1, nx + 1)], count, col, right,  tag, comm, &reqs[3]); // right

    MPI_Isend(&buf[IND(1, 1)],      count, row, top,    tag, comm, &reqs[4]); // top
    MPI_Isend(&buf[IND(ny, 1)],     count, row, bottom, tag, comm, &reqs[5]); // bottom
    MPI_Isend(&buf[IND(1, 1)],      count, col, left,   tag, comm, &reqs[6]); // left
    MPI_Isend(&buf[IND(1, nx)],     count, col, right,  tag, comm, &reqs[7]); // right
}

void redundancy_exchange(redundancy_cookie_t *cookie)
{
    int tag = 1;
    int cnt = 0;

    double *send      = cookie->send;
    double **receive  = cookie->receive;
    int *dest         = cookie->dest;
    int ranks         = cookie->ranks;
    int count         = cookie->count;
    MPI_Comm comm     = cookie->comm;
    MPI_Request *reqs = cookie->reqs;
    
    for (int i = 0; i < ranks; i++)
    {
        double *buf = receive[i];
        MPI_Irecv(buf, count, MPI_DOUBLE, dest[i], tag, comm, &reqs[cnt++]);
    }

    for (int i = 0; i < ranks; i++)
    {
        MPI_Isend(send, count, MPI_DOUBLE, dest[i], tag, comm, &reqs[cnt++]);
    }
}

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
    const int r_size = grid_task_redundancy_ranks_get(grid_task, rank, r_ranks);

    // Get redundancy local_grid and local_newgrid
    for (int i = 0; i < r_size; i++)
    {
        r_local_grid[i]    = grid_task_redundancy_local_grid_get(grid_task, r_ranks[i]);
        r_local_newgrid[i] = grid_task_redundancy_local_newgrid_get(grid_task, r_ranks[i]);
    }

    const int local_grid_size = (ny + 2) * (nx + 2);

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
            double x           = dx * (sj + j - 1);
            int ind            = IND(0, j);
            local_newgrid[ind] = local_grid[ind] = sin(PI * x);
        }
    }

    if (ranky == py - 1)
    {
        // Initialize bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
        for (int j = 1; j <= nx; j++)
        {
            // Translate col index to x coord in [0, 1]
            double x           = dx * (sj + j - 1);
            int ind            = IND(ny + 1, j);
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

    MPI_Request r_reqs[2 * r_size];
    MPI_Request reqs[8];
    double thalo   = 0;
    double treduce = 0;

    int niters = 0;

    while (1)
    {
        niters++;

        // Step 1: Update interior points
        update_interior_points(local_newgrid, local_grid, ny, nx);

        // Step 2: Check termination condition
        double maxdiff = check_termination_condition(local_newgrid, local_grid, ny, nx);

        // Step 3: Swap grids (after termination local_grid will contain result)
        double *p     = local_grid;
        local_grid    = local_newgrid;
        local_newgrid = p;

        treduce -= MPI_Wtime();
        // Step 4: All reduce (may fail)
        MPI_Allreduce(MPI_IN_PLACE, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        treduce += MPI_Wtime();

        if (maxdiff < EPS)
        {
            break;
        }

        halo_cookie_t halo_cookie = {
            .buf    = local_grid,
            .count  = 1,
            .row    = row,
            .col    = col,
            .top    = top,
            .bottom = bottom,
            .left   = left,
            .right  = right,
            .comm   = cartcomm,
            .reqs   = reqs,
            .nx     = nx,
            .ny     = ny 
        };

        redundancy_cookie_t redundancy_cookie = {
            .send    = local_grid,
            .receive = r_local_grid,
            .dest    = r_ranks,
            .ranks   = r_size,
            .count   = local_grid_size,
            .comm    = cartcomm,
            .reqs    = r_reqs
        };

        thalo -= MPI_Wtime();

        // Step 5: Halo exchange: T = 4 * (a + b * (rows / py)) + 4 * (a + b * (cols / px))
        halo_exchange(&halo_cookie);

        // Step 6: Wait all (may fail)
        MPI_Waitall(8, reqs, MPI_STATUS_IGNORE);

        // Step 7: Redundancy exchnage
        redundancy_exchange(&redundancy_cookie);

        // Step 8: Wait all (may fail)
        MPI_Waitall(2 * r_size, r_reqs, MPI_STATUS_IGNORE);

        thalo += MPI_Wtime();
    }

    MPI_Type_free(&row);
    MPI_Type_free(&col);

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

    double prof[3] = { ttotal, treduce, thalo };
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

    grid_task_free(grid_task);

    MPI_Finalize();
    return 0;
}

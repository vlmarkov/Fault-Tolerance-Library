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


#define EPS       0.001
#define PI        3.14159265358979323846
#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))
#define IND(i, j) ((i) * (nx + 2) + (j))

void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);

    if (!p)
    {
        fprintf(stderr, "No enough memory\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    return p;
}

int get_block_size(int n, int rank, int nprocs)
{
    int s = n / nprocs;

    if (n % nprocs > rank)
    {
        s++;
    }

    return s;
}

int get_sum_of_prev_blocks(int n, int rank, int nprocs)
{
    int rem = n % nprocs;
    return n / nprocs * rank + ((rank >= rem) ? rem : rank);
}

/*****************************************************************************************/
typedef enum
{
    DEAD_PROC = -2,
    BORDER    = -1,
    SHOW_REAL,
    SHOW_SHADOW,
} grid_task_e;

typedef struct
{
    int     rank;          // Assigned proc rank
    int     top;           // Neighbor top
    int     bottom;        // Neighbor down
    int     left;          // Neighbor left
    int     right;         // Neighbor right

    int     redudancy[16]; // Assigned proc rank todo
    int     r_counter;

    double *local_grid;
    double *local_newgrid;
} task_t;

typedef struct
{
    int      grid_height;
    int      grid_width;
    int      task_height;
    int      task_width;
    int      grid_height_proc;
    int      grid_width_proc;
    int      proc_per_node;

    task_t **tasks; // Grid[grid_height][grid_width]
} grid_task_t;

void grid_task_free(grid_task_t *grid)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Can't free memory\n", __FUNCTION__);
        exit(-1);        
    }

    for (int i = 1; i < grid->grid_height_proc; i++)
    {
        free(grid->tasks[i]);
    }

    free(grid->tasks);
    free(grid);
}

grid_task_t *grid_task_allocate(const int cols,
                                const int rows,
                                const int nx,
                                const int ny,
                                const int commsize)
{
    const int n = sqrt(commsize);
    const int m = sqrt(commsize);
    const int p = sqrt(commsize) / 2;

    if ((n < 4) || (m < 4) || (p < 2))
    {
        fprintf(stderr, "<%s> Not supproted grid size\n", __FUNCTION__);
        exit(-1);
    }

    grid_task_t *grid = (grid_task_t *)malloc(sizeof(grid_task_t));
    if (!grid)
    {
        fprintf(stderr, "<%s> Fail to allocate memory for grid\n",
            __FUNCTION__);
        exit(-1);
    }

    grid->grid_height      = cols;
    grid->grid_width       = rows;

    grid->grid_height_proc = n;
    grid->grid_width_proc  = m;

    grid->proc_per_node    = p;

    grid->task_height      = nx;
    grid->task_width       = ny;

    grid->tasks = (task_t **)malloc(sizeof(task_t *) * n);
    if (!grid->tasks)
    {
        grid_task_free(grid);
        fprintf(stderr, "<%s> Fail to allocate memory for tasks (real)\n",
            __FUNCTION__);
        exit(-1);
    }

    for (int i = 0; i < n; i++)
    {
        task_t *new_row = (task_t *)malloc(sizeof(task_t) * m);
        if (!new_row)
        {
            grid_task_free(grid);
            fprintf(stderr, "<%s> Fail to allocate memory for task\n", __FUNCTION__);
            exit(-1);
        }
        memset(new_row, 0, sizeof(sizeof(task_t) * m));
        
        grid->tasks[i] = new_row;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            task_t *task = &(grid->tasks[i][j]);
            task->local_grid    = xcalloc((ny + 2) * (nx + 2), sizeof(*task->local_grid));
            task->local_newgrid = xcalloc((ny + 2) * (nx + 2), sizeof(*task->local_newgrid));
        }
    }

    return grid;
}

void grid_task_assing_redudancy_ranks(grid_task_t *grid, int proc_num, int start_x, int start_y)
{
    const int x_times = grid->grid_width_proc  / grid->proc_per_node;
    const int y_times = grid->grid_height_proc / grid->proc_per_node;

    int ri = start_x + grid->proc_per_node;
    for (int x = 0; x < x_times; x++)
    {
        if (ri > grid->grid_height_proc - 1)
        {
            ri = ri - grid->grid_height_proc;
        }

        int rj = start_y + grid->proc_per_node;

        for (int y = 0; y < y_times; y++)
        {
            if (rj > grid->grid_width_proc - 1)
            {
                rj = rj - grid->grid_width_proc;
            }
            grid->tasks[ri][rj].redudancy[grid->tasks[ri][rj].r_counter++] = proc_num;
            
            rj += grid->proc_per_node;
        }

        ri += grid->proc_per_node;
    }
}

void grid_task_assing_top_neighbor(grid_task_t *grid, int x, int y)
{
    if (x == 0)
    {
        grid->tasks[x][y].top = BORDER;
    }
    else if (x > 0)
    {
        grid->tasks[x][y].top = grid->tasks[x - 1][y].rank;
    }
}

void grid_task_assing_bottom_neighbor(grid_task_t *grid, int x, int y)
{
    if (x == (grid->grid_height_proc - 1))
    {
        grid->tasks[x][y].bottom = BORDER;
    }
    else if (x < (grid->grid_height - 1))
    {
        grid->tasks[x][y].bottom = grid->tasks[x + 1][y].rank;
    }
}

void grid_task_assing_left_neighbor(grid_task_t *grid, int x, int y)
{
    if (y == 0)
    {
        grid->tasks[x][y].left = BORDER;
    }
    else if (y > 0)
    {
        grid->tasks[x][y].left = grid->tasks[x][y - 1].rank;
    }
}

void grid_task_assing_right_neighbor(grid_task_t *grid, int x, int y)
{
    if (y == (grid->grid_width_proc - 1))
    {
        grid->tasks[x][y].right = BORDER;
    }
    else if (y < (grid->grid_width_proc - 1))
    {
        grid->tasks[x][y].right = grid->tasks[x][y + 1].rank;
    }
}

static int grid_task_init(grid_task_t *grid_task)
{
    int counter = 0;

    if (!grid_task)
    {
        fprintf(stderr, "<%s> Can't init grid_task\n", __FUNCTION__);
        return -1;
    }

    if (!grid_task->tasks)
    {
        fprintf(stderr, "<%s> Can't init grid_task (real)\n", __FUNCTION__);
        return -1;
    }

    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            // Each struct task contains redudancy list
            grid_task_assing_redudancy_ranks(grid_task, counter, i, j);
            grid_task->tasks[i][j].rank = counter++;
        }
    }

    // Neighbors rank assing
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            grid_task_assing_top_neighbor(grid_task, i, j);
            grid_task_assing_bottom_neighbor(grid_task, i, j);
            grid_task_assing_left_neighbor(grid_task, i, j);
            grid_task_assing_right_neighbor(grid_task, i, j);
        }
    }

    return 0;
}

static double *grid_task_local_grid_get(grid_task_t *grid_task, int rank)
{
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            if (rank == grid_task->tasks[i][j].rank)
            {
                return grid_task->tasks[i][j].local_grid;
            }
        }
    }

    fprintf(stderr, "<%s>\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return NULL;
}

static double *grid_task_local_newgrid_get(grid_task_t *grid_task, int rank)
{
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            if (rank == grid_task->tasks[i][j].rank)
            {
                return grid_task->tasks[i][j].local_newgrid;
            }
        }
    }
    
    fprintf(stderr, "<%s>\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return NULL;
}

static int grid_task_bottom_neighbor_get(grid_task_t *grid_task, int rank)
{
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            if (rank == grid_task->tasks[i][j].rank)
            {
                return grid_task->tasks[i][j].bottom;
            }
        }
    }

    fprintf(stderr, "<%s>\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return -1;
}

static int grid_task_top_neighbor_get(grid_task_t *grid_task, int rank)
{
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            if (rank == grid_task->tasks[i][j].rank)
            {
                return grid_task->tasks[i][j].top;
            }
        }
    }

    fprintf(stderr, "<%s>\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return -1;
}

static int grid_task_left_neighbor_get(grid_task_t *grid_task, int rank)
{
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            if (rank == grid_task->tasks[i][j].rank)
            {
                return grid_task->tasks[i][j].left;
            }
        }
    }

    fprintf(stderr, "<%s>\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return -1;
}

static int grid_task_right_neighbor_get(grid_task_t *grid_task, int rank)
{
    for (int i = 0; i < grid_task->grid_height_proc; i++)
    {
        for (int j = 0; j < grid_task->grid_width_proc; j++)
        {
            if (rank == grid_task->tasks[i][j].rank)
            {
                return grid_task->tasks[i][j].right;
            }
        }
    }

    fprintf(stderr, "<%s>\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return -1;
}
/*****************************************************************************************/

int main(int argc, char *argv[])
{
    int rank;
    int commsize;

    MPI_Init(&argc, &argv);

    double ttotal = -MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create 2D grid of processes: commsize = px * py
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

    // Broadcast command line arguments
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

    // Allocate memory for local 2D subgrids with halo cells [0..ny + 1][0..nx + 1]
    int ny = get_block_size(rows, ranky, py);
    int nx = get_block_size(cols, rankx, px);

    /*
     * Depricated:
     * double *local_grid = xcalloc((ny + 2) * (nx + 2), sizeof(*local_grid));
     * double *local_newgrid = xcalloc((ny + 2) * (nx + 2), sizeof(*local_newgrid));
     */

    if (rank == 0)
    {
        printf("cols = %d\n", cols);
        printf("rows = %d\n", rows);
        printf("nx   = %d\n", nx);
        printf("ny   = %d\n", ny);
    }

    grid_task_t *grid_task = grid_task_allocate(cols, rows, nx, ny, commsize); // TODO

    grid_task_init(grid_task);

    double *local_grid    = grid_task_local_grid_get(grid_task, rank);
    double *local_newgrid = grid_task_local_newgrid_get(grid_task, rank);

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

    /*
     * Depricated:
     * int left, right, top, bottom;
     * MPI_Cart_shift(cartcomm, 0, 1, &left, &right);
     * MPI_Cart_shift(cartcomm, 1, 1, &top, &bottom);
     */
    
    int bottom = grid_task_bottom_neighbor_get(grid_task, rank);
    int top    = grid_task_top_neighbor_get(grid_task, rank);
    int left   = grid_task_left_neighbor_get(grid_task, rank);
    int right  = grid_task_right_neighbor_get(grid_task, rank);
    if (rank == 0)
    {
        printf("left = %d\n", left);
        printf("right = %d\n", right);
        printf("top   = %d\n", top);
        printf("bottom   = %d\n", bottom);
    }
    if (rank == 0)
    {
        int left, right, top, bottom;
        MPI_Cart_shift(cartcomm, 0, 1, &left, &right);
        MPI_Cart_shift(cartcomm, 1, 1, &top, &bottom);
        printf("left = %d\n", left);
        printf("right = %d\n", right);
        printf("top   = %d\n", top);
        printf("bottom   = %d\n", bottom);
    }

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
        double *p = local_grid;
        local_grid = local_newgrid;
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
    }

    MPI_Type_free(&row);
    MPI_Type_free(&col);

    free(local_newgrid);
    free(local_grid);

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

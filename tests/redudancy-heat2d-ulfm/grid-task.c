#include "grid-task.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/*****************************************************************************/
/* Main destruction function                                                 */
/*****************************************************************************/
void grid_task_free(grid_task_t *grid)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 1; i < grid->rows_per_proc; i++)
    {
        free(grid->tasks[i]);
    }

    free(grid->tasks);
    free(grid);
}

/*****************************************************************************/
/* Main construction function                                                */
/*****************************************************************************/
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
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    grid_task_t *grid = (grid_task_t *)malloc(sizeof(grid_task_t));
    if (!grid)
    {
        fprintf(stderr, "<%s> Fail to allocate memory for grid\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    grid->cols          = cols;
    grid->rows          = rows;

    grid->cols_per_task = nx;
    grid->rows_per_task = ny;

    grid->cols_per_proc = n;
    grid->rows_per_proc = m;

    grid->proc_per_node = p;

    grid->tasks = (task_t **)malloc(sizeof(task_t *) * n);
    if (!grid->tasks)
    {
        grid_task_free(grid);
        fprintf(stderr, "<%s> Fail to allocate memory for tasks\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++)
    {
        task_t *new_row = (task_t *)malloc(sizeof(task_t) * m);
        if (!new_row)
        {
            grid_task_free(grid);
            fprintf(stderr, "<%s> Fail to allocate memory for task\n", __FUNCTION__);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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

/*****************************************************************************/
/* Main intialize function                                                   */
/*****************************************************************************/
void grid_task_init(grid_task_t *grid)
{
    int counter = 0;

    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            // Each struct task contains redudancy list
            grid_task_redundancy_ranks_set(grid, counter, i, j);
            grid->tasks[i][j].rank = counter++;
        }
    }

    // Neighbors rank set
    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            grid_task_neighbor_top_set(grid, i, j);
            grid_task_neighbor_bottom_set(grid, i, j);
            grid_task_neighbor_left_set(grid, i, j);
            grid_task_neighbor_right_set(grid, i, j);
        }
    }
}

/*****************************************************************************/
/* Redundancy ranks setter                                                   */
/*****************************************************************************/
void grid_task_redundancy_ranks_set(grid_task_t *grid,
                                    const int    rank,
                                    const int    row,
                                    const int    col)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const int x_times = grid->cols_per_proc / grid->proc_per_node;
    const int y_times = grid->rows_per_proc / grid->proc_per_node;

    int ri = row + grid->proc_per_node;
    for (int x = 0; x < x_times; x++)
    {
        if (ri > grid->cols_per_proc - 1)
        {
            ri = ri - grid->cols_per_proc;
        }

        int rj = col + grid->proc_per_node;
        for (int y = 0; y < y_times; y++)
        {
            if (rj > grid->rows_per_proc - 1)
            {
                rj = rj - grid->rows_per_proc;
            }
            grid->tasks[ri][rj].redundancy[grid->tasks[ri][rj].r_counter++] = rank;
            
            rj += grid->proc_per_node;
        }

        ri += grid->proc_per_node;
    }
}

/*****************************************************************************/
/* Local grid getter                                                         */
/*****************************************************************************/
double *grid_task_local_grid_get(const grid_task_t *grid, const int rank)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                return grid->tasks[i][j].local_grid;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find local-grid\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return NULL; // Never reached
}

/*****************************************************************************/
/* Local new-grid getter                                                     */
/*****************************************************************************/
double *grid_task_local_newgrid_get(const grid_task_t *grid, const int rank)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                return grid->tasks[i][j].local_newgrid;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find local-newgrid\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return NULL; // Never reached
}

/*****************************************************************************/
/* Neighbor setters                                                          */
/*****************************************************************************/

/*****************************************************************************/
/* Top neighbor setter                                                       */
/*****************************************************************************/
void grid_task_neighbor_top_set(grid_task_t *grid, const int x, const int y)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (x == 0)
    {
        grid->tasks[x][y].top = BORDER;
    }
    else if (x > 0)
    {
        grid->tasks[x][y].top = grid->tasks[x - 1][y].rank;
    }

    return;
}

/*****************************************************************************/
/* Bottom neighbor setter                                                    */
/*****************************************************************************/
void grid_task_neighbor_bottom_set(grid_task_t *grid, const int x, const int y)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (x == (grid->cols_per_proc - 1))
    {
        grid->tasks[x][y].bottom = BORDER;
    }
    else if (x < (grid->cols_per_proc - 1))
    {
        grid->tasks[x][y].bottom = grid->tasks[x + 1][y].rank;
    }

    return;
}

/*****************************************************************************/
/* Left neighbor setter                                                      */
/*****************************************************************************/
void grid_task_neighbor_left_set(grid_task_t *grid, const int x, const int y)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (y == 0)
    {
        grid->tasks[x][y].left = BORDER;
    }
    else if (y > 0)
    {
        grid->tasks[x][y].left = grid->tasks[x][y - 1].rank;
    }

    return;
}

/*****************************************************************************/
/* Right neighbor setter                                                     */
/*****************************************************************************/
void grid_task_neighbor_right_set(grid_task_t *grid, const int x, const int y)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (y == (grid->rows_per_proc - 1))
    {
        grid->tasks[x][y].right = BORDER;
    }
    else if (y < (grid->rows_per_proc - 1))
    {
        grid->tasks[x][y].right = grid->tasks[x][y + 1].rank;
    }

    return;
}

/*****************************************************************************/
/* Neighbor getters                                                          */
/*****************************************************************************/

/*****************************************************************************/
/* Top neighbor getter                                                       */
/*****************************************************************************/
int grid_task_neighbor_top_get(const grid_task_t *grid, const int rank)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                return grid->tasks[i][j].top;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find neighbor\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return -1; // Never reached
}

/*****************************************************************************/
/* Bottom neighbor getter                                                    */
/*****************************************************************************/
int grid_task_neighbor_bottom_get(const grid_task_t *grid, const int rank)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                return grid->tasks[i][j].bottom;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find neighbor\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return -1; // Never reached
}

/*****************************************************************************/
/* Left neighbor getter                                                      */
/*****************************************************************************/
int grid_task_neighbor_left_get(const grid_task_t *grid, const int rank)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                return grid->tasks[i][j].left;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find neighbor\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return -1; // Never reached
}

/*****************************************************************************/
/* Right neighbor getter                                                     */
/*****************************************************************************/
int grid_task_neighbor_right_get(const grid_task_t *grid, const int rank)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                return grid->tasks[i][j].right;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find neighbor\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return -1; //  Never reached
}

/*****************************************************************************/
/* All neighbors getter                                                      */
/*****************************************************************************/
void grid_task_neighbors_get(const grid_task_t *grid,
                             const int          rank,
                             int               *top,
                             int               *bottom,
                             int               *left,
                             int               *right)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Bad grid pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Bad tasks pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                *top    = grid->tasks[i][j].top;
                *bottom = grid->tasks[i][j].bottom;
                *left   = grid->tasks[i][j].left;
                *right  = grid->tasks[i][j].right;
                return;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find neighbor\n", __FUNCTION__);

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return; // Never reached
}

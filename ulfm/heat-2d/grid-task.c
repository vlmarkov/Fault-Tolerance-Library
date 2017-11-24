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

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            free(grid->tasks[i][j].local_grid);
            free(grid->tasks[i][j].local_newgrid);
            free(grid->tasks[i][j].redundancy_task);
            free(grid->tasks[i][j].real_task);
        }
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

    // Allocate memory for grid-tasks
    grid->tasks = (task_t **)malloc(sizeof(task_t *) * n);
    if (!grid->tasks)
    {
        grid_task_free(grid);
        fprintf(stderr, "<%s> Fail to allocate memory for tasks\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Allocate memory for each row in grid-tasks
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

    /*
     * Allocate memory for local_grid and
     * local_newgrid and redundancy_task in each task
     */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            task_t *task          = &(grid->tasks[i][j]);

            task->local_grid      = xcalloc((ny + 2) * (nx + 2), sizeof(*task->local_grid));

            task->local_newgrid   = xcalloc((ny + 2) * (nx + 2), sizeof(*task->local_newgrid));

            task->redundancy_task = (task_t **)malloc(sizeof(task_t *) * commsize);
            if (!task->redundancy_task)
            {
                grid_task_free(grid);
                fprintf(stderr, "<%s> Fail to allocate memory for redundancy-tasks\n", __FUNCTION__);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            memset(task->redundancy_task, 0, sizeof(sizeof(task_t *) * commsize));

            task->real_task = (task_t **)malloc(sizeof(task_t *) * commsize);
            if (!task->real_task)
            {
                grid_task_free(grid);
                fprintf(stderr, "<%s> Fail to allocate memory for redundancy-tasks\n", __FUNCTION__);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            memset(task->real_task, 0, sizeof(sizeof(task_t *) * commsize));
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

    int commsize = (grid->cols_per_proc * grid->rows_per_proc);

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            grid->tasks[i][j].redundancy_counter  = 0;
            grid->tasks[i][j].redundancy_capacity = commsize;
            grid->tasks[i][j].real_counter        = 1;
            grid->tasks[i][j].real_capacity       = commsize;
            grid->tasks[i][j].real_task[0]        = &(grid->tasks[i][j]);
            grid->tasks[i][j].rank                = counter++;
            grid->tasks[i][j].x                   = i;
            grid->tasks[i][j].y                   = j;
        }
    }

    counter = 0;

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            // Each struct task contains redudancy list
            grid_task_redundancy_task_set(grid, counter++, i, j);
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
/* Main kill process's task fucntion                                         */
/*****************************************************************************/
void grid_task_kill_rank(grid_task_t *grid, const int rank)
{
    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (grid->tasks[i][j].rank == rank)
            {
                grid->tasks[i][j].rank = DEAD_PROC;
            }
        }
    }
}

/*****************************************************************************/
/* Main repair function                                                      */
/*****************************************************************************/
void grid_task_repair(grid_task_t *grid)
{
    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (grid->tasks[i][j].rank == DEAD_PROC)
            {
                task_t *dead_task = &grid->tasks[i][j];
                for (int r = 0; r < dead_task->redundancy_counter; r++)
                {
                    task_t *repair_task = dead_task->redundancy_task[r];

                    if (repair_task->rank != DEAD_PROC)
                    {
                        const int counter  = repair_task->real_counter;
                        const int capacity = repair_task->real_capacity;
                        if (counter == capacity)
                        {
                            fprintf(stderr, "<%s> Repair failed\n", __FUNCTION__);
                            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                        }

                        dead_task->rank = repair_task->rank;
                        repair_task->real_task[counter] = dead_task;
                        repair_task->real_counter++;

                        break;
                    }
                }
            }
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
void grid_task_redundancy_task_set(grid_task_t *grid,
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

            const int capacity = grid->tasks[ri][rj].redundancy_capacity;
            const int idx      = grid->tasks[ri][rj].redundancy_counter;

            if (idx == capacity)
            {
                fprintf(stderr, "<%s> Bad redundancy-counter\n", __FUNCTION__);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            grid->tasks[ri][rj].redundancy_task[idx] = grid_task_get(grid, rank);
            grid->tasks[ri][rj].redundancy_counter  += 1;

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

/*****************************************************************************/
/* Show/print functions                                                      */
/*****************************************************************************/

void grid_task_redundancy_ranks_send_show(const task_t *my_task)
{
    if (!my_task)
    {
        fprintf(stderr, "<%s> Bad my-task pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int rc = 0;
    char buf[256] = { 0 };
    int size = my_task->redundancy_counter;
    for (int r = 0; r < size; r++)
    {
        // Ignore self
        if (((my_task->x    == my_task->redundancy_task[r]->x)  &&
             (my_task->y    == my_task->redundancy_task[r]->y)) ||
              my_task->rank == my_task->redundancy_task[r]->rank)
        {
            continue;
        }
        else
        {
            rc += sprintf(buf + rc, "%04d ", my_task->redundancy_task[r]->rank);
        }
    }

    printf("I'am rank %04d, redundancy ranks to send   : %s\n", my_task->rank, buf);
    return;
}

void grid_task_redundancy_ranks_receive_show(const task_t *my_task)
{
    if (!my_task)
    {
        fprintf(stderr, "<%s> Bad my-task pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int rc = 0;
    char buf[256] = { 0 };
    int size = my_task->redundancy_counter;
    for (int r = 0; r < size; r++)
    {
        // Ignore self
        if (((my_task->x    == my_task->redundancy_task[r]->x)  &&
             (my_task->y    == my_task->redundancy_task[r]->y)) ||
              my_task->rank == my_task->redundancy_task[r]->rank)
        {
            continue;
        }
        else
        {
            rc += sprintf(buf + rc, "%04d ", my_task->redundancy_task[r]->rank);
        }
    }

    printf("I'am rank %04d, redundancy ranks to receive: %s\n", my_task->rank, buf);
    return;
}

static void grid_task_show_(const grid_task_t *grid, const int redundancy)
{
    const int      n     = grid->cols_per_proc;
    const int      m     = grid->rows_per_proc;
    const int      p     = grid->proc_per_node;
    const task_t **tasks = (const task_t **)grid->tasks;

    for (int i = 0; i < n; i++)
    {
        if (i % p == 0)
        {
            printf("\n");
        }
        for (int j = 0; j < m; j++)
        {
            if (j % p == 0)
            {
                printf(" ");
            }

            if (!redundancy)
            {
                printf("%04d ", tasks[i][j].rank);
            }
            else
            {
                int size = tasks[i][j].redundancy_counter;
                printf("[%04d]:{", tasks[i][j].rank);
                for (int r = 0; r < size; r++)
                {
                    printf("%04d ", tasks[i][j].redundancy_task[r]->rank);
                }
                printf("}");
            }
        }
        printf("\n");
    }
    printf("\n");
}

void grid_task_show(const grid_task_t *grid)
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

    grid_task_show_(grid, 0);

    return;
}

void grid_task_redundancy_show(const grid_task_t *grid)
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

    grid_task_show_(grid, 1);

    return;
}

void grid_task_task_show(const grid_task_t *grid)
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

    const task_t **tasks = (const task_t **)grid->tasks;

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            printf("Task->coordinates   : %04d:%04d\n", i, j);
            printf("Task->rank          : %04d\n", tasks[i][j].rank);
            printf("Task->top           : %04d\n", tasks[i][j].top);
            printf("Task->bottom        : %04d\n", tasks[i][j].bottom);
            printf("Task->left          : %04d\n", tasks[i][j].left);
            printf("Task->right         : %04d\n", tasks[i][j].right);
            printf("Task->red_capacity  : %04d\n", tasks[i][j].redundancy_capacity);
            printf("Task->red_counter   : %04d\n", tasks[i][j].redundancy_counter);
            printf("Task->real_capacity : %04d\n", tasks[i][j].real_capacity);
            printf("Task->real_counter  : %04d\n", tasks[i][j].real_counter);

            printf("Task->red_task      : { ");
            for (int r = 0; r < tasks[i][j].redundancy_counter; r++)
                printf("%04d ", tasks[i][j].redundancy_task[r]->rank);
            printf("}\n");

            printf("Task->real_task     : { ");
            for (int r = 0; r < tasks[i][j].real_counter; r++)
                printf("%04d ", tasks[i][j].real_task[r]->rank);
            printf("}\n"); 

            printf("Task->real_task*    : { ");
            for (int r = 0; r < tasks[i][j].real_counter; r++)
                printf("%04d(%04d:%04d) ", tasks[i][j].real_task[r]->rank,
                    tasks[i][j].real_task[r]->x, tasks[i][j].real_task[r]->y);
            printf("}\n"); 
            printf("\n");
        }
        printf("------------------------\n");
    }
    printf("\n");
}

/*****************************************************************************/
/* Redundancy ranks getter                                                   */
/*****************************************************************************/
int grid_task_redundancy_ranks_get(const grid_task_t *grid,
                                   const int          rank,
                                   int               *redundancy)
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

    int cnt = 0;

    for (int i = 0; i < grid->cols_per_proc; i++)
    {
        for (int j = 0; j < grid->rows_per_proc; j++)
        {
            if (rank == grid->tasks[i][j].rank)
            {
                int size = grid->tasks[i][j].redundancy_counter;
                for (int r = 0; r < size; r++)
                {
                    // Not include self
                    if (rank != grid->tasks[i][j].redundancy_task[r]->rank)
                    {
                        redundancy[cnt] = grid->tasks[i][j].redundancy_task[r]->rank;
                        cnt++;
                    }
                }
                return cnt;
            }
        }
    }

    fprintf(stderr, "<%s> Can't find neighbor\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return -1; // Never reached
}

/*****************************************************************************/
/* Redundancy local-grid getter                                              */
/*****************************************************************************/
double *grid_task_redundancy_local_grid_get(const grid_task_t *grid, const int rank)
{
    return grid_task_local_grid_get(grid, rank);
}

/*****************************************************************************/
/* Redundancy local-newgrid getter                                           */
/*****************************************************************************/
double *grid_task_redundancy_local_newgrid_get(const grid_task_t *grid,
                                               const int          rank)
{
    return grid_task_local_newgrid_get(grid, rank);
}

/*****************************************************************************/
/* Task getters                                                              */
/*****************************************************************************/

/*****************************************************************************/
/* Process's task getter                                                     */
/*****************************************************************************/
task_t *grid_task_get(const grid_task_t *grid, const int rank)
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
                return &grid->tasks[i][j];
            }
        }
    }

    fprintf(stderr, "<%s> Can't find task\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return NULL; // Never reached
}

/*****************************************************************************/
/* Redundancy task getter                                                    */
/*****************************************************************************/
int grid_task_redundancy_task_get(const task_t *task,
                                  int *ranks,
                                  double **grid,
                                  double **newgrid)
{
    if (!task)
    {
        fprintf(stderr, "<%s> Bad task pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < task->redundancy_capacity; i++)
    {
        ranks[i]   = 0;
        grid[i]    = 0;
        newgrid[i] = 0;
    }

    int counter = 0;
    for (int i = 0; i < task->redundancy_counter; i++)
    {
        /*
         * TODO
        // Ignore self
        if (((task->x    == task->redundancy_task[i]->x)  &&
             (task->y    == task->redundancy_task[i]->y)) ||
              task->rank == task->redundancy_task[i]->rank)
        {
            continue;
        }
        else
        {
            ranks[counter]   = task->redundancy_task[i]->rank;
            grid[counter]    = task->redundancy_task[i]->local_grid;
            newgrid[counter] = task->redundancy_task[i]->local_newgrid;
            counter++;
        }
        */

        ranks[counter]   = task->redundancy_task[i]->rank;
        grid[counter]    = task->redundancy_task[i]->local_grid;
        newgrid[counter] = task->redundancy_task[i]->local_newgrid;
        counter++;
    }

    return counter;
}

int grid_task_redundancy_counter_get(const task_t *task)
{
    if (!task)
    {
        fprintf(stderr, "<%s> Bad task pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    return task->redundancy_counter;
}

/*****************************************************************************/
/* Real task getter                                                          */
/*****************************************************************************/
int grid_task_real_task_get(const task_t *my_task, task_t **tasks)
{
    if (!my_task)
    {
        fprintf(stderr, "<%s> Bad task pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int size = my_task->real_counter;

    for (int i = 0; i < size; i++)
    {
        tasks[i] = my_task->real_task[i];
    }

    return size;
}

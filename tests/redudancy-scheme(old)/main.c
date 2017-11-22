#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "colors.h"

typedef enum
{
    BORDER    = -2,
    DEAD_PROC = -1,
    SHOW_REAL,
    SHOW_SHADOW,
} grid_task_e;

typedef struct
{
    int     rank;          // Assigned proc rank
    int     top;           // Neighbor top
    int     down;          // Neighbor down
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
    int      proc_per_node;

    task_t **tasks; // Grid[grid_height][grid_width]
} grid_task_t;


void *xcalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);
    if (p == NULL) {
        fprintf(stderr, "No enough memory\n");
        //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(-1);
    }
    return p;
}

/*****************************************************************************/
/* Basic destructor                                                          */
/*****************************************************************************/
void grid_task_free(grid_task_t *grid)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Can't free memory\n", __FUNCTION__);
        exit(-1);        
    }

    for (int i = 1; i < grid->grid_height; i++)
    {
        free(grid->tasks[i]);
    }

    free(grid->tasks);
    free(grid);
}

/*****************************************************************************/
/* Basic constructor                                                         */
/*****************************************************************************/
grid_task_t *grid_task_allocate(const int n,
                                const int m,
                                const int nx,
                                const int ny,
                                const int p)
{
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

    grid->grid_height   = n;
    grid->grid_width    = m;
    grid->proc_per_node = p;

    grid->task_height   = nx; // TODO
    grid->task_width    = ny; // TODO

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
    const int x_times = grid->grid_height / grid->proc_per_node;
    const int y_times = grid->grid_width  / grid->proc_per_node;

    int ri = start_x + grid->proc_per_node;
    for (int x = 0; x < x_times; x++)
    {
        if (ri > grid->grid_height - 1)
        {
            ri = ri - grid->grid_height;
        }

        int rj = start_y + grid->proc_per_node;

        for (int y = 0; y < y_times; y++)
        {
            if (rj > grid->grid_width - 1)
            {
                rj = rj - grid->grid_width;
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

void grid_task_assing_down_neighbor(grid_task_t *grid, int x, int y)
{
    if (x == (grid->grid_height - 1))
    {
        grid->tasks[x][y].down = BORDER;
    }
    else if (x < (grid->grid_height - 1))
    {
        grid->tasks[x][y].down = grid->tasks[x + 1][y].rank;
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
    if (y == (grid->grid_width - 1))
    {
        grid->tasks[x][y].right = BORDER;
    }
    else if (y < (grid->grid_width - 1))
    {
        grid->tasks[x][y].right = grid->tasks[x][y + 1].rank;
    }
}


/*****************************************************************************/
/* Internal initialize function (real)                                       */
/*****************************************************************************/
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

    for (int i = 0; i < grid_task->grid_height; i++)
    {
        for (int j = 0; j < grid_task->grid_width; j++)
        {
            // Each struct task contains redudancy list
            grid_task_assing_redudancy_ranks(grid_task, counter, i, j);
            grid_task->tasks[i][j].rank = counter++;
        }
    }

    // Neighbors rank assing
    for (int i = 0; i < grid_task->grid_height; i++)
    {
        for (int j = 0; j < grid_task->grid_width; j++)
        {
            grid_task_assing_top_neighbor(grid_task, i, j);
            grid_task_assing_down_neighbor(grid_task, i, j);
            grid_task_assing_left_neighbor(grid_task, i, j);
            grid_task_assing_right_neighbor(grid_task, i, j);
        }
    }

/*
    init_redudancy_by_rank(grid_task, 0); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 1); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 2); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 3); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 4); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 5); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 6); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 7); // TODO mpi rank
*/
/*
    init_redudancy_by_rank(grid_task, 8, 0, 8); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 9, 0, 9); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 10, 0, 10); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 11, 0, 11); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 12, 0, 12); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 13, 0, 13); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 14, 0, 14); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 15, 0, 15); // TODO mpi rank

    init_redudancy_by_rank(grid_task, 144, 9, 0); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 145, 9, 1); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 146, 9, 2); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 147, 9, 3); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 148, 9, 4); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 149, 9, 5); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 150, 9, 6); // TODO mpi rank
    init_redudancy_by_rank(grid_task, 151, 9, 7); // TODO mpi rank
*/
/*
    int x = 15, y = 15;
    printf("[%04d:%04d].rank %d top %d left %d down %d right %d\n",
        x, y, grid_task->tasks[x][y].rank,
        grid_task->tasks[x][y].top,
        grid_task->tasks[x][y].left,
        grid_task->tasks[x][y].down,
        grid_task->tasks[x][y].right);
*/
    return 0;
}

/*****************************************************************************/
/* Internal show function                                                    */
/*****************************************************************************/
static void grid_task_show_(const grid_task_t *grid, const int colored, const int real)
{
    const int      n     = grid->grid_height;
    const int      m     = grid->grid_width;
    const int      proc  = grid->proc_per_node;
    const task_t **tasks = (const task_t **)grid->tasks;

    for (int i = 0; i < n; i++)
    {
        if (i % 8 == 0)
        {
            printf("\n");
        }
        for (int j = 0; j < m; j++)
        {
            if (j % 8 == 0)
            {
                printf(" ");
            }

            if (colored)
            {
                if (tasks[i][j].rank < proc)
                {
                    if (real)
                    {
                        printf(RED("%04d "), tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf(RED("%04d "), tasks[i][j].redudancy[r]);
                        }
                    }
                }
                else if ((tasks[i][j].rank >= proc) &&
                         (tasks[i][j].rank < (2 * proc)))
                {
                    if (real)
                    {
                        printf(GRN("%04d "), tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf(GRN("%04d "), tasks[i][j].redudancy[r]);
                        }
                    }
                }
                else if ((tasks[i][j].rank >= (2 * proc)) &&
                         (tasks[i][j].rank < (3 * proc)))
                {
                    if (real)
                    {
                        printf(YEL("%04d "), tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf(YEL("%04d "), tasks[i][j].redudancy[r]);
                        }                    
                    }
                }
                else if ((tasks[i][j].rank >= (3 * proc)) &&
                         (tasks[i][j].rank < (4 * proc)))
                {
                    if (real)
                    {
                        printf(BLU("%04d "), tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf(BLU("%04d "), tasks[i][j].redudancy[r]);
                        }
                    }
                }
                else if ((tasks[i][j].rank >= (4 * proc)) &&
                         (tasks[i][j].rank < (5 * proc)))
                {
                    if (real)
                    {
                        printf(MAG("%04d "), tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf(MAG("%04d "), tasks[i][j].redudancy[r]);
                        }
                    }
                }
                else if ((tasks[i][j].rank >= (5 * proc)) &&
                         (tasks[i][j].rank < (6 * proc)))
                {
                    if (real)
                    {
                        printf(CYN("%04d "), tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf(CYN("%04d "), tasks[i][j].redudancy[r]);
                        }
                    }
                }
                else
                {
                    if (real)
                    {
                        printf("%04d ", tasks[i][j].rank);
                    }
                    else
                    {
                        for (int r = 0; r < tasks[i][j].r_counter; r++)
                        {
                            printf("%04d ", tasks[i][j].redudancy[r]);
                        }
                    }
                }
            }
            else
            {
                if (real)
                {
                    printf("%04d ", tasks[i][j].rank);
                }
                else
                {
                    printf("{");
                    for (int r = 0; r < tasks[i][j].r_counter; r++)
                    {
                        printf("%04d ", tasks[i][j].redudancy[r]);
                    }
                    printf("}");
                }
            }
        }
        printf("\n");
    }
    printf("\n");
}

/*****************************************************************************/
/* Simple show function                                                      */
/*****************************************************************************/
int grid_task_show_real(grid_task_t *grid, int colored)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Can't show grid_task\n", __FUNCTION__);
        return -1;
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Can't show grid_task (real)\n", __FUNCTION__);
        return -1;
    }
    
    grid_task_show_(grid, colored, 1);

    return 0;
}

int grid_task_show_redudancy(grid_task_t *grid, int colored)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Can't show grid_task\n", __FUNCTION__);
        return -1;
    }

    if (!grid->tasks)
    {
        fprintf(stderr, "<%s> Can't show grid_task (real)\n", __FUNCTION__);
        return -1;
    }
    
    grid_task_show_(grid, colored, 0);

    return 0;
}


void test(const int row, const int col, const int nx, const int ny, const int commsize)
{
    printf("[INFO] Test: start!\n");

    grid_task_t *grid_task = grid_task_allocate(row, col, nx, ny, commsize);

    grid_task_init(grid_task);
    grid_task_show_real(grid_task, 1);
    grid_task_show_redudancy(grid_task, 0);
    grid_task_free(grid_task);

    printf("[INFO] Test: done!\n");
}

int main(int argc, char const *argv[])
{
    int x = 4, y = 4, n = 2;
    int nx = 4, ny = 4;

    test(x, y, nx, ny, n);

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "colors.h"

typedef enum
{
    DEAD_PROC = -1,
    SHOW_REAL,
    SHOW_SHADOW
} grid_task_e;

typedef struct
{
    bool    reconfigure; // TODO
    int     rank;        // Assigned proc rank
    int     number;      // Block number
    size_t  x;           // Size of block by x
    size_t  y;           // Size of block by y
    void   *redudancy;   // Need to reconfigure grid_task ???
} task_t;

typedef struct
{
    int      n;       // Grid size by x
    int      m;       // Grid size by x
    int      proc;    // Proc size by node
    task_t **tasks;   // Grid[n][m] (real)
    task_t **tasks_r; // Grid[n][m] (redudancy)
} grid_task_t;

/*****************************************************************************/
/* Basic constructor                                                         */
/*****************************************************************************/
grid_task_t *grid_task_allocate(const int n, const int m, const int p)
{
    int i;

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

    grid->n    = n;
    grid->m    = m;
    grid->proc = p;

    grid->tasks = (task_t **)malloc(sizeof(task_t *) * n); // real
    if (!grid->tasks)
    {
        free(grid);
        fprintf(stderr, "<%s> Fail to allocate memory for tasks (real)\n",
            __FUNCTION__);
        exit(-1);
    }

    grid->tasks_r = (task_t **)malloc(sizeof(task_t *) * n); // redudancy
    if (!grid->tasks_r)
    {
        free(grid);
        fprintf(stderr, "<%s> Fail to allocate memory for tasks (redudancy)\n",
            __FUNCTION__);
        exit(-1);
    }

    for (i = 0; i < n; i++)
    {
        task_t *new_row = (task_t *)malloc(sizeof(task_t) * m); // real
        if (!new_row)
        {
            // TODO: memory leak !!
            free(grid);
            fprintf(stderr, "<%s> Fail to allocate memory for task\n", __FUNCTION__);
            exit(-1);
        }
        
        task_t *new_row_r = (task_t *)malloc(sizeof(task_t) * m); // redudancy
        if (!new_row_r)
        {
            // TODO: memory leak !!
            free(grid);
            fprintf(stderr, "<%s> Fail to allocate memory for task\n", __FUNCTION__);
            exit(-1);
        }

        memset(new_row, 0, sizeof(sizeof(task_t) * m));
        memset(new_row_r, 0, sizeof(sizeof(task_t) * m));
        
        grid->tasks[i]   = new_row;   // real
        grid->tasks_r[i] = new_row_r; // redudancy
    }

    return grid;
}

/*****************************************************************************/
/* Basic destructor                                                          */
/*****************************************************************************/
void grid_task_free(grid_task_t *grid)
{
    int i;

    if(!grid)
    {
        fprintf(stderr, "<%s> Can't free memory\n", __FUNCTION__);
        exit(-1);        
    }

    for (i = 0; i < grid->n; i++)
    {
        free(grid->tasks[i]);   // real
        free(grid->tasks_r[i]); // redudancy
    }

    free(grid->tasks);   // real
    free(grid->tasks_r); // redudancy
    free(grid);
}

/*****************************************************************************/
/* Internal show function                                                    */
/*****************************************************************************/
static void grid_task_show_(const task_t **tasks,
                            const int      n,
                            const int      m,
                            const int      proc,
                            const int      color)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (color)
            {
                if (tasks[i][j].rank < proc)
                {
                    printf(RED("%04d "), tasks[i][j].rank);
                }
                else if ((tasks[i][j].rank >= proc) &&
                         (tasks[i][j].rank < (2 * proc)))
                {
                    printf(GRN("%04d "), tasks[i][j].rank);
                }
                else if ((tasks[i][j].rank >= (2 * proc)) &&
                         (tasks[i][j].rank < (3 * proc)))
                {
                    printf(YEL("%04d "), tasks[i][j].rank);
                }
                else if ((tasks[i][j].rank >= (3 * proc)) &&
                         (tasks[i][j].rank < (4 * proc)))
                {
                    printf(BLU("%04d "), tasks[i][j].rank);
                }
                else if ((tasks[i][j].rank >= (4 * proc)) &&
                         (tasks[i][j].rank < (5 * proc)))
                {
                    printf(MAG("%04d "), tasks[i][j].rank);
                }
                else if ((tasks[i][j].rank >= (5 * proc)) &&
                         (tasks[i][j].rank < (6 * proc)))
                {
                    printf(CYN("%04d "), tasks[i][j].rank);
                }
                else
                {
                    printf("%04d ", tasks[i][j].rank);
                }
            }
            else
            {
                printf("%04d ", tasks[i][j].rank);
            }
        }
        printf("\n");
    }
    printf("\n");
}

static void grid_task_show_reconfigure_(const task_t **tasks,
    const int n, const int m, const int proc)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%04d[%c] ", tasks[i][j].rank, (tasks[i][j].reconfigure ? 'F' : 'G'));
        }
        printf("\n");
    }
    printf("\n");
}

/*****************************************************************************/
/* Color show function                                                       */
/*****************************************************************************/
int grid_task_show_color(grid_task_t *grid)
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

    if (!grid->tasks_r)
    {
        fprintf(stderr, "<%s> Can't show grid_task (redudancy)\n", __FUNCTION__);
        return -1;
    }

    const task_t **real      = (const task_t **)grid->tasks;
    const task_t **redudancy = (const task_t **)grid->tasks_r;

    printf("<%s> %s grid task\n", __FUNCTION__, "Real");
    grid_task_show_(real, grid->n, grid->m, grid->proc, 1);

    printf("<%s> %s grid task\n", __FUNCTION__, "Redudancy");
    grid_task_show_(redudancy, grid->n, grid->m, grid->proc, 1);

    return 0;
}

/*****************************************************************************/
/* Simple show function                                                      */
/*****************************************************************************/
int grid_task_show(grid_task_t *grid)
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

    if (!grid->tasks_r)
    {
        fprintf(stderr, "<%s> Can't show grid_task (redudancy)\n", __FUNCTION__);
        return -1;
    }

    const task_t **real      = (const task_t **)grid->tasks;
    const task_t **redudancy = (const task_t **)grid->tasks_r;

    printf("<%s> %s grid task\n", __FUNCTION__, "Real");
    grid_task_show_(real, grid->n, grid->m, grid->proc, 0);

    printf("<%s> %s grid task\n", __FUNCTION__, "Redudancy");
    grid_task_show_(redudancy, grid->n, grid->m, grid->proc, 0);

    return 0;
}

/*****************************************************************************/
/* Internal initialize function (real)                                       */
/*****************************************************************************/
static int grid_task_init_(grid_task_t *grid_task)
{
    int i, j, counter = 0;

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

    for (i = 0; i < grid_task->n; i++)
    {
        for (j = 0; j < grid_task->m; j++)
        {
            grid_task->tasks[i][j].reconfigure = false;
            grid_task->tasks[i][j].redudancy   = NULL;
            grid_task->tasks[i][j].x           = -1; // TODO
            grid_task->tasks[i][j].y           = -1; // TODO
            grid_task->tasks[i][j].rank        = counter;
            grid_task->tasks[i][j].number      = counter++;
        }
    }

    return 0;
}

/*****************************************************************************/
/* Internal initialize function (redudancy)                                  */
/*****************************************************************************/
static int grid_task_init_redudancy_(grid_task_t *grid_task)
{
    int i;  // real i-th index
    int j;  // real j-th index
    int ri; // redudancy i-th index
    int rj; // redudancy j-th index

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

        if (!grid_task->tasks_r)
    {
        fprintf(stderr, "<%s> Can't init grid_task (redudancy)\n", __FUNCTION__);
        return -1;
    }

    for (i = 0; i < grid_task->n; i++)
    {
        ri = i + grid_task->proc;
        
        if (ri >= grid_task->n)
        {
            ri = abs(grid_task->n - ri);
        }

        for (j = 0, rj = 0; j < grid_task->m; j++, rj += 2)
        {
            if (j == 0)
            {
                rj = j + grid_task->proc;
            }

            if (rj >= grid_task->m)
            {
                rj = abs(grid_task->m - rj);
            }

            if (j == grid_task->proc)
            {
                rj += 1;
            }

            grid_task->tasks_r[rj][ri].rank   = grid_task->tasks[i][j].rank;
            grid_task->tasks_r[rj][ri].number = grid_task->tasks[i][j].number;

            /*
             * Link real cell to redudancy cell
             */
            grid_task->tasks[i][j].redudancy  = (void *)&(grid_task->tasks_r[i][j]);
        }
    }

    return 0;
}

/*****************************************************************************/
/* Basic initialize function                                                 */
/*****************************************************************************/
int grid_task_init(grid_task_t *grid_task) // TODO: args
{
    if (grid_task_init_(grid_task) != 0)
        return -1;

    if (grid_task_init_redudancy_(grid_task) != 0)
        return -1;

    return 0;
}

// TDDO
int grid_task_reconfigure_redudancy_(grid_task_t *grid_task, const int x, const int y)
{
    const int row = grid_task->m;
    const int col = grid_task->n;

    // Check UP neighbor
    if (x + 1 < row)
    {
        task_t *up = &(grid_task->tasks_r[x + 1][y]);
        if (up && (!up->reconfigure))
        {
            grid_task->tasks_r[x][y].rank = up->rank;
            return 0;
        }
    }

    // Check RIGHT neighbor
    if (y + 1 < col)
    {
        task_t *right = &(grid_task->tasks_r[x][y + 1]);
        if (right && (!right->reconfigure))
        {
            grid_task->tasks_r[x][y].rank = right->rank;
            return 0;
        }
    }

    // Check DOWN neighbor
    if (x - 1 > 0)
    {
        task_t *down = &(grid_task->tasks_r[x - 1][y]);
        if (down && (!down->reconfigure))
        {
            grid_task->tasks_r[x][y].rank = down->rank;
            return 0;
        }
    }

    // Check LEFT neighbor
    if (y - 1 > 0)
    {
        task_t *left = &(grid_task->tasks_r[x][y - 1]);
        if (left && (!left->reconfigure))
        {
            grid_task->tasks_r[x][y].rank = left->rank;
            return 0;
        }
    }

    return -1;
}


int grid_task_reconfigure_(task_t *fail_task)
{
    task_t *redudancy_task = fail_task->redudancy;
    fail_task->rank = redudancy_task->rank;

    return 0;
}

int grid_task_reconfigure(grid_task_t *grid_task)
{
    printf("<%s> Invoke\n", __FUNCTION__);

    int i, j;

    if (!grid_task)
    {
        fprintf(stderr, "<%s> Can't kill fail rank\n", __FUNCTION__);
        return -1;
    }

    if (!grid_task->tasks)
    {
        fprintf(stderr, "<%s> Can't kill fail rank\n", __FUNCTION__);
        return -1;
    }

    for (i = 0; i < grid_task->n; i++)
    {
        for (j = 0; j < grid_task->m; j++)
        {
            if (grid_task->tasks[i][j].rank == DEAD_PROC)
            {
                /*
                 * Reconfigure real grid
                 */
                if (grid_task_reconfigure_(&grid_task->tasks[i][j]) != 0)
                {
                    return -1;
                }
            }
        }
    }

    printf("<%s> Done\n", __FUNCTION__);
    return 0;
}


int grid_task_reconfigure_redudancy(grid_task_t *grid_task)
{
    int i, j;

    if (!grid_task)
    {
        fprintf(stderr, "<%s> Can't kill fail rank\n", __FUNCTION__);
        return -1;
    }

    if (!grid_task->tasks)
    {
        fprintf(stderr, "<%s> Can't kill fail rank\n", __FUNCTION__);
        return -1;
    }

    for (i = 0; i < grid_task->n; i++)
    {
        for (j = 0; j < grid_task->m; j++)
        {
            task_t *redudancy_task = grid_task->tasks[i][j].redudancy;

            if(redudancy_task->reconfigure)
            {
                /*
                 * Reconfigure redudancy grid
                 */
                if (grid_task_reconfigure_redudancy_(grid_task, i, j) != 0)
                {
                    return -1;
                }

                redudancy_task->reconfigure = false;
            }
        }   
    }

    return 0;
}

int grid_task_kill_proc(grid_task_t *grid_task, int fail_rank)
{
    printf("<%s> Invoke\n", __FUNCTION__);

    int i, j;

    if (!grid_task)
    {
        fprintf(stderr, "<%s> Can't kill fail rank\n", __FUNCTION__);
        return -1;
    }

    if (!grid_task->tasks)
    {
        fprintf(stderr, "<%s> Can't kill fail rank\n", __FUNCTION__);
        return -1;
    }

    printf("Killing %d rank...\n", fail_rank);

    for (i = 0; i < grid_task->n; i++)
    {
        for (j = 0; j < grid_task->m; j++)
        {
            if (grid_task->tasks[i][j].rank == fail_rank)
            {
                grid_task->tasks[i][j].rank = DEAD_PROC;
            }

            if (grid_task->tasks_r[i][j].rank == fail_rank)
            {
                grid_task->tasks_r[i][j].reconfigure = true;
            }
        }   
    }

    printf("<%s> Done\n", __FUNCTION__);
    return 0;
}

int grid_task_show_reconfigure(grid_task_t *grid)
{
    printf("<%s> Invoke\n", __FUNCTION__);

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

    if (!grid->tasks_r)
    {
        fprintf(stderr, "<%s> Can't show grid_task (redudancy)\n", __FUNCTION__);
        return -1;
    }

    const task_t **real      = (const task_t **)grid->tasks;
    const task_t **redudancy = (const task_t **)grid->tasks_r;

    printf("<%s> %s grid task\n", __FUNCTION__, "Real");
    grid_task_show_reconfigure_(real, grid->n, grid->m, grid->proc);

    printf("<%s> %s grid task\n", __FUNCTION__, "Redudancy");
    grid_task_show_reconfigure_(redudancy, grid->n, grid->m, grid->proc);

    printf("<%s> Done\n", __FUNCTION__);
    return 0;
}

void test(const int row, const int col, const int proc)
{
    printf("[INFO] Test: start!\n");

    grid_task_t *grid_task = grid_task_allocate(row, col, proc);

    grid_task_init(grid_task);
    //grid_task_show(grid_task);
    grid_task_show_color(grid_task);

    grid_task_kill_proc(grid_task, 0);
    grid_task_kill_proc(grid_task, 1);
    grid_task_kill_proc(grid_task, 2);
    grid_task_kill_proc(grid_task, 3);
    grid_task_kill_proc(grid_task, 4);
    grid_task_kill_proc(grid_task, 5);
    grid_task_kill_proc(grid_task, 6);
    grid_task_kill_proc(grid_task, 7);
    grid_task_show(grid_task);
    grid_task_show_reconfigure(grid_task);

    grid_task_reconfigure(grid_task);
    grid_task_show(grid_task);
    grid_task_reconfigure_redudancy(grid_task);
    grid_task_show(grid_task);

    grid_task_show_reconfigure(grid_task);

    grid_task_show_color(grid_task);

    printf("\n*******************\n");
/* 
    grid_task_kill_proc(grid_task, 3);
    grid_task_show(grid_task);
    grid_task_show_reconfigure(grid_task);

    grid_task_reconfigure(grid_task);
    grid_task_show(grid_task);
    grid_task_reconfigure_redudancy(grid_task);
    grid_task_show(grid_task);

    grid_task_show_reconfigure(grid_task);

    grid_task_show_color(grid_task);
*/
/*
    grid_task_kill_proc(grid_task, 38);
    grid_task_show(grid_task);

    grid_task_reconfigure(grid_task);
    grid_task_show(grid_task);

    grid_task_reconfigure_redudancy(grid_task);
    grid_task_show(grid_task);
*/
    grid_task_free(grid_task);

    printf("[INFO] Test: done!\n");
}

int main(int argc, char const *argv[])
{
    int x = 16, y = 16, n = 8;

    test(x, y, n);

    return 0;
}

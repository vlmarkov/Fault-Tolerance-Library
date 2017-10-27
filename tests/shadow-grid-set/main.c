#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define RST  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define DEAD -1

typedef struct
{
    int rank;
    int number;
    size_t x;
    size_t y;
} task_t;

typedef struct
{
    int x; // Grid size by x
    int y; // Grid size by x
    int n; // Proc size by node
    task_t **tasks;
} grid_task_t;

grid_task_t *grid_task_allocate(int x, int y, int proc)
{
    int i;

    grid_task_t *grid = (grid_task_t *)malloc(sizeof(grid_task_t));
    if (!grid)
    {
        fprintf(stderr, "<%s> Fail to allocate memory for grid\n", __FUNCTION__);
        exit(-1);
    }

    grid->tasks = (task_t **)malloc(sizeof(task_t *) * x);
    if (!grid->tasks)
    {
        free(grid);
        fprintf(stderr, "<%s> Fail to allocate memory for tasks\n", __FUNCTION__);
        exit(-1);
    }

    grid->x = x;
    grid->y = y;
    grid->n = proc;

    for (i = 0; i < y; i++)
    {
        task_t *new_row = (task_t *)malloc(sizeof(task_t) * x);
        if (!new_row)
        {
            // TODO: memory leak !!
            free(grid);
            fprintf(stderr, "<%s> Fail to allocate memory for task\n", __FUNCTION__);
            exit(-1);
        }
        grid->tasks[i] = new_row;
    }

    return grid;
}

void grid_task_free(grid_task_t *grid)
{
    int i;

    if(!grid)
    {
        fprintf(stderr, "<%s> Can't free memory\n", __FUNCTION__);
        exit(-1);        
    }

    for (i = 0; i < grid->x; i++)
    {
        free(grid->tasks[i]);
    }

    free(grid->tasks);
    free(grid);
}

void grid_task_show(grid_task_t *grid)
{
    if (!grid)
    {
        fprintf(stderr, "<%s> Can't show grid_task\n", __FUNCTION__);
        exit(-1);
    }

    int i, j;

    for (i = 0; i < grid->y; i++)
    {
        for (j = 0; j < grid->x; j++)
        {
            if (grid->tasks[i][j].rank < 8)
            {
                printf(RED"%04d "RST, grid->tasks[i][j].rank);
            }
            else if ((grid->tasks[i][j].rank > 7) &&
                     (grid->tasks[i][j].rank < 16))
            {
                printf(GRN"%04d "RST, grid->tasks[i][j].rank);
            }
            else
            {
                printf("%04d ", grid->tasks[i][j].rank);
            }
        }
        printf("\n");
    }
    printf("\n");
}

void grid_task_init_real(grid_task_t *rgrid)
{
    if (!rgrid)
    {
        fprintf(stderr, "<%s> Can't init grid_task\n", __FUNCTION__);
        exit(-1);
    }

    if (!rgrid->tasks)
    {
        fprintf(stderr, "<%s> Can't init grid_tasks\n", __FUNCTION__);
        exit(-1);
    }

    int i, j, counter = 0;

    for (i = 0; i < rgrid->y; i++)
    {
        for (j = 0; j < rgrid->x; j++)
        {
            rgrid->tasks[i][j].rank   = counter;
            rgrid->tasks[i][j].number = counter++;
        }
    }
}

void grid_task_init_shadow(const grid_task_t *rgrid, grid_task_t *sgrid)
{
    if ((!sgrid) || (!rgrid))
    {
        fprintf(stderr, "<%s> Can't init grid_task\n", __FUNCTION__);
        exit(-1);
    }

    if ((!sgrid->tasks) || (!rgrid->tasks))
    {
        fprintf(stderr, "<%s> Can't init grid_tasks\n", __FUNCTION__);
        exit(-1);
    }

    int i, j;
    for (i = 0; i < sgrid->y; i++)
    {
        int si = i + sgrid->n;
        
        if (si >= sgrid->y)
        {
            si = abs(sgrid->y - si);
        }

        int sj = 0;
        for (j = 0; j < sgrid->x; j++)
        {
            if (j == 0)
            {
                sj = j + sgrid->n;
            }

            if (sj >= sgrid->x)
            {
                sj = abs(sgrid->x - sj);
            }

            if (j == sgrid->n)
            {
                sj += 1;
            }

            //printf("(%d:%d) ", sj, si);

            sgrid->tasks[sj][si].rank   = rgrid->tasks[i][j].rank;
            sgrid->tasks[sj][si].number = rgrid->tasks[i][j].number;

            sj += 2;
        }
        //printf("\n");
        //break;
    }
}


void test(int row, int col, int proc)
{
    printf("[INFO] Test: start!\n");

    grid_task_t *gt_real   = grid_task_allocate(row, col, proc);
    grid_task_t *gt_shadow = grid_task_allocate(row, col, proc);

    grid_task_init_real(gt_real);
    grid_task_show(gt_real);
    
    grid_task_init_shadow(gt_real, gt_shadow);
    grid_task_show(gt_shadow);

    grid_task_free(gt_real);
    grid_task_free(gt_shadow);

    printf("[INFO] Test: done!\n");
}

int main(int argc, char const *argv[])
{
    int x = 16, y = 16, n = 8;

    test(x, y, n);

    return 0;
}
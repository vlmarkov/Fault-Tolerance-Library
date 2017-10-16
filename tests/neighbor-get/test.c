#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define DEAD -1

typedef struct
{
    int rank;
    int number;
    size_t x;
    size_t y;
} task_t;

task_t **allocate_grid_task(int x, int y)
{
    int i;

    task_t **grid = (task_t **)malloc(sizeof(task_t *) * x);
    if (!grid)
    {
        fprintf(stderr, "%s\n", "error");
        exit(-1);
    }

    for (i = 0; i < x; i++)
    {
        task_t *new_row = (task_t *)malloc(sizeof(task_t) * y);
        if (!new_row)
        {
            fprintf(stderr, "%s\n", "error");
            exit(-1);
        }
        grid[i] = new_row;
    }

    return grid;
}

void init_grid_task(task_t **grid, int x, int y)
{
    int i, j, counter = 0;

    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            grid[i][j].rank   = counter;
            grid[i][j].number = counter++;
        }
    }
}

void show_grid_task(task_t **grid, int x, int y)
{
    int i, j;

    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            printf("%d\t", grid[i][j].rank);
        }
        printf("\n");
    }
    printf("\n");
}

void free_grid_task(task_t **grid, int x, int y)
{
    int i;

    for (i = 0; i < x; i++)
        free(grid[i]);

    free(grid);
}

bool try_to_reassing_rank(task_t **grid,
                          const int row, const int col,
                          const int x, const int y)
{
    // Check UP neighbor
    if (x + 1 < row)
    {
        task_t *up = &(grid[x + 1][y]);
        if (up && (up->rank != DEAD))
        {
            grid[x][y].rank = up->rank;
            return true;
        }
    }

    // Check RIGHT neighbor
    if (y + 1 < col)
    {
        task_t *right = &(grid[x][y + 1]);
        if (right && (right->rank != DEAD))
        {
            grid[x][y].rank = right->rank;
            return true;
        }
    }

    // Check DOWN neighbor
    if (x - 1 > 0)
    {
        task_t *down = &(grid[x - 1][y]);
        if (down && (down->rank != DEAD))
        {
            grid[x][y].rank = down->rank;
            return true;
        }
    }

    // Check LEFT neighbor
    if (y - 1 > 0)
    {
        task_t *left = &(grid[x][y - 1]);
        if (left && (left->rank != DEAD))
        {
            grid[x][y].rank = left->rank;
            return true;
        }
    }

    return false;
}

void repair_grid_task(task_t **grid, int x, int y)
{
    int i, j;

    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            if (grid[i][j].rank == DEAD)
            {
                if (try_to_reassing_rank(grid, x, y, i, j))
                    fprintf(stderr, "[DEBUG] repair - sucsees\n");
                else
                    fprintf(stderr, "[DEBUG] repair - fail\n");
            }
        }
    }
}

void test(int row, int col, int dead_x, int dead_y)
{
    printf("[INFO] Test: start!\n");

    task_t **grid_task = allocate_grid_task(row, col);

    init_grid_task(grid_task, row, col);

    show_grid_task(grid_task, row, col);

    grid_task[dead_x][dead_y].rank = DEAD;

    show_grid_task(grid_task, row, col);

    repair_grid_task(grid_task, row, col);

    show_grid_task(grid_task, row, col);

    free_grid_task(grid_task, row, col);

    printf("[INFO] Test: done!\n\n");
}

int main(int argc, char const *argv[])
{
    int i, j;
    int x = 4, y = 4;

    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            test(x, y, i, j);
        }
    }

    return 0;
}
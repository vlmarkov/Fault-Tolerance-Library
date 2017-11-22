#include "grid-task.h"

#include <stdio.h>

void test(const int row, const int col, const int nx, const int ny, const int commsize)
{
    printf("[INFO] Test: start!\n");

    grid_task_t *grid_task = grid_task_allocate(row, col, nx, ny, commsize);

    grid_task_init(grid_task);
    grid_task_show(grid_task);
    grid_task_redundancy_show(grid_task);


    for (int i = 0; i < commsize; i++)
    {
    	grid_task_redundancy_ranks_send_show(grid_task, i);
    	grid_task_redundancy_ranks_receive_show(grid_task, i);
    }

    grid_task_free(grid_task);

    printf("[INFO] Test: done!\n");
}

int main(int argc, char const *argv[])
{
	int x = 1024, y = 1024, n = 16;
    int nx = 256, ny = 256;

    test(x, y, nx, ny, n);
	return 0;
}
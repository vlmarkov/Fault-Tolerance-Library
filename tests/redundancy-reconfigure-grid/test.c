#include "utils.h"
#include "grid-task.h"

#include <stdio.h>

void test(const int row, const int col, const int nx, const int ny, const int commsize)
{
    printf("<%s> Inovke!\n", __FUNCTION__);

    printf("<%s> row      : %d\n", __FUNCTION__, row);
    printf("<%s> col      : %d\n", __FUNCTION__, col);
    printf("<%s> ny       : %d\n", __FUNCTION__, nx);
    printf("<%s> nx       : %d\n", __FUNCTION__, ny);
    printf("<%s> commsize : %d\n", __FUNCTION__, commsize);

    grid_task_t *grid_task = grid_task_allocate(row, col, nx, ny, commsize);

    grid_task_init(grid_task);

    grid_task_show(grid_task);

    task_t *my_task_0 = grid_task_get(grid_task, 0);
    task_t *my_task_1 = grid_task_get(grid_task, 1);

    grid_task_kill_rank(grid_task, 2);
    grid_task_kill_rank(grid_task, 3);
    grid_task_kill_rank(grid_task, 6);
    grid_task_kill_rank(grid_task, 7);
    grid_task_kill_rank(grid_task, 8);
    grid_task_kill_rank(grid_task, 9);
    grid_task_kill_rank(grid_task, 10);
    grid_task_kill_rank(grid_task, 11);
    grid_task_kill_rank(grid_task, 12);
    grid_task_kill_rank(grid_task, 13);
    grid_task_kill_rank(grid_task, 14);
    grid_task_kill_rank(grid_task, 15);

    grid_task_show(grid_task);

    grid_task_repair(grid_task);

    grid_task_task_show(grid_task);

    grid_task_show(grid_task);



    grid_task_redundancy_ranks_send_show(my_task_0);
    grid_task_redundancy_ranks_receive_show(my_task_0);
    grid_task_redundancy_ranks_send_show(my_task_1);
    grid_task_redundancy_ranks_receive_show(my_task_1);

    
    //grid_task_redundancy_show(grid_task);

    printf("<%s> Done!\n", __FUNCTION__);
}

int main(int argc, char const *argv[])
{
    const int x  = 1024;
    const int y  = 1024;
    const int n  = 16;
    const int nx = 256;
    const int ny = 256;

    test(x, y, nx, ny, n);

    return 0;
}

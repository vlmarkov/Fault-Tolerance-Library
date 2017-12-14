/*****************************************************************************/
/*Grid-task functionality unit test                                          */
/*****************************************************************************/

#include "utils.h"
#include "grid-task.h"

int main(int argc, char const *argv[])
{
    GridTask gt(1024, 1024, 256, 256, 16);

    gt.init(GRID_TASK_COMPUTE_REDUNDANCY);
/*
    gt.kill(15);

    gt.repair();

    gt.kill(5);

    gt.repair();

    gt.kill(7);

    gt.repair();

    gt.kill(13);

    gt.repair();
*/
    gt.show();

    return 0;
}

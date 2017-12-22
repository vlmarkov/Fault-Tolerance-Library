/*****************************************************************************/
/* Grid-task functionality unit test                                         */
/*****************************************************************************/

#include "src/utils.h"
#include "src/grid-task.h"

#include <iostream>

int main(int argc, char const *argv[])
{

    GridTask gt(1024, 1024, 256, 256, 16);

    gt.init(GRID_TASK_COMPUTE_REDUNDANCY);
/*
    gt.kill(15);
    gt.repair();
    gt.kill(14);
    gt.repair();
    gt.kill(13);
    gt.repair();
    gt.kill(12);
    gt.repair();
    gt.kill(11);
    gt.repair();
    gt.kill(10);
    gt.repair();
    gt.kill(9);
    gt.repair();
    gt.kill(8);
    gt.repair();
    gt.kill(7);
    gt.repair();
    gt.kill(6);
    gt.repair();
    gt.kill(3);
    gt.repair();
    gt.kill(2);
    gt.repair();
*/
    gt.show();

    gt.printGrid();

/* 
    std::vector<int> v;

    v.push_back(16);
    //v.push_back(64);
    //v.push_back(32);
    //v.push_back(48);


    for (int i = 0; i < (int)v.size(); i++)
    {
        GridTask test(1024, 1024, 256, 256, v[i]);
        test.init(GRID_TASK_COMPUTE_REDUNDANCY);
        test.show();
        test.printGridRedundancy();
        test.printGrid();
        std::cout << std::endl;
    }
*/
    return 0;
}

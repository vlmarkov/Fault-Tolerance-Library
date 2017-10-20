#include "fault-tollerance.h"

#include <mpi.h>

TaskCell::TaskCell(int _p, int _b, int _x, int _y) : proc(_p), block(_b), x(_x), y(_y) {}
TaskCell::~TaskCell() {}
        

TaskGrid::TaskGrid(int px, int py, int nx, int ny)
{
    int proc = 0;
    int block = 0;
    
    for (int i = 0; i < py; ++i)
    {
        std::vector<TaskCell> v;

        for (int j = 0; j < px; ++j)
        {
            v.push_back(TaskCell(proc++, block++, nx, ny));
        }

        taskGrid_.push_back(v);
    }
}

TaskGrid::~TaskGrid() {}

void TaskGrid::show()
{
    for (int i = 0; i < (int)taskGrid_.size(); ++i)
    {
        for (int j = 0; j < (int)taskGrid_[i].size(); ++j)
        {
            TaskCell c = taskGrid_[i][j];
            fprintf(stderr, "[p(%04i) b(%04i) x(%04i) y(%04i)] ",
                    c.proc, c.block, c.x, c.y);
        }
        fprintf(stderr, "\n");
    }
}

void TaskGrid::repair(int dead)
{
    printf("<%s> %d\n", __FUNCTION__, dead);

    for (int i = 0; i < (int)taskGrid_.size(); ++i)
    {
        for (int j = 0; j < (int)taskGrid_[i].size(); ++j)
        {
            if (taskGrid_[i][j].proc == dead)
            {
                taskGrid_[i][j].proc = -1;
            }
        }
    }
}

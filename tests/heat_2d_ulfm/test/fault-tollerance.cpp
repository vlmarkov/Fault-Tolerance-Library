#include "fault-tollerance.h"

#include <mpi.h>

#define DEAD -1

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
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < (int)taskGrid_.size(); ++i)
    {
        fprintf(stderr, "[DEBUG] Rank %d ", rank);
        for (int j = 0; j < (int)taskGrid_[i].size(); ++j)
        {
            TaskCell c = taskGrid_[i][j];

            //fprintf(stderr, "[p(%04i) b(%04i) x(%04i) y(%04i)] ", c.proc, c.block, c.x, c.y);

            fprintf(stderr, "[p(%04i)] ", c.proc);
        }
        fprintf(stderr, "\n");
    }
}

bool TaskGrid::tryToReassingRank(const int x, const int y)
{
    const int row = this->taskGrid_.size();
    const int col = this->taskGrid_[0].size();

    // Check UP neighbor
    if (x + 1 < row)
    {
        TaskCell *up = &(this->taskGrid_[x + 1][y]);
        if (up && (up->proc != DEAD))
        {
            this->taskGrid_[x][y].proc = up->proc;
            return true;
        }
    }

    // Check RIGHT neighbor
    if (y + 1 < col)
    {
        TaskCell *right = &(this->taskGrid_[x][y + 1]);
        if (right && (right->proc != DEAD))
        {
            this->taskGrid_[x][y].proc = right->proc;
            return true;
        }
    }

    // Check DOWN neighbor
    if (x - 1 > 0)
    {
        TaskCell *down = &(this->taskGrid_[x - 1][y]);
        if (down && (down->proc != DEAD))
        {
            this->taskGrid_[x][y].proc = down->proc;
            return true;
        }
    }

    // Check LEFT neighbor
    if (y - 1 > 0)
    {
        TaskCell *left = &(this->taskGrid_[x][y - 1]);
        if (left && (left->proc != DEAD))
        {
            this->taskGrid_[x][y].proc = left->proc;
            return true;
        }
    }

    return false;
}

void TaskGrid::repair()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < (int)taskGrid_.size(); ++i)
    {
        for (int j = 0; j < (int)taskGrid_[i].size(); ++j)
        {
            if (taskGrid_[i][j].proc == DEAD)
            {
                if (tryToReassingRank(i, j))
                    fprintf(stderr, "[DEBUG] Rank %d repair - sucsees\n", rank);
                else
                    fprintf(stderr, "[DEBUG] Rank %d repair - fail\n", rank);
            }
        }
    }
}

void TaskGrid::markDeadProc(int fail)
{
    for (int i = 0; i < (int)taskGrid_.size(); ++i)
    {
        for (int j = 0; j < (int)taskGrid_[i].size(); ++j)
        {
            if (taskGrid_[i][j].proc == fail)
            {
                taskGrid_[i][j].proc = DEAD;
            }
        }
    }
}

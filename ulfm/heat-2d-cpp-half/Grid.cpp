#include <iostream>

#include "Grid.h"

/*****************************************************************************/
/* This grid represent compute redundancy for half grid                      */
/*                                                                           */
/* +---+       +---+  - Each procces will migrate to i + gird / 2 (tasks)    */
/* | 0 | . . . | 1 |  - May occur undefined behavior if px and py not even   */
/* +---+       +---+  - need to test it                                      */
/*   .           .                                                           */
/*   .           .                                                           */
/*   .           .                                                           */
/* +---+       +---+                                                         */
/* | 2 | . . . | 3 |                                                         */
/* +---+       +---+                                                         */
/*   Pic 1 - Original                                                        */
/*                                                                           */
/* +---+       +---+                                                         */
/* |0/2| . . . |1/3|                                                         */
/* +---+       +---+                                                         */
/*   .           .                                                           */
/*   .           .                                                           */
/*   .           .                                                           */
/* +---+       +---+                                                         */
/* |2/0| . . . |3/1|                                                         */
/* +---+       +---+                                                         */
/*   Pic 2 - Redundancy                                                      */
/*****************************************************************************/

/**
 * Main constructor
 * @input: columns, rows, 
 *         cols by x, rows by y, 
 *         processes by x, processes by y
 */
Grid::Grid(int cols, int rows, int nx, int ny, int px, int py) : 
    cols_(cols), rows_(rows), nx_(nx), ny_(ny), px_(px), py_(py), alive_(py * px)
{
    // Allocate memory for grid-tasks
    for (int i = 0; i < this->py_; i++)
    {
        this->tasks_.push_back(std::vector<Task>(this->px_));
    }

    int rank   = 0;
    int tags   = 1;
    int repair = 1;
    int layer  = 0;

    // Init each task
    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            Task task(i, j, this->nx_, this->ny_, repair);

            this->setNeighbors_(task, i, j);

            this->setMpiRank_(task, rank++);

            this->setTags_(task, tags, layer);

            this->tasks_[i][j] = task;

            this->linkRanksTasks_(&this->tasks_[i][j], i, j);
        }
    }

    layer++;

    // Aditional MPI-tags assignment
    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            this->setTags_(this->tasks_[i][j], tags, layer);
        }
    }
}

/**
 * Destructor
 */
Grid::~Grid()
{
    // TODO: memory free
}

/*****************************************************************************/
/* Public methods                                                            */
/*****************************************************************************/

/**
 * Get task by MPI-rank
 * @input: MPI-rank
 * @return: pointer to Task
 */
Task* Grid::getTask(int rank)
{
    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            if (rank == *this->tasks_[i][j].getMpiRankPtr())
            {
                return &this->tasks_[i][j];
            }
        }
    }

    throw std::string("Can't find task by MPI-rank");
    return NULL; // depricated, change to nullptr
}

/**
 * Mark task as 'dead'
 * @input: MPI-rank
 */
void Grid::kill(int rank)
{
    this->alive_--; // Reduce alive processes

    if (this->alive_ < ((this->px_ * this->py_) * 0.5))
    {
        throw std::string("Reached the limit of reducibility"); // TODO: ?
    }

    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            if (rank == *this->tasks_[i][j].getMpiRankPtr())
            {
                this->tasks_[i][j].setStatus(DEAD_TASK);
            }
        }
    }
}

/**
 * Repair grid-tasks
 */
void Grid::repair()
{
    // TODO: ?

    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            if (this->tasks_[i][j].getStatus() == DEAD_TASK)
            {
                this->tasks_[i][j].repair();
            }
        }
    }
}

/**
 * Show grid-tasks
 */
void Grid::print()
{
    for (int i = 0; i < this->py_; i++)
    {
        for (int j = 0; j < this->px_; j++)
        {
            this->tasks_[i][j].print();
        }
    }
}

/*****************************************************************************/
/* Private methods                                                           */
/*****************************************************************************/

/**
 * Set neighbors for task
 * @input: reference to task, coordinates 'i' and 'j'
 */
void Grid::setNeighbors_(Task& task, int i, int j)
{
    // Safe assing neighbors

    // Set up neighbor
    if (i - 1 >= 0)
    {
        task.setUpNeighbor(&this->tasks_[i - 1][j]);
    }
    else
    {
        task.setUpNeighbor(NULL);
    }

    // Set left neighbor
    if (j - 1 >= 0)
    {
        task.setLeftNeighbor(&this->tasks_[i][j - 1]);
    }
    else
    {
        task.setLeftNeighbor(NULL);
    }

    // Set down neighbor
    if (i + 1 < this->py_)
    {
        task.setDownNeighbor(&this->tasks_[i + 1][j]);
    }
    else
    {
        task.setDownNeighbor(NULL);
    }

    // Set right neighbor
    if (j + 1 < this->px_)
    {
        task.setRightNeighbor(&this->tasks_[i][j + 1]);
    }
    else
    {
        task.setRightNeighbor(NULL);
    }
}

/**
 * Set MPI-rank for task
 * @input: MPI-rank
 */
void Grid::setMpiRank_(Task& task, int rank)
{
    task.setMpiRank(rank);
}

/**
 * Set MPI-tags for task
 * @input: reference to task, 
 *         reference to tag,
 *         redundancy layer
 */
void Grid::setTags_(Task& task, int& tag, int layer)
{
    // Safe assing tags

    Task* left = task.getLeftNeighbor();
    if (left)
    {
        int leftTag = left->getRightTag(layer);
        if (leftTag == -1)
        {
            task.addLeftTag(tag++);
        }
        else
        {
            task.addLeftTag(leftTag);
        }
    }
    else
    {
        task.addLeftTag(0);
    }

    Task* up = task.getUpNeighbor();
    if (up)
    {
        int upTag = up->getDownTag(layer);
        if (upTag == -1)
        {
            task.addUpTag(tag++);
        }
        else
        {
            task.addUpTag(upTag);
        }
    }
    else
    {
        task.addUpTag(0);
    }

    Task* right = task.getRightNeighbor();
    if (right)
    {
        int rightTag = right->getLeftTag(layer);
        if (rightTag == -1)
        {
            task.addRightTag(tag++);
        }
        else
        {
            task.addRightTag(rightTag);
        }
    }
    else
    {
        task.addRightTag(0);
    }

    Task* down = task.getDownNeighbor();
    if (down)
    {
        int downTag = down->getUpTag(layer);
        if (downTag == -1)
        {
            task.addDownTag(tag++);
        }
        else
        {
            task.addDownTag(downTag);
        }
    }
    else
    {
        task.addDownTag(0);
    }
}

/**
 * Compute next cordiantes for task (where task will be migrate)
 * @input: coordinates 'i' and 'j'
 */
void Grid::computeNextCoordinates_(int& i, int& j)
{
    int halfGrid = (this->py_ * this->px_ / 2);
    int cnt      = 0;

    for (; i < this->py_;)
    {
        for (; j < this->px_; ++j)
        {
            
            if (cnt == halfGrid)
            {
                return;
            }
            cnt++;
        }

        j = 0;

        i++;

        if (i == this->py_)
        {
            i = 0;
        }
    } 
}

/**
 * Link MPI-ranks, migrating tasks with task
 * @input: pointer to task, coordinates 'i' and 'j'
 */
void Grid::linkRanksTasks_(Task* task, int i, int j)
{
    // Safe assing

    this->computeNextCoordinates_(i, j);

    int* self = task->getMpiRankPtr();
    if (!self)
    {
        throw std::string("Can't get self MPI-rank");
    }

    int* redundancy = this->tasks_[i][j].getMpiRankPtr();
    if (!redundancy)
    {
        throw std::string("Can't get redundancy MPI-rank");
    }

    task->addRrank(self);
    task->addRrank(redundancy);

    task->addRtask(task);
    task->addRtask(&this->tasks_[i][j]);
}

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
Grid::Grid(int cols, int rows, int nx, int ny, int px, int py)
{
    this->cols_ = cols;
    this->rows_ = rows;
    this->nx_   = nx;
    this->ny_   = ny;
    this->px_   = px;
    this->py_   = py;

    // Allocate memory for grid-tasks
    for (int i = 0; i < this->py_; i++)
    {
        this->tasks_.push_back(std::vector<Task>(this->px_));
    }

    int rank = 0;
    int tags = 0;

    // Init each task
    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            Task t(i, j, this->nx_, this->ny_);

            this->setNeighbors(t, i, j);

            this->setMpiRank(t, rank++);

            this->setTags(t, tags++);

            this->tasks_[i][j] = t;

            this->linkRanksTasks(&this->tasks_[i][j], i, j);
        }
    }
}

Grid::~Grid()
{
    // TODO: memory free
}

void Grid::setNeighbors(Task& task, int i, int j)
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

void Grid::setMpiRank(Task& task, int rank)
{
    task.setMpiRank(rank);
}

void Grid::linkRanksTasks(Task* task, int i, int j)
{
    // Safe assing ranks

    this->computeNextCoordinates_(i, j);

    int* self = task->getMpiRankPtr();
    if (!self)
    {
        throw std::string("Can't get self mpi rank");
    }

    int* redundancy = this->tasks_[i][j].getMpiRankPtr();
    if (!redundancy)
    {
        throw std::string("Can't get redundancy mpi rank");
    }

    task->addRrank(self);
    task->addRrank(redundancy);

    task->addRtask(task);
    task->addRtask(&this->tasks_[i][j]);
}

void Grid::setTags(Task& task, int tag)
{
    // Safe assing tags

    task.addTag(tag);
    task.addTag(tag + (this->py_ * this->px_));
}

void setRedundancy(Task& task)
{

}

void Grid::kill(int rank)
{
    for (int i = 0; i < this->py_; ++i)
    {
        for (int j = 0; j < this->px_; ++j)
        {
            if (rank == *this->tasks_[i][j].getMpiRankPtr())
            {
                this->tasks_[i][j].setMpiRank(-1);
                return;
            }
        }
    }
}

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
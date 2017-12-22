#include "grid-task.h"

#include <iostream>
#include <algorithm>

Redundancy::Redundancy() { }
Redundancy::~Redundancy() { }

task_t *Redundancy::getSelfTask()
{
    if (!this->real.empty())
    {
        return this->real[0];
    }
    return NULL;
}

void Redundancy::addReal(task_t *t)
{
    this->real.push_back(t);
}

void Redundancy::addReplace(task_t *t)
{
    this->replace.push_back(t);
}

void Redundancy::addRedundancy(task_t *t)
{
    this->redundancy.push_back(t);
}

void Redundancy::realReorder()
{
    std::vector<task_t *> newReal;

    newReal.push_back(this->real[0]);

    for (int i = (int)this->real.size() - 1; i > 0; --i)
    {
        newReal.push_back(this->real[i]);
    }

    this->real = newReal;
}

void Redundancy::redundancyReorder()
{
    std::vector<task_t *> newRedundancy;

    newRedundancy.push_back(this->redundancy[0]);

    for (int i = (int)this->redundancy.size() - 1; i > 0; --i)
    {
        newRedundancy.push_back(this->redundancy[i]);
    }

    this->redundancy = newRedundancy;
}

void Redundancy::realSwapLast()
{
    int last = this->real.size() - 1;
    int prev = last - 1;

    std::swap(this->real[last], this->real[prev]);
}

void Redundancy::redundancySwapLast()
{
    int last = this->redundancy.size() - 1;
    int prev = last - 1;

    std::swap(this->redundancy[last], this->redundancy[prev]);
}

bool rsfunctor (task_t *a, task_t *b)
{
    return (a->rank < b->rank);
}

void Redundancy::redundancySort()
{
    std::sort(this->redundancy.begin(), this->redundancy.end(), rsfunctor);
}

std::vector<task_t *> Redundancy::getReal()
{
    return this->real;
}

std::vector<task_t *> Redundancy::getReplace()
{
    return this->replace;
}

std::vector<task_t *> Redundancy::getRedundancy()
{
    return this->redundancy;
}

int Redundancy::getRealSize()
{
    return this->real.size();
}

int Redundancy::getReplaceSize()
{
    return this->replace.size();
}

int Redundancy::getRedundancySize()
{
    return this->redundancy.size();
}

void Redundancy::printRealRank()
{
    for (int i = 0; i < (int)this->real.size(); ++i)
    {
        std::cout << this->real[i]->rank << " ";
    }
}

void Redundancy::printRealRankDetail()
{
    for (int i = 0; i < (int)this->real.size(); ++i)
    {
        std::cout << this->real[i]->rank << "("
                  << this->real[i]->x << ":"
                  << this->real[i]->y << ")"
                  << " ";
    }  
}

void Redundancy::printRedundancyRank()
{
    for (int i = 0; i < (int)this->redundancy.size(); ++i)
    {
        std::cout << this->redundancy[i]->rank << " ";
    }
}

void Redundancy::printReplaceRankDetail()
{
    for (int i = 0; i < (int)this->replace.size(); ++i)
    {
        std::cout << this->replace[i]->rank << "("
                  << this->replace[i]->x << ":"
                  << this->replace[i]->y << ")"
                  << " ";
    }
}

int Redundancy::repair(int row, int col)
{
/*
    printf("<%s> this->redundancy.size(%d)\n",
        __FUNCTION__, (int)this->redundancy.size());
*/
    if ((int)this->redundancy.size() == 0)
    {
        return GRID_TASK_BORDER; // null processor
    }

    //task_t *repairTask = this->redundancy[0];

    for (int i = 0; i < (int)this->redundancy.size(); i++)
    {
        task_t *findDead = this->redundancy[i];

        if ((findDead->x == row) && (findDead->y == col))
        {
/*
            printf("<%s> idex to erase(%d)\n",
                __FUNCTION__, i);
*/
            this->redundancy.erase(this->redundancy.begin() + i);
        }
    }
/*
    printf("<%s> this->redundancy.size(%d)\n",
        __FUNCTION__, (int)this->redundancy.size());
*/
    task_t *t = *this->redundancy.begin();

    return t->rank;
}

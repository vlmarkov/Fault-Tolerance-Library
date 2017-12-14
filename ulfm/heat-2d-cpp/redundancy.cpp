#include "grid-task.h"

#include <iostream>

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

void Redundancy::addRedundancy(task_t *t)
{
    this->redundancy.push_back(t);
}

std::vector<task_t *> Redundancy::getReal()
{
    return this->real;
}

std::vector<task_t *> Redundancy::getRedundancy()
{
    return this->redundancy;
}

int Redundancy::getRealSize()
{
    return this->real.size();
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

int Redundancy::repair()
{
    for (int i = 0; i < (int)this->redundancy.size(); ++i)
    {
        task_t *repairTask = this->redundancy[i];

        if (repairTask->rank != GRID_TASK_DEAD_PROC)
        {
            this->redundancy.erase(this->redundancy.begin() + i);
            return repairTask->rank;
        }
    }

    return GRID_TASK_BORDER; // null processor
}

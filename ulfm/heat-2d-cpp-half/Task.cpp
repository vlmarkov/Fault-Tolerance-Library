#include <iostream>
#include <cstring>
#include <string>

#include "Task.h"

// Default constructor
Task::Task() { }

// Main constructor
Task::Task(int i, int j, int nx, int ny)
{
    this->i_             = i;
    this->j_             = j;
    this->nx_            = nx;
    this->ny_            = ny;
    this->upNeighbor_    = NULL;
    this->downNeighbor_  = NULL;
    this->leftNeighbor_  = NULL;
    this->rightNeighbor_ = NULL;
    this->mpiRank_       = -1;

    this->grid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
    if (!this->grid_)
    {
        throw std::string("Can't allocate memory for grid");
    }

    this->newGrid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
    if (!this->newGrid_)
    {
        throw std::string("Can't allocate memory for new grid");
    }

    // TODO
}

// Copy constructor
Task::Task(const Task& rhs)
{
    this->i_             = rhs.i_;
    this->j_             = rhs.j_;
    this->nx_            = rhs.nx_;
    this->ny_            = rhs.ny_;
    this->upNeighbor_    = rhs.upNeighbor_;
    this->downNeighbor_  = rhs.downNeighbor_;
    this->leftNeighbor_  = rhs.leftNeighbor_;
    this->rightNeighbor_ = rhs.rightNeighbor_;
    this->tags_          = rhs.tags_;
    this->mpiRank_       = rhs.mpiRank_;
    this->rRanks_        = rhs.rRanks_;
    this->rTasks_        = rhs.rTasks_;

    this->grid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
    if (!this->grid_)
    {
        throw std::string("Can't allocate memory for grid");
    }

    this->newGrid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
    if (!this->newGrid_)
    {
        throw std::string("Can't allocate memory for new grid");
    }

    std::memcpy(this->grid_, rhs.grid_,
        sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));

    std::memcpy(this->newGrid_, rhs.newGrid_,
        sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));

    // TODO
}

// Destructor
Task::~Task()
{
    delete[] this->grid_;
    delete[] this->newGrid_;
}

// Assign operator
Task& Task::operator=(const Task& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    this->i_             = rhs.i_;
    this->j_             = rhs.j_;
    this->nx_            = rhs.nx_;
    this->ny_            = rhs.ny_;
    this->upNeighbor_    = rhs.upNeighbor_;
    this->downNeighbor_  = rhs.downNeighbor_;
    this->leftNeighbor_  = rhs.leftNeighbor_;
    this->rightNeighbor_ = rhs.rightNeighbor_;
    this->tags_          = rhs.tags_;
    this->mpiRank_       = rhs.mpiRank_;
    this->rRanks_        = rhs.rRanks_;
    this->rTasks_        = rhs.rTasks_;

    this->grid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
    if (!this->grid_)
    {
        throw std::string("Can't allocate memory for grid");
    }

    this->newGrid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
    if (!this->newGrid_)
    {
        throw std::string("Can't allocate memory for new grid");
    }

    std::memcpy(this->grid_, rhs.grid_,
        sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));

    std::memcpy(this->newGrid_, rhs.newGrid_,
        sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));

    // TODO

    return *this;
}

void Task::setUpNeighbor(Task* up)
{
    if (up)
    {
        this->upNeighbor_ = up;
    }
    else
    {
        this->upNeighbor_ = NULL;
    }
}

void Task::setDownNeighbor(Task* down)
{
    if (down)
    {
        this->downNeighbor_ = down;
    }
    else
    {
        this->downNeighbor_ = NULL;
    }
}

void Task::setLeftNeighbor(Task* left)
{
    if (left)
    {
        this->leftNeighbor_ = left;
    }
    else
    {
        this->leftNeighbor_ = NULL;
    }
}

void Task::setRightNeighbor(Task* right)
{
    if (right)
    {
        this->rightNeighbor_ = right;
    }
    else
    {
        this->rightNeighbor_ = NULL;
    }
}

void Task::setMpiRank(int rank)
{
    this->mpiRank_ = rank;
}

int* Task::getMpiRankPtr()
{
    return &this->mpiRank_;
}

void Task::addRrank(int* rank)
{
    this->rRanks_.push_back(rank);
}

void Task::addRtask(Task* task)
{
    this->rTasks_.push_back(task);
}

void Task::addTag(int tag)
{
    this->tags_.push_back(tag);
}

void Task::print()
{
    std::cout << "Task [ " << this->i_ << ", "
              << this->j_ << " ]" << std::endl;

    std::cout << "mpi rank       : " << this->mpiRank_ << std::endl;

    if (this->upNeighbor_)
    {
        std::cout << "Up neighbor    : [ "
                  << this->upNeighbor_->i_ << ", "
                  << this->upNeighbor_->j_ << "]"
                  << std::endl;
    }
    else
    {
        std::cout << "Up neighbor    : [ -1, -1 ]" << std::endl;
    }

    if (this->leftNeighbor_)
    {
        std::cout << "Left neighbor  : [ "
                  << this->leftNeighbor_->i_ << ", "
                  << this->leftNeighbor_->j_ << "]"
                  << std::endl;
    }
    else
    {
        std::cout << "Left neighbor  : [ -1, -1 ]" << std::endl;
    }

    if (this->downNeighbor_)
    {
        std::cout << "Down neighbor  : [ "
                  << this->downNeighbor_->i_ << ", "
                  << this->downNeighbor_->j_ << "]"
                  << std::endl;
    }
    else
    {
        std::cout << "Down neighbor  : [ -1, -1 ]" << std::endl;
    }

    if (this->rightNeighbor_)
    {
        std::cout << "Right neighbor : [ "
                  << this->rightNeighbor_->i_ << ", "
                  << this->rightNeighbor_->j_ << "]"
                  << std::endl;
    }
    else
    {
        std::cout << "Right neighbor : [ -1, -1 ]" << std::endl;
    }

    std::cout << "Ranks          : [ ";
    for (int i = 0; i < (int)this->rRanks_.size(); ++i)
    {
        std::cout << *this->rRanks_[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Tags           : [ ";
    for (int i = 0; i < (int)this->tags_.size(); ++i)
    {
        std::cout << this->tags_[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Tasks          : [ ";
    for (int i = 0; i < (int)this->rTasks_.size(); ++i)
    {
        std::cout << this->rTasks_[i]->i_ << ","
                  << this->rTasks_[i]->j_ << " ";
    }
    std::cout << "]" << std::endl;

    // TODO

    std::cout << std::endl;
}

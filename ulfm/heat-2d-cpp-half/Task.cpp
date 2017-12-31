#include <iostream>
#include <cstring>
#include <string>

#include "Task.h"

/**
 * Default constructor
 */
Task::Task() : 
    i_(-1), j_(-1), nx_(-1), ny_(-1),
    upNeighbor_(NULL), downNeighbor_(NULL),
    leftNeighbor_(NULL), rightNeighbor_(NULL),
    grid_(NULL), newGrid_(NULL),
    mpiRank_(-1), status_(UNKNOWN_TASK), repair_(-1)
{
    ;
}

/**
 * Main constructor
 */
Task::Task(int i, int j, int nx, int ny, int repair) : 
    i_(i), j_(j), nx_(nx), ny_(ny),
    upNeighbor_(NULL), downNeighbor_(NULL),
    leftNeighbor_(NULL), rightNeighbor_(NULL),
    grid_(NULL), newGrid_(NULL),
    mpiRank_(-1), status_(UNKNOWN_TASK), repair_(repair)
{
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

/**
 * Copy constructor
 */
Task::Task(const Task& rhs) :
    i_(rhs.i_), j_(rhs.j_), nx_(rhs.nx_), ny_(rhs.ny_),
    upNeighbor_(rhs.upNeighbor_), downNeighbor_(rhs.downNeighbor_),
    leftNeighbor_(rhs.leftNeighbor_), rightNeighbor_(rhs.rightNeighbor_),
    grid_(NULL), newGrid_(NULL),
    mpiRank_(rhs.mpiRank_), status_(rhs.status_), repair_(rhs.repair_)
{
    if (rhs.upNeighborTags_.size() > 0)
    {
        this->upNeighborTags_ = rhs.upNeighborTags_;
    }

    if (rhs.downNeighborTags_.size() > 0)
    {
        this->downNeighborTags_ = rhs.downNeighborTags_;
    }

    if (rhs.leftNeighborTags_.size() > 0)
    {
        this->leftNeighborTags_ = rhs.leftNeighborTags_;
    }

    if (rhs.rightNeighborTags_.size() > 0)
    {
        this->rightNeighborTags_ = rhs.rightNeighborTags_;
    }

    if (rhs.rRanks_.size() > 0)
    {
        this->rRanks_ = rhs.rRanks_;
    }

    if (rhs.rTasks_.size() > 0)
    {
        this->rTasks_  = rhs.rTasks_;
    }

    if (rhs.grid_)
    {
        this->grid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
        if (!this->grid_)
        {
            throw std::string("Can't allocate memory for grid");
        }
        std::memcpy(this->grid_, rhs.grid_,
            sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));
    }

    if (rhs.newGrid_)
    {
        this->newGrid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
        if (!this->newGrid_)
        {   
            throw std::string("Can't allocate memory for new grid");
        }

        std::memcpy(this->newGrid_, rhs.newGrid_,
            sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));
    }

    // TODO
}

/**
 * Destructor
 */
Task::~Task()
{
    delete[] this->grid_;
    delete[] this->newGrid_;
}

/**
 * Assign operator
 */
Task& Task::operator=(const Task& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    const_cast <int&> (this->i_)  = rhs.i_;
    const_cast <int&> (this->j_)  = rhs.j_;
    const_cast <int&> (this->nx_) = rhs.nx_;
    const_cast <int&> (this->ny_) = rhs.ny_;

    this->mpiRank_ = rhs.mpiRank_;
    this->status_  = rhs.status_;
    this->repair_  = rhs.repair_;

    this->upNeighbor_ = rhs.upNeighbor_;
    this->downNeighbor_ = rhs.downNeighbor_;
    this->leftNeighbor_ = rhs.leftNeighbor_;
    this->rightNeighbor_ = rhs.rightNeighbor_;

    if (rhs.upNeighborTags_.size() > 0)
    {
        this->upNeighborTags_ = rhs.upNeighborTags_;
    }

    if (rhs.downNeighborTags_.size() > 0)
    {
        this->downNeighborTags_ = rhs.downNeighborTags_;
    }

    if (rhs.leftNeighborTags_.size() > 0)
    {
        this->leftNeighborTags_ = rhs.leftNeighborTags_;
    }

    if (rhs.rightNeighborTags_.size() > 0)
    {
        this->rightNeighborTags_ = rhs.rightNeighborTags_;
    }

    if (rhs.rRanks_.size() > 0)
    {
        this->rRanks_ = rhs.rRanks_;
    }

    if (rhs.rTasks_.size() > 0)
    {
        this->rTasks_  = rhs.rTasks_;
    }

    if (rhs.grid_)
    {
        this->grid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
        if (!this->grid_)
        {
            throw std::string("Can't allocate memory for grid");
        }
        std::memcpy(this->grid_, rhs.grid_,
            sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));
    }

    if (rhs.newGrid_)
    {
        this->newGrid_ = new double [((this->ny_ + 2) * (this->nx_ + 2))];
        if (!this->newGrid_)
        {   
            throw std::string("Can't allocate memory for new grid");
        }

        std::memcpy(this->newGrid_, rhs.newGrid_,
            sizeof(double) * ((this->ny_ + 2) * (this->nx_ + 2)));
    }

    // TODO

    return *this;
}

/*****************************************************************************/
/* Public methods                                                            */
/*****************************************************************************/

/* */
void Task::setMpiRank(int rank)
{
    this->mpiRank_ = rank;
    this->status_  = ALIVE_TASK;
}

int* Task::getMpiRankPtr()
{
    return &this->mpiRank_;
}

int Task::getMpiRank()
{
    return this->mpiRank_;
}

/* */
int Task::getStatus()
{
    return this->status_;
}

void Task::setStatus(int status)
{
    if (status == DEAD_TASK)
    {
#ifdef MPI_SUPPORT
        this->mpiRank_ = MPI_PROC_NULL;
#else
        this->mpiRank_ = -1;
#endif /* MPI_SUPPORT */
    }

    this->status_ = status;
}

/* */
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

Task* Task::getUpNeighbor()
{
    return this->upNeighbor_;
}

Task* Task::getDownNeighbor()
{
    return this->downNeighbor_;
}

Task* Task::getLeftNeighbor()
{
    return this->leftNeighbor_;
}

Task* Task::getRightNeighbor()
{
    return this->rightNeighbor_;
}

int Task::getUpNeighborRank(int layer)
{
    Task* up = this->getUpNeighbor();
    if (up)
    {
        return up->getNextRank_(layer);
    }
    else
    {
#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */
    }
}

int Task::getDownNeighborRank(int layer)
{
    Task* down = this->getDownNeighbor();
    if (down)
    {
        return down->getNextRank_(layer);
    }
    else
    {
#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */
    }
}

int Task::getLeftNeighborRank(int layer)
{
    Task* left = this->getLeftNeighbor();
    if (left)
    {
        return left->getNextRank_(layer);
    }
    else
    {
#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */
    }
}

int Task::getRightNeighborRank(int layer)
{
    Task* right = this->getRightNeighbor();
    if (right)
    {
        return right->getNextRank_(layer);
    }
    else
    {
#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */
    }
}

/* */
void Task::setLocalGrid(double* grid)
{
    if (!grid)
    {
        throw std::string("Can't set local grid (bad pointer)");
    }

    this->grid_ = grid;
}

void Task::setLocalNewGrid(double* newGrid)
{
    if (!newGrid)
    {
        throw std::string("Can't set local new grid (bad pointer)");
    }

    this->newGrid_ = newGrid;
}

double* Task::getLocalGrid()
{
    return this->grid_;
}

double* Task::getLocalNewGrid()
{
    return this->newGrid_;
}

/* */
void Task::addRrank(int* rank)
{
    this->rRanks_.push_back(rank);
}

/* */
void Task::addRtask(Task* task)
{
    this->rTasks_.push_back(task);
}

/* */
void Task::addUpTag(int tag)
{
    this->upNeighborTags_.push_back(tag);
}

void Task::addDownTag(int tag)
{
    this->downNeighborTags_.push_back(tag);
}

void Task::addLeftTag(int tag)
{
    this->leftNeighborTags_.push_back(tag);
}

void Task::addRightTag(int tag)
{
    this->rightNeighborTags_.push_back(tag);
}

int Task::getUpTag(int layer)
{
    return this->getNextUpTag_(layer);
}

int Task::getDownTag(int layer)
{
    return this->getNextDownTag_(layer);
}

int Task::getLeftTag(int layer)
{
    return this->getNextLeftTag_(layer);
}

int Task::getRightTag(int layer)
{
    return this->getNextRightTag_(layer);
}

/**
 * Swap grid and new-grdi fields
 */
void Task::swapLocalGrids()
{
    std::swap(this->grid_, this->newGrid_);
}

/**
 * Repair task
 */
void Task::repair()
{
    if  (!this->repair_)
    {
        throw std::string("Can't repair task (reached repair limit)");
    }

    Task* task = this->getNextRepair_();

    // This automaticly set rank and status
    this->setMpiRank(task->mpiRank_);

    task->reduceRepairAbility_();

    this->reduceRepairAbility_();

    this->status_ = ALIVE_TASK;
}

/**
 * Show whole infomation about task
 */
void Task::print()
{
    std::cout << "Task [ " << this->i_ << ", "
              << this->j_ << " ]" << std::endl;

    std::cout << "MPI rank       : " << this->mpiRank_ << std::endl;

    switch(this->status_)
    {
        case UNKNOWN_TASK:
            std::cout << "Status         : UNKNOWN_TASK" << std::endl;
            break;
        case DEAD_TASK:
            std::cout << "Status         : DEAD_TASK" << std::endl;
            break;
        case ALIVE_TASK:
            std::cout << "Status         : ALIVE_TASK" << std::endl;
            break;
    }

    std::cout << "Repair         : " << this->repair_ << std::endl;

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

    std::cout << "Up tags        : [ ";
    for (int i = 0; i < (int)this->upNeighborTags_.size(); ++i)
    {
        std::cout << this->upNeighborTags_[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Down tags      : [ ";
    for (int i = 0; i < (int)this->downNeighborTags_.size(); ++i)
    {
        std::cout << this->downNeighborTags_[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Left tags      : [ ";
    for (int i = 0; i < (int)this->leftNeighborTags_.size(); ++i)
    {
        std::cout << this->leftNeighborTags_[i] << " ";
    }
    std::cout << "]" << std::endl;

    std::cout << "Right tags     : [ ";
    for (int i = 0; i < (int)this->rightNeighborTags_.size(); ++i)
    {
        std::cout << this->rightNeighborTags_[i] << " ";
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

/*****************************************************************************/
/* Private methods                                                           */
/*****************************************************************************/

/**
 *
 */
Task* Task::getNextRepair_()
{
    return this->rTasks_[1];
}

/**
 *
 */
int Task::getNextRank_(int layer)
{
    return *this->rRanks_[layer];
}

/**
 *
 */
int Task::getNextUpTag_(int layer)
{
    if (layer < (int)this->upNeighborTags_.size())
    {
        return this->upNeighborTags_[layer];
    }
    else
    {
        return -1;
    }
}

/**
 *
 */
int Task::getNextDownTag_(int layer)
{
    if (layer < (int)this->downNeighborTags_.size())
    {
        return this->downNeighborTags_[layer];
    }
    else
    {
        return -1;
    }
}

/**
 *
 */
int Task::getNextLeftTag_(int layer)
{
    if (layer < (int)this->leftNeighborTags_.size())
    {
        return this->leftNeighborTags_[layer];
    }
    else
    {
        return -1;
    }
}

/**
 *
 */
int Task::getNextRightTag_(int layer)
{
    if (layer < (int)this->rightNeighborTags_.size())
    {
        return this->rightNeighborTags_[layer];
    }
    else
    {
        return -1;
    }
}

/**
 *
 */
void Task::reduceRepairAbility_()
{
    this->repair_--;
}

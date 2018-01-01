#include <iostream>
#include <cstring>
#include <string>

#include "Task.h"

/**
 * Default constructor
 */
Task::Task() : i_(-1), j_(-1), nx_(-1), ny_(-1)
{
    this->upNeighbor_    = NULL;
    this->downNeighbor_  = NULL;
    this->leftNeighbor_  = NULL;
    this->rightNeighbor_ = NULL;
    this->grid_          = NULL;
    this->newGrid_       = NULL;
    this->mpiRank_       = -1;
    this->status_        = UNKNOWN_TASK;
    this->repair_        = -1;
}

/**
 * Main constructor
 * @input: coodinate 'i', coordinate 'j',
 *         size by x, size by y,
 *         repair counter
 */
Task::Task(int i, int j, int nx, int ny, int repair) : i_(i), j_(j), nx_(nx), ny_(ny)
{
    this->upNeighbor_    = NULL;
    this->downNeighbor_  = NULL;
    this->leftNeighbor_  = NULL;
    this->rightNeighbor_ = NULL;
    this->grid_          = NULL;
    this->newGrid_       = NULL;
    this->mpiRank_       = -1;
    this->status_        = UNKNOWN_TASK;
    this->repair_        = repair;

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
 * @input: task object
 */
Task::Task(const Task& rhs) : i_(rhs.i_), j_(rhs.j_), nx_(rhs.nx_), ny_(rhs.ny_)
{
    this->upNeighbor_    = rhs.upNeighbor_;
    this->downNeighbor_  = rhs.downNeighbor_;
    this->leftNeighbor_  = rhs.leftNeighbor_;
    this->rightNeighbor_ = rhs.rightNeighbor_;
    this->grid_          = NULL;
    this->newGrid_       = NULL;
    this->mpiRank_       = rhs.mpiRank_;
    this->status_        = rhs.status_;
    this->repair_        = rhs.repair_;

    this->upNeighborTags_    = rhs.upNeighborTags_;
    this->downNeighborTags_  = rhs.downNeighborTags_;
    this->leftNeighborTags_  = rhs.leftNeighborTags_;
    this->rightNeighborTags_ = rhs.rightNeighborTags_;
    this->rRanks_            = rhs.rRanks_;
    this->rTasks_            = rhs.rTasks_;
    this->replacements_      = rhs.replacements_;

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

    this->upNeighbor_    = rhs.upNeighbor_;
    this->downNeighbor_  = rhs.downNeighbor_;
    this->leftNeighbor_  = rhs.leftNeighbor_;
    this->rightNeighbor_ = rhs.rightNeighbor_;

    this->upNeighborTags_    = rhs.upNeighborTags_;
    this->downNeighborTags_  = rhs.downNeighborTags_;
    this->leftNeighborTags_  = rhs.leftNeighborTags_;
    this->rightNeighborTags_ = rhs.rightNeighborTags_;
    this->rRanks_            = rhs.rRanks_;
    this->rTasks_            = rhs.rTasks_;
    this->replacements_      = rhs.replacements_;

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

/**
 *
 */
void Task::setMpiRank(int rank)
{
    this->mpiRank_ = rank;
    this->status_  = ALIVE_TASK;
}

/**
 *
 */
int* Task::getMpiRankPtr()
{
    return &this->mpiRank_;
}

/**
 *
 */
int Task::getMpiRank()
{
    return this->mpiRank_;
}

/**
 * Get task status
 */
int Task::getStatus()
{
    return this->status_;
}

/**
 * Set task status
 * @input: status
 */
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

/**
 *
 */
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

/**
 *
 */
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

/**
 *
 */
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

/**
 *
 */
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

/**
 *
 */
Task* Task::getUpNeighbor()
{
    return this->upNeighbor_;
}

/**
 *
 */
Task* Task::getDownNeighbor()
{
    return this->downNeighbor_;
}

/**
 *
 */
Task* Task::getLeftNeighbor()
{
    return this->leftNeighbor_;
}

/**
 *
 */
Task* Task::getRightNeighbor()
{
    return this->rightNeighbor_;
}

/**
 *
 */
int Task::getUpNeighborRank(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        Task* up = this->rTasks_[layer]->getUpNeighbor();
        if (up)
        {
            return up->getNextRank_(layer);
        }
    }
    else
    {
        throw std::string("Can't get up neighbor");
    }

#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */

}

/**
 *
 */
int Task::getDownNeighborRank(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        Task* down = this->rTasks_[layer]->getDownNeighbor();
        if (down)
        {
            return down->getNextRank_(layer);
        }
    }
    else
    {
        throw std::string("Can't get down neighbor");
    }

#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */

}

/**
 *
 */
int Task::getLeftNeighborRank(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        Task* left = this->rTasks_[layer]->getLeftNeighbor();
        if (left)
        {
            return left->getNextRank_(layer);
        }
    }
    else
    {
        throw std::string("Can't get left neighbor");
    }

#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */

}

/**
 *
 */
int Task::getRightNeighborRank(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        Task* right = this->rTasks_[layer]->getRightNeighbor();
        if (right)
        {
            return right->getNextRank_(layer);
        }
    }
    else
    {
        throw std::string("Can't get right neighbor");
    }

#ifdef MPI_SUPPORT
        return MPI_PROC_NULL;
#else
        return -1;
#endif /* MPI_SUPPORT */

}

/**
 *
 */
void Task::setLocalGrid(double* grid)
{
    if (!grid)
    {
        throw std::string("Can't set local grid (bad pointer)");
    }

    this->grid_ = grid;
}

/**
 *
 */
void Task::setLocalNewGrid(double* newGrid)
{
    if (!newGrid)
    {
        throw std::string("Can't set local new grid (bad pointer)");
    }

    this->newGrid_ = newGrid;
}

/**
 * Get local grid
 * @input: redundancy layer
 * @return: pointer to grid
 */
double* Task::getLocalGrid(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        return this->rTasks_[layer]->grid_;
    }
    else
    {
        throw std::string("Can't get local grid");
    }

    return NULL;
}

/**
 * Get local new-grid
 * @input: redundancy layer
 * @return: pointer to new-grid
 */
double* Task::getLocalNewGrid(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        return this->rTasks_[layer]->newGrid_;
    }
    else
    {
        throw std::string("Can't get local new-grid");
    }

    return NULL;
}

/**
 *
 */
void Task::addRrank(int* rank)
{
    this->rRanks_.push_back(rank);
}

/**
 *
 */
void Task::addRtask(Task* task)
{
    this->rTasks_.push_back(task);
}

/**
 *
 */
void Task::addUpTag(int tag)
{
    this->upNeighborTags_.push_back(tag);
}

/**
 *
 */
void Task::addDownTag(int tag)
{
    this->downNeighborTags_.push_back(tag);
}

/**
 *
 */
void Task::addLeftTag(int tag)
{
    this->leftNeighborTags_.push_back(tag);
}

/**
 *
 */
void Task::addRightTag(int tag)
{
    this->rightNeighborTags_.push_back(tag);
}

/**
 *
 */
int Task::getUpTag(int layer)
{
    return this->getNextUpTag_(layer);
}

/**
 *
 */
int Task::getDownTag(int layer)
{
    return this->getNextDownTag_(layer);
}

/**
 *
 */
int Task::getLeftTag(int layer)
{
    return this->getNextLeftTag_(layer);
}

/**
 *
 */
int Task::getRightTag(int layer)
{
    return this->getNextRightTag_(layer);
}

/**
 * Swap grid and new-grid fields
 * @input: redundancy layer
 */
void Task::swapLocalGrids(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        std::swap(this->rTasks_[layer]->grid_,
                  this->rTasks_[layer]->newGrid_);
    }
    else
    {
        throw std::string("Can't swap grids");
    }
}

/**
 * Get number of redundancy layers
 */
int Task::getLayersNumber()
{
    return (int)this->rTasks_.size();
}

/**
 * Returns replacements vector
 */
std::vector<Task*> Task::getReplacements()
{
    return this->replacements_;
}

/**
 * Get x coordinates
 * @input: layer
 */
int Task::getX(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        return this->rTasks_[layer]->i_;
    }
    else
    {
        throw std::string("Can't get x");
    }
}

/**
 * Get y coordinates
 * @input: layer
 */
int Task::getY(int layer)
{
    if (layer < (int)this->rTasks_.size())
    {
        return this->rTasks_[layer]->j_;
    }
    else
    {
        throw std::string("Can't get y");
    }
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

    Task* replacement = this->getNextReplacement_();

    // This automaticly set rank and status
    this->setMpiRank(replacement->mpiRank_);

    replacement->reduceRepairAbility_();

    this->reduceRepairAbility_();

    replacement->addReplacement_(this);

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

    std::cout << "Replacements   : [ ";
    for (int i = 0; i < (int)this->replacements_.size(); ++i)
    {
        std::cout << this->replacements_[i]->mpiRank_ << " ("
                  << this->replacements_[i]->i_ << ", "
                  << this->replacements_[i]->j_ << ") ";
    }
    std::cout << "]" << std::endl;

    // TODO

    std::cout << std::endl;
}

/**
 * Show whole infomation about task
 * by redundancy layers
 */
void Task::printByLayers()
{
    std::cout << "Task [ " << this->i_ << ", "
              << this->j_ << " ]" << std::endl;

    std::cout << "MPI rank       : ";
    for (int i = 0; i < this->getLayersNumber(); ++i)
    {
        std::cout << this->rTasks_[i]->mpiRank_ << " ";
    }
    std::cout << std::endl;

    std::cout << "Up neighbor    : ";
    for (int i = 0; i < this->getLayersNumber(); ++i)
    {
        if (this->upNeighbor_)
        {
            std::cout << this->upNeighbor_->rTasks_[i]->mpiRank_ << " ";
        }
        else
        {
            std::cout << "-1" << " ";
        }
    }
    std::cout << std::endl;

    std::cout << "Left neighbor  : ";
    for (int i = 0; i < this->getLayersNumber(); ++i)
    {
        if (this->leftNeighbor_)
        {
            std::cout << this->leftNeighbor_->rTasks_[i]->mpiRank_ << " ";
        }
        else
        {
            std::cout << "-1" << " ";
        }
    }
    std::cout << std::endl;

    std::cout << "Down neighbor  : ";
    for (int i = 0; i < this->getLayersNumber(); ++i)
    {
        if (this->downNeighbor_)
        {
            std::cout << this->downNeighbor_->rTasks_[i]->mpiRank_ << " ";
        }
        else
        {
            std::cout << "-1" << " ";
        }
    }
    std::cout << std::endl;

    std::cout << "Right neighbor : ";
    for (int i = 0; i < this->getLayersNumber(); ++i)
    {
        if (this->rightNeighbor_)
        {
            std::cout << this->rightNeighbor_->rTasks_[i]->mpiRank_ << " ";
        }
        else
        {
            std::cout << "-1" << " ";
        }
    }
    std::cout << std::endl;

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
 * Add task to vector replacements
 * @input: task
 */
void Task::addReplacement_(Task* task)
{
    this->replacements_.push_back(task);
}

/**
 *
 */
Task* Task::getNextReplacement_()
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

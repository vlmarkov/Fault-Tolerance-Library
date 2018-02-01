#include "FailureSimulation.h"

#include <ctime>
#include <cstdlib>
#include <csignal>
#include <iostream>

FailureSimulation::FailureSimulation(const int size, const double rate, const int policy)
{
    this->policy_      = policy;
    this->alive_       = size - 1;
    this->procSize_    = size;
    this->failureRate_ = 1.0 / rate;
    this->invokeCnt_   = 0;

    srand((unsigned)time(0));
}

FailureSimulation::~FailureSimulation() { }

void FailureSimulation::generateFailure(const int mpiRank)
{
    int rankToKill = 0;

    switch (this->policy_)
    {
        case FAILURE_POLICY_RANDOM_HALF:
        {
            if (this->alive_ < this->procSize_ / 2)
            {
                std::cerr << "There is nothing to kill" << std::endl;
                return;
            }
            rankToKill = rand() % this->alive_--;
            break;
        }
        case FAILURE_POLICY_SERIAL_HALF_TAIL:
        {
            if (this->alive_ < this->procSize_ / 2)
            {
                std::cerr << "There is nothing to kill" << std::endl;
                return;
            }
            rankToKill = this->alive_--;
            break;
        }
        case FAILURE_POLICY_SERIAL_HALF_FRONT:
        {
            if (this->alive_ < this->procSize_ / 2)
            {
                std::cerr << "There is nothing to kill" << std::endl;
                return;
            }
            rankToKill = this->procSize_ - 1 - this->alive_--;
            break;
        }
        default:
        {
            throw std::string("Unknow failure policy");
            break;
        }
    }

    if (++this->invokeCnt_ == this->failureRate_)
    {
        if (rankToKill == mpiRank)
        {
            std::cerr << "Failure acquired for " << rankToKill << " rank" << std::endl;
            raise(SIGKILL);
        }
        this->invokeCnt_ = 0;
    }
}

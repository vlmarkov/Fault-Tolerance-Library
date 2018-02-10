/** @file */
#ifndef _FAILURE_SIMULATOR_H_
#define _FAILURE_SIMULATOR_H_

/**
 * @brief Available policies to failure simulation
 */
enum
{
    FAILURE_POLICY_RANDOM_HALF,       /**< Pure chaotic */
    FAILURE_POLICY_SERIAL_HALF_TAIL,  /**< From last to 1/2 */
    FAILURE_POLICY_SERIAL_HALF_FRONT  /**< From first to 1/2 */
};

/**
 * @brief Failure simulation abstraction
 * @details Generates failure of MPI-process via system call raise(SIGKILL) according failure rate and failure policy.
 * @author Markov V.A.
 * @date 2018
 * @version 0.1
 * @warning Work In Progress
 */
class FailureSimulator
{
    public:
        /**
         * @brief Main constructor
         * @param size an integer argument MPI-comm-world size
         * @param rate an real argument failure rate
         * @param policy an integer argument failure simulation policy
         * @details Creates failure simulation object with specific rate and policy
         */
        FailureSimulator(const int size, const double rate, const int policy);

        /**
         * @brief Destructor
         */
        ~FailureSimulator();

        /**
         * @brief Main failure generate method
         * @param mpiRank an integer argument self MPI-rank
         * @details Generates MPI-process failure
         * @code
         * if ((invokeCounter == failureRate) && (generateRank == mpiRank))
         * {
         *     raise(SIGKILL);
         * }
         * @endcode
         */
        void generateFailure(const int mpiRank);

    private:
        int    policy_;      /**< Failure policy */
        int    alive_;       /**< How many processes are alive */
        int    procSize_;    /**< How many processes we have */
        int    invokeCnt_;   /**< Invoke counter */
        double failureRate_; /**< How often does the MPI-process fail */
};

#endif /* _FAILURE_SIMULATOR_H_ */

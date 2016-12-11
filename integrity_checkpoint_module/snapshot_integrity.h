#ifndef _SNAPSHOT_INTEGRETY_
#define _SNAPSHOT_INTEGRETY_

#include <string>
#include <vector>

#define SNAPSHOT_DIR_NAME "snapshot"

/*****************************************************************************/
/* SnapshotDataFile class defenition                                         */
/*****************************************************************************/
class SnapshotDataFile
{
    public:
        int rank_;
        std::string fileName_;

        SnapshotDataFile(int index, std::string fileName);
        ~SnapshotDataFile();

        SnapshotDataFile & operator= (const SnapshotDataFile &rhs);
};

/*****************************************************************************/
/* SnapshotIntegrity class defenition                                        */
/*****************************************************************************/
class SnapshotIntegrity
{
    private:
        int commSize_;

        void addToSnapshotVector_( std::vector<std::vector<SnapshotDataFile>> &sVector, 
                                   int rank, 
                                   std::string fileName );

        static void sortSnapshotVector_(std::vector<SnapshotDataFile> &sVector);
        static bool sortSnapshotVectorValue_(SnapshotDataFile a, SnapshotDataFile b);

        void checkIntegity_(std::vector<std::vector<SnapshotDataFile>> &sVector);
        bool integritySnapshots_(std::vector<SnapshotDataFile> &sVector, int size);
        bool fileNameMatch_(std::string compareFileName, 
                            std::vector<SnapshotDataFile> sVector,
                            SnapshotDataFile &returnObject);

        static void printSnapshotVector_(const std::vector<SnapshotDataFile> sVector);

    public:
        SnapshotIntegrity(int commSize);
        ~SnapshotIntegrity();

        void getIntegritySnapshots();
    
};

#endif /* _SNAPSHOT_INTEGRETY_ */

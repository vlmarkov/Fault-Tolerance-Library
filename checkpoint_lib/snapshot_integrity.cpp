#include "snapshot_integrity.h"

#include <dirent.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

/*****************************************************************************/
/* SnapshotDataFile class implementation                                     */
/*****************************************************************************/
SnapshotDataFile::SnapshotDataFile(int rank, string fileName) \
                                        : rank_(rank), fileName_(fileName) {}

SnapshotDataFile::~SnapshotDataFile() {}

SnapshotDataFile& SnapshotDataFile::operator=(const SnapshotDataFile &rhs)
{
    this->rank_     = rhs.rank_;
    this->fileName_ = rhs.fileName_;

    return *this;
}


/*****************************************************************************/
/* SnapshotIntegrity class implementation                                    */
/*****************************************************************************/

/* Main constructor */
SnapshotIntegrity::SnapshotIntegrity(int commSize) : commSize_(commSize) {}

/* Destructor */
SnapshotIntegrity::~SnapshotIntegrity() {}

/* Main method */
void SnapshotIntegrity::getIntegritySnapshots()
{
    DIR *dir = NULL;
    vector<vector<SnapshotDataFile>> snapshotVector(this->commSize_);

    // Step one: gather snapshot files
    for (auto i = 0; i < this->commSize_; i++) {
        char snapshot_path[256] = { 0 };
        sprintf(snapshot_path, "./%s/%d/", SNAPSHOT_DIR_NAME, i);

        try {
            dir = opendir(snapshot_path);
            if (!dir)
                throw string("Can't open directory: " + string(snapshot_path));
        }
        catch (string err) {
            cerr << "[ERROR]" << err << endl;
            exit(1);
        }

        for (struct dirent *file = readdir(dir); file != NULL; file = readdir(dir)) {
            // If not a current direacory and not backwards directory
            if ((strcmp(".", file->d_name) == 0) || (strcmp("..", file->d_name) == 0))
                continue;

            addToSnapshotVector_(snapshotVector, i, string(file->d_name));    
        }
        closedir(dir);
    }
    // Step two: sort snapshot files
    for_each(snapshotVector.begin(), snapshotVector.end(), sortSnapshotVector_);
    
    //for_each(snapshotVector.begin(), snapshotVector.end(), printSnapshotVector_);

    // Step three: check integity of snapshot files
    checkIntegity_(snapshotVector);
}

void SnapshotIntegrity::addToSnapshotVector_(vector<vector<SnapshotDataFile>> &sVector, 
                                             int rank, 
                                             string fileName )
{
    try {
        if ((rank < 0) && (rank > (int)(sVector.size() - 1)))
            throw string("rank value is out of range");

        if (fileName.size() == 0)
            throw string("fileName value is empty");
    }

    catch (string err) {
        cerr << "[ERROR]" << err << endl;
        exit(1);
    }

    SnapshotDataFile newSnapshot(rank, fileName);
    sVector[rank].push_back(newSnapshot);
}

void SnapshotIntegrity::sortSnapshotVector_(vector<SnapshotDataFile> &sVector)
{
    std::sort(sVector.begin(), sVector.end(), sortSnapshotVectorValue_);
}

bool SnapshotIntegrity::sortSnapshotVectorValue_(SnapshotDataFile a, SnapshotDataFile b)
{
    if (a.fileName_.compare(b.fileName_) < 0) {
        return true;
    } else {
        return false;
    }
}

void SnapshotIntegrity::checkIntegity_(vector<vector<SnapshotDataFile>> &sVector)
{
    int globalSize = sVector.size();

    // Iterate 0 sub-vector
    for (auto i = sVector[0].rbegin(); i != sVector[0].rend(); i++) {
        vector<SnapshotDataFile> integrityVector;
       
        integrityVector.push_back(*i);

        for (auto j = 1; j < globalSize; j++) {
            SnapshotDataFile returnObject(-1, "NULL");
            if (fileNameMatch_(i->fileName_, sVector[j], returnObject)) {
                integrityVector.push_back(returnObject);
            } else {
                break;
            }
        }

        if (integritySnapshots_(integrityVector, globalSize)) {
            saveSnapshotVector_(integrityVector);
            return;
        }
    }
}

bool SnapshotIntegrity::integritySnapshots_(vector<SnapshotDataFile> &sVector, int size)
{
    auto counter = 0;

    for (auto i = 0; i < (int)sVector.size(); i++) {
        string line;
        string file;
        
        file = "./";
        file += SNAPSHOT_DIR_NAME;
        file += "/";
        file += to_string(sVector[i].rank_);
        file += "/";
        file += sVector[i].fileName_;

        ifstream checkSnapohot(file);

        try {
            if (checkSnapohot.is_open()) {
                while (getline(checkSnapohot, line)) {
                    cout << line << endl;
                    if (line.compare(INTEGRITY_SNAPSHOT) == 0) {
                        counter++;
                        break;
                    }
                }
                checkSnapohot.close();
            } else {
                throw string("Can't open file " + file);
            }
        }

        catch (string err) {
            cerr << "[ERROR] " << err << endl;
            exit(1);
        }

    }

    if (counter == size) {
        return true;
    } else {
        return false;
    }
}


bool SnapshotIntegrity::fileNameMatch_( string compareFileName, 
                                        vector<SnapshotDataFile> sVector,
                                        SnapshotDataFile &returnObject )
{
    for (auto i = sVector.rbegin(); i != sVector.rend(); i++) {
        if (compareFileName.compare(i->fileName_) == 0) {
            returnObject = *i;
            return true;
        }
    }
    return false;
}

void SnapshotIntegrity::printSnapshotVector_(const vector<SnapshotDataFile> sVector)
{
    for (auto i: sVector) {
        std::cout << "[" << i.rank_ << "=" << i.fileName_ << "]";
    }
    cout << endl;
}

void SnapshotIntegrity::saveSnapshotVector_(const vector<SnapshotDataFile> sVector)
{
    ofstream outFile;

    outFile.open(INTEGRITY_SNAPSHOT_FILE);
    try {
        if (outFile.is_open()) {
            for (auto i: sVector) {
                outFile << string(  to_string(i.rank_)    \
                                    + string("=")         \
                                    + string(i.fileName_) \
                                    + string("\n"));
            }
            outFile.close();
        } else {
            throw string("Can't open file " + string(INTEGRITY_SNAPSHOT_FILE));
        }
    }

    catch (string err) {
        cerr << "[ERROR] " << err << endl;
        exit(1);
    }
}

int main(int argc, char const *argv[])
{
    try {
        if (argc < 1) {
            throw "not enought arguments";
        }
        SnapshotIntegrity snapshotIntegrityObject(atoi(argv[1]));

        snapshotIntegrityObject.getIntegritySnapshots();
    }
    catch (string err) {
        cerr << "[ERROR] " << err << endl;
        exit(0); 
    }
    return 0;
}

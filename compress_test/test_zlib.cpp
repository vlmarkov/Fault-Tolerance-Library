#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#include "zlib.h"

using namespace std;

int main(int argc, char const *argv[])
{
    double *dest = NULL;
    double *source = NULL;

    uLongf *destLen = NULL;

    uLong sourceLen = 65536;
    uLong sourceLen_ = sourceLen + 100;
    
    destLen = &sourceLen_;

    
    source = new double[sourceLen];
    dest = new double[sourceLen];

    for (uLong i = 0; i < sourceLen; ++i) {
        source[i] = i;
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    int rc = compress((Bytef*)dest, destLen, (Bytef*)source, sourceLen);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    switch(rc)
    {
        case Z_OK: //в случае успеха
            cout << "Compress Succes" << endl;
            break;

        case Z_MEM_ERROR: //если недостаточно памяти
            cerr << "Not enought memory" << endl;
            break;

        case Z_BUF_ERROR: //если недостаточно места в буфере назначения
            cerr << "Not enought memory in destination buffer" << endl;
            break;

        default:
            cerr << "Unexpected error" << endl;
            break;
    }

    cout << "Size before compress : " << sourceLen << endl;
    cout << "Size atfter compress : " << *destLen << endl;
    cout << "Compress ratio       : " << sourceLen / *destLen << endl;

    //using namespace std::chrono::;

    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    cout << "It took me " << time_span.count() << " seconds." << endl;

    delete[] source;
    delete[] dest;

    return 0;
}
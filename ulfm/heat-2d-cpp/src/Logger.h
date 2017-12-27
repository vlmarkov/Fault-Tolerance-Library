#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <string>
#include <fstream>
#include <iostream>

class Logger
{
    public:
        Logger(const int rank);
        ~Logger();
        Logger& operator<<(std::string str);
        void trace(const char* func, const int line, std::string str);

    private:
        std::string   name_;
        std::ofstream file_;
};

#endif /* _LOGGER_H_ */

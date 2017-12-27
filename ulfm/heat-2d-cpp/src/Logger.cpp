#include "Logger.h"

Logger::Logger(const int rank)
{
    char tmp[128];
    sprintf(tmp, "log-%d-rank.txt", rank);

    this->name_ = std::string(tmp);

    this->file_.open(this->name_);

    if (!this->file_.is_open())
    {
        throw std::string("Can't open log file");
    }
}

Logger::~Logger()
{
    if (this->file_.is_open())
    {
        this->file_.close();
    }
}

Logger& Logger::operator<<(std::string str)
{
    if (this->file_.is_open())
    {
        this->file_ << str << std::endl;
    }
    else
    {
        throw std::string("Can't write to log file");
    }
    
    return *this;
}

void Logger::trace(const char *func, const int line, std::string str)
{
    if (this->file_.is_open())
    {
        this->file_ << "<" << func << ">"
                    << "<" << line << ">"
                    << " "
                    << str << std::endl;
    }
    else
    {
        throw std::string("Can't write to log file");
    }
}

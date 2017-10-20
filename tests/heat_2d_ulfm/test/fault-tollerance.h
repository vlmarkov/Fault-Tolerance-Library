#ifndef _FAULT_TOLLERANCE_H_
#define _FAULT_TOLLERANCE_H_

#include <vector>

class TaskCell
{
    public:
        TaskCell(int _p, int _b, int _x, int _y);
        ~TaskCell();
        
        int proc;
        int block;
        int x;
        int y;    
};

class TaskGrid
{
    public:
        TaskGrid(int px, int py, int nx, int ny);
        ~TaskGrid();

        void show();
        void markDeadProc(int fail);
        void repair();
        bool tryToReassingRank(const int x, const int y);
    
    private:
        std::vector<std::vector<TaskCell> > taskGrid_;
};

#endif /* _FAULT_TOLLERANCE_H_ */
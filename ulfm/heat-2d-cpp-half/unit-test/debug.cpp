#include <iostream>
#include <string>

#include "../Grid.h"
#include "../Task.h"

int main(int argc, char const *argv[])
{
    try
    {
        Grid testGrid(1024, 1024, 256, 256, 8, 8);

        testGrid.printPretyTags();
        testGrid.printPretyNeighbors();
/*
        testGrid.kill(3);
        testGrid.repair();

        testGrid.kill(2);
        testGrid.repair();

        testGrid.print();

        for (int i = 0; i < 72; i++)
        {
            testGrid.kill(143 - i);
            testGrid.repair();
        }
*/
    }
    catch (std::string err)
    {
        std::cerr << err << std::endl;
        return -1;
    }

    return 0;
}

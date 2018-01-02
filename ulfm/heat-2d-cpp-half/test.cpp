#include <iostream>
#include <string>

#include "Grid.h"
#include "Task.h"

int main(int argc, char const *argv[])
{
	try
	{
		Grid testGrid(1024, 1024, 256, 256, 2, 2);

		//testGrid.kill(15);

		//testGrid.repair();

		testGrid.kill(1);

		testGrid.repair();

		testGrid.print();

		//Task *t = testGrid.getTask(5);

		//std::cout << t->getUpNeighborRank(0) << std::endl;
		//std::cout << t->getUpNeighborRank(1) << std::endl;
		//std::cout << t->getDownNeighborRank(0) << std::endl;
		//std::cout << t->getDownNeighborRank(1) << std::endl;
	}
	catch (std::string err)
	{
		std::cerr << err << std::endl;
		return -1;
	}

	return 0;
}

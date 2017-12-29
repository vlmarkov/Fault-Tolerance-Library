#include <iostream>
#include <string>

#include "Grid.h"
#include "Task.h"

int main(int argc, char const *argv[])
{
	try
	{
		Grid testGrid(1024, 1024, 256, 256, 4, 4);

		testGrid.kill(15);

		testGrid.print();
	}
	catch (std::string err)
	{
		std::cerr << err << std::endl;
		return -1;
	}

	return 0;
}

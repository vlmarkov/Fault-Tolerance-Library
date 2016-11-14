#include "../checkpoint_lib/checkpoint_lib.h"

#include <signal.h>

void time_handler(int sig)
{
    printf("%s():%d\n",__func__, __LINE__);
}

int main(int argc, char const *argv[])
{
	CHECKPOINT_LIB_INIT(3, 5);

	CHECKPOINT_ASSIGN(&&point1);
	CHECKPOINT_ASSIGN(&&point2);
	CHECKPOINT_ASSIGN(&&point3);

	CHECKPOINT_HANDLER(time_handler);
	CHECKPOINT_TIMER();

	printf("%s():%d\n",__func__, __LINE__);

	CHECKPOINT_GOTO(2);

	printf("%s():%d\n",__func__, __LINE__);

	CHECKPOINT_SET(point1);
	printf("%s():%d\n",__func__, __LINE__);
	
	CHECKPOINT_SET(point2);
	printf("%s():%d\n",__func__, __LINE__);
	
	CHECKPOINT_SET(point3);
	printf("%s():%d\n",__func__, __LINE__);

	while(1) {

	}
	
	return 0;
}
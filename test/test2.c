
#include <stdio.h>
#include <signal.h>
#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>

void **CPL_GLOBAL_JUMP_TABLE;

void **init_table_(int size);
int  get_checkpoint_idx_by_name_(void **table, int size, void *name);

int get_checkpoint_idx_by_name_(void **table, int size, void *name)
{
    int i;
    for (i = 0; i < size; i++) {
        if (table[i] == name) {
            break;
        }
    }
    return i;
}


void **init_table_(int size)
{
    void **jump_table = (void**) malloc (sizeof(void*) * size);
    if (!jump_table) {
        fprintf(stderr, "[ERROR] can't allocate memory for CPL_GLOBAL_JUMP_TABLE\n");
        exit(1);
    }
    return jump_table;
}


int main(int argc, char const *argv[])
{
    void **table = init_table_(3);
    int cnt = 0;

    table[cnt++] = &&name1;
    table[cnt++] = &&name2;
    table[cnt++] = &&name3;



    name1 :

    cnt = 5;

    name2 :

    cnt = 6;

    name3 :

    printf("idx = %d\n", get_checkpoint_idx_by_name_(table, 3, &&name1));
    printf("idx = %d\n", get_checkpoint_idx_by_name_(table, 3, &&name2));
    printf("idx = %d\n", get_checkpoint_idx_by_name_(table, 3, &&name3));

    printf("addr = %x\n", &&name1);
    printf("addr = %x\n", &&name2);
    printf("addr = %x\n", &&name3);


    return 0;
}
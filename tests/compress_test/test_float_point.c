#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <inttypes.h>
/*
double double_trunc(double x, int zerobits)
{
    // mask is e.g. 0xffffffffffff0000 for zerobits==16
    uint64_t mask      = -(1LL << zerobits);
    uint64_t floatbits = (*((uint64_t*)(&x)));
    floatbits          &= mask;
    
    return *((double*)(&floatbits));
}
*/
double double_trunc(double x, int zerobits)
{
    uint64_t mask      = -(1LL << 11);
    uint64_t floatbits = (*((uint64_t*)(&x)));
    //uint64_t mantisa   = floatbits & 0xFFFFFFFFFFFFF;
    floatbits          &= floatbits ^ mask;

    //printf("%x\n", mask);
    //printf("%x\n", mantisa);
    //printf("%x\n", mantisa ^ mask);
    
    return *((double*)(&floatbits));
}

void print_double_raw_bit(void *double_number)
{
    char *a  = (char *)double_number;
    int size = sizeof(double) / sizeof(char); 

    for (int i = 0; i < size; i++) {
        printf("[%x]\n", *a++);
    }
    printf("\n");
}

void save_double_trunc(int fd, void *double_number)
{
    char *a  = (char *)double_number;
    int size = sizeof(double) / sizeof(char);
    char trunc_raw_data[size - 3];
    char new_line = '\n';

    a += 3;

    for (int i = 3; i < size; i++) {
        trunc_raw_data[i] = *a++;
        write(fd, &trunc_raw_data[i], 1);
    }
    write(fd, &new_line, 1);

    double *c = NULL;

    c = (double *)trunc_raw_data;

    //printf("%f\n", *c);

}

int main(int argc, char const *argv[])
{
    double a = 1.23456;

    double b = double_trunc(a, 30);

    print_double_raw_bit(&a);
    print_double_raw_bit(&b);

    FILE *file_original = fopen("file_original.txt", "wb");
    //FILE *file_trunc    = open("file_trunc.txt", "wb");
    int file_trunc    = open("file_trunc.txt", O_CREAT | O_TRUNC | O_WRONLY, S_IRWXU);

    for (int i = 0; i < 1000; i++) {
        fprintf(file_original, "%f\n", a);
        save_double_trunc(file_trunc, &b);
    }

    fclose(file_original);
    close(file_trunc);

    printf("->%f\n", a);
    printf("->%f\n", b);


    unsigned int n1 = 0x0;
    unsigned int n2 = 0xf;

    n1 = n2 ^ n1;

    printf("%x\n", n1);

    return 0;
}
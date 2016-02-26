#include <stdio.h>

int main() {
    unsigned int number_of_lines = 0;
    FILE *infile = fopen("pat00001_N_1.fq", "r");  /* #test.txt", "r"); */
    int ch;

    while (EOF != (ch=getc(infile)))
        if ('\n' == ch)
            ++number_of_lines;
    printf("%u\n", number_of_lines);
    return 0;
}

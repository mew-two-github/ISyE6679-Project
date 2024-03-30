#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main()
{
    char str[128];
    int result;
    FILE* f = fopen("data/Case 14 Branch.csv", "r");

    /*...*/

    do {
        result = fscanf(f, "%127[^,\n]", str);

        if(result == 0)
        {
            result = fscanf(f, "%*c");
        }
        else
        {
            //Put here whatever you want to do with your value.
            printf("%s\n", str);
        }

    } while(result != EOF);

    return 0;
}
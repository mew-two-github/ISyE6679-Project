#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* getfield(char* line, int num)
{
    char* tok;
    for (tok = strtok(line, ",");
            tok && *tok;
            tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}

int main()
{
    FILE* stream = fopen("data/Case 14 Branch.csv", "r");
    char** From_Bus,To_Bus,Resistance,Reactance,Susceptance;
    double RateA,RateB,RateC,Ratio,Angle,Status,MinAngleDiff,MaxAngleDiff,G,B;
    char line[1024];
    int i = 0;
    while (fgets(line, 1024, stream))
    {
        
        char* tmp = strdup(line);
        // printf("Field 3 would be %s\n", getfield(tmp, 3));
        // NOTE strtok clobbers tmp
        *(From_Bus[i]) = *(getfield(tmp, 1));
        printf("%s\n", From_Bus[i]);
        free(tmp);
        i = i + 1;
    }
    for ( int j = 0; j < i; j++)
    {
        printf( "%f\n", From_Bus[j] );
    }
}
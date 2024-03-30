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
    char* From_Bus[100],To_Bus,Resistance,Reactance,Susceptance;
    double RateA,RateB,RateC,Ratio,Angle,Status,MinAngleDiff,MaxAngleDiff,G,B;
    char line[1024];
    int i = 0;
    while (fgets(line, 1024, stream))
    {
        
        char* tmp = strdup(line);
        // printf("Field 3 would be %s\n", getfield(tmp, 3))
        char * token = strtok(tmp, ",");
        // loop through the string to extract all other tokens
        char* split[100];
        int j =0;
        // printf( "%s\n", token );
        while( token != NULL ) 
        {
            // printf( "%s\n", token ); //printing each token
            if(j == 0)
            {
                
                From_Bus[i] = strdup(token);
                //printf( "%s\n", token );
            }
            token = strtok(NULL, ",");
            j = j+1;
        }
        free(token);
        free(tmp);
        i = i + 1;
    }
    printf("%i\n",i);

    for ( int j = 0; j < i; j++)
    {
        printf( "%s\n", From_Bus[j] );
    }
    double lol = 1;
    for ( int j = 1; j < i; j++)
    {
        char* hmm = From_Bus[j];
        lol = strtod(hmm,NULL );
        // printf( "%s\n", hmm );
        printf( "%f \n", lol );
    }

}
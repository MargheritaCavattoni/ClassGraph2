#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {


    // OPENING THE FILE
    char const* const fileName = argv[1];
    //char fileName[] = "/Users/francescoandreace/Downloads/S1_12.paf";
    FILE* file = fopen(fileName,"r");
    if (file == NULL)
    {
        printf("Sorry, file doesn't exist.");
        return 0;
    }

    // VARIABLES INITIALIZATION
    char line[2000]; // line of the file
    char *saver[9]; // array in which strings are stored
    int countLine = 0; // number of the line being read
    int countEl = 0; // number of elements opened

    // READING THE LINES AND SAVING THE TOKENS
    while (fgets(line, sizeof(line), file)){

        saver[countEl] = strtok(line, "\t");
        while( saver[countEl] != NULL ) {
            if (countEl == 4){
                if (strcmp(saver[countEl],"+") == 0){
                    saver[countEl] = "0";
                }
                else{
                    saver[countEl] = "1";
                }
            }

            saver[++countEl] = strtok(NULL, "\t");

            // BREAK IF READING ELEMENT NUMBER 9
            if (countEl == 9) {
                countEl = 0;
                break;
            }
        }

        //PRINTING
        if (++countLine == 1) {
            printf("VT\t0\t%s\n",saver[1]);
        }
        printf("ED\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\n",saver[0],saver[5],saver[2],saver[3],saver[1],saver[7],saver[8],saver[6],saver[4]);
    }

    //CLOSING

    fclose(file);
    return 0;
}

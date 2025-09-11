//=======================
// To compile this file: 
// g++ -c helloWorld.C helloWorld
// g++ -o helloWorld helloWorld.o
//=======================

#include <stdio.h>

int main(int argc, char* argv[])
{
	printf("Hello World!\n");

    printf("Total Argument Counts (argc) =%d \n", argc);
	for( int i=0; i<argc; i++ )
	{
	    printf("  argv[%d] = %s\n", i, argv[i]);
	}

    printf("\n");

    return 0;
}

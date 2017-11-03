/*** Allard Group jun.allard@uci.edu                    ***/

void binaryConversion();


/*******************************************************************************/
//  GLOBAL VARIABLES for binary conversion
/*******************************************************************************/

int remainder,base,num,no_of_1s;

/*******************************************************************************/

void binaryConversion()
{
    
    /*
     * C program to accept a decimal number and convert it to binary
     * and count the number of 1's in the binary number
     */
    // Adapted from: http://www.sanfoundry.com/c-program-decimal-binary-count-1-binary/
    
    for (i=0;i<sizeOfRateMatrix;i++)
    {
        base     = 1;
        binary   = 0;
        no_of_1s = 0;
        num      = i;
        while (num > 0)
        {
            remainder = num % 2;
            /*  To count no.of 1s */
            if (remainder == 1)
            {
                no_of_1s++;
            }
            binary = binary + remainder * base;
            num = num / 2;
            base = base * 10;
        }
        totalBound[i] = no_of_1s;
     }
    
    //debugging
    if(1)
    {
        // print matrix of number bound
        for (i=0;i<sizeOfRateMatrix;i++)
        {
            printf("Number bound: %lf \n", totalBound[i]);
        }
    }
    
    
}




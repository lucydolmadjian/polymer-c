/*** Allard Group jun.allard@uci.edu                    ***/

void runGillespie();

void runGillespie()
{
    char *line = malloc(1000);
    const char delim[2] = " ";
    char *token;
    
    size_t count;
    
    //read transitionMatrix from file
    ratesFile = fopen(matrixName, "r");
    
    if (ratesFile == NULL )
    {
        printf("Could not open file.");
        exit(0);
    }
    else
    {

        while(getline(&line, &count, ratesFile) != -1)
        {
            printf("This is the next row: %s\n", line);
            printf("This is the count: %zu\n",count);
            printf("This is i: %d\n",i);
            
            token = strtok(line,delim);
            
            while ( token != NULL )
            {
                transitionMatrix[i][j]=atof(token);
                
                token = strtok(NULL, delim);
                printf("%.8f\t",transitionMatrix[i][j]);
                j++;
            }
//            
//            for (j=0;j<64;j++)
//            {
//                sscanf(line,"%lf",&transitionMatrix[i][j]);
//                printf("%lf\t",transitionMatrix[i][j]);
//            }
            printf("\n");
            i++;
            
        }
        
//        //debugging
//        for (i=0;i<64;i++)
//        {
//            for (j=0; j<64;j++)
//            {
//                printf("%f\t",transitionMatrix[i][j]);
//            }
//            printf("\n");
//        }

        
//        for (i=0;i<64;i++)
//        {
//            for (j=0;j<64;j++)
//            {
//                if(!fscanf(ratesFile,"%lf", &transitionMatrix[i][j]))
//                      break;
//                printf("%lf\n",transitionMatrix[i][j]);
//            }
//        }
        
    }
    
    
    
//    char line[200];
//    iy=0;
//    
//    while (fgets(line, sizeof(line), ratesFile))
//    {
//        iSite[iy]=atoi(line);
//        iy++;
//    }
//    
//    fclose(iSiteList);
//    
//    iSiteTotal=iy;
    
  

    
    
    
    
    
}




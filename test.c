#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
//walter malone
//7/13/2015






void quicksort(float arr[], int left, int right)
{
    // Base case
    if (left >= right)
    {
        return;
    }

    // Choose pivot to be the last element in the subarray
    float pivot = arr[right];

    // Index indicating the "split" between elements smaller than pivot and 
    // elements greater than pivot
    int counter = left;
    // Traverse through array from l to r
    for (int i = left; i <= right; i++)
    {
        // If an element less than or equal to the pivot is found...
        if (arr[i] <= pivot)
        {
            // is to the left of all elements greater than pivot
                   float  temporary = arr[counter];
                   arr[counter] = arr[i];
                   arr[i] = temporary;

            //swap(&arr[cnt], &arr[i]);

            // Make sure to increment cnt so we can keep track of what to swap
            // arr[i] with
            counter++;
        }
    }
    quicksort(arr, left, counter-2); // Recursively sort the left side of pivot
    quicksort(arr, counter, right);   // Recursively sort the right side of pivot
}



int main(int argc, char * argv[] )
{
  //defines some basic variables and pointers
  //"coverg1" is a percent error convergence criteria you may modify
  //right now it is set at 10% error
  //"coverge2" is an absolute value error cpnvergence criteria you may modify
  //right now it is set 0.2.  It covers the cases our output data is close to zero
  //the rest of the variables are for storing data, positions in arrays, or are flags
  FILE * output;
  FILE * testoutput;
  int hasFailed = 0;
  float converg1 = 10;
  float converg2 = 0.2;
  int halt = 0;
  int position1 = 0;
  int position2  = 0;
  float out2[1000000];
  float out1[1000000];



  //opens the files "output" and "testoutput"
  //output comes with code.  testoutput is generated with the make file
   if((output = fopen ("output", "r"))== NULL)                              //opens appropiate files
     {
         printf("%s", "output does not exist.  Try again\n");
         return -1;
     }

   if((testoutput = fopen ("testoutput", "r"))== NULL)                      
     {
         printf("%s", "The test failed to write an output.  Try again\n");
         return -1;
     }

   //we store the data using these two pointers, garbage1 and garbage2.   
    char ** garbage1 = malloc(sizeof(char * )*1000000);                     //points to 100k char *, can adjust size as needed
    char ** garbage2 = malloc(sizeof(char * )*1000000);
    char ** back1  = garbage1;                                             // points to the begining of our list of input
    char ** back2  = garbage2;


   //some counters
   int o = 0;                                                              
   int t = 0;
   int count = 0;
   //initials some more memory for us to store our data


   while(count < 1000000)
   {
       *( garbage1+count) = malloc(sizeof(char)*80);                       // each string will be able to hold 80 characters. 
       *( garbage2+count) = malloc(sizeof(char)*80);
       count++;
   }




    //TestResults is a file that simply tells us if we failed the test
    //DetailedTestResults is a file that tells us how we failed
    FILE * results;
    FILE * detailedResults;
    detailedResults = fopen("DetailedTestResults","a");                                                       
    results = fopen("TestResults","a");

   //another counter, and a print line to tell us what test we are on.
   //the name of the test is taken from, user input
   int WillRead = 0;                                                             //WillRead tells us when to store data
   fprintf(results,"Current Test: %s \n",argv[1]  );
   fprintf(detailedResults,"Current Test: %s \n",argv[1]  );


   //we scan the first file, "output"
   while(fscanf(output,"%s",*(garbage1+o)) != EOF)                        //reading our first input file "output", the test file
   {
      //num tells us if we scanned a number
      int isNum = 0;

      //I needed some keywords to scan from the file.  This was the most challenging aspect 
      //if the file contains the word "Structure" then I have found we are always reading the line "Final Strucutre" 
      //if we come across "Finial Strucutre" we store the input from that point on, those are our results
      if(strcmp(*(garbage1+o), "Structure") == 0)
      {
         halt = 1;
         WillRead =1;
      }


      //This line looks for a few keywords to begin storing data if we never come across "Final Strucutre"
      if( (strcmp(*(garbage1+o), "Frequencies") == 0) || (strcmp(*(garbage1+o), "Forces") == 0) || (strcmp(*(garbage1+o), "Energies") == 0) || (strcmp(*(garbage1+o), "Total") == 0))                          
      {
         WillRead=1;

      }

     //This line looks for a few keywords to stop storing data
      if( (strcmp(*(garbage1+o), "time") == 0 ) || (strcmp(*(garbage1+o), "\n") == 0)  || (strcmp(*(garbage1+o), "Transition") == 0  || (strcmp(*(garbage1+o), "guesses")) == 0))                            
      {
         WillRead=0;
      }

      //if  WillRead=1 then we have come across an area of interest and can test if it a floating point number
      if(WillRead == 1)                                                          
      {
         int c = 0;
         while(*((*(garbage1+o))+c) != '\0' && c < 80)                    //loops through each string 
         {
            if (*((*(garbage1+o))+c) == '.')                              //determines if the string is a float
            {                                                             //There is a predetermined function for this but it fails to read "-"
                isNum = 1;
            }
            c++;

         }  
         c=0;

         //I need to look through the input again to look for false positives, (A.U.) ect.
         while(*((*(garbage1+o))+c) != '\0' && c < 80)
         {
            if (*((*(garbage1+o))+c) == 'A'|| *((*(garbage1+o))+c) == 'M')
            {
                isNum = 0;
            }
            c++;

         }

      }

      if (isNum == 1 && WillRead ==1)                                                  //if the string is a float and we are in the right part of the output we store it 
      {
          //If we ran into the word "Structure we must reord the location of the next float
          if( halt ==1 )
            {
               position1 = o;
               halt = 0;
            }                                                                   // by iterating the pointer
          o++;
      }

    }


   WillRead = 0;  



   //here we repeat the process for the next file
   //should probably make a function for this
   while(fscanf(testoutput,"%s",*(garbage2+t)) != EOF)                   
   {
      int isNum = 0;

         if(strcmp(*(garbage2+t), "Structure") == 0)
         {
            halt = 1;
            WillRead =1;
         }


      if( (strcmp(*(garbage2+t), "Frequencies") == 0) || (strcmp(*(garbage2+t), "Forces") == 0) || (strcmp(*(garbage2+t), "Energies") == 0) || (strcmp(*(garbage2+t), "Total") == 0))
      {
         WillRead=1;


      }

      if( (strcmp(*(garbage2+t), "time") == 0) || (strcmp(*(garbage2+t), "\n")) == 0  ||  (strcmp(*(garbage2+t), "Transition")) == 0  || (strcmp(*(garbage2+t), "guesses")) == 0    )
      {
         WillRead=0;
      }


      if(WillRead == 1)
      {
         int c = 0;
         while(*((*(garbage2+t))+c) != '\0' && c < 80)
         {
            if (*((*(garbage2+t))+c) == '.')
            {
                isNum = 1;
            }
            c++;

         }
         c=0;
         while(*((*(garbage2+t))+c) != '\0' && c < 80)
         {
            if (*((*(garbage2+t))+c) == 'A'|| *((*(garbage2+t))+c) == 'M')
            {
                isNum = 0;
            }
            c++;

         }



      }


      if (isNum == 1 && WillRead ==1)
      {
            if( halt ==1 )
            {
               position2 = t;
               halt = 0;
            }                                                 
             t++;
      }
    }




   //resets our pointers
   garbage1 = back1;                                                          
   garbage2 = back2;

   int ofinal = (o-position1);
   int tfinal = (t-position2);
   fprintf(detailedResults, "The total number of significant  paramters are %d and  %d \n", ofinal,tfinal);
   if ( (o != t) && (ofinal) != (tfinal))                                                               //if the quaniity of our outputs differ we have failed the test
   {
     hasFailed = 1;
     printf("Test failed.  The outpout, %d, and testoutput, %d, files contain a different number of results\n.", ofinal, tfinal);
     fprintf(detailedResults, "Test failed.  The outpout, %d, and testoutput, %d, files contain a different number of results\n.", ofinal, tfinal); 
     goto end;
   }
   int k = 0;
   while(k<ofinal)                                                                 //tests to see if our outputs are the same
   {

      out1[k] = atof(*(garbage1+k+position1));
      out2[k] = atof(*(garbage2+k+position2));
      k++;
   }
   quicksort(out1, 0, ofinal-1);

   //printf("The first number is: %f\n", out1[1]);
   quicksort(out2, 0, tfinal-1);


   k = 0;
   fprintf(detailedResults, " output value,  testoutput value,  difference, error\n ");
   while(k<ofinal)                                                                 //tests to see if our outputs are the same
   {
      float x = out1[k];
      float y = out2[k];
      float convergence1 = 0;
      if( x != 0 && y != 0)
      {
         convergence1 = ((x-y)/x)*100;
      }
      float convergence2 = x-y;
      fprintf(detailedResults, "%f %f %f %f \n",x, y,  convergence2, convergence1);
      if( (convergence1 > converg1 || convergence1 < -converg1) && (convergence2 > converg2 || convergence2 < -converg2))
      {
         hasFailed = 1;
         fprintf(detailedResults, "Test failed.  %f and %f are not within the covergence criteria\n", x, y);
      }
      k++;

   }

end:                                                              
   //determines if we failed the test and appends the output accordingly
   printf(" Test complete: cleaning\n");
   if(hasFailed == 0)
   {
   fprintf(results, "The test was SUCCESSFUL.\n\n");
   fprintf(detailedResults, "The test was SUCCESSFUL. \n\n");
   printf( "The test was successful. \n");
   }
   else
   {
   fprintf(results, "The test was UNSUCCESSFUL.  Better check %s \n \n",argv[1] );
   fprintf(detailedResults, "The test was UNSUCCESSFUL.  Better check %s \n \n",argv[1]);
   printf("The test was unsuccessful.  Better check %s \n", argv[1]);
   }

    //frees our memory.
    count = 0;
    garbage1 = back1;
    garbage2 = back2;
    while(count < 1000000)
    {
       free(*(garbage1+count));
       free(*(garbage2+count));
       count++;
    }
    garbage1 = back1;
    garbage2 = back2;

    free (garbage1);
    free (garbage2);

    //closes our files
    fclose(results);
    fclose(output);
    fclose(testoutput);
    fclose(detailedResults);
    return 0 ;
}

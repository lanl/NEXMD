#include <stdlib.h>
#include <stdio.h>
#include "KFc.h"
// Note: see KFc.h for a short description of the public functions,
// and KFReader.c for the function definitions.

double *gradient_unordered;

void read_adf_results_(int *natoms, double *energy, double *gradient, int *use_dftb, double *dipxyz, char *keyfile, int *do_grad )
{
    KFFile kf;
    double *gradients;
    int *ordering;
    int internal_order[2][*natoms];
    
    // Note: opening a file is a heavy-weight operation;
    // open it once and close it when done reading data.
    
    // Open keyfile; will exit if the file DNE or fails to open
    if(*use_dftb==0)
    {
      if (openKFFile(&kf, keyfile) < 0) exit(1);
      gradient_unordered = malloc(*natoms*3*sizeof(double));
      ordering = malloc(*natoms*2*sizeof(int));
    }
    else //Doing DFTB
    {
      if (openKFFile(&kf, keyfile) < 0) exit(1);
    }
    
    // Allocate memory...
    gradients = malloc(*natoms*3*sizeof(double));
    
    // Get data...
    if(*use_dftb==0) // Reorder atoms if we are not using DFTB
    {
      if (do_grad)
        getKFData(&kf, "GeoOpt%Gradients_CART", (void*) gradient_unordered);
      getKFData(&kf, "Energy%Bond Energy", energy);
      getKFData(&kf, "Geometry%atom order index", (void*) ordering);
      getKFData(&kf, "Properties%Dipole", (void*) dipxyz);

      //Get ordering data;
      int i,j,counter=0;
      for (i = 0; i < 2; i++)
          for (j = 0; j < *natoms; j++)
          internal_order[i][j]=ordering[counter++];
      i=0;
      //Re-order the gradients into gradients_ordered
      for (j = 0; j < *natoms*3; j+=3)
      {
         gradient[j]=gradient_unordered[3*(internal_order[0][i]-1)];
         gradient[j+1]=gradient_unordered[3*(internal_order[0][i]-1)+1];
         gradient[j+2]=gradient_unordered[3*(internal_order[0][i]-1)+2];
         i++;
      }

    }
    else //We are doing DFTB
    {
      getKFData(&kf, "finalresults%.dftb.gradient", (void*) gradient);
      getKFData(&kf, "finalresults%.dftb.energy", (void*) energy);
    }
   
    closeKFFile(&kf);
    return;
}

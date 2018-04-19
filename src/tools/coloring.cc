// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <string.h>
#include <cassert>


/* C version of greedy multicoloring algorithm adapted from SPARSKIT (LGPL) */
void greedyMC(int n, int *ja, int *ia, int *num_colors, int *iord, int *colors)
{
   int i, j, k, maxcol, mycol, ncols;
   
   memset(colors, -1, n*sizeof(int));
//for(i=0; i<n; i++)
//   colors[i] = -1;
   
   maxcol = n;
   ncols = 0;   
   int *tmp = new int[n];
   memset(tmp, 0, n*sizeof(int));
//for(i=0; i<n; i++)
//   tmp[i] = 0;
   
/* scan all nodes */   
   for(i=0; i<maxcol; i++)
   {
   /* look at adjacent nodes for previously assigned colors */
      int mcol = -1;
      int ii = iord[i];
      int k1 = ia[ii];
      int k2 = ia[ii+1];
      for(k=k1; k<k2; k++)
      {
         j = ja[k];
         int icol = colors[j];
         if(icol != -1)
         {
            mcol = mcol > icol ? mcol : icol;
            /* mark neighbor's color */
            tmp[icol] = -1; 
         }
      }

      /* scan marked colors to assign color */
      for(mycol=0; mycol<=mcol; mycol++)
         if(tmp[mycol] != -1) break;
      /* reset tmp array */
      memset(tmp, 0, (mcol+1)*sizeof(int));
     
      /* assign color and update num_colors */
      colors[ii] = mycol;
      ncols = ncols > mycol ? ncols : mycol;
   }
   /* done with coloring. Assign values */
   *num_colors = ncols+1;
   
   /* check result */
#ifndef NDEBUG
   for(i=0; i<n; i++)
   {
      int k1 = ia[i];
      int k2 = ia[i+1];
      for(k=k1; k<k2; k++)
      {
         j = ja[k];
         if(i == j) continue;
         assert(colors[i] != colors[j]);
      }
   }
#endif

   delete [] tmp;
   return;
} 


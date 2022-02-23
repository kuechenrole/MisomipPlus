/* 
 * fort_part_v2.c
 *
 */

/* 

 a Fortran interface to METIS-5.1 and simple block/range partitioners,
 adapted by Natalja Rakowsky from fort_part.c of FoSSI, 
  the Family of Simplified Solver Interfaces by
  Stephan Frickenhaus, Computing-Center, AWI-Bremerhaven, Germany
  www.awi-bremerhaven.de
 
*/

/* SF 06/2001 */
/* NR 07/2014 */
/* The partinioner Interface to METIS-5.1 */
/* n: dimension of system (number of equations=number of variables) */
/*    compressed matrix on input: */
/* ptr : rowptr or colptr */
/* adj : colind or rowind */
/* wgt: weights at nodes (FESOM/FVOM: number of levels for each 2D node) */
/* np  : partitioning for np processors */
/* part : partitioning vector */

#include "metis.h"
#include <stdio.h>




#ifdef UNDER_
void partit_block_
#endif
#ifdef UPCASE
void PARTIT_BLOCK
#endif
#ifdef MATCH
void partit_block
#endif
(int *n, int *ptr, int *adj, int *np, int *part)
{
int i ;

if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}
for (i=0;i<*n;i++) {
 part[i]=*np * i / *n; 
 if (part[i]==*np) part[i]=*np-1;
 }
}



#ifdef UNDER_
void partit3_
#endif
#ifdef UPCASE
void PARTIT3
#endif
#ifdef MATCH
void partit3
#endif
(int *n, int *ptr, int *adj, int *wgt, int *np, int *part)


{

int i ;
//idx_t *null;
idx_t opt[METIS_NOPTIONS];
idx_t ec, ncon;
idx_t *wgt_2d3d;


if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}

 if (sizeof(idx_t) != sizeof(int)){
   printf("Problem in partit: Size of int (%d) and Metis' idx_t (%d) differ!",sizeof(idx_t),sizeof(int));
   return;
 }


 METIS_SetDefaultOptions(opt);
 opt[METIS_OPTION_NUMBERING]=1; /* Fortran numbering */

printf("Partit %d %d %d %d %d %d\n",*n,*np,ptr[0],ptr[1],adj[0],adj[1]);


 ncon=2;
 wgt_2d3d = (idx_t*) malloc(2* *n *sizeof(idx_t));
 for (i=0; i<*n; i++){
   wgt_2d3d[2*i]   = 1;
   wgt_2d3d[2*i+1] = wgt[i]; 
 }

 METIS_PartGraphKway(n,&ncon,ptr,adj,wgt_2d3d,NULL,NULL,np,NULL,NULL,opt,&ec,part);


printf("Edgecut %d\n", ec);
 for(i=0;i<*n;i++) part[i]--;

 
 free(wgt_2d3d);
}




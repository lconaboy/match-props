#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define MINPART         50000

#define MAXSTRING       4096

#define pow2(x)         ((x)*(x))
#define pow3(x)         ((x)*(x)*(x))

// Structs
typedef struct HALO *HALOptr;
typedef struct HALO
{   
  long   npart;
  
  /* read_halos() will convert lengths into grid units [0,LGRID-1]! 
     
     LC 1/1/21 I don't think it does!
   */
  long   haloID;
  double Xc;
  double Yc;
  double Zc;
  double Vx;
  double Vy;
  double Vz;
  double Rvir;
  double Mvir;
  double sigV;
  
  // LC TODO we also want fmhires and gas/stellar props, or maybe just
  // output a list of IDs and extract the props later
  
  /* feel free to add whatever halo property you fancy... */
  /* ...but for the time being we are only concerned with position and radius */
  
}HALO;

// Globals
HALOptr halo;
long    nhalos_unc;
double  Xhalo, Yhalo, Zhalo, BoxSize;


void read_halos     (char *halofile);

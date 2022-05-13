/*==============================================================================

 *  FindHalo:   read *.AHF_halos file and finds a certain halo based upon user supplied (X,Y,Z)
 *
 *
 *  input:    - 1x *.AHF_halos file
 *            - (X,Y,Z)
 *
 *  output:   - ID (and position) of halo closest to (X,Y,Z)
 * 
 * LC 1/1/21 Ripped from ahf/analysis/ahfFindHalo.c
 */

#include "read_ahf.h"

//#define DEBUG


// Declarations
// void read_halos     (char *halofile);
// LC: do this in read_ahf.h

/*==============================================================================
 * read in initial halo positions from AHF_halos
 *==============================================================================*/
void read_halos(char *infile)
{
  FILE   *fpin;
  char    dummyline[MAXSTRING];
  long    ihalo, haloID, hostHalo, numSubStruct;
  long    nbins;
  long    nhalos_tot = 0;  // total number of halos in file
  // long    nhalos_unc = 0;  // total number of uncontaminated halos // global
  double  UNCONTAM = 0.9999;
  double  Xc, Yc, Zc, VXc, VYc, VZc, npart, Mvir, Rvir, Rmax, r2, mbp_offset;
  double  com_offset, Vmax, v_esc, sigV, lambda, lambdaE, Lx, Ly, Lz, b, c, Eax, Eay;
  double  Eaz, Ebx, Eby, Ebz, Ecx, Ecy, Ecz, ovdens, fMhires;
  double  Ekin, Epot, SurfP, Phi0, cNFW, n_gas, M_gas, lambda_gas, lambdaE_gas;
  double  Lx_gas, Ly_gas, Lz_gas, b_gas, c_gas, Eax_gas, Eay_gas, Eaz_gas;
  double  Ebx_gas, Eby_gas, Ebz_gas, Ecx_gas, Ecy_gas, Ecz_gas, Ekin_gas;
  double  Epot_gas, n_star, M_star, M;

  
  
  printf("-- reading _halos file %s\n",infile);
  
  /* open AHF_halos file */
  if((fpin=fopen(infile,"r"))==NULL)
   {
    printf("I cannot open %s\n", infile);
    exit(1);
   }
  
  
  /* overread header line */
  fgets(dummyline,MAXSTRING,fpin);
  
  /* how many halos are there? */
  ihalo = 0;
  while(fgets(dummyline,MAXSTRING,fpin) != NULL) ihalo++;
  // This is the total number of halos in that file
  nhalos_tot = ihalo;

  /* rewind file back to start and overread header */
  rewind(fpin);
  fgets(dummyline,MAXSTRING,fpin);

  // Reset global variabls
  nhalos_unc = 0;
  
  // Now count number of uncontaminated haloes
  for(ihalo=0; ihalo<nhalos_tot; ihalo++) {
     /* read info from file */
     fgets(dummyline,MAXSTRING,fpin);
    
     /* extract information from last read dummyline */
     //ss/canf(dummyline,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&npart,&nvpart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Mvir,&Rvir);
     sscanf(dummyline,"%ld %ld %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %ld %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir,&Rmax,&r2,&mbp_offset,&com_offset,&Vmax,&v_esc,&sigV,&lambda,&lambdaE,&Lx,&Ly,&Lz,&b,&c,&Eax,&Eay,&Eaz,&Ebx,&Eby,&Ebz,&Ecx,&Ecy,&Ecz,&ovdens,&nbins,&fMhires);

     // hostHalo is 0 or -1 if it is not a subhalo
     if (hostHalo < 1) {
       if (fMhires > UNCONTAM) nhalos_unc++;
     }
   }

  // Assign number of uncontaminated haloes to nhalos
  printf("---- %ld uncontaminated halos out of %ld\n", nhalos_unc, nhalos_tot);
  
  /* allocate memory for halos */
  halo = (HALOptr) calloc(nhalos_unc, sizeof(HALO));
  long _ihalo = 0;
  /* eventually read halos file */
  /* rewind file back to start and overread header */
  rewind(fpin);
  fgets(dummyline,MAXSTRING,fpin);

  for(ihalo=0; ihalo<nhalos_tot; ihalo++)
   {
     /* read info from file */
     fgets(dummyline,MAXSTRING,fpin);
    
     /* extract information from last read dummyline */
     //ss/canf(dummyline,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&npart,&nvpart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Mvir,&Rvir);
     sscanf(dummyline,"%ld %ld %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir,&Rmax,&r2,&mbp_offset,&com_offset,&Vmax,&v_esc,&sigV,&lambda,&lambdaE,&Lx,&Ly,&Lz,&b,&c,&Eax,&Eay,&Eaz,&Ebx,&Eby,&Ebz,&Ecx,&Ecy,&Ecz,&ovdens,&nbins,&fMhires,&Ekin,&Epot,&SurfP,&Phi0,&cNFW,&n_gas,&M_gas,&lambda_gas,&lambdaE_gas,&Lx_gas,&Ly_gas,&Lz_gas,&b_gas,&c_gas,&Eax_gas,&Eay_gas,&Eaz_gas,&Ebx_gas,&Eby_gas,&Ebz_gas,&Ecx_gas,&Ecy_gas,&Ecz_gas,&Ekin_gas,&Epot_gas,&n_star,&M_star);
     // sscanf(dummyline,"%ld %ld %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %ld %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir,&Rmax,&r2,&mbp_offset,&com_offset,&Vmax,&v_esc,&sigV,&lambda,&lambdaE,&Lx,&Ly,&Lz,&b,&c,&Eax,&Eay,&Eaz,&Ebx,&Eby,&Ebz,&Ecx,&Ecy,&Ecz,&ovdens,&nbins,&fMhires);

          // hostHalo is 0 or -1 if it is not a subhalo
     if (hostHalo < 1) {
       if (fMhires > UNCONTAM) {
	 M = Mvir - M_gas - M_star;  // DM mass

	 /* transfer info to halo structure */
	 halo[_ihalo].npart  = (long) npart;
	 halo[_ihalo].haloID = haloID;  // this is already long
	 halo[_ihalo].Xc     = Xc;
	 halo[_ihalo].Yc     = Yc;
	 halo[_ihalo].Zc     = Zc;
	 halo[_ihalo].Rvir   = Rvir;
	 halo[_ihalo].Mvir   = M; 
	 halo[_ihalo].sigV   = sigV; 
	 halo[_ihalo].Vx     = VXc;
	 halo[_ihalo].Vy     = VYc;
	 halo[_ihalo].Vz     = VZc;
	 _ihalo++;
       }
     }
   }  
  fclose(fpin);
  
  printf("-- done with %s\n", infile);
}

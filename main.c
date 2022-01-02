#include "read_ahf.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* TODO

     match<run><set>

     Working?
 */

long match(int iout, char* run1, char* run2);

long *match1, *match2;

int main(int argc, char *argv[]) {
  char *run1 = argv[1];
  char *run2 = argv[2];
  char *run3 = argv[3];
  char *p;
  long match_tot1, match_tot2;
  
  long _iout = strtol(argv[4], &p, 10);  // Convert the iout arg to long
  int iout = _iout;
  printf("%d\n", iout);
  
  match_tot1 = match(iout, run1, run2);
  long match11[match_tot1], match21[match_tot1];
  for (long i=0; i<match_tot1; i++) {
    match11[i] = match1[i];
    match21[i] = match2[i];
  }

  match_tot2 = match(iout, run1, run3);
  long match12[match_tot2], match32[match_tot2];
  for (long i=0; i<match_tot2; i++) {
    match12[i] = match1[i];
    match32[i] = match2[i];
  }


  // Do final swap so smallest is 1
  /*
  if (match_tot1 > match_tot2) {
    long _match_tot1;
    long _match11[match_tot2], _match21[match_tot2];

    for (long i=0; i<match_tot2; i++) {
      _match11[i] = match12[i];
      _match21[i] = match32[i];
    }
  }
  */
  
  // Now copmpare the two lists
  long match_tot_min = MIN(match_tot1, match_tot2);
  long match_tot_max = MAX(match_tot1, match_tot2);
  long match_tot = 0;
  long* match1_f;
  long* match2_f;
  long* match3_f;
  
  int bcomp[match_tot_min];
  long icomp[match_tot_min];
  
  if (match_tot1 < match_tot2) {
    // This branch means there are fewer elements in set 1 (i.e. run1 and run2), so we loop over match_tot1
    for (long i=0; i<match_tot1; i++) {
      bcomp[i] = 0;
      icomp[i] = -1;
      for (long j=0; j<match_tot2; j++) {
	if (match11[i] == match12[j]) {
	  bcomp[i] = 1;
	  icomp[i] = j;
	  match_tot++;
	  break; // No need to check further
	}
      }
    }
    // long match1_f[match_tot], match2_f[match_tot], match3_f[match_tot];
    match1_f = (long*) calloc(match_tot, sizeof(long));
    match2_f = (long*) calloc(match_tot, sizeof(long));
    match3_f = (long*) calloc(match_tot, sizeof(long));
    long j=0;
    // Now final iteration
    for (long i=0; i<match_tot_min; i++) {
      if (bcomp[i]) {
	match1_f[j] = match11[i];
	match2_f[j] = match21[i];
	match3_f[j] = match32[icomp[i]];
	j++;
      }
    }
    assert(j == match_tot);

  } else {
    // This branch means there are fewer elements in the second set of
    // runs tested (i.e. run1 and run3, so we loop over match_tot2
    // first)
    for (long i=0; i<match_tot2; i++) {
      bcomp[i] = 0;
      icomp[i] = -1;
      for (long j=0; j<match_tot1; j++) {
	if (match12[i] == match11[j]) {
	  bcomp[i] = 1;
	  icomp[i] = j;
	  match_tot++;
	  break; // No need to check further
	}
      }
    }
    //long match1_f[match_tot], match2_f[match_tot], match3_f[match_tot];
    match1_f = (long*) calloc(match_tot, sizeof(long));
    match2_f = (long*) calloc(match_tot, sizeof(long));
    match3_f = (long*) calloc(match_tot, sizeof(long));

    long j=0;
    // Now final iteration
    for (long i=0; i<match_tot_min; i++) {
      if (bcomp[i]) {
	match1_f[j] = match12[i];
	match3_f[j] = match32[i];
	match2_f[j] = match21[icomp[i]];
	j++;
      }
    }
    assert(j == match_tot);
  }


  FILE *fp;
  char fnout[128];
  sprintf(fnout, "match_%s_%s_%s_%03d.txt", run1, run2, run3, iout);
  fp = fopen(fnout, "w");
  for (int i=0; i<match_tot; i++) {
    fprintf(fp, "%ld %ld %ld\n", match1_f[i], match2_f[i], match3_f[i]);	
  } 
  fclose(fp);

  
  // match1, match3 = match(iout, run1, run3);
}

long match(int iout, char* run1, char* run2){
  HALOptr halos1, halos2;
  long nhalos1, nhalos2;
  // double dx2, dv2, dm, sigx2, sigv2, sigm;
  double sigx2, sigv2, sigm;
  double dx2_tot, dv2_tot, dm_tot;
  // int iout;
  int bswap = 0;
  

  printf("-- comparing %s and %s for output %d\n", run1, run2, iout);
  
  // Read first set of haloes
  char fn1[128];
  sprintf(fn1, "%s/all_%03d.AHF_halos", run1, iout);
  read_halos(fn1);
  halos1 = halo;
  nhalos1 = nhalos_unc;
  printf("---- %ld haloes read for halos1\n", nhalos1);

  // Read second set of haloes
  char fn2[128];
  sprintf(fn2,"%s/all_%03d.AHF_halos", run2, iout);
  read_halos(fn2);
  halos2 = halo;
  nhalos2 = nhalos_unc;
  printf("---- %ld haloes read for halos2\n", nhalos2);

  // Make sure nhalos1 is longest
  if (nhalos1 > nhalos2) {
    long _nhalos1;
    HALOptr _halos1;
    
    printf("-- swapping halos1 and halos2 so halos2 is largest\n");

    _nhalos1 = nhalos2;
    _halos1 = halos2;
    nhalos2 = nhalos1;
    halos2 = halos1;
    nhalos1 = _nhalos1;
    halos1 = _halos1;
    /* run1 = argv[2]; */
    /* run2 = argv[1]; */
    bswap = 1;
  } else {
    printf("-- halos1 already smallest, not swapping\n");
  }

  // TODO make ji an array
  long ji[nhalos1];
  double dx2[nhalos1], dv2[nhalos1], dm[nhalos1];
  
  // Loop over first list
  for (long i=0; i<nhalos1; i++) {
    double dx2_tot=0, dv2_tot=0, dm_tot=0;
    // Reset comparison variables
    dx2[i] = 100.;  // 10 * Rvir
    dv2[i] = 25.;  // 5 * |V|
    dm[i] = 2.;  // 2 * dMvir
    ji[i] = -1;

    // Normalisation vars
    sigx2 = pow2(halos1[i].Rvir);
    // sigv2 = pow2(halos1[i].Vx) + pow2(halos1[i].Vy) + pow2(halos1[i].Vz)
    sigv2 = pow2(halos1[i].sigV);
    // sigv2 = pow2(halos1[i].Vx) + pow2(halos1[i].Vy) + pow2(halos1[i].Vz);
    sigm = halos1[i].Mvir;
    
    // Now loop over the second list
    for (long j=0; j<nhalos2; j++) {
      double _dx2=0, _dv2=0, _dm;

      _dx2 += pow2(halos2[j].Xc - halos1[i].Xc);
      _dx2 += pow2(halos2[j].Yc - halos1[i].Yc);
      _dx2 += pow2(halos2[j].Zc - halos1[i].Zc);
      _dx2 /= sigx2;

      _dv2 += pow2(halos2[j].Vx - halos1[i].Vx);
      _dv2 += pow2(halos2[j].Vy - halos1[i].Vy);
      _dv2 += pow2(halos2[j].Vz - halos1[i].Vz);
      _dv2 /= sigv2;

      _dm = fabs(halos2[j].Mvir - halos1[i].Mvir) / sigm;

      dx2_tot += _dx2;
      dv2_tot += _dv2;
      dm_tot += _dm;
      
      if (_dm < dm[i]) {
	// printf("closer in mass %lf %lf\n", dm, _dm);
	dm[i] = _dm;
	if (_dx2 < dx2[i]) {
	  // printf("closer in pos %lf %lf\n", dx2, _dx2);
	  dx2[i] = _dx2;
	  if (_dv2 < dv2[i]) {
	    // printf("closer in vel %lf %lf\n", dv2, _dv2);
	    // printf("closer than all \n");
	    dv2[i] = _dv2;
	    ji[i] = j; // If this is a closer match, then store it	    
	  }
	}
      }
      
    }
    
    /* if (ji > -1) { */
    /*   printf("MATCH i %ld ji %ld M1 %lf M2 %lf\n", i, ji, halos1[i].Mvir, halos2[ji].Mvir); */
    /*   printf("min dx2 %lf dv2 %lf dm %lf\n", dx2, dv2, dm); */
    /*   printf("mean dx2 %lf dv2 %lf dm %lf\n", dx2_tot/nhalos2, dv2_tot/nhalos2, dm_tot/nhalos2); */
    /* } */
  }

  // NB this was me playing with pointers when I don't properly
  // understand them aka it's all crap
  // Now check for duplicates in ji
  /* int *ip, *jp; */
  /* for (ip = &ji[0]; ip < &ji[nhalos2]; ip++){ */
  /*   int match=0; */
    
  /*   for (jp = &ji[0]; jp < &ji[nhalos2]; jp++) { */
  /*     if (ip == jp) match++; */
  /*     printf("%ld %ld\n", ip, jp); */
  /*   } */
  /*   // TODO might have to fix this if there are multiple matches */
  /*   if (match > 1) printf("%ln has %d matches\n", jp, match); */
  /* } */

  // Now check for duplicates in ji
  long match_tot=0;
  long dup[nhalos1];
  for (long i=0; i<nhalos1; i++) {
    int match=0;
    for (long j=0; j<nhalos1; j++) {
      dup[j] = -1;
      if (i == ji[j]) {
	match++;
	dup[j] = 1;
      }
    }
  
    // TODO might have to fix this if there are multiple matches
    if (match > 0) match_tot++;
    if (match > 1) {
      printf("---- %ld has %d matches, removing\n", i, match);
      for (long j=0; j<nhalos1; j++) {
	if (dup[j] > -1) {
	  ji[j] = -1;
	}
      }
      match_tot--;
    }
  }

  printf("-- total matches after cleaning %d\n", match_tot);
  

  // Now extract matching IDs
  // long match1[match_tot], match2[match_tot], imatch=0;
  long imatch=0;
  match1 = (long*) calloc(match_tot, sizeof(long));
  match2 = (long*) calloc(match_tot, sizeof(long));
  for (long i=0; i<nhalos1; i++) {
    if (ji[i] > -1) {
      match1[imatch] = halos1[i].haloID;
      match2[imatch] = halos2[ji[i]].haloID;
      imatch++;
    }
  }

  printf("%ld %ld\n", imatch, match_tot);
  assert(imatch == match_tot);

  // We might've swapped halos1 and halos2, if so swap back
  if (bswap) {
    long _match1[match_tot];
    printf("-- swapping halos1 and halos2 back\n");
    for (long i=0; i<match_tot; i++) {
      _match1[i] = match2[i];
      match2[i] = match1[i];
      match1[i] = _match1[i];
    }
  }
  
  FILE *fp;
  char fnout[128];
  sprintf(fnout, "match_%s_%s_%03d.txt", run1, run2, iout);
  fp = fopen(fnout, "w");
  for (int i=0; i<match_tot; i++) {
    fprintf(fp, "%ld %ld\n", match1[i], match2[i]);	
  } 
  fclose(fp);

  return match_tot;
}

  



  
  

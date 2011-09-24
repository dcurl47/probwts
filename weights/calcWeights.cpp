#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "kdtrees.h"



using namespace std;


void find_NR_NMAG(const char *filename, int *NR, int *NMAG);
void read_data(const char *filename, int NR, int NMAG, double *zp, double *zs, double *sed, double *mags);
void find_NR_NMAG2(const char *filename, int *NR, int *NMAG);

int main(int argc, char **argv) {
  
  
  if (argc < 4) {
    fprintf(stderr, "Usage: %s train_tbl photo_tbl NNEAR\n", argv[0]);
    exit(1);
  }
 
  
  
  int NNEAR = 5;
  if (argc > 3) NNEAR = atoi(argv[3]);
  
  int tNR, tNMAG;
  int pNR, pNMAG;
  
  find_NR_NMAG(argv[1], &tNR, &tNMAG);
  find_NR_NMAG(argv[2], &pNR, &pNMAG);

 
  vector<string> tid(tNR);
  double *tzp = (double*)malloc(sizeof(double)*tNR);
  double *tzs = (double*)malloc(sizeof(double)*tNR);
  double *tze = (double*)malloc(sizeof(double)*tNR);
  double *tzw = (double*)malloc(sizeof(double)*tNR);
  double *tra = (double*)malloc(sizeof(double)*tNR);
  double *tde = (double*)malloc(sizeof(double)*tNR);
  double *tm = (double*)malloc(sizeof(double)*tNR*tNMAG);
  double *te = (double*)malloc(sizeof(double)*tNR*tNMAG);
  
  vector<string> pid(pNR);
  double *pzp = (double*)malloc(sizeof(double)*pNR);
  double *pzs = (double*)malloc(sizeof(double)*pNR);
  double *pze = (double*)malloc(sizeof(double)*pNR);
  double *pzw = (double*)malloc(sizeof(double)*pNR);
  double *pra = (double*)malloc(sizeof(double)*pNR);
  double *pde = (double*)malloc(sizeof(double)*pNR);
  double *pm = (double*)malloc(sizeof(double)*pNR*pNMAG);
  double *pe = (double*)malloc(sizeof(double)*pNR*pNMAG);
  
  //============================================================
  int *pnum = (int*)malloc(sizeof(double)*pNR);
  int *prank = (int*)malloc(sizeof(double)*pNR);
  pnum=0;
  prank=0;
  //============================================================
  
  if (tNMAG != pNMAG) {
    fprintf(stderr, "number of magnitudes does not match\n");
    exit(3);
  }
  
  read_data(argv[1], tNR, tNMAG, tzp, tzs, tzw, tm);
  read_data(argv[2], pNR, pNMAG, pzp, pzs, pzw, pm);
  
  vector<double> dt(tNR);
  vector<double> wt(tNR);
  vector<double> dp(pNR);
  vector<double> dact(tNR);
  vector<double> ds(NNEAR);
  vector<int> indx(NNEAR);
  vector<double> zptemp(NNEAR);
  vector<double> zstemp(NNEAR);
  vector<double> dd(NNEAR);
  int t10 = pNR / 10;
 
  double totalWeights = 0.0;
  
  //=============================================================================
  /////////////////////////////////////
  //Construct KD-tree for training set/
  /////////////////////////////////////
 
 
  const int dim=5;

  //Copy tm into vector of points:
	
  vector< Point<dim> > mptt(tNR);  
  for (int i=0;i<tNR;i++){
    vector<double> mmag(dim);
    for(int j=0;j<dim;j++){
      mmag[j]=tm[j+dim*i];
    }
    Point<dim> pt;

    for(int j=0;j<dim;j++){ 
      pt.x[j]=mmag[j];
    }
    mptt[i]=pt;
  }
 

//Construct kd-tree
  kdtree<dim> kdt(mptt);
 
 //===========================================================================


 ////////////////////////////////////////
  //Construct KD-tree for photometric set/
  ////////////////////////////////////////
 //Copy pm into vector of points:
  
  vector< Point<dim> > mptp(pNR); 
  
  for (int i=0;i<pNR;i++){
    vector<double> mmag(dim);
    for(int j=0;j<dim;j++){
      mmag[j]=pm[j+dim*i];
    }
    Point<dim> pt;
  
    for(int j=0;j<dim;j++){ 
      pt.x[j]=mmag[j];
    }
    mptp[i]=pt;
  }
 
	
//Construct kd-tree
  kdtree<dim> kdp(mptp);
  


  //==========================================================================
  /////////////////////////////////////////////////////////
  //Find Nth nearest neighbor (nnt)for each object (Ot) in/ 
  //training set and its distance d to Ot                 /
  /////////////////////////////////////////////////////////
  const int nmax =100000;
  int *nnt=new int[NNEAR];
  double *ndt=new double[NNEAR];
  int *listp=new int[nmax];
  int *nnp;
  double *ndp;
  vector<double> nNNEAR(tNR);
  vector<double> distNNEAR(tNR);
   vector<int> num(pNR);

  for (int i=0;i<pNR;i++) num[i]=0;
  for (int i=0;i<nmax;i++) listp[i]=0;

  vector<double> nNN(tNR);
  vector<double> distNN(tNR);
  vector<int> nret(tNR);

//======================================================
//timing code
//======================================================
  timeval *t1 = new timeval();
  timeval *t2 = new timeval();

	
  gettimeofday(t1,NULL);

    for(int i=0;i<tNR;i++){

    if (i % 1000 == 0) {
      gettimeofday(t2, NULL);
      double ms = 1000.*(t2->tv_sec-t1->tv_sec) + 0.001*(t2->tv_usec-t1->tv_usec);
      ms = ms / 1000.0;
      printf("i = %d, time per obj = %lf ms\n", i, ms);
      fflush(NULL);
      gettimeofday(t1, NULL);
    }

//=========================================================

	kdt.nneigh(NNEAR,ndt,nnt,i);
	
    nNN[i]=nnt[0];
    distNN[i]=ndt[0];

    Point<dim> pt;
    for(int j=0;j<dim;j++){ 
      pt.x[j]=mptt[i].x[j];
    }

 ///////////////////////////////////////////////////////////////////
 //Find all photometric set objects within distance distNN[i] of pt/
 ///////////////////////////////////////////////////////////////////


    nret[i]= kdp.pointsr(distNN[i],mptt[i],listp);
    //nret[i] Contains number of neighbors within distNN[i] of the pt.
    for(int j=0;j<nret[i];j++){
    

      num[listp[j]]++; //counts how many times each photometric set object
                       //was used.

    }
 
    wt[i] = (double)nret[i] / ((double)NNEAR);
    totalWeights += wt[i];
    
    }
//=========================================================
//Output stuff
//=========================================================


  FILE *outfile = fopen("nnweight.prw", "w");
  
  for(int i=0; i<tNR; ++i) {
    fprintf(outfile, "%10.7f %10.7f %16.12e ",	    
	    tzs[i], tzp[i], wt[i]/totalWeights);
    for(int j=0; j<tNMAG; ++j) {
      fprintf(outfile, "%10.7f ", tm[i*tNMAG+j]);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);


 
  //=============================================================================
 
  FILE *outfile2 = fopen("phot.num", "w");
  
  for(int i=0; i<pNR; ++i) {
    fprintf(outfile2, "%10.7f %10.7f %d ",pzs[i], pzp[i], num[i]);
    fprintf(outfile2, "\n");
  }
  fclose(outfile2);

 
  return 0;
}




void find_NR_NMAG(const char *filename, int *NR, int *NMAG) {
  //This function is used to find out the number of rows and number of magnitudes in a table file
  //Find out the number of rows
  FILE *in = fopen(filename, "r");
  char line[256];
  char *token = NULL;
  fgets(line, 256, in);
  token = strtok(line, " \n");
  *NMAG = 0;
  int NS = 3;
  while (token != NULL) {
    (*NMAG)++;
    token = strtok(NULL, " \n");
  }
  *NMAG = (*NMAG - NS);
  //Find out the number of rows
  *NR = 0;
  while (!feof(in)) {
    (*NR)++;
    fgets(line, 256, in);
  }
  printf("NR = %d, NMAG = %d\n", *NR, *NMAG);
}

void read_data(const char *filename, int NR, int NMAG, double *zp, double *zs, double *sed, double *mags) {
  FILE *in = fopen(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filename);
    exit(1);
  }

  int i, j;
  for(i=0; i<NR; ++i) {
    fscanf(in, "%lf %lf %lf", &zs[i], &zp[i], &sed[i]);
    for(j=0; j<NMAG; ++j) {
      fscanf(in, "%lf", &mags[i*NMAG + j]);
      // mags[i*NMAG+j] = (mags[i*NMAG+j]-10.0)/22.0;
    }
  }
  fclose(in);
}

void find_NR_NMAG2(const char *filename, int *NR, int *NMAG) {
  //This function is used to find out the number of rows and number of magnitudes in the input photometric set
  //The column format for this input should be:
  // 1: zphot
  // 2-5: g, r, i, z

  //Find out the number of rows
  FILE *in = fopen(filename, "r");
  char line[256];
  char *token = NULL;
  fgets(line, 256, in);
  token = strtok(line, " ");
  *NMAG = 0;
  int NS = 1;
  while (token != NULL) {
    (*NMAG)++;
    token = strtok(NULL, " \n");
  }
  *NMAG = *NMAG - NS;
  //Find out the number of rows
  *NR = 0;
  while (!feof(in)) {
    (*NR)++;
    fgets(line, 256, in);
  }
  printf("NR = %d, NMAG = %d\n", *NR, *NMAG);
}

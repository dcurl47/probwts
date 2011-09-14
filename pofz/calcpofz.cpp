#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/time.h>


#include "kdtrees.h"

void sort2(vector<double> &arr, vector<double> &brr);
void find_NR_NMAG(const char *filename, int *NR, int *NMAG);
void read_data(const char *filename, int NR, int NMAG, double *zs, double *zp, double *sed, double *mags);
int smuthistwei(vector<double> &hx,vector<double> &hy, vector<double> &rhos, vector<double> &wei, int res);
void swap(double& a,double& b);
 double rmin=0.;
 double rmax=1.1;



int main(int argc, char **argv) {

  if (argc < 6) {
    fprintf(stderr, "Usage: %s train_tbl photo_tbl NNEAR grid res\n", argv[0]);
    exit(1);
  }

  int NNEAR = atoi(argv[3]);
  int grid=atoi(argv[4]);
  vector<double> hx(grid);
  vector<double> Nz1(grid);
  
  if (argc>6) rmin=atof(argv[6]);
  if (argc>7) rmax=atof(argv[7]);
  cout<<"rmin = "<<rmin<<"  rmax = "<<rmax<<endl;

  int res = atoi(argv[5]);
  int tNR, tNMAG;
  int pNR, pNMAG;
  FILE *out = fopen("output.etbl", "w");
  if (out == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", "output.etbl");
    exit(1);
  }
 
  FILE *out3 = fopen("output.grid", "w");
  if (out3 == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", "output.etbl");
    exit(1);
  }
  find_NR_NMAG(argv[1], &tNR, &tNMAG);
  find_NR_NMAG(argv[2], &pNR, &pNMAG);

  double *tzp1 = (double*)malloc(sizeof(double)*tNR);
  double *tzs1 = (double*)malloc(sizeof(double)*tNR);
  double *twei = (double*)malloc(sizeof(double)*tNR);

  double *tm = (double*)malloc(sizeof(double)*tNR*tNMAG);

  double *uzp1 = (double*)malloc(sizeof(double)*tNR);
  double *uzs1 = (double*)malloc(sizeof(double)*tNR);
  double *uwei = (double*)malloc(sizeof(double)*tNR);

  double *um = (double*)malloc(sizeof(double)*tNR*tNMAG);

  double *pzp1 = (double*)malloc(sizeof(double)*pNR);
  double *pzs1 = (double*)malloc(sizeof(double)*pNR);
  double *pwei = (double*)malloc(sizeof(double)*pNR);
  double *pm = (double*)malloc(sizeof(double)*pNR*pNMAG);

  if (tNMAG != pNMAG) {
    fprintf(stderr, "number of magnitudes does not match\n");
    exit(3);
  }

  read_data(argv[1], tNR, tNMAG, tzs1, tzp1, twei, tm);
  read_data(argv[2], pNR, pNMAG, pzs1, pzp1, pwei, pm);
  
  vector<double> ds(NNEAR);
  vector<int> indx(NNEAR);
  vector<double> zp1temp(NNEAR);
  vector<double> zs1temp(NNEAR);
  vector<double> weitemp1(NNEAR);
  
 
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
    pt.x[0]=mmag[0];
    pt.x[1]=mmag[1];
    pt.x[2]=mmag[2];
    pt.x[3]=mmag[3];
    if (dim>4) pt.x[4]=mmag[4];
    if (dim>5) pt.x[5]=mmag[5];
    if (dim>6) pt.x[6]=mmag[6];
    if (dim>7) pt.x[7]=mmag[7];
    if (dim>8) pt.x[8]=mmag[8];
    if (dim>9) pt.x[9]=mmag[9];


   mptt[i]=pt;
  }
  

//Construct kd-tree
  kdtree<dim> kdt(mptt);
 
 //===========================================================================

 //Copy pm into vector of points:
  vector< Point<dim> > mptp(pNR);  
  for (int i=0;i<pNR;i++){
    vector<double> mmag(dim);
    for(int j=0;j<dim;j++){
      mmag[j]=pm[j+dim*i];
    }
    Point<dim> pt;
    pt.x[0]=mmag[0];
    pt.x[1]=mmag[1];
    pt.x[2]=mmag[2];
    pt.x[3]=mmag[3];
    if (dim>4) pt.x[4]=mmag[4];
    if (dim>5) pt.x[5]=mmag[5];
    if (dim>6) pt.x[6]=mmag[6];
    if (dim>7) pt.x[7]=mmag[7];
    if (dim>8) pt.x[8]=mmag[8];
    if (dim>9) pt.x[9]=mmag[9];


    mptp[i]=pt;
  }
 

  /////////////////////////////////////////////////////////
  //Find Nth nearest training set neighbors (nnt)for      / 
  //each object (Ot) in photometric set                   /
  /////////////////////////////////////////////////////////
  const int nmax =100000;
  int *nnt=new int[NNEAR];
  double *ndt=new double[NNEAR];
  int *listp=new int[nmax];
  int *nnp;//=new int[NNEAR];
  double *ndp;//=new double[NNEAR];
  vector<double> nNNEAR(tNR);
  vector<double> distNNEAR(tNR);
  for (int i=0;i<nmax;i++) listp[i]=0;

  vector<double> nNN(tNR);
  vector<double> distNN(tNR);
  vector<int> nret(tNR);


  int t10 = pNR / 10;
  timeval *t1 = new timeval();
  timeval *t2 = new timeval();
  gettimeofday(t1, NULL);

 

  for(int i=0; i<pNR; ++i) {
		//cout << i << endl;   
	  if (i % t10 == 0) {
      gettimeofday(t2, NULL);
      printf(" %4.2f, time taken per step = %4.2f\n", (float)(i/t10)*10.0, (1000.*(t2->tv_sec-t1->tv_sec) + 0.001*(t2->tv_usec - t1->tv_usec))/t10);
      
	  gettimeofday(t1, NULL);
      fflush(NULL);
    }

	 
  kdt.nneigh(NNEAR,ndt,nnt,mptp[i]);
	
  for(int j=0;j<NNEAR;j++){
    zs1temp[j]=tzs1[nnt[j]];
    zp1temp[j]=tzp1[nnt[j]];
    weitemp1[j]=twei[nnt[j]];
  }
  


  smuthistwei(hx,Nz1,zs1temp,weitemp1,res);
  fprintf(out, "%lf %lf  ", pzs1[i], pzp1[i]);
   

  for(int j=0; j<grid; ++j) {
    if (Nz1[j] < 1e-20)
      fprintf(out, " %1.0lf", Nz1[j]);
    else
      fprintf(out, " %lf", Nz1[j]);
    
  }
  fprintf(out, "\n");
  
  }
  
  for(int j=0; j<grid; ++j) {
    fprintf(out3,"%lf \n", hx[j]);
  }
  
  
  fclose(out);
  fclose(out3);
  
  
  
  return 0;
}



void sort2(vector<double> &arr, vector<double> &brr)
{
	const int M=7,NSTACK=50;
	int i,ir,j,k,jstack=-1,l=0;
	double a,b;
	vector<int> istack(NSTACK);

	int n=arr.size();
	ir=n-1;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (jstack < 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			swap(arr[k],arr[l+1]);
			swap(brr[k],brr[l+1]);
			if (arr[l] > arr[ir]) {
				swap(arr[l],arr[ir]);
				swap(brr[l],brr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
				swap(arr[l+1],arr[ir]);
				swap(brr[l+1],brr[ir]);
			}
			if (arr[l] > arr[l+1]) {
				swap(arr[l],arr[l+1]);
				swap(brr[l],brr[l+1]);
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			b=brr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				swap(arr[i],arr[j]);
				swap(brr[i],brr[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			brr[l+1]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack >= NSTACK){ printf("NSTACK too small in sort2.");
			exit(1);
			}
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
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
  *NMAG = *NMAG - NS;
  //Find out the number of rows
  *NR = 0;
  while (!feof(in)) {
    (*NR)++;
    fgets(line, 256, in);
  }
  printf("NR = %d, NMAG = %d\n", *NR, *NMAG);
}


void read_data(const char *filename, int NR, int NMAG, double *zs, double *zp, double *sed, double *mags) {
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

int smuthistwei(vector<double> &hx,vector<double> &hy, vector<double> &rhos, vector<double> &wei, int res) { 
  //first, find the min and the max of the data                 
  //int option = 1;                            
  int size=rhos.size();
  int ndiv=hy.size();
  //res=size/res;

  double hfwth, rho_min=rmin, rho_max=rmax; 
  vector<double> xmin(ndiv),xmax(ndiv);  
  double count=0;                                               
  double drho = (rho_max - rho_min)/(double)res;              
  double step=(rho_max - rho_min - drho)/((double)ndiv-1.);
  for (int i = 0; i < ndiv; ++i) {                             
    xmin[i] = rho_min + step*(double)i;                    
    xmax[i] = xmin[i] + drho;                
    hy[i] = 0.;  
    hx[i]=(xmin[i]+xmax[i])/2.;     
  }  
  double sumwei=0.;
  for (int i = 0; i < size; ++i) { 
    sumwei+=wei[i];
  }


  for (int j = 0; j < ndiv; ++j) {                             
    if (j == ndiv-1){                                          
      xmax[j]=xmax[j]+0.0001;                          
    }                                                           
    for (int i = 0; i < size; ++i) {                               
      if ((rhos[i] >= xmin[j]) && (rhos[i] < xmax[j])){ 
	if (sumwei != 0){
	hy[j]+=wei[i]/sumwei;                              
	} else  hy[j] =0.;
      } 
    }                                                           
    count +=hy[j];                                 
  }                                                           
                                                                 
  return 0;                                                     
}                         



double calcSigma(vector<double> &zs, vector<double> &zp) {
  if (zs.size() < 2) {
    return 0.0;
  }
  double r = 0.0;
  for(int i=0; i<zs.size(); ++i) {
    r += (zs[i]-zp[i])*(zs[i]-zp[i]);
  }
  return sqrt(r/zs.size());
}

double calcSigma68(vector<double> &zs, vector<double> &zp, vector<double> &wei) {
  if (zs.size() < 2) {
    return 0.0;
  }
  double ws=0.,ws1=0.;
  vector<double> d(zs.size()), wfrac(zs.size());
  for(int i=0; i<zs.size(); ++i) {
    d[i] = fabs(zs[i]-zp[i]);
    ws+=wei[i];
  }
  sort2(d,wei);
  int s68=0,k=0;
  for(int i=0; i<zs.size(); ++i) {
    ws1+=wei[i];
    //wfrac[i]=ws1/ws;
    k=i;
    if ((ws1/ws) > 0.68) break;
  }
  s68=k;

  return d[s68];
}
void swap(double& a,double& b) //swaps values of a and b
{
	double temp = a;
	a = b;
	b = temp;
} 

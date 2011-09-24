#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/time.h>


#include "kdtrees.h"

void sort2(vector<double> &arr, vector<double> &brr);
void find_NR_NMAG(const char *filename, int *NR, int *NMAG);
void read_data(const char *filename, int NR, int NMAG, double *zs, double *zp, double *sed, double *mags);
//void read_data(const char *filename, int NR, int NMAG, double *zs, double *weight, double *ra, double *dec, int *id, double *mags);
int smuthistwei(vector<double> &hx,vector<double> &hy, vector<double> &rhos, vector<double> &wei, int res);
void swap(double& a,double& b);
void gaussj(vector<vector<double> > &a, vector<vector<double> > &b);
void covsrt(vector<vector<double> > &covar, vector<bool> &ia, const int mfit);
void lfit(vector<vector<double> > &x, vector<double> &y, vector<double> &sig, vector<double> &a,vector<bool> &ia, vector<vector<double> > &covar, double &chisq, void funcs(vector<double>, vector<double> &));
void func(vector<double > var, vector<double> &a);
double evalfit(vector<double>  var, vector<double> &afunc, vector<double> &a);
int factorial(int number);
 double rmin=0.;
 double rmax=1.35;

int order=1;

int main(int argc, char **argv) {

  if (argc < 4) {
    fprintf(stderr, "Usage: %s train_tbl photo_tbl NNEAR \n", argv[0]);
    exit(1);
  }
 
  int NNEAR = 100;
  if (argc>3) NNEAR=atoi(argv[3]);

  if (argc>4) order=atoi(argv[4]);
  if (argc>5) rmin=atof(argv[4]);
  if (argc>6) rmax=atof(argv[5]);
//   cout<<"rmin = "<<rmin<<"  rmax = "<<rmax<<endl;

  //int res = atoi(argv[5]);
  int tNR, tNMAG;
  int pNR, pNMAG;
  FILE *out = fopen("output.etbl", "w");
  if (out == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", "output.etbl");
    exit(1);
  }
 

  find_NR_NMAG(argv[1], &tNR, &tNMAG);
  find_NR_NMAG(argv[2], &pNR, &pNMAG);
  int  *tid=(int*)malloc(sizeof(double)*tNR);
  double *tzp1 = (double*)malloc(sizeof(double)*tNR);
  double *tzs1 = (double*)malloc(sizeof(double)*tNR);
  double *twei = (double*)malloc(sizeof(double)*tNR);
  double *tm = (double*)malloc(sizeof(double)*tNR*tNMAG);
  double *tra = (double*)malloc(sizeof(double)*tNR);
  double *tdec = (double*)malloc(sizeof(double)*tNR);


  int  *pid=(int*)malloc(sizeof(double)*pNR);
  double *pzp1 = (double*)malloc(sizeof(double)*pNR);
  double *pzs1 = (double*)malloc(sizeof(double)*pNR);
  double *pwei = (double*)malloc(sizeof(double)*pNR);
  double *pm = (double*)malloc(sizeof(double)*pNR*pNMAG);
  double *pra = (double*)malloc(sizeof(double)*pNR);
  double *pdec = (double*)malloc(sizeof(double)*pNR);

  //dmatrix **blah= blah(NNEAR,NNEAR);
  vector<vector<double> > blah ( 6, vector<double> ( 8 ) );
  


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
  vector<vector<double> > mtemp (NNEAR, vector <double> (tNMAG));
  vector<double> sig(NNEAR);
  int ncoeffs;
  ncoeffs=factorial(tNMAG+order)/(factorial(tNMAG)*factorial(order));

  if (order<1 || order >3) {
    cout<<"The order of the polynomial you speficied is not allowed!!!"<<endl;
    exit(1);
  }
  
  vector<double> a(ncoeffs);
  double chisq;
  vector<bool>ia(ncoeffs);
  vector<vector<double> > covar (ncoeffs, vector <double> (ncoeffs));
  for (int i=0;i<ncoeffs;i++){
  ia[i]=1;
  }
  for (int i=0;i<NNEAR;i++){
  sig[i]=1.;
  }

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

    for(int k=0;k<tNMAG;k++){
      mtemp[j][k]=tm[k+nnt[j]*tNMAG];
    }
  }


  //==============================================
  //Fitting part
  //==============================================
  //find best-fitting coefficients

  lfit(mtemp,zs1temp,sig,a,ia,covar,chisq,func);
  //evaluate fit at the magnitudes of the photometric set object
  vector<double> afunc(ncoeffs);
  vector<double> mmagtemp(dim);
  for(int j=0;j<dim;j++){
    mmagtemp[j]=pm[j+dim*i];
  }
  double yyy=evalfit(mmagtemp,afunc,a);


 
  fprintf(out, "%lf %lf  \n", pzs1[i], yyy);
   
  }
  
  
  fclose(out);
   
  
  
  return 0;
  }



void sort2(vector<double> &arr, vector<double> &brr){
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
 
// void read_data(const char *filename, int NR, int NMAG, double *zs, double *weight, double *ra, double *dec, int *id, double *mags) {
//   FILE *in = fopen(filename, "r");
//   if (in == NULL) {
//     fprintf(stderr, "Error: failed to open file %s\n", filename);
//     exit(1);
//   }

//   int i, j;
//   for(i=0; i<NR; ++i) {
//     fscanf(in, "%lf %lf %d %lf %lf", &ra[i],&dec[i],&id[i],&zs[i], &weight[i]);
//     for(j=0; j<NMAG; ++j) {
//       fscanf(in, "%lf", &mags[i*NMAG + j]);
//       // mags[i*NMAG+j] = (mags[i*NMAG+j]-10.0)/22.0;
//     }
 
//   }
//   fclose(in);
// }

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
void covsrt(vector<vector<double> > &covar, vector<bool> &ia, const int mfit)
{
	int i,j,k;

	int ma=ia.size();
	for (i=mfit;i<ma;i++)
		for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit-1;
	for (j=ma-1;j>=0;j--) {
		if (ia[j]) {
			for (i=0;i<ma;i++) swap(covar[i][k],covar[i][j]);
			for (i=0;i<ma;i++) swap(covar[k][i],covar[j][i]);
			k--;
		}
	}
}
void gaussj(vector<vector<double> > &a, vector<vector<double> > &b)
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;

	int n=a.size();
	int m=b[0].size();
	vector<int> indxc(n),indxr(n),ipiv(n);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) swap(a[irow][l],a[icol][l]);
			for (l=0;l<m;l++) swap(b[irow][l],b[icol][l]);
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) cerr<<"gaussj: Singular Matrix"<<endl;
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				swap(a[k][indxr[l]],a[k][indxc[l]]);
	}
}


void lfit(vector<vector<double> > &x, vector<double> &y, vector<double> &sig, vector<double> &a,
	vector<bool> &ia, vector<vector<double> > &covar, double &chisq,
	  void funcs(vector<double>, vector<double> &))
{

  for (int j=0;j<a.size();j++){
    a[j]=0.;
  }
  int i,j,k,l,m,mfit=0;
  double ym,wt,sum,sig2i;
  
  int ndat=x.size();
  int nmag=x[0].size();
  
  int ma=a.size();
  vector<double> afunc(ma);
  vector<vector<double> > beta(ma,vector<double>(1));
  for (j=0;j<ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) cerr<<"lfit: no parameters to be fitted"<<endl;
  for (j=0;j<mfit;j++) {
    for (k=0;k<mfit;k++) covar[j][k]=0.0;
    beta[j][0]=0.0;
  }
  //cout<<"mfit = "<<mfit<<endl;
  for (i=0;i<ndat;i++) {
    vector<double> tempmag(nmag);
    for (int ii=0;ii<nmag;ii++){
      tempmag[ii]=x[i][ii];
    }
    funcs(tempmag,afunc);
    ym=y[i];
    if (mfit < ma) {
      for (j=0;j<ma;j++)
	if (!ia[j]) ym -= a[j]*afunc[j];
    }
    sig2i=1.0/(sig[i]*sig[i]);
    
    for (j=0,l=0;l<ma;l++) {
      if (ia[l]) {
	wt=afunc[l]*sig2i;
	for (k=0,m=0;m<=l;m++)
	  if (ia[m]) covar[j][k++] += wt*afunc[m];
	beta[j++][0] += ym*wt;
      }
    }
  }
  for (j=1;j<mfit;j++)
    for (k=0;k<j;k++)
      covar[k][j]=covar[j][k];
  vector<vector<double> > temp (mfit,vector<double>(mfit));
  for (j=0;j<mfit;j++)
    for (k=0;k<mfit;k++)
      temp[j][k]=covar[j][k];
  gaussj(temp,beta);
  for (j=0;j<mfit;j++)
    for (k=0;k<mfit;k++)
      covar[j][k]=temp[j][k];
  for (j=0,l=0;l<ma;l++)
    if (ia[l]) a[l]=beta[j++][0];
  chisq=0.0;
  for (i=0;i<ndat;i++) {
    vector<double> tempmag(nmag);
    for (int ii=0;ii<nmag;ii++){
      tempmag[ii]=x[i][ii];
    }
    funcs(tempmag,afunc);
    sum=0.0;
    for (j=0;j<ma;j++) sum += a[j]*afunc[j];
    chisq += ((y[i]-sum)/sig[i])*((y[i]-sum)/sig[i]);
  }
  covsrt(covar,ia,mfit);
}
int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}

void func(vector<double>  var, vector<double> &a){
  
  int size=a.size();
  int sizevar=var.size();
  int cor=0;
  int tot=1;
  a[0]=1.0;
  
  for (int i=0;i<sizevar;i++){
    a[i+1]=var[i];
  }
  if (size>sizevar +1){
    for (int i=0;i<sizevar;i++){
      cor+=i;
      
      for (int j=i;j<sizevar;j++){
	a[i*sizevar +j+1+sizevar-cor]=var[i]*var[j];

      }
    }
  }
  int part=sizevar*(sizevar+1)/2+sizevar +1;
  if (size>sizevar*(sizevar+1)/2+sizevar +1){
    int count=0;
    for (int i=0;i<sizevar;i++){

      
      for (int j=i;j<sizevar;j++){
	for (int k=j;k<sizevar;k++){
	  
	  a[count+part]=var[i]*var[j]*var[k];
	  count++;
	  
	  
	}
      }
    }
  }
  
}

double evalfit(vector<double>  var, vector<double> &afunc, vector<double> &a){

  double y=0.;
  func(var,afunc);
  int size=a.size();
  int size2=afunc.size();
  if(size!=size2)
    { cout<<"Sizes don't match in evalfit"<<endl;
      exit(1);
    }
  for (int i=0;i<size;i++){
    y+=a[i]*afunc[i];
  }
  return y;
}

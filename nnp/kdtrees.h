#include <math.h>
#include <vector>
#include <iostream>

using namespace std;
//const int dim = 5;  //number of dimensions
const double big = 150; //upper bound on the coordinates. Lower bound set to -50.
  

//A point in magnitude space
template<int dim> struct Point
{
	
	double x[dim];
	double z;
	Point operator = (Point y)
	{
		for(int q = 0; q < dim; q++)
				x[q] = y.x[q];
			z = y.z;
			return *this;
	}
	Point(double y[dim])
	{
		for(int i = 0; i < dim; i++)
			x[i] = y[i];
		z=0;
	}
	Point()
	{
		z=0;
		for(int i = 0; i < dim; i++)
			x[i] = 0;
	}
	Point(double y[],double zin)
	{
		z=zin;
		for(int i=0;i < dim;i++)
			x[i] = y[i];
	}
	bool operator == (Point y)
	{
		
		for(int i = 0; i < dim; i++)
			if(x[i] != y.x[i])
				return false;
		if(z != y.z)
			return false;
		return true;
	}
};

//A k-dim hypercube 
template<int dim> struct Hcube
{
	
	Point<dim> low,  high;   //The two corners with lowest and highest values
	Hcube(){}
	Hcube(Point<dim> &ilo, Point<dim> &ihi): low(ilo), high(ihi){}
};

//euclidean distance between hyperbox and point 
//equivalent to the infimum of all distances to points in the box
template<int dim> double d(Hcube<dim> &a, Point<dim> &b)
{
	double sum = 0;
	
	for(int y=0; y < dim; y++)
	{
	
		if(b.x[y] < a.low.x[y])
			sum += (a.low.x[y] - b.x[y])*(a.low.x[y] - b.x[y]);
		if(b.x[y] > a.high.x[y])
			sum += (b.x[y] - a.high.x[y])*(b.x[y] - a.high.x[y]);
	}
	return sqrt(sum);
}

//euclidean distance between two points in magnitude space
template<int dim> double d(Point<dim> &a, Point<dim> &b)
{
	double sum = 0;
	for(int y = 0; y < dim; y++)
		sum+= (a.x[y]-b.x[y])*(a.x[y]-b.x[y]);
	return sqrt(sum);
}

//if Point b is is the hypercube a 
template<int dim> bool inCube(Hcube<dim> &a, Point<dim> &b)
{
	if(d(a,b)==0)
		return true;
	return false;
}

//Node of the binary tree
template<int dim> struct treenode: Hcube<dim>
{
	//int currentdim; //Dimension of node 
	int schild;  //Left node or all values with greater mag in currentdim
	int bchild; //right node
	int parent; //parent node 
	int lowind; // lowest number in array of indicies
	int highind;  // highesst number in array of indices
	treenode(Point<dim> hi, Point<dim> lo, int lowi, int highi, int pari, int left, int right): Hcube<dim>(lo,hi)
	{
		lowind = lowi;
		highind = highi;
		parent = pari;
		schild = left;
		bchild = right;
	}
	treenode() {}
};



// SWAPS THE VALUES OF a AND b
void swap(int& a,int& b) //swaps values of a and b
{
	int temp = a;
	a = b;
	b = temp;
} 

//partition elements by quick-sort algorithm using an array of indexes
int partition(int k, double *mags,int *ind, int n)
{
	int parl, parh, up, down, med, vali;
	parh = n-1;	
	parl = 0;
	double val;
	while(true)
	{
		if(parh <= parl + 1)
		{
			if(parh == (parl + 1) && mags[ind[parh]] < mags[ind[parl]])
				swap(ind[parl],ind[parh]);
			return ind[k];
		}
		else
		{
			med = (parh + parl)/2;
			swap(ind[med],ind[parl+1]);
			if(mags[ind[parl]] > mags[ind[parh]])
				swap(ind[parl], ind[parh]);
			if(mags[ind[parl+1]] > mags[ind[parh]])
				swap(ind[parl+1],ind[parh]);
			if(mags[ind[parl]] > mags[ind[parl+1]])
				swap(ind[parl],ind[parl+1]);
			up = parl + 1;
			down = parh;
			vali = ind[parl+1];
			val = mags[vali];
			while(true)
			{
				do up++;
				while(mags[ind[up]] < val);
				do down--;
				while(mags[ind[down]] > val);
				if(down < up) break;
				swap(ind[up], ind[down]);
			}
			ind[parl + 1] = ind[down];
			ind[down] = vali;
			if(down >= k) parh = down - 1;
			if(down <= k) parl = up;
		}

	}
	return 1; 
}



template<int dim> struct kdtree
{
	
	treenode<dim> *nodes;
	vector< Point<dim> > &pts;
	double *cords;
	vector<int> ptind;
	int *revind;
	int npts;
	//CONSTRUCTOR takes in a vector of Points
	kdtree(vector< Point<dim> > &mags): pts(mags)
	{
		
		
		int cur, tpar, tdim, curlo, curhi, np, k, boxind;  //Task loop variables
		int boxnum, m;  //Number of box variables
		int pars[50], dims[50];
		int n = mags.size();  //Number of galaxies
		npts = n;
		int *indp;
		double *corp;
		
		for(int i = 0; i < n; i++) //Initializing index array
			ptind.push_back(i);
		
		cords = new double[n*dim];  //Creating an array of all the first magnitudes followed by second, third, fourth, and fifth
		
		for(int i = 0; i < n; i++)
			for(int j = 0; j < dim; j++)
				cords[i + n*j] = pts[i].x[j];
		
		//Creating Upper and Lower Boundaries
		double temp[dim];
		for(int i = 0; i < dim; i++)
			temp[i] = -big;
		Point<dim> lo(temp);
		for(int i = 0; i < dim; i++)
			temp[i] = big;
		Point<dim> hi(temp);

		//Calculating the number of box neededs: boxnum
		m = 1;
		while(m < n)
			m*= 2;
		boxnum = m -1;
		if(2*n - m/2 -1 < boxnum)
			boxnum= 2*n - m/2 -1;
		
		//initializing array of nodes
		nodes = new treenode<dim>[boxnum];
		nodes[0] =  treenode<dim>(hi,lo,0,n-1,0,0,0);
		
		//Task LOOP
		pars[1] = 0;
		dims[1] = 0;
		cur = 1;
		boxind = 0;
		while(cur)
		{
			tpar = pars[cur];
			tdim = dims[cur];
			cur--;
			curlo = nodes[tpar].lowind;
			curhi = nodes[tpar].highind;
			indp = &ptind[curlo];
			corp = &cords[tdim*n];
			np	 = curhi- curlo + 1; //points 
			k = (np-1)/2; 
			(void) partition(k,corp,indp,np);
			hi = nodes[tpar].high;
			lo = nodes[tpar].low;
			hi.x[tdim] = cords[tdim*n + indp[k]];
			lo.x[tdim] = cords[tdim*n + indp[k]];
			boxind++;
			nodes[boxind] = treenode<dim>(hi,nodes[tpar].low,curlo,curlo+k,tpar,0,0);  //Creates the smalles daughter box
			boxind++;
			nodes[boxind] = treenode<dim>(nodes[tpar].high,lo,curlo+k+1,curhi,tpar,0,0); // Creates the larger daughter box
			nodes[tpar].schild = boxind - 1;  //sets the children
			nodes[tpar].bchild = boxind;
			if(k > 1) //if left box needs to be subdivided
			{
				cur++;
				pars[cur] = boxind-1;
				dims[cur] = (tdim+1)%dim; // Increments the dimension. sets back to 0 if tdim = dim
			}
			if(np-k > 3) //if right box needs subdivisions
			{
				cur++;
				pars[cur] = boxind;
				dims[cur] = (tdim+1)%dim; 
			}
		}	
		revind = new int[n];

		for(int j = 0; j < npts; j++)
			revind[ptind[j]] = j; 
		
	} 
	int findCube(Point<dim> & p) //returns the index of the cube containig Point p
	{
		int num = 0;
		int curdim = 0;
		int ldau;
		while(nodes[num].schild != 0) //if the node isn't a leaf
		{
			ldau = nodes[num].schild;
			if(p.x[curdim] <= nodes[ldau].high.x[curdim]) 
				num = ldau;
			else
				num = nodes[num].bchild;
			curdim = (curdim+1)%dim;
		}
	
		return num;
	}
	int findBox(int p)
	{
		int num = 0;
		int ind = revind[p];
		int ldau;
		int curdim = 0;
		while(nodes[num].schild > 0)
		{
			ldau = nodes[num].schild;
			if(ind <= nodes[ldau].highind)
				num = ldau;
			else 
				num = nodes[num].bchild;
			curdim = (curdim+1)%dim;
		}
		return num;
	}
	double d2(Point<dim>& a, Point<dim>& b) //returns a large distance if the points are the same
	{
		if(a==b)
			return big*dim;
		return d(a,b);
	}
	double d2(int i, int j)
	{
		return d2(pts[i],pts[j]);
	}
	
	//n-nearest neughbours to Point p
	//ds distances with ds[0] biggest, ns corresponding indicies
	void nneigh(int n, double *ds, int *ns, int j)
	{
		int boxi;
		double dcur;
		if(n > npts -1)
			throw("Not Enough Points");
		for(int i = 0; i < n; i++)
			ds[i] = big*dim;
		boxi = nodes[findBox(j)].parent;
		while(nodes[boxi].highind - nodes[boxi].lowind < n)  
			boxi = nodes[boxi].parent;
		for(int i = nodes[boxi].lowind; i <= nodes[boxi].highind; i++)
		{
			if(j == ptind[i]) continue;
			dcur = d2(j,ptind[i]);
			if(dcur < ds[0])
			{
				ds[0] = dcur;
				ns[0] = ptind[i];
				if(n>0) fixheap(ds,ns,n);
			}
		}

		int cur = 1;
		int task[100];
		task[1] = 0;
		int curbox;
		while(cur)//cur > 0
		{
			curbox = task[cur];
			cur--;
			if(boxi == curbox) continue;
			if( d(nodes[curbox],pts[j]) < ds[0])
			{
				if(nodes[curbox].schild > 0)
				{	cur++;
					task[cur] = nodes[curbox].schild;
					cur++;
					task[cur] = nodes[curbox].bchild;
				}
				else
				{
					for(int i = nodes[curbox].lowind; i <= nodes[curbox].highind; i++)
					{
						dcur = d2(ptind[i],j);
						if(dcur < ds[0])
						{
							ds[0] = dcur;
							ns[0] = ptind[i];
							if(n>1) fixheap(ds,ns,n);
						}
					}
				}
			}
		}
		return;
	}

	
	void nneigh(int n, double *ds, int *ns, Point<dim> &p)  
	{
		
		int boxi;
		double dcur;
		
		if(n > npts-1)                            // Not enough points to return the n neighbours
			throw("Not enough points");
		for(int i = 0; i < n; i++)
			ds[i] = big*dim;
		
		//find a box containing n-Points and given point
		boxi = nodes[findCube(p)].parent;
		while(nodes[boxi].highind - nodes[boxi].lowind < n)  
			boxi = nodes[boxi].parent;
		
		//Keep the n closest Points in this box
		for(int i = nodes[boxi].lowind; i <= nodes[boxi].highind; i++)
		{
			if(p == pts[ptind[i]]) 
				continue;
			dcur = d2(pts[ptind[i]],p);
			if(dcur < ds[0])
			{
				ds[0] = dcur;
				ns[0] = ptind[i];
				if(n > 1) fixheap(ds,ns,n);
			}
		}
		
		int cur = 1;
		int task[100];
		task[1] = 0;
		int curbox;
		while(cur)//cur > 0
		{
			curbox = task[cur];
			cur--;
			if(boxi == curbox) continue;
			if( d(nodes[curbox],p) < ds[0])
			{
				if(nodes[curbox].schild > 0)
				{	cur++;
					task[cur] = nodes[curbox].schild;
					cur++;
					task[cur] = nodes[curbox].bchild;
				}
				else
				{
					for(int i = nodes[curbox].lowind; i <= nodes[curbox].highind; i++)
					{
						dcur = d2(pts[ptind[i]],p);
						if(dcur < ds[0])
						{
							ds[0] = dcur;
							ns[0] = ptind[i];
							if(n>1) fixheap(ds,ns,n);
						}
					}
				}
			}
		}
			
			
		
	
		
	}
	
	void fixheap(double *ds, int *ns, int ni) //fixes a heap where the 0-th element is possibly out of place
	{
		int n = ni -1;
		double v = ds[0];
		int vind = ns[0];
		int jhi = 0;
		int jlo = 1;
		while(jlo <= n)
		{
			if(jlo < n && ds[jlo] < ds[jlo+1])  //If the right node is larger
				jlo++;
			if(v >= ds[jlo]) //it forms a heap already
				break;
			ds[jhi] = ds[jlo]; //promotes the bigger of the branches
			ns[jhi] = ns[jlo];
			jhi = jlo; //move down the heap
			jlo = 2*jlo + 1; // calculates position of left branch
		}
		ds[jhi] = v; //places v, vind at correct position in heap
		ns[jhi] = vind;
	
	}

	//Returns the distances and indices of all points within radius r of Point p
	int pointsr(double r, Point<dim> & p,  int *inds)
	{
		int box = 0;
		int curdim = 0;
		int oldbox,dau1,dau2;
		int nret = 0;
		while(nodes[box].schild > 0)
		{
			oldbox = box;
			dau1 = nodes[box].schild;
			dau2 = nodes[box].bchild;
			if(p.x[curdim] + r <= nodes[dau1].high.x[curdim]) box = dau1;
			else if(p.x[curdim] - r >= nodes[dau2].low.x[curdim]) box = dau2;
			curdim++;
			curdim%dim;
			if(box == oldbox) break;

		}
		int task[100];
		int cur = 1;
		task[1] = box;
		int curbox;
		while(cur)
		{
			box = task[cur];
			cur--;
			if(nodes[box].schild != 0)
			{
				if(d(nodes[nodes[box].schild],p) <= r)
				{
					cur++;
					task[cur] = nodes[box].schild;
				}
				if(d(nodes[nodes[box].bchild],p) <= r)
				{
					cur++;
					task[cur] = nodes[box].bchild;
				}
			}
			else
			{
				for(int j= nodes[box].lowind; j <= nodes[box].highind; j++)
				if(d(pts[ptind[j]],p) <= r)
				{
					
					inds[nret] = ptind[j];
					nret++;
				}

			}
		}




		return nret;
	}

	
};

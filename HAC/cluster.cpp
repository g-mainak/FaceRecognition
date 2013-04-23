#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<algorithm>
#include <ctime>

using namespace std;

typedef struct
{
	double centroid[8];
	double prev_mean[8];
    double curr_mean[8];
    int n[400];
    unsigned int count;
} cluster;

int kmeans();
int find_cluster();
int dec_split( vector<cluster> & );
int split( vector<cluster> & );
int merge(vector<cluster> & );
int c_merge(vector<cluster> &, int,int);
int kmeans_part(vector<cluster> &);
unsigned int i;
unsigned int j;
unsigned int k;
//unsigned int numFeatures=1, numPictures=4;
vector <cluster> CL;
unsigned int N = 271;
unsigned int nI = N*4;
unsigned int nA = 4;
double **img, threshold;
int fileCounter = 0;
double avg_count = 0 , avg_total = 0;
int cc[60] = {0};
int total = 0;

int sz = CL.size();
    double STD[ 300 ][ 6 ];
    double std_max[ 300 ];
    double sum = 0 , diff = 0;
    double AVG[ 300 ];
    int decide[ 300 ] = {0};
    int pos_max[ 300 ]={0};
    int Ni = N * 3;

double euclid_dist(double *d1, double *d2)
{
	double dist = 0.0;
	for(unsigned int i=0; i<nA; i++)
	{
		dist += (d1[i]-d2[i])*(d1[i]-d2[i]);
	}
	return sqrt(dist);
}

double man_dist(double *d1, double *d2)
{
	double dist = 0.0;
	for(unsigned int i=0; i<nA; i++)
	{
		dist += abs(d1[i]-d2[i]);
	}
	return (dist / nA);
}

string int2string(int n)
{
       int i,j;
       string s="";
       while( n > 0 )
       {
              char c = '0' + n%10;
              s = c + s;
              n = n / 10;
       }
       return s;
}

void combine( int first, int second, vector<vector<double> > *distMat, vector <cluster> *clusters)
{
	//int m = clusters->at(first).count, n = clusters->at(second).count;
	int m=1 , n=1;
    for(unsigned int i=0; i<nA; i++)
	{
		clusters->at(first).centroid[i] = (m*clusters->at(first).centroid[i] + n*clusters->at(second).centroid[i])/(m+n);
	}
	for(unsigned int i=0; i<clusters->at(second).count; i++)
		clusters->at(first).n[i+ clusters->at(first).count] = clusters->at(second).n[i];
    clusters->at(first).count += clusters->at(second).count;
	vector<cluster>::iterator i = clusters->begin();
	i += second;
	clusters->erase(i);
}

void printCluster(vector<cluster> clusters , int fileC)
{
	string s = "cluster_" ; 
	s = s + int2string(fileC) + ".txt";
    ofstream out(s.c_str());
	//for(vector<cluster>::iterator i = clusters.begin(); i != clusters.end(); i++)
		//cout<< "Cluster centre:\t" << i->centroid[0]<<"\t"<<i->centroid[1]<<"\t"<<i->centroid[2]<<"\t"<<i->centroid[3]<<endl;
	for(unsigned int j = 0 ; j < clusters.size() ; j++ )
	{
		//out<<"Cluster No. "<<j<<" "<<" No of elements in Cluster "<<CL[j].count<<endl;
		out<<j<<" "<<clusters[j].count<<" ";
        for(unsigned int k=0 ; k<nA ; k++ )
			//out<<CL[j].centroid[k]<<" ";
		sort(clusters[j].n,clusters[j].n+clusters[j].count);
		for(unsigned int  k =0 ; k < clusters[j].count ; k++ )
			out<<clusters[j].n[k]<<" ";
		out<<endl;
	}
	out.close();
}

void HAC( vector<cluster> *clusters)
{
	int first, second, row, col;
	double leastDist=1e300, f;
	for(vector<cluster>::iterator i = clusters->begin(); i!= clusters->end(); ++i)
	{
		for(vector<cluster>::iterator j = clusters->begin(); j!=clusters->end(); ++j)
		{
			double temp = euclid_dist(i->centroid, j->centroid);
			if (temp!=0)
				leastDist = (temp<leastDist)?temp:leastDist;
		}
	}
	//cout<<"LeastDist: "<<leastDist<<endl;
	f=leastDist;
	do
	{
		leastDist=1e300;
		vector<vector<double> > distMat (clusters->size(), vector<double> (clusters->size()));
		for(vector<double>::iterator i = distMat[0].begin(); i!= distMat[0].end(); i++)
		{
			row = i - distMat[0].begin();
			for(vector<double>::iterator j = distMat.at(row).begin(); row > (j - distMat.at(row).begin()) ; j++)
			{
				col = j - distMat.at(row).begin();
				distMat.at(row).at(col) = euclid_dist( clusters->at(row).centroid , clusters->at(col).centroid);
				if(leastDist>*j)
				{
					leastDist=*j;
					first= row;
					second= col;
				}
			}
		}
		combine(second, first, &distMat, clusters);
	} while( clusters->size()>23);
	//printCluster(*clusters,99);
}


int main(int argc, char** argv)
{
	clock_t start, end;
	double runTime;
	start = clock();
    img = new double *[nI];
	for( i=0 ; i<nI ; i++ )
	{
         img[i] = new double [8];
    }
    ifstream in1;
	in1.open("dist4.txt");
	for( i=0 ;  i<nI ; i++ )
    {
         for( j=0 ; j<4 ; j++ )
         {
              in1>>img[i][j];
         }
         //swap(img[i][0],img[i][4]);
         //swap(img[i][1],img[i][5]);
         
    }
    in1.close();
    kmeans();
    
    //HAC(&CL);
    cout<<CL.size()<<endl;
    while(1)
    {
            dec_split(CL);
            merge(CL);
            if( (CL.size() < 140) )
                break;
    }
    find_cluster();
    printCluster(CL,2);
    cout<<CL.size()<<endl;
    cout<<avg_count / N <<endl; 
    for( i=0 ; i<60 ; i++ )
         cout<< (cc[i] *100 ) / N <<" "<<((23 + i*10)*100)/(N*3)<<endl;
    end = clock();
    	runTime = (end - start) / (double) CLOCKS_PER_SEC ;
    	printf ("Run time is %g seconds\n", runTime);
    system("pause");

    return 0;
}



int find_cluster()
{
     ofstream out("fd.txt");
     ofstream ot("rank.txt");
     i = 3;
     for( i=3 ; i<nI ; i=i+4 )
     {
            
            double dist[ CL.size() ] , *test_img;
            test_img = img[i];
            int pos;
            for( j=0 ; j < CL.size() ; j++ )
                 dist[j] = 0;
     
             for( j =0 ; j < CL.size() ; j++ )
             {
                  dist[j] = euclid_dist( CL[j].centroid , test_img );
                  //dist[j] = man_dist( CL[j].centroid , test_img );
             }
             ot<<"Image No. "<<i<<"\t";
             out<<i<<" ";
             /*
             for( j =0 ; j < CL.size() ; j++ )
                  out<<dist[j]<<" ";
             out<<endl;
             */
             
             int maxi = 1000;
             int flag = 0; 
             int pos_f = -1;
             int count_img = 0;
             int crank[CL.size()];
             int exact_rank = -1;
             vector <cluster> newCL;
             count_img = 0 ;
             for( j = 0 ; j < CL.size() ; j++ ) // 5 best Ranks
             {
                  pos = min_element( dist , dist + CL.size() ) - dist;
                  dist[ pos ] = maxi;
                  out<<pos<<" ";
                  if( flag==0  )
                      count_img += CL[pos].count;   
                  for(int g=0; g<CL[pos].count; g++)
               	      if(CL[pos].n[g]<=(i-1) && CL[pos].n[g]>=(i-3) && flag ==0 )
               	      {
                            pos_f = j;
                            flag = 1;
                      }
                  //cout<<endl;

             }
             
             avg_count += count_img;
             avg_total += nI;
             for( j=0 ; j<60 ; j++ )
               if( count_img <= 23 + j*10 )
               {
                   cc[j]++;
               }
             
             ot<<pos_f<<endl;
             out<<pos_f<<endl;
             printCluster( newCL , 102 );
             
     }
     
     out.close();
     ot.close();
     return 0;
}



int kmeans()
{
  
  
  for( i=0 ; i < nI ; i+=4 ) // formation of centroid by 1st image of each person
    {
         cluster c1;
         c1.count = 0;
         c1.n[ c1.count ] = i;
         for( j=0 ; j<nA ; j++ )
         {
                  c1.centroid[ j ] = img[ c1.n[c1.count] ][j];
         }
         c1.count++;
         CL.push_back( c1 );
         
    }
    
    double dist[ CL.size() ];
    for( i=0 ; i<nI ; i++ )
    {
        if( (i % 4 == 0) || (i % 4 == 3 ) ) // escaping the 0th and 3rd image of every person
            continue;
             int pos;
             for( j=0 ; j < CL.size() ; j++ )
                dist[j] = 0.0;
             for( j =0 ; j < CL.size() ; j++ )
             {
                 dist[j] = euclid_dist( CL[j].centroid , img[i] );
                 //dist[j] = man_dist( CL[j].centroid , img[i] );
             }
             pos = min_element( dist , dist + CL.size() ) - dist;
             CL[ pos ].n[ CL[ pos ].count ] = i ;

             CL[ pos ].count++;

    }
    for( i=0 ; i<CL.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
                 for( k=1 ; k<CL[i].count ; k++ )
                 {
                      CL[i].centroid[j] += img[ CL[i].n[k] ][j] ;
                 }
                 CL[i].centroid[j] /= CL[i].count;
            }
       }
  
    
  for( int x = 0 ; x < 100000 ; x++ )
  { 
       //printCluster(CL,x);
       
       for( i=0 ; i<CL.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
              CL[i].prev_mean[j] = CL[i].centroid[j];
              CL[i].centroid[j] = 0;
            }
            CL[i].count = 0;
            
       }
       for( i=0 ; i<nI ; i++ )
       {
             if( (i % 4 == 3 ) ) // escaping the 3rd image of every person
                continue;
             int pos;
             for( j=0 ; j < CL.size() ; j++ )
                dist[j] = 0.0;
             for( j =0 ; j < CL.size() ; j++ )
             {
                dist[j] = euclid_dist( CL[j].prev_mean , img[i] );
                 //dist[j] = man_dist( CL[j].prev_mean , img[i] );
             }
             pos = min_element( dist , dist + CL.size() ) - dist;
             CL[ pos ].n[   CL[ pos ].count   ] = i;

             CL[ pos ].count++;

       }
       //cout<<CL.size();
       for( i=0 ; i<CL.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
                 for( k=0 ; k<CL[i].count ; k++ )
                 {
                      CL[i].centroid[j] += img[ CL[i].n[k] ][j] ;
                 }
                 CL[i].centroid[j] /= CL[i].count;
            }
       }
       
       int flag = 0;    
       for( i=0 ; i<CL.size() ; i++ )
       {
             if( CL[i].count == 0 )
             {
                 CL.erase(CL.begin() + i);
                 continue;
             }
             for( j=0 ; j < nA ; j++ )
             {
                  if( abs(CL[i].centroid[j] - CL[i].prev_mean[j]) > 0.0001 )
                      flag = 1;
             }
       }
       if( flag == 0) 
           break;
              
       
  }
  
  
}


    

int dec_split( vector<cluster> &clusters )
{
    int i,j,k;
    int sz = clusters.size();
    double sum = 0 , diff = 0;
    double AVG[ sz ];
    //int decide[ sz ] = {0};
    //int pos_max[sz]={0};
    //int Ni = 208 * 3;
    double avg_thresh = 0;
    double std_thresh = 0;
    double no_thresh = 10;
    for( i = 0 ; i<sz ; i++ )
    {
           diff = 0;
           for( k=0 ; k<nA ; k++ )
           {
                STD[i][k] = 0;
                for( j=0 ; j<clusters[i].count ; j++ )
                {
                     diff += pow( ( img[ clusters[i].n[j] ][k] - clusters[i].centroid[k] ) , 2 );
                }
                STD[i][k] = sqrt( diff / clusters[i].count );
           }           
    }
    
    for( i = 0 ; i<sz ; i++ )
    {
           diff = 0;
           int  pos = max_element( STD[i] , STD[i] + nA ) - STD[i] ; 
           std_max[i] = STD[i][pos];
           pos_max[i]=pos;
    }
    
    for( i = 0 ; i<sz ; i++ )
    {
           diff = 0;
           for( j=0 ; j<clusters[i].count ; j++ )
           {
                diff += euclid_dist( img[ clusters[i].n[j] ], clusters[i].centroid ); 
           }          
           AVG[ i ] = diff / clusters[i].count ;
           avg_thresh += diff;
           //cout<<i<<" "<<AVG[i]<<endl;
           
    }
    avg_thresh /= (double)Ni;
    //cout<<avg_thresh<<endl;
    for( i=0 ; i<sz ; i++ )
    {
         if( ( std_max[i] > std_thresh ) && ( AVG[i] > avg_thresh ) && (clusters[i].count > no_thresh) )
         {
             //cout<<i<<endl;
             decide [i] = 1;
         }
    }
    split(clusters);   
    
}
    
int split(vector<cluster> &clusters )
{
    cluster c1;
    double alpha = 0.2;
    int sz = clusters.size();
    for( i=0; i<sz; i++)
    {
         if (decide[i]==1)
         {
              //cout<<i<<" "<<pos_max[i]<<" "<<std_max[i]<<endl;
              c1.count=0;
              for (j=0;j<nA;j++)
              {    
                  c1.centroid[j]=clusters[i].centroid[j];
              }
              c1.centroid[pos_max[i]]-=alpha*std_max[i];
              clusters[i].count=0;
              clusters[i].centroid[pos_max[i]] +=alpha*std_max[i];
              clusters.push_back(c1);
         }
         
    }
    double dist[clusters.size()];
  for( int x = 0 ; x < 100000 ; x++ )
  { 
       //printCluster(CL,x);
       
       for( i=0 ; i<clusters.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
              clusters[i].prev_mean[j] = clusters[i].centroid[j];
              clusters[i].centroid[j] = 0;
            }
            clusters[i].count = 0;
            
       }
       for( i=0 ; i<nI ; i++ )
       {
            if( (i % 4 == 3 ) ) // escaping the 0th and 3rd image of every person
                continue;
            int pos;
             for( j=0 ; j < clusters.size() ; j++ )
                dist[j] = 0.0;
             for( j =0 ; j < clusters.size() ; j++ )
             {
                 dist[j] = euclid_dist( clusters[j].prev_mean , img[i] );
                 //dist[j] = man_dist( CL[j].prev_mean , img[i] );
             }
             pos = min_element( dist , dist + clusters.size() ) - dist;
             clusters[ pos ].n[   clusters[ pos ].count   ] = i;

             clusters[ pos ].count++;

       }
       //cout<<CL.size();
       for( i=0 ; i<clusters.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
                 for( k=0 ; k<clusters[i].count ; k++ )
                 {
                      clusters[i].centroid[j] += img[ clusters[i].n[k] ][j] ;
                 }
                 clusters[i].centroid[j] /= clusters[i].count;
            }
       }
       
       int flag = 0;    
       for( i=0 ; i<clusters.size() ; i++ )
       {
             /*if( clusters[i].count == 0 )
             {
                 clusters.erase(clusters.begin() + i);
                 continue;
             }*/
             for( j=0 ; j < nA ; j++ )
             {
                  if( abs(clusters[i].centroid[j] - clusters[i].prev_mean[j]) > 0.0001 )
                      flag = 1;
             }
       }
       
       if( flag == 0) 
           break;
   }   
   //printCluster(clusters , 111);         
}    


int merge(vector<cluster> & clusters)
{
    int sz = clusters.size();
    double dist[sz][sz];
    double theta = 0.05;
    int a1 , a2;
    int cnt=0;
    int P = clusters.size() / 2;
    for( i=0 ; i<sz ; i++ )
    {
         for( j=0 ; j<sz ; j++ )
         {
              if( i < j )
                  dist[i][j] = euclid_dist( clusters[i].centroid , clusters[j].centroid );
              else
                  dist[i][j] = 1e300;
         }
    }
    ofstream out("distmatrix.txt");
    for( i=0 ; i<sz ; i++ )
    {
         for( j=0 ; j<sz ; j++ )
         {
              out<<i*sz+j<<" "<<dist[i][j]<<endl;
              
         }
         out<<endl;
    }
    for( i=0 ; i<sz ; i++ )
    {
         for( j=0 ; j<sz ; j++ )
         {
              if (cnt<P)
              {  int pos = min_element( dist[0] , dist[0] + sz*sz ) - dist[0];
                 a1 = pos / sz;
                 a2 = pos % sz;
                 if (dist[a1][a2]<theta)
                 {       dist[ a1 ][ a2 ] =1e300;
                         for( k=a1 + 1 ; k<sz ; k++ )
                              dist[ a1 ][ k ] = 1e300;
                         for( k = 0  ; k < a1 ; k++ )
                              dist[ k ][ a1 ] = 1e300;
                         for( k=a2 + 1 ; k<sz ; k++ )
                              dist[ a2 ][ k ] = 1e300;
                         for( k = 0  ; k < a2 ; k++ )
                              dist[ k ][ a2 ] = 1e300;
                 }
                 c_merge(clusters,a1,a2);
                 cnt++;
              }
         }
    }
    kmeans_part(clusters);     
    //printCluster(clusters , 112); 

}
                  
    
int c_merge(vector<cluster> & clusters,int a1,int a2)
{
    cluster c1;
    c1.count=0;
    for (i=0;i<nA;i++)
    {
         c1.centroid[i]= clusters[a1].count*clusters[a1].centroid[i]+clusters[a2].count*clusters[a2].centroid[i];
         c1.centroid[i]/=(clusters[a1].count+clusters[a2].count);
    }
    clusters.erase(clusters.begin()+a1);
    clusters.erase(clusters.begin()+a2);
    clusters.push_back(c1);
}   

int kmeans_part(vector<cluster> & clusters)
{
     double dist[clusters.size()];
  for( int x = 0 ; x < 100000 ; x++ )
  { 
       //printCluster(CL,x);
       
       for( i=0 ; i<clusters.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
              clusters[i].prev_mean[j] = clusters[i].centroid[j];
              clusters[i].centroid[j] = 0;
            }
            clusters[i].count = 0;
            
       }
       for( i=0 ; i<nI ; i++ )
       {
            if( (i % 4 == 3 ) ) // escaping the 0th and 3rd image of every person
                continue;
            int pos;
             for( j=0 ; j < clusters.size() ; j++ )
                dist[j] = 0.0;
             for( j =0 ; j < clusters.size() ; j++ )
             {
                 dist[j] = euclid_dist( clusters[j].prev_mean , img[i] );
                 //dist[j] = man_dist( CL[j].prev_mean , img[i] );
             }
             pos = min_element( dist , dist + clusters.size() ) - dist;
             clusters[ pos ].n[   clusters[ pos ].count   ] = i;

             clusters[ pos ].count++;

       }
       //cout<<CL.size();
       for( i=0 ; i<clusters.size() ; i++ )
       {
            for( j=0 ; j<nA ; j++ )
            {
                 for( k=0 ; k<clusters[i].count ; k++ )
                 {
                      clusters[i].centroid[j] += img[ clusters[i].n[k] ][j] ;
                 }
                 clusters[i].centroid[j] /= clusters[i].count;
            }
       }
       
       int flag = 0;    
       for( i=0 ; i<clusters.size() ; i++ )
       {
             if( clusters[i].count == 0 )
             {
                 clusters.erase(clusters.begin() + i);
                 continue;
             }
             for( j=0 ; j < nA ; j++ )
             {
                  if( abs(clusters[i].centroid[j] - clusters[i].prev_mean[j]) > 0.0001 )
                      flag = 1;
             }
       }
       
       if( flag == 0) 
           break;
   }   
}

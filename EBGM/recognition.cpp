/*
 * recognition.cpp
 *
 *  Created on: Jun 2, 2010
 *      Author: Mainak
 */
#include "EBGMperipherals.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>
using namespace std;


double DEPredictiveStep(Node* j1, Node* j2, double *tdx, double *tdy, Mask* masks){
    //BEST!!!
	double j12 = 0;
    double j11 = 0;
    double j22 = 0;
    int i;
    double sim = 0.0;
    double dx  = 0.0,
          dy  = 0.0;
    double Gxx, Gxy, Gyx, Gyy;
    double Px, Py;
    assert(j1 &&  j1->length && j2 && j1->length==j2->length);

    Gxx = 0;
    Gxy = 0;
    Gyx = 0;
    Gyy = 0;

    Px  = 0;
    Py  = 0;

    for(i = 0; i < j1->length; i++){
        double ang = j1->ang[i] - j2->ang[i];

        while(ang >  PI){ ang -= 2*PI; }
        while(ang < -PI){ ang += 2*PI; }

        Gxx += j1->mag[i]*j2->mag[i]*masks->kx[2*i]*masks->kx[2*i];
        Gxy += j1->mag[i]*j2->mag[i]*masks->kx[2*i]*masks->ky[2*i];
        Gyx += j1->mag[i]*j2->mag[i]*masks->kx[2*i]*masks->ky[2*i];
        Gyy += j1->mag[i]*j2->mag[i]*masks->ky[2*i]*masks->ky[2*i];

        Px  += j1->mag[i]*j2->mag[i]*masks->kx[2*i]*(ang);
        Py  += j1->mag[i]*j2->mag[i]*masks->ky[2*i]*(ang);
    }

    if(Gxx*Gyy-Gxy*Gyx != 0.0){
        dx = (Gyy*Px-Gyx*Py)/(Gxx*Gyy-Gxy*Gyx);
        dy = (-Gxy*Px+Gxx*Py)/(Gxx*Gyy-Gxy*Gyx);
    } else {
        dx = 0.0;
        dy = 0.0;
    }

    j12 = 0.0;
    j11 = 0.0;
    j22 = 0.0;
    for(i = 0; i < j1->length; i++){
        j12 += j1->mag[i]*j2->mag[i]*cos(j1->ang[i] - j2->ang[i] - (dx * masks->kx[2*i] + dy * masks->ky[2*i]));
        j11 += j1->mag[i]*j1->mag[i];
        j22 += j2->mag[i]*j2->mag[i];
    }
    sim = j12/sqrt(j11*j22);

    *tdx = dx;
    *tdy = dy;

    return sim;
}

double JetSimilarity(Node* j1, Node* j2, Mask *masks){
    double dx = 0.0, dy = 0.0;
    return DEPredictiveStep(j1, j2, &dx, &dy, masks);
}

double fgSimMagnitude(FaceGraph f1, FaceGraph f2, Mask *masks){
    double totalSim = 0.0;
    int i;

    for(i = 0; i < f1.numNodes; i++){
        totalSim += JetSimilarity( &f1.Nodes[i], &f2.Nodes[i], masks );
    }
    totalSim = totalSim / f1.numNodes;

    return -totalSim;
}

void quickSort(double **arr, int left, int right)
{
      int i = left, j = right;
      double tmp;
      double pivot = arr[(left + right)/2][1];

      do
      {
            while (arr[i][1] < pivot)
                  i++;
            while (arr[j][1] > pivot)
                  j--;
            if (i <= j)
            {
                  tmp = arr[i][1];
                  arr[i][1] = arr[j][1];
                  arr[j][1] = tmp;
                  tmp = arr[i][0];
                  arr[i][0] = arr[j][0];
                  arr[j][0] = tmp;
                  i++;
                  j--;
            }
      }while(i<=j);

      if (left < j)
            quickSort(arr, left, j);
      if (i < right)
            quickSort(arr, i, right);
}

void printArr(int size, double **arr)
{
	for(int i=0; i<size; i++)
	{
		cout<<arr[i][0]<<" ";
	}
	cout<<endl;
}

void count(int k, double** arr)
{
	ofstream a_file ( "count.txt", ios::app );
	for(int i=0; i< 100; i++)
	{
		if(arr[i][0]== k)
		{
			a_file<<i+1<<endl;
			return;
		}
	}
	a_file.close();
}

void recognize(FaceGraph graph, BunchGraph *bg, Mask *masks, int k)
{
	double bestSim = 1e300;
	int bestFace = 0;
    double **arr = new double*[bg->numFaces];
    for(int i=0; i<bg->numFaces; i++)
    {
    	arr[i] = new double[2];
    	arr[i][0] = bg->faces[i].ID;
    }
	for(int i=0; i< bg->numFaces; i++)
	{
		double sim = fgSimMagnitude(graph, bg->faces[i], masks);
		arr[i][1] = sim;
		if (sim<bestSim)
		{
			bestSim = sim;
			bestFace = i;
		}
	}
	quickSort(arr, 0, (bg->numFaces)-1);
	printArr(bg->numFaces, arr);
	count(k, arr);
}

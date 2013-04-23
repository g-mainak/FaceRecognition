/*
 * EBGMCalc.cpp
 *
 *  Created on: May 25, 2010
 *      Author: Mainak
 */
#include "EBGMperipherals.h"
#include "EBGMCalc.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define RANDOM   (fabs(((double)rand())/RAND_MAX))

using namespace std;
const int NUMNODES = 15;
const int numberOfTemplates = 15;

double SIM_DISPLACE(Node* j1, Node*  j2, double dx, double dy, Mask* masks) {
    double j12 = 0.0;
    double j11 = 0.0;
    double j22 = 0.0;
    int i;
    for(i = 0; i < j1->length; i++){
        j12 += j1->mag[i]*j2->mag[i]*
	       cos(j1->ang[i] - j2->ang[i] -
	       (dx * masks->kx[2*i] + dy * masks->ky[2*i]));
        j11 += j1->mag[i]*j1->mag[i];
        j22 += j2->mag[i]*j2->mag[i];
    }
    return j12/sqrt(j11*j22);
}

double dispEst(Node* j1, Node* j2, double *tdx, double *tdy, Mask *masks){
    double sim  = 0.0;
    double dx   = 0.0,
          dy   = 0.0;
    int change = 1,
        iter   = 50;
    double nextx, nexty;
    double bestsim = 0.0;
    double tol  = .2;
    double step = 1.0;
//NARROWING SEARCH
    assert(j1 && j1->length  && j2 && j1->length==j2->length);

    sim = SIM_DISPLACE(j1, j2, dx,dy, masks);
    bestsim = sim;
    *tdx = 0.0;
    *tdy = 0.0;

    change = 1;
    sim = SIM_DISPLACE(j1, j2, dx,dy, masks);
    bestsim = sim;

    nextx = dx;
    nexty = dy;
    while(change && iter--){
        dx = nextx;
        dy = nexty;
        change = 0;

        sim = SIM_DISPLACE(j1, j2, dx+step,dy+step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx+step;
            nexty = dy+step;
            change = 1;
        }

        sim = SIM_DISPLACE(j1, j2, dx+step,dy-step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx+step;
            nexty = dy-step;
            change = 1;
        }

        sim = SIM_DISPLACE(j1, j2, dx-step,dy-step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx-step;
            nexty = dy-step;
            change = 1;
        }

        sim = SIM_DISPLACE(j1, j2, dx-step,dy+step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx-step;
            nexty = dy+step;
            change = 1;
        }
        dx = nextx;
        dy = nexty;

        if(change == 0 && step > tol){
            change = 1;
            step = step * 0.5;
        }
    }

    *tdx = dx;
    *tdy = dy;

    return bestsim;
}


void getJetsFromImage(Node *n, Mask *masks, Image image)
{
	if( masks->numMasks )
	{
		/* This algrithm aproximates the gabor parameters at a location in continuous space
	    By computing the gabor convoluation at an integer location close to the contiuous
		point and then shifting the phase parameters the apropreate amount. It is much faster
		then using linear interpolation and it is probably more accurate. */
		int i;
	    float rx = n->x;
	    float ry = n->y;
	    float dx = n->x - rx;
	    float dy = n->y - ry;
	    n->length= masks->numMasks/2;
	    n->ang = new double[masks->numMasks/2];
	    n->imagPart = new double[masks->numMasks/2];
	    n->mag = new double[masks->numMasks/2];
	    n->realPart = new double[masks->numMasks/2];
	    /* Compute the jet coeffecents */
	    for( i = 0; i < n->length; i++)
	    {
	    	n->realPart[i] = convolvePoint(rx, ry, 0, &image, &masks->masks[2*i]);
	        n->imagPart[i] = convolvePoint(rx, ry, 0, &image, &masks->masks[2*i+1]);
	        //FINITE(jet->realPart[i]);
	        //FINITE(jet->imagPart[i]);
	    }
	    computePolar(n);
	    /* Adjust phase coordinates */
	    for( i = 0; i < n->length; i++)
	    {
	    	n->ang[i] +=  dx * masks->kx[2*i] + dy * masks->ky[2*i];
	    }
	    /* Recompute coeffecents based on phase change */
	    for( i = 0; i < n->length; i++)
	    {
	    	n->realPart[i] = n->mag[i] * cos(n->ang[i]);
	        n->imagPart[i] = n->mag[i] * sin(n->ang[i]);
	    }
	}
	else
	{
	  	n->ang = new double[0];
	   	n->imagPart = new double[0];
	    n->mag = new double[0];
	    n->realPart = new double[0];
	}
}

Node findClosestTemplate(BunchGraph *bg, int nodeNum, Mask* masks, Image image)
{

	double scalex = image.width;
	double scaley= image.height;
	int best = 0;
	double dx = 0.0, dy = 0.0;
    double bestsim = -1.0e300;
	Node *temp = new Node;
	for(int i=0; i<numberOfTemplates; i++)	//Runs for each template
	{
		temp->x = bg->faces[i].Nodes[nodeNum].x*scalex; //Multiplies to get the scale
		temp->y = bg->faces[i].Nodes[nodeNum].y*scaley;
		getJetsFromImage(temp, masks, image);	//Extracts jet from coordinates
		double sim = dispEst( &bg->faces[i].Nodes[nodeNum], temp, &dx, &dy, masks); //calculates similarity between probe image and ith template
		if(sim > bestsim)
        {
            bestsim = sim;
            best = i;	//Stores best template
        }
	}
	dispEst( &bg->faces[best].Nodes[nodeNum], temp, &dx, &dy, masks); //Fine tuning the coordinates by looking at the neighbouring pixels
    temp->x += dx;
    temp->y += dy;
    getJetsFromImage(temp, masks, image); //Final Coordinates of this fiducial point
    return *temp;
}

FaceGraph putGraph(BunchGraph *bg, Mask* masks, Image image)
{
	FaceGraph fg;
	fg.numNodes = NUMNODES;
	fg.Nodes = new Node[NUMNODES];
	for(int i=0; i<NUMNODES; i++) //Runs for each fiducial point.
	{
		fg.Nodes[i] = findClosestTemplate(bg, i, masks, image);
	}
	return fg;
}

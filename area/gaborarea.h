/*
 * gaborarea.h
 *
 *  Created on: Jun 24, 2010
 *      Author: Mainak
 */

#ifndef GABORAREA_H_
#define GABORAREA_H_

#include <assert.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#define IE( img , x , y, c ) 	( (img)->data[(x)+ img->width*(y)] )
#define ie( img , x , y, c ) ((x < img->width && y < img->height && x >= 0 && y >=  0 ) ? IE(img, x, y, c) : 0.0 )
#define EQUAL_ZERO(v,tol)  	( ABS(v) < tol )
#define TRUNC(v)    		( (int) (v) )
#define INT_FLOOR(a) ((int)(a))
#define INT_CEIL(a)  ((int)(a)+1)
#define ABS(v)      		( ((v) < 0)     ? -(v) : (v)  )
#define SQR(v)      		( (v) * (v) )

const double PI = 3.14159;
double sum = 0;
double cmp_area[4];
std::ofstream out1;
typedef struct Image
{
	int width, height, maxVal;
	double *data;
};

typedef struct Mask
{
	int numMasks;
	Image *masks;
    int length;
    double* wavelength;
    double* angle;
    double* phase;
    double* aspect;
    double* radius;
    int*    size;
    /* Precomputed parameters for displacement estimation */
    double* kx;
    double* ky;
};

typedef struct Node
{
	int length;
	double x, y;
	double *realPart, *imagPart, *mag, *ang;
};

typedef struct FaceGraph
{
	int numNodes;
	Node *Nodes;
};

typedef struct BunchGraph
{
	int numFaces;
	FaceGraph* faces;
};

struct Point
{
	double x, y;
};
#endif /* GABORAREA_H_ */

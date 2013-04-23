/*
 * EBGM.h
 *
 *  Created on: May 20, 2010
 *      Author: Mainak
 */

#ifndef EBGM_H_
#define EBGM_H_

#define IE( img , x , y, c ) 	( (img)->data[(x)+ img->width*(y)] )
#define ie( img , x , y, c ) ((x < img->width && y < img->height && x >= 0 && y >=  0 ) ? IE(img, x, y, c) : 0.0 )
#define EQUAL_ZERO(v,tol)  	( ABS(v) < tol )
#define TRUNC(v)    		( (int) (v) )
#define INT_FLOOR(a) ((int)(a))
#define INT_CEIL(a)  ((int)(a)+1)
#define ABS(v)      		( ((v) < 0)     ? -(v) : (v)  )
#define SQR(v)      		( (v) * (v) )

const double PI = 3.14159;
struct Image
{
	int width, height, maxVal;
	double *data;
};

struct Mask
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

struct Node
{
	int length;
	double x, y;
	double *realPart, *imagPart, *mag, *ang;
};

struct FaceGraph
{
	int numNodes;
	int ID;
	Node *Nodes;
};

struct BunchGraph
{
	int numFaces;
	FaceGraph* faces;
};

Image readImage(char*);
void writeImage(Image, const char*);
FaceGraph readGraph(char*);
void writeGraph(char*, FaceGraph);
void computePolar(Node*);
double convolvePoint(double x, double y, int c, Image *im, Image *mask);
double interpLinear(Image *img, double x, double y, int c);
void readBunchGraph();
void readMask();
Image makeGaborMask(double lambda, double theta, double phi,
                    double gamma, double sigma,int size);
void ZeroMeanUnitLength(Image *im);
BunchGraph* addFacetoBunch(BunchGraph* bg, FaceGraph *face);
void writeBunchGraph(char* fileName, BunchGraph *bg);
BunchGraph* readBunchGraph(char* fileName);
Mask* readMask(char* filename);
void writePoints(char*);
Image rotate(Image, double);
Image crop(Image, int);
Image flip(Image);
Image join(Image, Image);

#endif /* EBGM_H_ */

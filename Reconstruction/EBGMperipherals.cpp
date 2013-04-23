/*
 * EBGM.cpp
 *
 *  Created on: May 20, 2010
 *      Author: Mainak
 */
#include "EBGMperipherals.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include <cstring>
#include <cstdlib>

using namespace std;

Image readImage(char* filename){
    int  width, height, max, x, y, min;
    int  val;
    char fchar;
    char line[100];
    char ftype[16];
    FILE *infile;
    Image im;
    /* Read first line, and test if it is propoer PNM type */
    infile = (FILE*) fopen(filename, "rb");
    assert(infile);
    fgets(line,100,infile);
    sscanf(line," %s",ftype);
    if (! (strcmp(ftype,"P5") == 0)) {
        printf("Error <%s,%d>: Currenlty only binary pgm files, type P5, supported",__FILE__,__LINE__);fflush(stdout);
        exit(1);
    }

    /* Read lines, ignoring those starting with Comment Character, until the
        Image Dimensions are read. */
    fchar = '#';
    while (fchar == '#') {
        fgets(line,100,infile);
        sscanf(line, " %c", &fchar);
    }
    sscanf(line, " %d %d", &width, &height);
    /* Read lines, ignoring those starting with Comment Character, until the
        maximum pixel value is read. */
    fchar = '#';
    while (fchar == '#') {
        fgets(line,100,infile);
        sscanf(line, " %c", &fchar);
    }
    sscanf(line, "%d", &max);

    if (! (max == 255)) {
        fprintf(stdout,"readImagePGM: Warning, max value %d for pixels in image %s is not 255\n", max, filename); fflush(stdout);
    }
    im.height = height;
    im.width = width;
    im.maxVal = max;
    im.data = new double[width*height];
    min=256;
    max=0;
    val = fgetc(infile);
    for(y = 0; y < height; y++)
    {
        for(x = 0; x < width; x++)
        {
            im.data[y*width + x] = val;
            //if(im.data[y*im.width + x] >max) max= im.data[y*im.width + x];
            //if(im.data[y*im.width + x] <min) min= im.data[y*im.width + x];
            val = fgetc(infile);
        }
    }
    /*double scale = 256.0/(max-min);
    for(int i=0; i<im.height; i++)
    	for(int j=0; j<im.width; j++)
       		im.data[i*im.width + j] = (im.data[i*im.width + j]-min)*scale;*/
    fclose(infile);
    return( im );
}

void writeImage(Image im, const char* filename)
{
	FILE* outfile = fopen(filename,"wb");
    int x, y;
    if(!outfile)
    {
        printf("could not open %s for writing.\n",filename);
        exit(1);
    }
    fprintf(outfile,"P5\n");
    fprintf(outfile,"%d %d\n",im.width, im.height);
    fprintf(outfile,"%d\n",255);
    for(y = 0; y < im.height; y++)
        for(x = 0; x < im.width; x++)
            fputc((int)im.data[y*im.width + x],outfile);
    fclose(outfile);
}

FaceGraph readGraph(char* fileName)
{
	ifstream in;
	FaceGraph *graph = new FaceGraph;
	in.open(fileName);
	in>>graph->numNodes;
	in>>graph->ID;
	graph->Nodes = new Node[graph->numNodes];
	for(int i=0; i< graph->numNodes; i++)
	{
		in>>graph->Nodes[i].x>>graph->Nodes[i].y>>graph->Nodes[i].length;
		graph->Nodes[i].ang = new double[graph->Nodes[i].length];
		graph->Nodes[i].imagPart = new double[graph->Nodes[i].length];
		graph->Nodes[i].mag = new double[graph->Nodes[i].length];
		graph->Nodes[i].realPart = new double[graph->Nodes[i].length];
		for(int j=0; j<graph->Nodes[i].length; j++)
			in>>graph->Nodes[i].ang[j]>>graph->Nodes[i].imagPart[j]>>graph->Nodes[i].mag[j]>>graph->Nodes[i].realPart[j];
	}
	in.close();
	return *graph;
}

void writeGraph(char* fileName, FaceGraph graph)
{
	ofstream out;
	out.open(fileName);
	out<<graph.numNodes<<endl;
	out<<graph.ID<<endl;
	for(int i=0; i< graph.numNodes; i++)
	{
		out<<graph.Nodes[i].x<<" "<<graph.Nodes[i].y<<" "<<graph.Nodes[i].length<<endl;
		for(int j=0; j<graph.Nodes[i].length; j++)
			out<<graph.Nodes[i].ang[j]<<" "<<graph.Nodes[i].imagPart[j]<<" "<<graph.Nodes[i].mag[j]<<" "<<graph.Nodes[i].realPart[j]<<endl;
	}
	out.close();
}

/*
This is a utility function that computes the polar coordinage transform
for the values in a gabor jet.
*/
void computePolar(Node *jet)
{
    int i;

    for( i = 0; i < jet->length; i++){
        /* Find the magnitude of the imaginary points */
        jet->mag[i] = sqrt((jet->realPart[i])*(jet->realPart[i])+(jet->imagPart[i])*(jet->imagPart[i]));
        /* Find the angle of the imaginary points */
        if(jet->realPart[i] != 0){
            /* determine which quadrint the point lies in */
            if(jet->realPart[i] >= 0){
                jet->ang[i] = atan(jet->imagPart[i]/jet->realPart[i]);
            }
            else{
                jet->ang[i] = PI + atan(jet->imagPart[i]/jet->realPart[i]);
            }
        }
        else{
	        if(jet->imagPart[i] >= 0){
                jet->ang[i] = PI/2;
            }
            else{
                jet->ang[i] = -PI/2;
            }
        }
    }
}

double convolvePoint(double x, double y, int c, Image *im, Image *mask){
    int i, j;
    double mysum = 0;
    double offsetx = x - mask->width/2.0  - 0.5;
    double offsety = y - mask->height/2.0 - 0.5;

    if(EQUAL_ZERO(offsetx - TRUNC(offsetx),.01) && EQUAL_ZERO(offsetx - TRUNC(offsetx),.01))
    {
        mysum = 0;
        for(j = 0; j < mask->height; j++)
        {
            for(i = 0; i < mask->width; i++)
            {
                mysum += ie(im,i+TRUNC(offsetx),j+TRUNC(offsety),c) * IE(mask,i,j,c);
            }
        }
    }
    else
    {
        mysum = 0;
        for(j = 0; j < mask->height; j++)
        {
            for(i = 0; i < mask->width; i++)
            {
               mysum += interpLinear(im,i+offsetx,j+offsety,c) * IE(mask,i,j,c);
            }
        }
    }
    return mysum;
}

double interpLinear(Image *img, double x, double y, int c)
{
    double xfrac = (x - floor(x));
    double yfrac = (y - floor(y));
    int xLower = INT_FLOOR(x);
    int xUpper = INT_CEIL(x);
    int yLower = INT_FLOOR(y);
    int yUpper = INT_CEIL(y);
    double valUpper, valLower;
    valUpper = ie(img,xLower,yUpper,c)*(1.0-xfrac) + ie(img,xUpper,yUpper,c)*(xfrac);
    valLower = ie(img,xLower,yLower,c)*(1.0-xfrac) + ie(img,xUpper,yLower,c)*(xfrac);
    return valLower*(1.0-yfrac) + valUpper*(yfrac);
}

Image makeGaborMask(double lambda, double theta, double phi,
                    double gamma, double sigma,int size)
{
    Image *mask = new Image;
    mask->height = size;
    mask->width = size;
    mask->data = new double[size*size];
    int i, j;
    for(j = 0; j < size; j++){
        for(i = 0; i < size; i++){
            double x = size/2.0 - size + i;
            double y = size/2.0 - size + j;
            double xp =  x*cos(theta)+y*sin(theta);
            double yp = -x*sin(theta)+y*cos(theta);
            double tmp1 = -(xp*xp+gamma*gamma*yp*yp)/(2*sigma*sigma);
            double tmp2 = (2*PI*xp/lambda)+phi;
            mask->data[j*mask->width + i] = exp(tmp1)*cos(tmp2);
            //cout<<x<<" "<<y<<" "<<tmp1<<" "<<tmp2<<" "<<exp(tmp1)*(cos(tmp2))<<endl;
        }
    }
    /* Rescale the pixel values to have a standard total length */
    ZeroMeanUnitLength(mask);
    return *mask;
}

void ZeroMeanUnitLength(Image *im){
    int i, j;
    double mean = 0; double sqrsum = 0; double invlength = 0;
    for(j = 0; j < im->height; j++){
        for(i = 0; i < im->width; i++){
           {
                mean += IE(im,i,j,1);
            }
        }
    }
    mean = mean / ( im->height*im->width );
    /* printf("mean: %f\n",mean); */
    for(j = 0; j < im->height; j++){
        for(i = 0; i < im->width; i++){
        	{
                IE(im,i,j,1) = IE(im,i,j,1)-mean;
                sqrsum += SQR(IE(im,i,j,1));
            }
        }
    }
    /* printf("sqrsum: %f\n",sqrsum); */
    if(sqrsum != 0){
        invlength = 1.0/sqrt(sqrsum);
    } else {
        invlength = 1.0;
    }
    for(j = 0; j < im->height; j++){
        for(i = 0; i < im->width; i++){
        	{
                IE(im,i,j,1) *= invlength;
            }
        }
    }
}

BunchGraph* addFacetoBunch(BunchGraph* bg, FaceGraph *face)
{
	BunchGraph *newBunch = new BunchGraph;
	newBunch->numFaces = bg->numFaces+1;
	newBunch->faces = new FaceGraph[newBunch->numFaces];
	for(int i=0; i<bg->numFaces; i++)
	{
		newBunch->faces[i].Nodes = new Node[bg->faces[i].numNodes];
		newBunch->faces[i].numNodes = bg->faces[i].numNodes;
		newBunch->faces[i].ID = bg->faces[i].ID;
		for(int j=0; j<bg->faces[i].numNodes; j++)
		{
			newBunch->faces[i].Nodes[j].length = bg->faces[i].Nodes[j].length;
			newBunch->faces[i].Nodes[j].x = bg->faces[i].Nodes[j].x;
			newBunch->faces[i].Nodes[j].y = bg->faces[i].Nodes[j].y;
			newBunch->faces[i].Nodes[j].ang = new double[bg->faces[i].Nodes[j].length];
			newBunch->faces[i].Nodes[j].imagPart = new double[bg->faces[i].Nodes[j].length];
			newBunch->faces[i].Nodes[j].mag = new double[bg->faces[i].Nodes[j].length];
			newBunch->faces[i].Nodes[j].realPart = new double[bg->faces[i].Nodes[j].length];
			for(int k=0; k<bg->faces[i].Nodes[j].length; k++)
			{
				newBunch->faces[i].Nodes[j].ang[k] = bg->faces[i].Nodes[j].ang[k];
				newBunch->faces[i].Nodes[j].imagPart[k] = bg->faces[i].Nodes[j].imagPart[k];
				newBunch->faces[i].Nodes[j].mag[k] = bg->faces[i].Nodes[j].mag[k];
				newBunch->faces[i].Nodes[j].realPart[k] = bg->faces[i].Nodes[j].realPart[k];
			}
		}
	}
	newBunch->faces[bg->numFaces].Nodes = new Node[face->numNodes];
	newBunch->faces[newBunch->numFaces-1].numNodes = face->numNodes;
	newBunch->faces[newBunch->numFaces-1].ID = face->ID;
	for(int i=0; i<face->numNodes; i++)
	{
		newBunch->faces[newBunch->numFaces-1].Nodes[i].length = face->Nodes[i].length;
		newBunch->faces[newBunch->numFaces-1].Nodes[i].x = face->Nodes[i].x;
		newBunch->faces[newBunch->numFaces-1].Nodes[i].y = face->Nodes[i].y;
		newBunch->faces[newBunch->numFaces-1].Nodes[i].ang = new double[face->Nodes[i].length];
		newBunch->faces[newBunch->numFaces-1].Nodes[i].imagPart = new double[face->Nodes[i].length];
		newBunch->faces[newBunch->numFaces-1].Nodes[i].mag = new double[face->Nodes[i].length];
		newBunch->faces[newBunch->numFaces-1].Nodes[i].realPart = new double[face->Nodes[i].length];
		for(int k=0; k<face->Nodes[i].length; k++)
		{
			newBunch->faces[newBunch->numFaces-1].Nodes[i].ang[k] = face->Nodes[i].ang[k];
			newBunch->faces[newBunch->numFaces-1].Nodes[i].imagPart[k] = face->Nodes[i].imagPart[k];
			newBunch->faces[newBunch->numFaces-1].Nodes[i].mag[k] = face->Nodes[i].mag[k];
			newBunch->faces[newBunch->numFaces-1].Nodes[i].realPart[k] = face->Nodes[i].realPart[k];
		}
	}
	return newBunch;
}

void writeBunchGraph(char* fileName, BunchGraph *bg)
{
	ofstream out;
	out.open(fileName);
	out<<bg->numFaces<<endl;
	for(int i=0; i<bg->numFaces; i++)
	{
		out<<bg->faces[i].numNodes<<endl;
		out<<bg->faces[i].ID<<endl;
		for(int j=0; j< bg->faces[i].numNodes; j++)
		{
			out<<bg->faces[i].Nodes[j].x<<" "<<bg->faces[i].Nodes[j].y<<" "<<bg->faces[i].Nodes[j].length<<endl;
			for(int k=0; k<bg->faces[i].Nodes[j].length; k++)
				out<<bg->faces[i].Nodes[j].ang[k]<<" "<<bg->faces[i].Nodes[j].imagPart[k]<<" "<<bg->faces[i].Nodes[j].mag[k]<<" "<<bg->faces[i].Nodes[j].realPart[k]<<endl;
		}
	}
	out.close();
}

BunchGraph* readBunchGraph(char* fileName)
{
	ifstream in;
	BunchGraph *bg = new BunchGraph;
	in.open(fileName);
	in>>bg->numFaces;
	bg->faces = new FaceGraph[bg->numFaces];
	for(int i=0; i<bg->numFaces; i++)
	{
		in>>bg->faces[i].numNodes;
		in>>bg->faces[i].ID;
		bg->faces[i].Nodes = new Node[bg->faces[i].numNodes];
		for(int j=0; j< bg->faces[i].numNodes; j++)
		{
			in>>bg->faces[i].Nodes[j].x;
			in>>bg->faces[i].Nodes[j].y;
			in>>bg->faces[i].Nodes[j].length;
			bg->faces[i].Nodes[j].ang = new double[bg->faces[i].Nodes[j].length];
			bg->faces[i].Nodes[j].imagPart = new double[bg->faces[i].Nodes[j].length];
			bg->faces[i].Nodes[j].mag = new double[bg->faces[i].Nodes[j].length];
			bg->faces[i].Nodes[j].realPart = new double[bg->faces[i].Nodes[j].length];
			for(int k=0; k<bg->faces[i].Nodes[j].length; k++)
				in>>bg->faces[i].Nodes[j].ang[k]>>bg->faces[i].Nodes[j].imagPart[k]>>bg->faces[i].Nodes[j].mag[k]>>bg->faces[i].Nodes[j].realPart[k];
		}
	}
	in.close();
	return bg;
}

Mask* readMask(char* filename)
{
    int maskCount;
    double lambda, angle, phase, gama, sigma; int maskSize;
    int i;
    FILE* file = fopen(filename,"r");
    if(!file){
        printf("Error opening file: %s",filename);
        exit(1);
    }
    /* Read in the number of Masks. */
    fscanf(file, "%d",&maskCount);
    Mask *masks = new Mask;
    masks->numMasks = maskCount;
    masks->masks = new Image[maskCount];
    masks->wavelength = new double[maskCount];
    masks->angle = new double[maskCount];
    masks->phase = new double[maskCount];
    masks->aspect = new double[maskCount];
    masks->radius = new double[maskCount];
    masks->kx = new double[maskCount];
    masks->ky = new double[maskCount];
    masks->size = new int[maskCount];
    /* Read in mask parameters and create masks. */
    for(i = 0; i < maskCount; i++)
    {
        if(fscanf(file, "%lf %lf %lf %lf %lf %d", &lambda, &angle, &phase, &gama, &sigma, &maskSize) != 6)
        {
            printf("Error reading mask file: %s line: %d", filename, i);
            exit(1);
        }
        masks->wavelength[i] = lambda;
        masks->angle[i]      = angle;
        masks->phase[i]      = phase;
        masks->aspect[i]     = gama;
        masks->radius[i]     = sigma;
        masks->kx[i]         = 2.0*PI*cos(angle)/lambda;
        masks->ky[i]         = 2.0*PI*sin(angle)/lambda;
        masks->size[i]       = maskSize;
        masks->masks[i] = makeGaborMask(lambda, angle, phase, gama, sigma, maskSize);
    }
    fclose( file );
    return masks;
}

void writePoints(char* fileName)
{
	int len = strlen(fileName);
	fileName[len-3]= 't';
	fileName[len-2]= 'x';
	fileName[len-1]= 't';
	ofstream out(fileName);
	ifstream in("novelGraph.gph");
	int num, garbage;
	double ignore;
	in>>num;
	out<<num<<endl;
	in>>garbage;
	for(int i=0; i< num; i++)
	{
		double x, y;
		in>>x;
		in>>y;
		out<<x<<" "<<y<<endl;
		int n;
		in>>n;
		for(int j=0; j<n; j++)
		{
			in>>ignore;
			in>>ignore;
			in>>ignore;
			in>>ignore;
		}
	}
	fileName[len-3]= 'p';
	fileName[len-2]= 'g';
	fileName[len-1]= 'm';
	in.close();
}

Image rotate(Image im, double angle)
{
	Image newImage;
	double angle2 = PI/2 - angle;
	double cosine = cos(angle);
	double sine = sin(angle);
	int new_cols = im.height*fabs(cos(angle2)) + im.width * fabs(cosine);
	int new_rows = im.height * fabs(cosine) + im.width * fabs(cos(angle2));
	int centrex = im.width/2;
	int centrey= im.height/2;
	double floatx, floaty;
	int intx, inty;
	int xdiff = (new_cols - im.width)/2;
	int ydiff = (new_rows - im.height)/2;
	newImage.height = new_rows;
	newImage.width = new_cols;
	newImage.data = new double[new_rows*new_cols];
	newImage.maxVal = 255;
	int index = 0, offset=0;
	unsigned char* line_buf = (unsigned char*)malloc(new_cols);
	for(int y = -ydiff; y<new_rows - ydiff; y++)
	{
		index = 0;
		for(int x=-xdiff; x< new_cols - xdiff; x++)
		{
			floatx = (x-centrex) * cosine + (y - centrey)* sine;
			floatx += centrex;
			floaty = (y - centrey)* cosine - (x-centrex) * sine;
			floaty += centrey;
			intx = (int)floatx;
			inty = (int)floaty;
			if((intx<0) ||(intx>=im.width-1) || (inty<0) || (inty >=im.height-1))
				line_buf[index++] = 0;
			else
			{
				float EWweight = floatx - (float)intx;
				float NSweight = floaty - (float)inty;
				unsigned long source_add = (unsigned long) inty * im.width + intx;
				float NW = (float) im.data[source_add];
				float NE = (float) im.data[source_add + 1];
				float SW = (float) im.data[source_add+ im.width];
				float SE = (float) im.data[source_add+ im.width + 1];
				float top = NW + EWweight*(NE-NW);
				float bottom = SW + EWweight*(SE-SW);
				line_buf[index++] = (char)(top + NSweight * (bottom-top));
			}
		}
		for(int i=0; i<new_cols; i++)
			newImage.data[offset + i]=line_buf[i];
		offset+=new_cols;
	}
	return newImage;
}

Image crop(Image im, int xPos)
{
	Image newImage;
	newImage.height = im.height;
	newImage.width = im.width - xPos;
	newImage.maxVal = 255;
	newImage.data = new double[newImage.height*newImage.width];
	for(int i=0; i<im.height; i++)
		memcpy(&newImage.data[i*newImage.width], &im.data[i*im.width + xPos], newImage.width*sizeof(double));
	return newImage;
}

Image flip(Image im)
{
	Image newImage;
	newImage.height = im.height;
	newImage.width = im.width;
	newImage.maxVal = 255;
	newImage.data = new double[newImage.height*newImage.width];
	for(int i=0; i<im.height; i++)
		for(int j=0; j<im.width; j++)
			newImage.data[i*im.width + j] = im.data[i*im.width+ (im.width-j-1)];
	return newImage;
}

Image join(Image im1, Image im2)
{
	Image newImage;
	newImage.height = im1.height;
	newImage.width = im1.width + im2.width;
	newImage.maxVal = 255;
	newImage.data = new double[newImage.width * newImage.height];
	memset(newImage.data, 0, newImage.width * newImage.height*sizeof(double));
	for(int i=0; i<newImage.height; i++)
		{
			memcpy(&newImage.data[i*newImage.width], &im1.data[i*im1.width], im1.width*sizeof(double));
			memcpy(&newImage.data[i*newImage.width + im1.width], &im2.data[i*im2.width], im2.width*sizeof(double));
		}
	return newImage;
}

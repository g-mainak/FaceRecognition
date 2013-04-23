/*
 * main.cpp
 *
 *  Created on: 20-May-2011
 *      Author: Mainak
 */

#include "EBGM.h"
#include "EBGMperipherals.h"
#include "EBGMCalc.h"
#include <math.h>
#include "regression.h"
#include <iostream>
#include <stdlib.h>
using namespace std;
int main(int argc, char *argv[])
{
	/*
	 * Read the IMage.
	 * Find the points.
	 * Find the Nose Bridge, Nose Tip and Centre of lip.
	 * Draw regression line.
	 * Rotate whole picture to make regression line 90 degree.
	 * cut away half picture
	 * duplicate half picture
	 * crop
	 */
	BunchGraph *bg = readBunchGraph((char*)"bunchGraph.bg");
	Mask *mask = readMask((char*)"gbm.wavelet");
	Image image = readImage(argv[1]);
	FaceGraph face = putGraph(bg, mask, image);
	double a = regA(&face.Nodes[2], &face.Nodes[3], &face.Nodes[7]);
	double b = regB(&face.Nodes[2], &face.Nodes[3], &face.Nodes[7]);
	double angle = atan(a);
	cout<<a<<"\t"<<b<<"\t"<<angle*180/PI<<endl;
	Image rotatedImage;
	if(angle>0)
		rotatedImage = rotate(image, -PI/2+angle);
	else
		rotatedImage = rotate(image, PI/2+angle);
	Image croppedImage = crop(rotatedImage, (face.Nodes[2].x+face.Nodes[3].x+face.Nodes[7].x)/3);
	Image flipped = flip(croppedImage);
	writeImage(join( flipped, croppedImage), "a.pgm");
	//writeImage(rotatedImage, "r.pgm");
}

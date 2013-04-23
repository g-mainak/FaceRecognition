#include <fstream>
#include <iostream>
#include <string.h>
#include "EBGMperipherals.h"
using namespace std;
void interpret(char* images)
{
	Image im = readImage(images);
	int n = strlen(images);
	images[n-4]='p';
	images[n-3]='.';
	images[n-2]='p';
	images[n-1]='g';
	images[n]='m';
	images[n+1]= '\0';
	ifstream in;
	in.open("novelGraph.gph");
	int num, len, ix, iy, garbage;
	double x, y;
	double junk;
	in>>num;
	//in>>garbage;
	for(int i=0; i<num ;i++)
	{
		in>>x>>y>>len;
		ix = x;
		iy = y;
		if(x<im.width && y<im.height)
			im.data[(iy*im.width + ix)] = 255;
		for(int j=0; j<4*len; j++)
			in>>junk;
	}
	writeImage(im, images);
	in.close();
}

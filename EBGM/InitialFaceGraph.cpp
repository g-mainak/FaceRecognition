#include <fstream>
#include "EBGMperipherals.h"
#include "EBGMCalc.h"
#include <iostream>
#include <string>
using namespace std;

void initial(char *images, char *txt, int iter)
{
	BunchGraph *bg;
	if(iter == 0)
	{
		bg = new BunchGraph;
		bg->numFaces = 0;
	}
	else
		bg = readBunchGraph((char*)"bunchGraph.bg");
	Mask *masks = readMask((char*)"gbm.wavelet");
	Image image = readImage(images);
	ifstream in;
	in.open(txt);
	FaceGraph *face = new FaceGraph;
	in>>face->numNodes;
	face->Nodes = new Node[face->numNodes];
	for(int i=0; i<face->numNodes; i++)
	{
		in>>face->Nodes[i].x>>face->Nodes[i].y;
		getJetsFromImage(&face->Nodes[i], masks, image);
	}
	bg = addFacetoBunch(bg, face);
	writeBunchGraph((char*)"bunchGraph.bg", bg);
}

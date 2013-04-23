/*
 * start.cpp
 *
 *  Created on: May 20, 2010
 *      Author: Mainak
 */
#include "EBGMperipherals.h"
#include "EBGMCalc.h"
#include <iostream>
using namespace std;

FaceGraph novelFitting(BunchGraph *bg, Mask *mask, Image image)
{
	FaceGraph graph = putGraph(bg, mask, image);
	return graph;
}

void novelFitting2(BunchGraph *bg, Mask *mask, Image image)
{
	FaceGraph graph = putGraph(bg, mask, image);
	writeGraph("novelGraph.gph", graph);
	delete[] graph.Nodes;
}

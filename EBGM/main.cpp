/*
 * main.cpp
 *
 *  Created on: May 25, 2010
 *      Author: Mainak
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include "EBGM.h"
#include "EBGMPeripherals.h"
#include "EBGMCalc.h"
#include <ctime>
#include <cstdlib>
using namespace std;


int main()
{
	int choice, begin, finish, num,  correct=0, incorrect =0, count =0;;
	cout<<"Do you want to train (1) or test (2) or Automate Marking(3)?"<<endl;
	cin>> choice;
	clock_t start, end;
	double runTime;
	start = clock();
	if(choice == 1 )
	{
		ifstream in("train.txt");
		in >> num;
		char garbage[2];
		in.getline(garbage, 1);
		BunchGraph *bg = new BunchGraph;
		Mask *masks = readMask((char*)"gbm.wavelet");
		bg->numFaces = 0;
		for(int i=0; i<num; i++)
		{
			char imageName[80], landMarkName[80];
			in.getline(imageName, 80);	//read ImageName
			in.getline(landMarkName, 80);	//read file that contains the coordinates
			cout<<imageName<<" "<<landMarkName<<endl;
			Image image = readImage(imageName);	//read the image
			ifstream in;
			in.open(landMarkName);
			FaceGraph *face = new FaceGraph;
			in>>face->numNodes;
			face->ID = i+1;	//set the face ID. Currently, it is the name of the folder. for the actual application, some other way to set the ID must be found.
			face->Nodes = new Node[face->numNodes];
			for(int i=0; i<face->numNodes; i++)
			{
				in>>face->Nodes[i].x>>face->Nodes[i].y; //stores the coordinates of the fiducial points.
				getJetsFromImage(&face->Nodes[i], masks, image);	//extracts information from these fiducial points
				face->Nodes[i].x/=image.width;	//Scaling
				face->Nodes[i].y/=image.height;	//Scaling
			}
			bg = addFacetoBunch(bg, face);	//Adds new face to existing bunch
		}
		writeBunchGraph((char*)"bunchGraph.bg", bg); //Writes to file.
	}
	else if (choice == 2)
	{
		cout<<"enter the starting folder:";
		cin >> begin;
		cout<<"enter the ending folder:";
		cin >> finish;
		BunchGraph *bg = readBunchGraph((char*)"bunchGraph.bg");
		Mask *mask = readMask((char*)"gbm.wavelet");
		for(int i=begin; i<=finish; i++)
		{
			for(int j=2; j<10; j++)
			{
				string imageName = "F:/workspace/Pix/"; //This code generates the image path.
				char folder[5], picture[5];
				itoa(i, folder, 10);
				itoa(j, picture, 10);
				imageName += folder ;
				imageName += "/" ;
				imageName += (picture);
				imageName += ".pgm";
				cout<<imageName<<endl;
				FILE* f = fopen( imageName.c_str(), "rb" ); //Check to see if that image exists or not.
				if( f)
				{
					count++;
					Image image = readImage((char*)imageName.c_str()); //read Image from file.
					FaceGraph face = novelFitting(bg, mask, image);	//Generate the FaceGraph
					face.ID = i;	//Get Face ID. This ID does not affect result. It is used only to verify the rank.
					recognize(face, bg, mask, face.ID); //Recognition.
					delete[] face.Nodes;
				}
				else
				{
					cout<<imageName<<" does not exist.\n";
					continue;
				}
			}
		}
	}
	else if (choice==3)
	{//Basically automatic marking also generates a face graph. Just that inistead of doing recognition, it makes a text file.
		ifstream in("mark.txt");
		in >> num;
		char garbage[2], imageName[80];
		Image image;
		string graphName;
		in.getline(garbage, 1);
		BunchGraph *bg = readBunchGraph((char*)"bunchGraph.bg");
		Mask *mask = readMask((char*)"gbm.wavelet");
		for(int i=0; i<num; i++)
		{
			in.getline(imageName, 80);
			cout<<imageName<<endl;
			FILE* f = fopen( imageName, "rb" );
			if(f)
			{
				image = readImage(imageName);
				novelFitting2(bg, mask, image); //same as novel fitting, but it generates a file.
				//writePoints(imageName);	//Makes a text file containing the coordinates of the fiducial points.
				interpret(imageName);	//Generates an image showing teh location of the fiducial points.
			}
		}
	}
	cout<<"Correct: "<<correct<<" Incorrect: "<<incorrect<<endl;
	end = clock();
	runTime = (end - start) / (double) CLOCKS_PER_SEC ;
	printf ("Run time is %g seconds\n", runTime);
	if(count) printf ("Run time per image is %g seconds\n", runTime/count);
	return 0;
}

/*
 * face_recognition.cpp
 *
 *  Created on: 10-May-2010
 *      Author: Sawan Das
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "cv.h"
#include "cvaux.h"
#include "highgui.h"
//" Global variables
int nTrainFaces = 0;     //number of training images
int nEigens             =  0;// number of eigenvalues
IplImage  ** faceImgArr  =  0; //array of face images
CvMat *personNumTruthMat=   0;	// array of person numbers
IplImage *pAvgTrainImg = 0; //the average image
IplImage **eigenVectArr = 0; //eigenvectors
CvMat *eigenValMat = 0;   //eigenvalues
CvMat  *projectedTrainFaceMat = 0;//projected training faces
///Function prototypes
void learn();
void recognize();
void doPCA();
void storeTrainingData();
int loadTrainingData(CvMat **pTrainPersonNumberMat);
int findNearestNeighbor(float *projectedTestFaces);
int loadFaceImgArray(char *filename);
void printUsage();
int main( int argc, char ** argv)
{
	int scanVal;
	printf("Enter value (1 -- learning 2 -- detection) = ");
	scanf("%d", &scanVal);
	if(scanVal == 1)
	{
		learn();
	}
	else
	{
		clock_t start, stop;
		double t = 0.0;
		/* Start timer */
		assert((start = clock())!=-1);
		recognize();
		stop = clock();
		t = (double) (stop-start)/CLOCKS_PER_SEC;
		printf("Run time: %f\n", t);
	}
	/*if(argc != 2)
	{
		printUsage();
		return 1;
	}
	if(!strcmp(argv[1], "train"))
		learn();
	else if(!strcmp(argv[1], "test"))
		recognize();
	else
	{
		printf("Unknown command: %s\n", argv[1]);
		printUsage();
	}*/
	return 0;
}



void printUsage()
{
	printf("Usage: eigenface<conmnd>\n"
			"Valid comnands are \n"
			"train\n"
			"test \n");
}

void learn()
{
	int i;
	nTrainFaces = loadFaceImgArray("train.txt");
	if( nTrainFaces < 2)
	{
		printf("loadFaceImgArray error\n");
		return;
	}
	//do PCA on the training faces
	doPCA() ;
	projectedTrainFaceMat = cvCreateMat(nTrainFaces, nEigens, CV_32FC1);
	for(i = 0; i<nTrainFaces; i++)
	{
		cvEigenDecomposite(
				faceImgArr[i],
				nEigens,
				eigenVectArr,
				0, 0,
				pAvgTrainImg,
				projectedTrainFaceMat->data.fl + i * nEigens);
	}
	// store the recogn:t .on oast,, a: an x:.1. f i:k-
	storeTrainingData();
}

int loadFaceImgArray(char *filename)
{
	FILE *imgListFile = 0;
	char imgFilename [512];
	int iFace, nFaces=0;
	// open the input file
	imgListFile = fopen (filename, "r");
	// count the number of faces
	while( fgets(imgFilename, 512, imgListFile))++nFaces;
	rewind(imgListFile);
	//allocate the face-image array and per&on number matrix
	faceImgArr = (IplImage**)cvAlloc(nFaces *sizeof(IplImage*));
	personNumTruthMat = cvCreateMat(1, nFaces, CV_32SC1);
	// store the face images in an array
	for(iFace = 0; iFace < nFaces; iFace++)
	{
		// read person number and name of image file
		fscanf(imgListFile, "%d %s", personNumTruthMat->data.i + iFace, imgFilename);
		// load the face image
		faceImgArr[iFace] = cvLoadImage(imgFilename, CV_LOAD_IMAGE_GRAYSCALE);
	}
	fclose(imgListFile);
	return nFaces;
}

void doPCA()
{
	int i;
	CvTermCriteria calcLimit;
	CvSize faceImgSize;
	//set the number of eigenvalue& to u.e
	nEigens = nTrainFaces - 1;
	// allocate the eigenvectoz image.
	faceImgSize.width = faceImgArr[0]->width;
	faceImgSize.height = faceImgArr[0]->height;
	eigenVectArr = (IplImage**)cvAlloc(sizeof(IplImage*)* nEigens);
	for(i = 0; i < nEigens; i++)
		eigenVectArr[i] = cvCreateImage(faceImgSize, IPL_DEPTH_32F, 1);
	//allocate the eigenvalue array
	eigenValMat = cvCreateMat(1, nEigens, CV_32FC1);
	// allocate the averaged image
	pAvgTrainImg = cvCreateImage(faceImgSize, IPL_DEPTH_32F, 1);
	// yet the PCA termination criterion
	calcLimit = cvTermCriteria( CV_TERMCRIT_ITER, nEigens,   1);
	///compute average image, eigenvaiuec, and eigenvector3
	cvCalcEigenObjects(
			nTrainFaces,
			(void*)faceImgArr,
			(void*)eigenVectArr,
			CV_EIGOBJ_NO_CALLBACK,
			0,
			0,
			&calcLimit,
			pAvgTrainImg,
			eigenValMat->data.fl);
}


void storeTrainingData()
{
	CvFileStorage * fileStorage;
	int i;
	// create a file-storage interface
	fileStorage = cvOpenFileStorage("facedata.xml", 0, CV_STORAGE_WRITE);
	//store all the data
	cvWriteInt( fileStorage, "nEigens", nEigens);
	cvWriteInt( fileStorage, "nTrainFaces", nTrainFaces) ;
	cvWrite(fileStorage, "trainPersonNumMat", personNumTruthMat, cvAttrList(0, 0));
	cvWrite(fileStorage, "eigenValMat", eigenValMat, cvAttrList(0,0));
	cvWrite(fileStorage, "projectedTrainFaceMat", projectedTrainFaceMat, cvAttrList(0,0));
	cvWrite(fileStorage, "avgTrainImg", pAvgTrainImg, cvAttrList(0,0));
	for(i = 0; i < nEigens; i++)
	{
		char varname[200];
		sprintf(varname, "eigenVect_%d", i);
		cvWrite(fileStorage, varname, eigenVectArr[i], cvAttrList(0,0));
	}
	// release the file-storage interface
	cvReleaseFileStorage( &fileStorage);
}

void recognize()
{

	int i, nTestFaces = 0;   //the number of test images
	CvMat *trainPersonNumMat = 0;// the person numbers during training
	float *projectedTestFace = 0;
	// load test images and ground truth for person number
	nTestFaces = loadFaceImgArray("test.txt");
	printf("%d test faces loaded\n", nTestFaces);
	// load the saved training data
	if(!loadTrainingData(&trainPersonNumMat)) return;
	// project the test images onto the PCA subspace
	projectedTestFace  = (float*)cvAlloc( nEigens * sizeof(float));
	int incorrect=0;
	for(i = 0; i<nTestFaces; i++)
	{
		int iNearest, nearest, truth;
		// project the test image onto the PCA subspace
		cvEigenDecomposite(
				faceImgArr[i],
				nEigens,
				eigenVectArr,
				0, 0,
				pAvgTrainImg,
				projectedTestFace);
		iNearest = findNearestNeighbor(projectedTestFace);
		truth = personNumTruthMat->data.i[i];
		nearest = trainPersonNumMat->data.i[iNearest];
		if (nearest!=truth) incorrect++;
		printf("nearest = %d, Truth = %d\n", nearest, truth);
	}
	printf("%d incorrect Matches!\n", incorrect);
}

int loadTrainingData(CvMat **pTrainPersonNumMat)
{
	CvFileStorage *fileStorage;
	int i;
	//create a tile-storage interface
	fileStorage = cvOpenFileStorage("facedata.xml", 0, CV_STORAGE_READ);
	if(!fileStorage)
	{
		fprintf(stderr, "Can't open facedata.xml\n");
		return 0;
	}
	nEigens = cvReadIntByName(fileStorage, 0, "nEigens", 0);
	nTrainFaces = cvReadIntByName(fileStorage, 0, "nTrainFaces", 0);
	*pTrainPersonNumMat = (CvMat*)cvReadByName(fileStorage, 0, "trainPersonNumMat", 0);

	eigenValMat = (CvMat *)cvReadByName(fileStorage, 0, "eigenValMat", 0);
	projectedTrainFaceMat = (CvMat *)cvReadByName(fileStorage, 0, "projectedTrainFaceMat", 0);
	pAvgTrainImg =(IplImage*)cvReadByName(fileStorage, 0, "avgTrainImg", 0);
	eigenVectArr = (IplImage **)cvAlloc(nTrainFaces * sizeof(IplImage*));
	for(i = 0; i<nEigens; i++)
	{
		char varname[200];
		sprintf( varname, "eigenVect_%d", i);
		eigenVectArr[i] = (IplImage *)cvReadByName(fileStorage, 0, varname, 0);
	}
	//release the file-storage interface
	cvReleaseFileStorage(&fileStorage);
	return 1;


}

int findNearestNeighbor(float *projectedTestFaces)
{

	double leastDistSq = DBL_MAX;
	int i, iTrain, iNearest=0;
	for (iTrain = 0; iTrain < nTrainFaces; iTrain++)
	{
		double distSq = 0;
		for (i = 0; i < nEigens; i++)
		{
			float d_i =
					projectedTestFaces[i]-projectedTrainFaceMat->data.fl[iTrain * nEigens + i];
			distSq += d_i * d_i;
		}
		if(distSq < leastDistSq)
		{
			leastDistSq = distSq;
			iNearest = iTrain;
		}
	}
	return iNearest;
}



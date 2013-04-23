/*
 * Normalize.cpp
 *
 *  Created on: May 28, 2010
 *      Author: Mainak
 */
#include <iostream>

using namespace std;
struct Image
{
	int width, height, maxVal;
	double *data;
};

Image readImage(char* filename){
    int  width, height, max, x, y;
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
    val = fgetc(infile);
    for(y = 0; y < height; y++)
    {
        for(x = 0; x < width; x++)
        {
            im.data[y*width + x] = val;
            val = fgetc(infile);
        }
    }
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
}

int main(int argc, char** argv)
{
	double max = 0, min=256;
	Image im = readImage(argv[1]);
	for(int i=0; i<im.height; i++)
		for(int j=0; j<im.width; j++)
		{
			if(im.data[i*im.width + j] >max) max= im.data[i*im.width + j];
			if(im.data[i*im.width + j] <min) min= im.data[i*im.width + j];
		}
	double scale = 256/(max-min);
	for(int i=0; i<im.height; i++)
		for(int j=0; j<im.width; j++)
			im.data[i*im.width + j] = (im.data[i*im.width + j]-min)*scale;
	writeImage(im, argv[1]);
	return 0;
}

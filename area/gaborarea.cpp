/*
 * gaborarea.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: Mainak
 */
#include "gaborarea.h"


using namespace std;
const int NUMNODES = 15;

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

Image readImage(  char* filename){
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


void interpret(Point* p, char* images)
{
	Image im = readImage(images);
	int n = strlen(images);
	images[n-4]='p';
	images[n-3]='.';
	images[n-2]='p';
	images[n-1]='g';
	images[n]='m';
	images[n+1]= '\0';
	for(int i=0; i<NUMNODES ;i++)
	{
		im.data[(int)(p[i].y)*im.width + (int)(p[i].x)] = 255;
	}
	writeImage(im, images);
}

void writeGraph(char* fileName, FaceGraph graph)
{
	ofstream out;
	out.open(fileName);
	out<<graph.numNodes<<endl;
	for(int i=0; i< graph.numNodes; i++)
	{
		out<<graph.Nodes[i].x<<" "<<graph.Nodes[i].y<<" "<<graph.Nodes[i].length<<endl;
		for(int j=0; j<graph.Nodes[i].length; j++)
			out<<graph.Nodes[i].ang[j]<<" "<<graph.Nodes[i].imagPart[j]<<" "<<graph.Nodes[i].mag[j]<<" "<<graph.Nodes[i].realPart[j]<<endl;
	}
	out.close();

	/*ofstream o;
	o.open("sid.txt", ios::app);
	for(int i=0; i< graph.numNodes; i++)
	{
		for(int j=0; j<graph.Nodes[i].length; j++)
			o<<graph.Nodes[i].imagPart[j]<<" "<<graph.Nodes[i].realPart[j]<<" ";
	}
	o<<endl<<endl;
	o.close();*/
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

double SIM_DISPLACE(Node* j1, Node*  j2, double dx, double dy, Mask* masks) {
    double j12 = 0.0;
    double j11 = 0.0;
    double j22 = 0.0;
    int i;
    for(i = 0; i < j1->length; i++){
        j12 += j1->mag[i]*j2->mag[i]*
	       cos(j1->ang[i] - j2->ang[i] -
	       (dx * masks->kx[2*i] + dy * masks->ky[2*i]));
        j11 += j1->mag[i]*j1->mag[i];
        j22 += j2->mag[i]*j2->mag[i];
    }
    return j12/sqrt(j11*j22);
}

double dispEst(Node* j1, Node* j2, double *tdx, double *tdy, Mask *masks){
    double sim  = 0.0;
    double dx   = 0.0,
          dy   = 0.0;
    int change = 1,
        iter   = 50;
    double nextx, nexty;
    double bestsim = 0.0;
    double tol  = .2;
    double step = 1.0;
//NARROWING SEARCH
    assert(j1 && j1->length  && j2 && j1->length==j2->length);

    sim = SIM_DISPLACE(j1, j2, dx,dy, masks);
    bestsim = sim;
    *tdx = 0.0;
    *tdy = 0.0;

    change = 1;
    sim = SIM_DISPLACE(j1, j2, dx,dy, masks);
    bestsim = sim;

    nextx = dx;
    nexty = dy;
    while(change && iter--){
        dx = nextx;
        dy = nexty;
        change = 0;

        sim = SIM_DISPLACE(j1, j2, dx+step,dy+step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx+step;
            nexty = dy+step;
            change = 1;
        }

        sim = SIM_DISPLACE(j1, j2, dx+step,dy-step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx+step;
            nexty = dy-step;
            change = 1;
        }

        sim = SIM_DISPLACE(j1, j2, dx-step,dy-step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx-step;
            nexty = dy-step;
            change = 1;
        }

        sim = SIM_DISPLACE(j1, j2, dx-step,dy+step, masks);
        if( bestsim < sim ){
            bestsim = sim;
            nextx = dx-step;
            nexty = dy+step;
            change = 1;
        }
        dx = nextx;
        dy = nexty;

        if(change == 0 && step > tol){
            change = 1;
            step = step * 0.5;
        }
    }

    *tdx = dx;
    *tdy = dy;

    return bestsim;
}


void getJetsFromImage(Node *n, Mask *masks, Image image)
{
	if( masks->numMasks )
	{
		/* This algrithm aproximates the gabor parameters at a location in continuous space
	    By computing the gabor convoluation at an integer location close to the contiuous
		point and then shifting the phase parameters the apropreate amount. It is much faster
		then using linear interpolation and it is probably more accurate. */
		int i;
	    float rx = n->x;
	    float ry = n->y;
	    float dx = n->x - rx;
	    float dy = n->y - ry;
	    n->length= masks->numMasks/2;
	    n->ang = new double[masks->numMasks/2];
	    n->imagPart = new double[masks->numMasks/2];
	    n->mag = new double[masks->numMasks/2];
	    n->realPart = new double[masks->numMasks/2];
	    /* Compute the jet coeffecents */
	    for( i = 0; i < n->length; i++)
	    {
	    	n->realPart[i] = convolvePoint(rx, ry, 0, &image, &masks->masks[2*i]);
	        n->imagPart[i] = convolvePoint(rx, ry, 0, &image, &masks->masks[2*i+1]);
	        //FINITE(jet->realPart[i]);
	        //FINITE(jet->imagPart[i]);
	    }
	    computePolar(n);
	    /* Adjust phase coordinates */
	    for( i = 0; i < n->length; i++)
	    {
	    	n->ang[i] +=  dx * masks->kx[2*i] + dy * masks->ky[2*i];
	    }
	    /* Recompute coeffecents based on phase change */
	    for( i = 0; i < n->length; i++)
	    {
	    	n->realPart[i] = n->mag[i] * cos(n->ang[i]);
	        n->imagPart[i] = n->mag[i] * sin(n->ang[i]);
	    }
	}
	else
	{
	  	n->ang = new double[0];
	   	n->imagPart = new double[0];
	    n->mag = new double[0];
	    n->realPart = new double[0];
	}
}


Node findClosestTemplate(BunchGraph *bg, int nodeNum, Mask* masks, Image image)
{
	double scalex = image.width;
	double scaley= image.height;
	int best = 0;
	double dx = 0.0, dy = 0.0;
    double bestsim = -1.0e300;
	Node *temp = new Node;
	for(int i=0; i<5; i++)
	{
		temp->x = bg->faces[i].Nodes[nodeNum].x*scalex;
		temp->y = bg->faces[i].Nodes[nodeNum].y*scaley;
		getJetsFromImage(temp, masks, image);
		double sim = dispEst( &bg->faces[i].Nodes[nodeNum], temp, &dx, &dy, masks);
		if(sim > bestsim)
        {
            bestsim = sim;
            best = i;
        }
	}
	dispEst( &bg->faces[best].Nodes[nodeNum], temp, &dx, &dy, masks);
    temp->x += dx;
    temp->y += dy;
    getJetsFromImage(temp, masks, image);
    return *temp;
}

FaceGraph putGraph(BunchGraph *bg, Mask* masks, Image image)
{
	FaceGraph fg;
	fg.numNodes = NUMNODES;
	fg.Nodes = new Node[NUMNODES];
	for(int i=0; i<NUMNODES; i++)
	{
		fg.Nodes[i] = findClosestTemplate(bg, i, masks, image);
	}
	return fg;
}

Point* novelFitting(BunchGraph *bg, Image image, Mask *mask)
{
	FaceGraph graph = putGraph(bg, mask, image);
	Point *p = new Point[NUMNODES];
	for(int i=0; i<graph.numNodes; i++)
	{
		p[i].x = graph.Nodes[i].x;
		p[i].y = graph.Nodes[i].y;
	}
	delete[] graph.Nodes;
	return p;
}


void findArea(Point* p)
{
	FILE* out1 = fopen("Areas.txt", "a");
	double eyes = SQR(p[0].x - p[1].x) + SQR(p[0].y - p[1].y);
	double nose = SQR(p[2].x - p[3].x) + SQR(p[2].y - p[3].y);
	double leyen = SQR(p[0].x - p[2].x) + SQR(p[0].y - p[2].y);
	double reyen = SQR(p[1].x - p[2].x) + SQR(p[1].y - p[2].y);
	double noseTip = SQR(p[4].x - p[5].x) + SQR(p[4].y - p[5].y);
	double lentip = SQR(p[0].x - p[3].x) + SQR(p[0].y - p[3].y);
	double rentip = SQR(p[1].x - p[3].x) + SQR(p[1].y - p[3].y);
	double ntip2mouth = SQR(( p[8].x - p[6].x )*( p[ 3 ].y - p[6].y ) - ( p[8].y - p[6].y )*( p[ 3 ].x - p[6].x )) / ( SQR ( p[8].x - p[6].x) + SQR( p[8].y - p[6].y ) );
	double totalNorm = eyes + leyen + reyen;
	double total =  nose + ntip2mouth + noseTip;
	double totalEN = lentip + rentip + nose + eyes;
	fprintf(out1, "%f %f %f %f %f %f %f %f\n", nose/total, noseTip/total, ntip2mouth/total, eyes/totalNorm, leyen/totalNorm, reyen/totalNorm, lentip/totalEN, rentip/totalEN);
}

int main(int argc, char** argv)
{
	Mask *mask = readMask("gbm.wavelet");
	BunchGraph *bg = readBunchGraph("bunchGraph.bg");
	for( int i = 44 ; i <= (200) ; i++ )
	{
		char cc[10];
		itoa(i,cc,10);
		string str = cc;
		int j;
		for(j=1; j>0; j++)
		{
			string c1;
			char jay[6];
			itoa(j, jay, 10);
			c1 = "F:/workspace/Pix/" + str + "/" + jay + ".pgm";
			FILE* f = fopen( c1.c_str(), "r" );
			if(!f)
				break;
			cout<<c1.c_str()<<endl;
			Image image = readImage((char*)c1.c_str());
			Point *p = novelFitting(bg, image, mask);
			interpret(p, (char*)c1.c_str());
			findArea(p);
			delete p;
		}
	}
	out1.close();
	return 0;
}

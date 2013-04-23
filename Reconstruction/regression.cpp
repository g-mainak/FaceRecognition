/*
 * regression.c
 *
 *  Created on: 20-May-2011
 *      Author: Mainak
 */
#include "regression.h"

double regA(Node* a, Node* b, Node* c)
{
	double sigmaxy = a->x*a->y +b->x*b->y +c->x*c->y;
	double sigmax = a->x +b->x +c->x;
	double sigmay = a->y +b->y +c->y;
	double sigmaxx = a->x*a->x +b->x*b->x +c->x*c->x;
	double sigmayy = a->y*a->y +b->y*b->y +c->y*c->y;
	return (3*sigmaxy - sigmax*sigmay)/(3*sigmaxx - sigmax*sigmax);
}

double regB(Node* a, Node* b, Node* c)
{
	double sigmaxy = a->x*a->y +b->x*b->y +c->x*c->y;
	double sigmax = a->x +b->x +c->x;
	double sigmay = a->y +b->y +c->y;
	double sigmaxx = a->x*a->x +b->x*b->x +c->x*c->x;
	double sigmayy = a->y*a->y +b->y*b->y +c->y*c->y;
	return (sigmay*sigmax*sigmax - sigmax*sigmaxy)/(3*sigmaxx - sigmax*sigmax);
}

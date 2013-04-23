/*
 * NovelFit.h
 *
 *  Created on: May 25, 2010
 *      Author: Mainak
 */

#ifndef NOVELFIT_H_
#define NOVELFIT_H_
#include "EBGMperipherals.h"

FaceGraph novelFitting(BunchGraph* , Mask* , Image );
void novelFitting2(BunchGraph* , Mask* , Image );
void initial(char*, char*, int);
void interpret(char*);
void recognize( FaceGraph, BunchGraph*, Mask*, int );

#endif /* NOVELFIT_H_ */

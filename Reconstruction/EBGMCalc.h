/*
 * EBGMCalc.h
 *
 *  Created on: May 25, 2010
 *      Author: Mainak
 */

#ifndef EBGMCALC_H_
#define EBGMCALC_H_


double dispEst(Node*, Node*, double, double, Mask* );
void getJetsFromImage(Node*, Mask*, Image);
void getFinalLocation(BunchGraph*, Node*, int, Image, Mask*);
void guessLocation(BunchGraph*, Node*, int );
Node getAverageNode(BunchGraph*, int);
FaceGraph putGraph(BunchGraph*, Mask*, Image);

#endif /* EBGMCALC_H_ */

/*
 * train.cpp
 *
 *  Created on: Jun 7, 2010
 *      Author: Mainak
 */
#include <iostream>
#include <fstream>
using namespace std;
int num = 40;
int main(int argc, char** argv)
{
	ofstream out("F://workspace//ebgm//debug//mark.txt");
	for(int i=1; i<=200; i++)
		for(int j=1; j<=10; j++)
		{
			out<<"F:/workspace/Pix/" << i<< "/" <<j << ".pgm"<<endl;
			//out<<"F:/workspace/Pix/" << i<< "/" <<j << ".txt"<<endl;
		}
	return 0;
}

/*
 * test.cpp
 *
 *  Created on: Jun 8, 2010
 *      Author: Mainak
 */
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
	ifstream in;
	in.open("F:\\workspace\\EBGM\\Debug\\count.txt");
	ofstream out ;
	out.open("out.txt");
	int arr[200];
	for(int i=0; i<200 ;i++)
		arr[i]=0;
	int num=0;
	for(int i=0; i<329; i++)
	{
		in>> num;
		arr[num]++;
		num =0;
	}
	for(int i=0; i<200 ;i++)
	{
		cout<<"rank "<<(i)<<" = "<<arr[i]<<endl;
		out<<arr[i]<<endl;
	}
	return 0;
}

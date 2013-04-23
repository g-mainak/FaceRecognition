/*
 * php.cpp
 *
 *  Created on: 11-Jun-2011
 *      Author: Mainak
 */
#include <fstream>
using namespace std;
int main()
{
	ofstream out ;
	out.open("php.txt");
	char* words[] = {"school", "entryDate", "leavingDate", "board", "degree", "distinction", "subject", "percent", "cpi", "maxCpi", "pass", "organization", "position", "work", "joinDate", "exitDate", "lastPay", "experience"};
	for(int i=0; i<18; i++){
		for(int j=1; j<=5; j++)
		out<<words[i]<<j<<"= '$"<<words[i]<<j<<"', ";
		out<<endl;
	}
}

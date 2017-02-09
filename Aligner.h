#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <stdlib.h>


using namespace std;


class SubMtrx{
	map <string,int> mapping;
	string aa;
	float avg_rate = 0;
	
public:
	SubMtrx(string file="BLOSUM62");
	
	float score(string key);
	
	float score(char c1, char c2);
	
	float get_rate();

	char best_match(char c);

	char best_match(char c1, float r1, char c2, float r2);
};


float inferAln(vector<string> &aln, string S0, float d1, string S1, float d2, SubMtrx m, float idrate);
float pairAln(vector<string> &aln, string S0, string S1, SubMtrx m, float idrate);

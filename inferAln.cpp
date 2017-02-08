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
	SubMtrx(string file){
		ifstream fin;
		fin.open(file);
		string line;
		getline(fin,line);
		for (int i=0; i<line.length(); i++){
			if (line[i] != ' ')
				aa.push_back(line[i]);
		}
		avg_rate = 0;
		int count = 0;
		while (!fin.eof()){
			char a;
			fin >> a;
			int score;
			for (int i =0; i<aa.length(); i++){
				fin >> score;
				string key = "xx";
				key[0] = a;
				key[1] = aa[i];
				mapping[key] = score;
				avg_rate += pow(float(2),score);
				count++;
			}
		}
		avg_rate /= count;
	}
	float score(string key){
		return mapping[key];
	}

	float score(char c1, char c2){
		string key = "xx";
		key[0] = c1;
		key[1] = c2;
		return mapping[key];
	}
	
	float get_rate(){
		return this->avg_rate;
	}
	char best_match(char c){
		char best = this->aa[0];
		float best_score = this->score(c,best);

		for (int i = 1; i < this->aa.length(); i++){
			float curr_score = this->score(c,aa[i]);
			if (curr_score > best_score || (curr_score == best_score && aa[i] == c)){
				best_score = curr_score;
				best = aa[i];
			}
		}
		return best;
	}
	char best_match(char c1, float r1, char c2, float r2){
		char best = this->aa[0];
		float best_score = r1*this->score(c1,best) + r2*this->score(c2,best);
		
		for (int i = 1; i < this->aa.length(); i++){
			float curr_score = r1*this->score(c1,aa[i]) + r2*this->score(c2,aa[i]);
			if (curr_score > best_score || (curr_score == best_score && (aa[i] == c1|| aa[i] == c2))){
				best_score = curr_score;
				best = aa[i];
			}
		}
		return best;
	}
};

float inferAln(string* &aln, string S0, float d1, string S1, float d2, SubMtrx m, float idrate){
	int alpha = floor(log2(m.get_rate()*idrate));
	float d = d1+d2;
	float r1 = d/d1, r2 = d/d2;
	
	int **scoring = new int*[S1.length()+1];
	char **backtrack = new char*[S1.length()+1];
	char **infer = new char*[S1.length()+1];

	for (int i = 0; i < S1.length()+1; i++){
		scoring[i] = new int[S0.length()+1];
		backtrack[i] = new char[S0.length()+1];
		infer[i] = new char[S0.length()+1];
	}

	scoring[0][0] = 0;
	backtrack[0][0] = '$';
	infer[0][0] = '$';
	
	for (int j = 1; j < S0.length()+1; j++){
		char x = m.best_match(S0[j-1]);
		infer[0][j] = x;
		scoring[0][j] = r1*m.score(S0[j-1],x) + r2*alpha + alpha;
		backtrack[0][j] = 'L';
	}

	for (int i = 1; i < S1.length()+1; i++){
		char x = m.best_match(S1[i-1]);
		infer[i][0] = x;
		scoring[i][0] = r1*alpha + r2*m.score(S1[i-1],x) + alpha;
		backtrack[i][0] = 'U';
	}
	
	for (int i = 1; i < S1.length()+1; i++){
		for (int j = 1; j < S0.length()+1; j++){
		// match
			char x_match = m.best_match(S0[j-1],r1,S1[i-1],r2);
			float score_match = scoring[i-1][j-1] + r1*m.score(S0[j-1],x_match) + r2*m.score(S1[i-1],x_match) + m.score(S0[j-1],S1[i-1]);
		// gap in S0
			char x_g0;
			float score_g0;
			char x = m.best_match(S1[i-1]);
			float sm = r1*alpha + r2*m.score(x,S1[i-1]) + alpha;
			float sg = r2*alpha + alpha;
			if (sm > sg){
				x_g0 = x;
				score_g0 = scoring[i-1][j] + sm;
			} else {
				x_g0 = '-';
				score_g0 = scoring[i-1][j] + sg;
			}
				
		// gap in S1				
			char x_g1;
			float score_g1;
			x = m.best_match(S0[j-1]);
			sm = r2*alpha + r1*m.score(x,S0[j-1]) + alpha;
			sg = r1*alpha + alpha;
			if (sm > sg){
				x_g1 = x;
				score_g1 = scoring[i][j-1] + sm;
			} else {
				x_g1 = '-';
				score_g1 = scoring[i][j-1] + sg;
			}
		
		// find the best of all
			if (score_match > score_g0 && score_match > score_g1){
				scoring[i][j] = score_match;
				backtrack[i][j] = 'D'; // diagonal
				infer[i][j] = x_match;
			}
			else if (score_g0 > score_match && score_g0 > score_g1){
				scoring[i][j] = score_g0;
				backtrack[i][j] = 'U'; // up
				infer[i][j] = x_g0;
			}
			else {
				scoring[i][j] = score_g1;
				backtrack[i][j] = 'L'; // left
				infer[i][j] = x_g1;
			}
		}
	}
	// backtrack
	int i = S1.length(), j = S0.length();
	aln = new string[3];
	
	while (i > 0 || j > 0){
		aln[2] = infer[i][j] + aln[2];
		switch (backtrack[i][j]){
			case 'D':
				aln[0] = S0[j-1] + aln[0];
				aln[1] = S1[i-1] + aln[1];
				i--;
				j--;
				break;
			case 'L':
				aln[0] = S0[j-1] + aln[0];
				aln[1] = '-' + aln[1];
				j--;
				break;
			case 'U':
				aln[0] = '-' + aln[0];
				aln[1] = S1[i-1] + aln[1];
				i--;
				break;
		}	
	}
	return scoring[S1.length()][S0.length()];
}

int main(int argc, char* argv[]){
	SubMtrx M("BLOSUM62");
	//cout << M.score('B','B');
	string S0 = "PMIDAVPEGNDMTAKNKDLYIRSILGDKFDQQAALSKLDDAYLLECMYIANIKTVSILDNGLPSVTQPTYVSADTRRLKIRIAEWVDWNMSEDLIIYICYRLPSDFTVSAIKSIIAIAYSNRQADLAEAQAVLTKMCRGRFVFEGVKAKENSAQELAARTRKLKCNGRVEYAATFFMGVALKAVTSAQPGQIRGAKLD";
	string S1 = "PLSDARPEINEMTAKNKEVYWRSILGDKFDQQAALSKLDEAYLIEGMYIANIETQHVSILDNGLPIVVQPIYVDADDRRLKIRIASWIDWYISEDLIIYICVRLPSDFNISAIKSMAQKQLVLAKMCRGRFVFEGVKDKENVAQELPARTRDLECNGRNEYPTTFFLGVAIKAITTAQPTQMKGAKLD";
	
	string *aln;
	cout << inferAln(aln,S0,2,S1,1,M,0.01) << endl;
	
	cout << S0 << endl << S1 << endl << endl;
	cout << aln[2] << endl << aln[0] << endl << aln[1] << endl;

	delete [] aln;

	return 0;
}


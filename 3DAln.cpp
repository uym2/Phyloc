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
	float avg_rate = 0;
	
public:
	SubMtrx(string file){
		ifstream fin;
		fin.open(file);
		string line;
		string aa;
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
		//cout << avg_rate << endl;
	}
	int score(string key){
		return mapping[key];
	}
	int sub_rate(){
		return this->avg_rate;
	}
};


string* Align(string S1, string S2, string S3, SubMtrx m,float indel_rate){
	int alpha = floor(log2(m.sub_rate()*indel_rate));
	//cout << alpha << endl;
	string *aln = new string[3];
	int ***scoring = new int**[S1.length()+1];
	int ***backtrack = new int**[S1.length()+1];

	for (int i = 0; i < S1.length()+1; i++){
		scoring[i] = new int* [S2.length()+1];
		backtrack[i] = new int* [S2.length()+1];
		for (int j = 0; j < S2.length()+1; j++){
			scoring[i][j] = new int [S3.length()+1];
			backtrack[i][j] = new int [S3.length()+1];
			for (int k = 0; k < S3.length()+1; k++){
				scoring[i][j][k] = 0;
				backtrack[i][j][k] = -1;
			}
		}
	}

	for (int i = 0; i < S1.length()+1; i++){
		for (int j = 0; j < S2.length()+1; j++){
			for (int k = 0; k < S3.length()+1; k++){
				vector<int> aln_scores;
				int aln_all = -10000, aln_12 = -10000,aln_13 = -10000,aln_23 = -10000,
					aln_1 = -10000, aln_2 = -10000,aln_3 = -10000;
				string str="";
				if (i>0 && j>0 && k>0){
					aln_all = scoring[i-1][j-1][k-1] + m.score(str+S1[i-1]+S2[j-1])
							+ m.score(str+S1[i-1]+S3[k-1]) + m.score(str+S2[j-1]+S3[k-1]);
					aln_scores.push_back(aln_all);
				}
				if (i>0 && j>0){
					aln_12 =  scoring[i-1][j-1][k] + 2*alpha + m.score(str+S1[i-1]+S2[j-1]);
					aln_scores.push_back(aln_12);
				} 
				if (i>0 && k>0){
					aln_13 =  scoring[i-1][j][k-1] + 2*alpha + m.score(str+S1[i-1]+S3[k-1]);
					aln_scores.push_back(aln_13);
				}
				if (j>0 && k>0){
					aln_23 =  scoring[i][j-1][k-1] + 2*alpha + m.score(str+S2[j-1]+S3[k-1]);
					aln_scores.push_back(aln_23);
				}
				if (i>0){
					aln_1 =  scoring[i-1][j][k] + 2*alpha;// + m.score(str+S2[j]+S3[k]);
					aln_scores.push_back(aln_1);
				} 
				if (j>0){
					aln_2 =  scoring[i][j-1][k] + 2*alpha;// + m.score(str+S1[i]+S3[k]);
					aln_scores.push_back(aln_2);
				}
				if (k>0){
					aln_3 =  scoring[i][j][k-1] + 2*alpha;// + m.score(str+S1[i]+S2[j]);
					aln_scores.push_back(aln_3);
				}
				if (aln_scores.size()>0){
					int max = aln_scores[0];
					for (int p = 1; p < aln_scores.size(); p++){
						if (max < aln_scores[p])
							max = aln_scores[p];
					}
					scoring[i][j][k] = max;
					
					if (max == aln_all)
							backtrack[i][j][k] = 0;	
					else if (max == aln_1)
							backtrack[i][j][k] = 1;
					else if (max == aln_2)
							backtrack[i][j][k] = 2;
					else if (max == aln_3)
							backtrack[i][j][k] = 3;
					else if (max == aln_12)
							backtrack[i][j][k] = 12;
					else if (max == aln_13)
							backtrack[i][j][k] = 13;
					else if (max == aln_23)
							backtrack[i][j][k] = 23;
				}
			}
		}
	}
	
	//cout << scoring[S1.length()][S2.length()][S3.length()] << endl;

	// backtracking ...
	string str;
	int i = S1.length(), j = S2.length(), k = S3.length();
	while (i>0 || j>0 || k>0){
		switch (backtrack[i][j][k]){
			case 0:
				//cout << 0 << endl;
				aln[0] = S1[i-1]+aln[0];
				aln[1] = S2[j-1]+aln[1];
				aln[2] = S3[k-1]+aln[2];
				//cout << S1[i-1] << " " << S2[j-1] << " " << m.score(str+S1[i-1]+S2[j-1]) << endl;
				//cout << S1[i-1] << " " << S3[k-1] << " " << m.score(str+S1[i-1]+S3[k-1]) << endl;
				//cout << S2[j-1] << " " << S3[k-1] << " " << m.score(str+S2[j-1]+S3[k-1]) << endl;
				i--;
				j--;
				k--;
				break;
			case 23:
				aln[0] = '-'+aln[0];
				aln[1] = S2[j-1]+aln[1];
				aln[2] = S3[k-1]+aln[2];
				//cout << "23" << endl; 
				//cout << S2[j-1] << " " << S3[k-1] << " " << m.score(str+S2[j-1]+S3[k-1]) << endl;
				j--;
				k--;
				break;
			case 13:
				//cout << 13 << endl;
				//cout << S1[i-1] << " " << S3[k-1] << " " << m.score(str+S1[i-1]+S3[k-1]) << endl;
				aln[1] = '-'+aln[1];
				aln[0] = S1[i-1]+aln[0];
				aln[2] = S3[k-1]+aln[2]; 
				i--;
				k--;
				break;
			case 12:
				//cout << 12 << endl;
				//cout << S1[i-1] << " " << S2[j-1] << " " << m.score(str+S1[i-1]+S2[j-1]) << endl;
				aln[2] = '-'+aln[2];
				aln[1] = S2[j-1]+aln[1];
				aln[0] = S1[i-1]+aln[0]; 
				i--;
				j--;
				break;
			case 3:
				//cout << 3 << endl;
				aln[0] = '-'+aln[0];				
				aln[1] = '-'+aln[1];
				aln[2] = S3[k-1]+aln[2];
				k--;
				break;
			case 2:
				//cout << 2 << endl;
				aln[0] = '-'+aln[0];
				aln[2] = '-'+aln[2];
				aln[1] = S2[j-1]+aln[1];
				j--;
				break;
			case 1:
				//cout << 1 << endl;
				aln[1] = '-'+aln[1];
				aln[2] = '-'+aln[2];
				aln[0] = S1[i-1]+aln[0];
				i--;
				break;
		}
		//cout << aln[0] << endl << aln[1] << endl << aln[2] << endl ;
		
		//cout << scoring[i][j][k] << endl << endl;
	}
	

	//aln[0] = S1;
	//aln[1] = S2;
	//aln[2] = S3;

	for (int i = 0; i < S1.length()+1; i++){
		for (int j = 0; j < S2.length()+1; j++)
			delete[] scoring[i][j];
		delete [] scoring[i];
	}
	delete [] scoring;

	return aln;
}  
int main(int argc, char* argv[]){
	cout << "Reading input ..." << endl;
	string fasta = argv[1];
	SubMtrx M("BLOSUM62");
	ifstream fin;
	fin.open(fasta);
	string name1,name2,name3;
	string seq1,seq2,seq3;
	fin >> name1 >> seq1 >> name2 >> seq2 >> name3 >> seq3;
	fin.close();

	ofstream fout;
	string *aln;
	float indel_rate = stof(argv[3]);
	cout << "Aligning sequences ... " << endl;
	aln = Align(seq1,seq2,seq3,M,indel_rate);
	
	cout << "Writing output ..." << endl;
	fout.open(argv[2]);
	fout << name1 << "\n" << aln[0] << "\n" << name2 << "\n" << aln[1] << "\n" << name3 << "\n" << aln[2] << endl;
	fout.close();
	
	cout << "Alignment written to " << argv[2] << endl;
	delete [] aln;

	return 0;
}


#include "Aligner.h"

SubMtrx::SubMtrx(string file){
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

float SubMtrx::score(string key){
		return mapping[key];
	}

float SubMtrx::score(char c1, char c2){
		string key = "xx";
		key[0] = c1;
		key[1] = c2;
		return mapping[key];
	}
	
float SubMtrx::get_rate(){
		return this->avg_rate;
	}

char SubMtrx::best_match(char c){
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

char SubMtrx::best_match(char c1, float r1, char c2, float r2){
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


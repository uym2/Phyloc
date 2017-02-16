#include "AlnTree.h"
#include <iostream>
#include <fstream>
#include <map>
#include <string>

using namespace std;

bool AlnTree::mapSeq2Tree(string seqfile){
	map<string,string> seq_map;
	string str;
	ifstream fin;
	fin.open(seqfile);
	while (fin >> str){
		if (str[0] == '>'){
			string seq;
			fin >> seq;
			seq_map[str.substr(1)] = seq;
			//cout << str << endl;
			//cout << seq << endl << endl;
		}
	}
	fin.close();

	return this->root->mapSeq(seq_map);
}

bool AlnTree::btmupAln(){
	this->root->btmupAln(this->m,this->idrate);
	return true;
}

bool AlnNode::mapSeq(map<string,string> seq_map){
	if (this->is_leaf()){
		if (!seq_map.count(this->get_label())){
			cout << this->get_label() << " not found!" << endl;
			return false;
		}
		this->mySeq = seq_map[this->get_label()];
	} else {
		for (int i = 0; i<this->nchild(); i++){
			if (!this->get_child(i)->mapSeq(seq_map))
				return false;
		}
	}
	return true;
}

bool AlnNode::set_seq(){
	string seq;
	for (int i = 0; i < this->Aln.back().length(); i++){
		if (Aln.back()[i] != '-')
			seq.push_back(Aln.back()[i]);
	}
	this->mySeq = seq;
	return true;
}

string AlnNode::get_seq(){
	return this->mySeq;
}

float AlnNode::btmupAln(SubMtrx m,float idrate){
	float score = 0;
	if (this->is_leaf()){
		this->Aln.push_back(this->get_seq());
		this->leafLabel.push_back(this->get_label());
	}
	else{ // go down until we reach a leaf
		for (int i = 0; i < this->nchild(); i++)
			this->get_child(i)->btmupAln(m,idrate);

		// now we come back 
		if (this->nchild() == 1){
			this->Aln.push_back(this->get_child(0)->get_seq());
			this->Aln.push_back(this->get_child(0)->get_seq());
		}
		else {
			string S0 = this->get_child(0)->get_seq();
			string S1 = this->get_child(1)->get_seq();

			float d1 = this->get_child(0)->get_edge();
			float d2 = this->get_child(1)->get_edge();
			
			score = inferAln(this->Aln,S0,d1,S1,d2,m,idrate);	
				
			this->set_seq();
			
			//cout << this->get_child(0)->getAln().size() << endl;
			this->leafLabel = this->get_child(0)->getLeafLabel();
			this->leafLabel.insert(this->leafLabel.end(),this->get_child(1)->getLeafLabel().begin(),this->get_child(1)->getLeafLabel().end());
			if (!this->get_child(0)->is_leaf()){
				this->get_child(0)->transitAln(this->Aln,0);
				this->Aln.erase(this->Aln.begin());
				this->Aln.insert(this->Aln.begin(),this->get_child(0)->getAln().begin(),this->get_child(0)->getAln().end()-1);
				//cout << this->get_child(0)->getLeafLabel()[0] << endl;
			}
			//cout << this->get_child(0)->getAln().size() << endl;

			if (!this->get_child(1)->is_leaf()){
				this->get_child(1)->transitAln(this->Aln,this->Aln.size()-2);
				this->Aln.erase(this->Aln.end()-2);
				this->Aln.insert(this->Aln.end()-1,this->get_child(1)->getAln().begin(),this->get_child(1)->getAln().end()-1);
			}
			
			//cout << this->leafLabel.size() << endl;
			//cout << this->Aln.size() << endl;

			if (this->nchild() > 2)
				this->polyAln(m,idrate);		

			//this->printAln();
			//cout << endl << endl;
			
		}
	}
	return score;
}

bool AlnNode::polyAln(SubMtrx m, float idrate){
	if (this->nchild() <= 2)
		return false;
	for (int i = 2; i < this->nchild(); i++){
		vector<string> newAln;
		pairAln(newAln,this->get_child(i)->get_seq(),this->get_seq(),m,idrate); // notice: order matters here
		this->transitAln(newAln,1);
		this->add_aln(newAln[0]);			
		this->get_child(i)->transitAln(this->Aln,this->Aln.size()-2);

		if (!this->get_child(i)->is_leaf()){
			this->Aln.erase(this->Aln.end()-2);
			this->Aln.insert(this->Aln.end()-1,this->get_child(i)->getAln().begin(),this->get_child(i)->getAln().end()-1);
		}
		this->leafLabel.insert(this->leafLabel.end(),this->get_child(i)->getLeafLabel().begin(),this->get_child(i)->getLeafLabel().end()); 
	}
	return true;
}

bool AlnNode::add_gaps(int pos){
	for (int i = 0; i < this->Aln.size(); i++)
		this->Aln[i].insert(this->Aln[i].begin()+pos,'-');
	return true;
}

bool insert_gaps(vector<string> &aln, int pos){
	for (int i = 0; i < aln.size(); i++)
		aln[i].insert(aln[i].begin()+pos,'-');
	return true;
}

bool AlnNode::transitAln(vector<string> &refAln,int shareIdx){
	for (int i = 0; i < refAln[shareIdx].size(); i++){
		if (this->Aln.back()[i] != refAln[shareIdx][i]){
			if (this->Aln.back()[i] == '-'){
				insert_gaps(refAln,i);
			}
			else if (refAln[shareIdx][i] == '-')
				this->add_gaps(i);
			else{
				cout << "reference sequence and inferred sequence are not matched" << endl;
				return false;
			}
		}
	}
	return true;
}

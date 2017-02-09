#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <fstream>
#include "Aligner.h"

using namespace std;

class Node{
public:
	string label;
	double edge;
	vector<Node*> children;
//public:
	//vector<Node*> children;
	Node();
	Node(string label);
	~Node();
	bool set_label(string s);
	string get_label();
	double get_edge();
	bool set_edge(double e);
	bool writeNewick(ofstream &fout);
	bool add_child(Node *ch);
	bool is_leaf(){ return this->children.empty(); }

	
	virtual bool add_aln(string newAln) { return false; };
	
	virtual bool mapSeq(map<string,string> seq_map) {return false;}
	
	virtual float btmupAln(SubMtrx m, float idrate) {return false;}
	
	virtual bool polyAln(SubMtrx m, float idrate) {return false;} // align polytomies
	
	virtual string get_seq() {return "";}
	
	virtual bool set_seq() {return false;}

	virtual bool transitAln(vector<string> &refAln, int shareIdx) {return false;}

	virtual bool add_gaps(int pos) {return false;}

	virtual bool printAln(ofstream &fout) {return false;}

	virtual vector<string>& getAln() { vector<string> dummy; return dummy;}
	
	virtual vector<string>& getLeafLabel() { vector<string> dummy; return dummy;}
};

class Tree{
public:
	Node *root;
//public:
	Tree();
	~Tree();
	bool readNewick(string treefile);
	bool writeNewick(string treefile);
	virtual Node* create_node();
};

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <fstream>
#include "Aligner.h"

using namespace std;

class Node{
	vector<Node*> children;
	string label;
	double edge;
public:
	Node();
	Node(string label);
	~Node();
	
	bool set_label(string s);
	string get_label();
	
	double get_edge();
	bool set_edge(double e);

	Node* get_child(int n);	
	int nchild();

	bool writeNewick(ofstream &fout);
	bool add_child(Node *ch);
	bool is_leaf(){ return this->children.empty(); }
	
	virtual bool add_aln(string newAln) = 0;
	
	virtual bool mapSeq(map<string,string> seq_map) = 0;

	virtual float btmupAln(SubMtrx m, float idrate) = 0;	

	virtual bool polyAln(SubMtrx m, float idrate) = 0 ; 	

	virtual string get_seq() = 0 ;
	
	virtual bool set_seq() = 0;

	virtual bool transitAln(vector<string> &refAln, int shareIdx) = 0;

	virtual bool add_gaps(int pos) = 0;

	virtual bool printAln(ofstream &fout) = 0;

	virtual vector<string>& getAln() = 0;
	
	virtual vector<string>& getLeafLabel() = 0;
};

class Tree{
	Node *root;
public:
	Tree();
	~Tree();
	bool readNewick(string treefile);
	bool writeNewick(string treefile);
	Node*& get_root() { return this->root; }
	virtual Node* create_node() = 0;
};

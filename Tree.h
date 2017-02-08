#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <fstream>
using namespace std;

class Node{
	string label;
	double edge;
	vector<Node*> children;
public:
	Node();
	Node(string label);
	~Node();
	bool set_label(string s);
	string get_label();
	double get_edge();
	bool set_edge(double e);
	bool writeNewick(ofstream &fout);
	bool add_child(Node *ch);
};

class Tree{
	Node *root;
public:
	Tree();
	~Tree();
	bool readNewick(string treefile);
	bool writeNewick(string treefile);
	Node* create_node();
};

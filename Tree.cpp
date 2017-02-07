#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <fstream>
using namespace std;

class Node{
	string label;
	int edge;
public:
	vector<Node*> children;
	Node(){
		this->label = "";
		this->edge = 0;
	}
	Node(string label){
		this->label = label;
		this->edge = 0;
	}
	~Node(){
		for (int i = 0; i < this->children.size(); i++)
			delete this->children[i];
	}
	string get_label(){
		return this->label;
	}
	bool set_label(string s){
		this->label = s;
		return true;
	}

	bool writeNewick(ofstream &fout){
		if (!this->children.empty())
		{
			fout << "(";
			this->children[0]->writeNewick(fout);
			for (int i = 1; i < this->children.size(); i++){
				fout << ",";
				this->children[i]->writeNewick(fout);
			}
			fout << ")";
		}
		fout << this->get_label();
	}
};

class Tree{
	Node *root;
public:
	Tree(){
		root = NULL;
	}
	~Tree(){
		delete root;
	}


	bool readNewick(string treefile){
		char c;
		stack<Node*> stk;

		string label = "";
		ifstream fin;
		fin.open(treefile);
		while (!fin.eof()){
			fin >> c;
			if (c == '('){
				stk.push(NULL);
				label = "";
			}
			else if (c == ',')  {
				if (label != ""){
					Node *p = new Node;
					p->set_label(label);
					stk.push(p);
					//cout << label << endl;
				}
				label = "";
			} else if (c == ')'){
				if (label != ""){
					Node *p = new Node;
					p->set_label(label);
					stk.push(p);
					//cout << label << endl;
				}
				Node *p = new Node;
				Node *q;
				while (1){
					q = stk.top();
					stk.pop();
					if (q)
					{
						//cout << "I pushed back something!" << endl;
						p->children.push_back(q);
					}
					else
						break;
				}
				stk.push(p);
				label = "";
			} else if (c == ';'){
				this->root = stk.top();
				//this->root->set_label("root");
				stk.pop();
				break;
			}
			else {
				label += c;
			}
		}
		fin.close();
		return true;
	}
	
	bool writeNewick(string treefile){
		ofstream fout;
		fout.open(treefile);
		//this->root->set_label("root");
		//fout << this->root->get_label();
		this->root->writeNewick(fout);
		fout << ";";
		fout.close();
	}
};

int main(){
	Tree a_tree;
	a_tree.readNewick("avian_topo.tre");
	a_tree.writeNewick("out.tre");
}

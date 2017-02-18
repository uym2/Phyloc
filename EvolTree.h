#include "Tree.h"
#include <boost/expm.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
using namespace boost::numeric::ublas;

class EvolNode: public Node{
	string mySeq;
public:
	string get_seq(){ return this->mySeq; }
	bool set_seq(string seq) { this->mySeq = seq; return true; }
	

// nonsense functions that but cannot remove ...
	bool add_aln(string newAln) { return false; }
	bool mapSeq(map<string,string> seq_map) { return false; }
	float btmupAln(SubMtrx m, float idrate) { return -1; }
	bool polyAln(SubMtrx m, float idrate) { return false; }
	bool transitAln(std::vector<string> &refAln, int shareIdx) { return false; }
	bool add_gaps(int pos) { return false; }
	bool printAln(ofstream &fout) { return false; }
	bool set_seq() { return false; }
	std::vector<string>& getAln() { std::vector<string> foo; return foo; }
	std::vector<string>& getLeafLabel() { std::vector<string> foo; return foo; }
};


class EvolTree: public Tree{
	matrix<double> R;
public:
	EvolTree(){ this->R.resize(4,4); }

	bool set_R(std::vector<double> pi, std::vector<double> rate);
	bool set_R(string param_file);

	Node* create_node(){
		//cout << "Created an EvolNode" << endl;
		Node* p = new EvolNode;
		return p;
	}

	bool print_R(){
		cout << this->R << endl;
		return true;
	}
};

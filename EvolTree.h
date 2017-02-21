#include "Tree.h"
#include <boost/expm.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
using namespace std;
using namespace boost::numeric::ublas;

class EvolNode: public Node{
	string mySeq;
public:
	string get_seq(){ return this->mySeq; }
	bool set_seq(string seq) { 
		for (int i = 0; i < seq.length(); i++)
			if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T')
				return false;
		this->mySeq = seq; 
		return true; 
	}

	bool gen_seqs(const matrix<double> &R);	
	char randit(const std::vector<double> &dtrb);

};


class EvolTree: public Tree{
	matrix<double> R;
public:
	EvolTree(){ this->R.resize(4,4);}

	bool set_R(const std::vector<double> &pi, const std::vector<double> &rate);
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

	bool seq_root(const std::vector<double> &pi,int N);
	bool seq_root(string seq) { return this->get_root()->set_seq(seq); }
	bool gen_seqs() { return this->get_root()->gen_seqs(this->R); }
};

#include "EvolTree.h"
#include <cmath>
#include <random>

#define EPSILON 0.000001

bool EvolNode::gen_seqs(const matrix<double> & R){
	if (this->is_leaf())
		return true;
	for (int i = 0; i < this->nchild(); i++){
		this->get_child(i)->gen_seqs(R);
	}
}


char EvolNode::randit(const std::vector<double> &dtrb){
	std::random_device rd;
	std::mt19937 gen(rd());
    	std::uniform_real_distribution<> dis(0, 1);
	double r = dis(gen);
	int i = 0;
	double s = dtrb[0];

	while (s < r)
		s += dtrb[++i];

	switch (i)
	{
		case 0:
			this->mySeq.push_back('A');
			return 'A';
		case 1:
			this->mySeq.push_back('C');
			return 'C';
		case 2:
			this->mySeq.push_back('G');
			return 'G';
		case 3:
			this->mySeq.push_back('T');
			return 'T';
	}
}

bool EvolTree::seq_root(const std::vector<double> &pi, int N){
	for (int i = 0; i < N; i++)
		this->get_root()->randit(pi);
	return true;
}

bool EvolTree::set_R(string param_file){
	ifstream fin;
	fin.open(param_file);
	if (!fin){
		cout << "Cannot open param_file for input" << endl;
		return false;
	}
	
	std::vector<double> pi, rate;
	double x;
	for (int i = 0; i<4; i++){
		fin >> x;
		pi.push_back(x);
	}

	for (int i = 0; i<6; i++){
		fin >> x;
		rate.push_back(x);
	}

	fin.close(); 

	return this->set_R(pi,rate);
}
bool EvolTree::set_R(const std::vector<double> &pi, const std::vector<double> &rate){
	// rate order: CT AT GT AC CG AG

	if (pi.size() != 4 || rate.size() != 6){
		cout << "invalid parameters!" << endl;
		return false;
	}

	double sum = 0;
	for (int i = 0; i < pi.size(); i++)
		sum += pi[i];
	if (abs(sum-1.0) > EPSILON){
		cout << "invalid parameters!" << endl;
		return false;
	}

	double alpha,beta,gamma,delta,eta,epsilon;
	alpha = rate[5]/pi[2];
	beta = rate[3]/pi[1];
	gamma = rate[1]/pi[3];
	delta = rate[4]/pi[2];
	eta = rate[0]/pi[3];
	epsilon = rate[2]/pi[3];

	double GA,CA,TA,GC,TC,TG,AA,CC,GG,TT;
	GA = alpha*pi[0];
	CA = beta*pi[0];
	TA = gamma*pi[0];
	GC = delta*pi[1];
	TC = eta*pi[1];
	TG = epsilon*pi[2];

	AA = -(alpha*pi[2]+beta*pi[1]+gamma*pi[3]);
	CC = -(beta*pi[0]+delta*pi[2]+eta*pi[3]);
	GG = -(alpha*pi[0]+delta*pi[1]+eta*pi[3]);
	TT = -(gamma*pi[0]+epsilon*pi[2]+eta*pi[1]);

	double r = -(pi[0]*AA + pi[1]*CC + pi[2]*GG + pi[3]*TT);

	this->R(0,0) = AA/r;
	this->R(1,1) = CC/r;
	this->R(2,2) = GG/r;
	this->R(3,3) = TT/r;

	this->R(0,1) = rate[3]/r;
	this->R(0,2) = rate[5]/r;
	this->R(0,3) = rate[1]/r;

	this->R(1,0) = CA/r;
	this->R(1,2) = rate[4]/r;
	this->R(1,3) = rate[0]/r;

	this->R(2,0) = GA/r;
	this->R(2,1) = GC/r;
	this->R(2,3) = rate[2]/r;

	this->R(3,0) = TA/r;
	this->R(3,1) = TC/r;
	this->R(3,2) = TG/r;


	return true;
}	

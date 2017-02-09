#include "Tree.h"

class AlnNode: public Node{
//	vector<string> Aln;
	string mySeq;
public:
	//vector<AlnNode*> children;
	vector<string> leafLabel;
	vector<string> Aln;
	bool add_aln(string newAln){
		this->Aln.insert(Aln.end()-1,newAln);
		return true;}
	
	bool mapSeq(map<string,string> seq_map);
	
	float btmupAln(SubMtrx m, float idrate);
	
	bool polyAln(SubMtrx m, float idrate); // align polytomies
	
	string get_seq();
	
	bool set_seq();

	bool transitAln(vector<string> &refAln,int shareIdx);

	bool add_gaps(int pos);

	bool printAln(ofstream &fout){
		for (int i = 0; i < this->Aln.size()-1; i++){
			fout << ">" << this->leafLabel[i] << endl;
			fout << this->Aln[i] << endl;
		}
		return true;
	}

	vector<string>& getAln() { return this->Aln;}
	vector<string>& getLeafLabel() {return this->leafLabel;}
};

class AlnTree: public Tree{
	SubMtrx m;
	float idrate;
public:
	//AlnNode* root;
	AlnTree(SubMtrx m, float idrate)
		{ this->m=m; this->idrate=idrate; }

	Node* create_node(){
		//cout << "Created an AlnNode" << endl;
		Node* p = new AlnNode;
		return p;
	}

	bool mapSeq2Tree(string seqfile);

	bool btmupAln();

	bool printTopAln(string file){ 
		ofstream fout;
		fout.open(file);
		if (!fout)
			return false;
		this->root->printAln(fout);
		fout.close();
		return true;
	}
};

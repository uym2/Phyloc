#include "Tree.h"

class AlnNode: public Node{
	string mySeq;
	std::vector<string> leafLabel;
	std::vector<string> Aln;
public:
	bool add_aln(string newAln){
		this->Aln.insert(Aln.end()-1,newAln);
		return true;}
	
	bool mapSeq(map<string,string> seq_map);
	
	float btmupAln(SubMtrx m, float idrate);
	
	bool polyAln(SubMtrx m, float idrate); // align polytomies
	
	string get_seq();
	
	bool set_seq();

	bool transitAln(std::vector<string> &refAln,int shareIdx);

	bool add_gaps(int pos);

	bool printAln(ofstream &fout){
		for (int i = 0; i < this->Aln.size()-1; i++){
			fout << ">" << this->leafLabel[i] << endl;
			fout << this->Aln[i] << endl;
		}
		return true;
	}

	std::vector<string>& getAln() { return this->Aln;}
	std::vector<string>& getLeafLabel() {return this->leafLabel;}
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
		this->get_root()->printAln(fout);
		fout.close();
		return true;
	}
};

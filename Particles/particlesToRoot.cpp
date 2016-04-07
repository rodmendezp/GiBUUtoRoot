#include <fstream>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

TString fName;
TString Metal;
TString tLine;
TString rootfName;

int main(int argc,  char **argv){
	// Arg 1 is the filename
	rootfName = "part.root";
	if(argc > 1)
		fName = (TString) argv[1];
	if(argc > 2)
		Metal = (TString) argv[2];
		
	TFile *f = new TFile(rootfName, "RECREATE");
	TNtuple *ntuple_elec;
	TNtuple *ntuple_part;
	
	ntuple_elec = new TNtuple("gibuu_elec", "gibuu_elec", "TargType:Q2:Nu:Xb:Xf");
	ntuple_part = new TNtuple("gibuu_parts", "gibuu_part", "PID:TargType:Q2:Nu:Xb:W:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:P:T4:Betta");
		
	ifstream in;
	in.open(fName);
	if(in.is_open())
		std::cout << "File " << fName << " is open" << std::endl;
	else
	{
		std::cout << "Error opening " << fName << std::endl;
		return -1;	
	}
	
	Int_t nParts;
	Int_t partID;
	Int_t eventType;
	Float_t weight;
	Float_t unused;
	Float_t px, py, pz;
	Float_t e;
	Float_t m;
	Float_t spin;
	Float_t nu;
	Float_t q2;
	Float_t eps;
	Float_t phiL;
			
	while(!in.eof()){
		in >> tLine;
		if(tLine == "<event>"){
			in >> nParts >> unused >> weight >> unused >> unused >> unused;
			for(Int_t i = 0; i < nParts; i++){
				in >> partID >> unused >> unused >> unused >> unused >> unused >> px >> py >> pz >> e >> m >> unused >> spin;
				// Calculate all the missing variables
				//std::cout << partID << "\t" << px << "\t" << py << "\t" << pz << "\t" << e << "\t" << m << "\t" << spin << std::endl;
			}
			in >> tLine >> unused >> nu >> q2 >> eps >> phiL >> eventType;
			//std::cout << tLine << "\t" << nu << "\t" << q2 << "\t" << eps << "\t" << phiL << "\t" << eventType << std::endl;
			//ntuple_elec->Fill(1., q2, nu,)
		}
	}
	
	f->cd();
	ntuple_elec->Write();
	ntuple_part->Write();
	
	
	return 0;
}

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
	rootfName = "part.root"
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
			
	while(!in.eof()){
		in >> tLine;
		if(tLine == "<event>")		
			std::cout << tLine << std::endl;
	}
	
	f->cd();
	ntuple_elec->Write();
	ntuple_part->Write();
	
	
	return 0;
}

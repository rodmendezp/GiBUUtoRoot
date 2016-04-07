#include <fstream>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"

TString fName;
TString Metal;
TString tLine;
TString rootfName;

Float_t kMassP = 0.938272;
Float_t kMassPi = 0.139570;
Float_t kEbeam = 5.0;

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
	Int_t *partID;
	Int_t eventType;
	Float_t weight;
	Float_t unused;
	Float_t *px;
	Float_t *py;
	Float_t *pz;
	Float_t *e;
	Float_t *m;
	Float_t *spin;
	Float_t nu;
	Float_t q2;
	Float_t eps;
	Float_t phiL;
	
	Float_t xb;
	
	Float_t *w;
	Float_t *theta;
	Float_t *phi;
	Float_t *zh;
	Float_t *pt;
	Float_t *w2p;
	Float_t *xf;
	Float_t *p;
	Float_t *t4;
	Float_t *betta;
	Float_t *plcm;
	
			
	while(!in.eof()){
		in >> tLine;
		if(tLine == "<event>"){
			in >> nParts >> unused >> weight >> unused >> unused >> unused;
			partID = new Int_t[nParts];
			px = new Float_t[nParts];
			py = new Float_t[nParts];
			pz = new Float_t[nParts];
			e = new Float_t[nParts];
			m = new Float_t[nParts];
			spin = new Float_t[nParts];
			w = new Float_t[nParts];
			theta = new Float_t[nParts];
			phi = new Float_t[nParts];
			zh = new Float_t[nParts];
			pt = new Float_t[nParts];
			w2p = new Float_t[nParts];
			xf = new Float_t[nParts];
			p = new Float_t[nParts];
			t4 = new Float_t[nParts];
			betta = new Float_t[nParts];
			plcm = new Float_t[nParts];
			for(Int_t i = 0; i < nParts; i++){
				in >> partID[i] >> unused >> unused >> unused >> unused >> unused;
				in >> px[i] >> py[i] >> pz[i] >> e[i] >> m[i] >> unused >> spin[i];
				//std::cout << partID[i] << "\t" << px[i] << "\t" << py[i] << "\t" << pz[i];
				//std::cout << "\t" << e[i] << "\t" << m[i] << "\t" << spin[i] << std::endl;
			}
			in >> tLine >> unused >> nu >> q2 >> eps >> phiL >> eventType;
			xb = q2/(2*nu*kMassP);
			ntuple_elec->Fill(1, q2, nu, xb, weight);
			for(Int_t i = 0; i < nParts; i++){
				p[i] = TMath::Sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
				w[i] = TMath::Sqrt(kMassP * kMassP + 2 * kMassP * nu - q2); // Check
				plcm[i] = (nu + kMassP) * (TMath::Sqrt(1.) - TMath::Sqrt(q2 + nu * nu) * zh[i] * nu / (nu + kMassP)) / w[i];
				theta[i] = 1.;
				phi[i] = 1.;
				zh[i] = e[i]/nu; // Check
				//pt[i] = TMath::Sqrt(p[i]*p[i]*(1 - TMath::Cos(theta[i])*TMath::Cos(theta[i]))); // Check
				w2p[i] = 1.;
				xf[i] = 1.;
				t4[i] = 1.;
				betta[i] = TMath::Sqrt(p[i]*p[i]/(p[i]*p[i] + m[i]*m[i])) ; // Check
			}
			delete[] partID;
			delete[] px;
			delete[] py;
			delete[] pz;
			delete[] e;
			delete[] m;
			delete[] spin;
			delete[] w;
			delete[] theta;
			delete[] phi;
			delete[] zh;
			delete[] pt;
			delete[] w2p;
			delete[] xf;
			delete[] p;
			delete[] t4;
			delete[] betta;
			//std::cout << tLine << "\t" << nu << "\t" << q2 << "\t" << eps << "\t" << phiL << "\t" << eventType << std::endl;
			//ntuple_elec->Fill(1., q2, nu,)
		}
	}
	
	f->cd();
	ntuple_elec->Write();
	ntuple_part->Write();
	
	
	return 0;
}

#include <fstream>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"

TString fName;
TString Metal;
TString tLine;
TString rootfName;

Float_t kMassP = 0.938272;
Float_t kMassPi = 0.139570;
Float_t kEbeam = 5.0;
Float_t kMntr = 0.939565;

void MomElectron(Float_t &pex, Float_t &pey, Float_t &pez, Float_t q2, Float_t nu, Float_t phie);
void MomPhoton(Float_t &pfx, Float_t &pfy, Float_t &pfz, Float_t pex, Float_t pey, Float_t pez, Float_t nu);
Float_t ThetaPQ(Float_t pfx, Float_t pfy, Float_t pfz, Float_t px, Float_t py, Float_t pz);
Float_t PhiPQ(Float_t pfx, Float_t pfy, Float_t pfz, Float_t px, Float_t py, Float_t pz);
Float_t PmaxCM(Float_t W);

int main(int argc,  char **argv){
	rootfName = "part.root";
	if(argc > 1)
		fName = (TString) argv[1];
	else if(argc > 2)
		Metal = (TString) argv[2];
	else{
		std::cout << "Please add the name of the input file" << std::endl;
		std::cout << "./particlesToRoot input.file" << std::endl;
		return 0;
	}
		
	TFile *f = new TFile(rootfName, "RECREATE");
	TNtuple *ntuple_elec;
	TNtuple *ntuple_part;
	
	ntuple_elec = new TNtuple("gibuu_elec", "gibuu_elec", "TargType:Q2:Nu:Xb:Xf:Weight");
	ntuple_part = new TNtuple("gibuu_parts", "gibuu_part", "PID:TargType:Q2:Nu:Xb:W:ThetaPQ:PhiPQ:Zh:Pt:Xf:P:T4:Betta:Weight");
		
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
	Float_t *betta;
	Float_t *plcm;
	
	Float_t pex, pey, pez;
	Float_t pfx, pfy, pfz;
			
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
			betta = new Float_t[nParts];
			plcm = new Float_t[nParts];
			for(Int_t i = 0; i < nParts; i++){
				in >> partID[i] >> unused >> unused >> unused >> unused >> unused;
				in >> px[i] >> py[i] >> pz[i] >> e[i] >> m[i] >> unused >> spin[i];
			}
			in >> tLine >> unused >> nu >> q2 >> eps >> phiL >> eventType;
			xb = q2/(2*nu*kMassP);
			ntuple_elec->Fill(1, q2, nu, xb, weight);
			for(Int_t i = 0; i < nParts; i++){
				p[i] = TMath::Sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
				w[i] = TMath::Sqrt(kMassP * kMassP + 2 * kMassP * nu - q2); // Check
				plcm[i] = (nu + kMassP) * (TMath::Sqrt(1.) - TMath::Sqrt(q2 + nu * nu) * zh[i] * nu / (nu + kMassP)) / w[i];
				MomElectron(pex, pey, pez, q2, nu, phiL);
				MomPhoton(pfx, pfy, pfz, pex, pey, pez, nu);
				theta[i] = ThetaPQ(pfx,  pfy,  pfz,  p[i],  p[i],  p[i]);
				phi[i] = PhiPQ(pfx,  pfy,  pfz,  p[i],  p[i],  p[i]);
				zh[i] = e[i]/nu; // Check
				pt[i] = TMath::Sqrt(p[i]*p[i]*(1 - TMath::Cos(theta[i])*TMath::Cos(theta[i]))); // Check
				xf[i] = plcm[i]/PmaxCM(w[i]);
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
			delete[] betta;
		}
	}
	
	f->cd();
	ntuple_elec->Write();
	ntuple_part->Write();
	
	return 0;
}

void MomElectron(Float_t &pex, Float_t &pey, Float_t &pez, Float_t q2, Float_t nu, Float_t phi){
	Float_t costheta;
	Float_t kEscat;
	Float_t pexy;
	
	kEscat = kEbeam - nu;
	costheta = (2*kEbeam*kEscat - q2)/(2*kEbeam*kEscat);
	pez = (kEbeam - nu)*costheta;
	pexy = (kEbeam - nu)*TMath::Sqrt(1 - costheta*costheta);
	pex = pexy*TMath::Cos(phi);
	pey = pexy*TMath::Sin(phi);
		
	return;
}

void MomPhoton(Float_t &pfx, Float_t &pfy, Float_t &pfz, Float_t pex, Float_t pey, Float_t pez, Float_t nu){

	pfx = -pex;
	pfy = -pey;
	pfz = TMath::Sqrt(nu*nu-pfx*pfx-pfy*pfy);
	
	return;
}

Float_t ThetaPQ(Float_t pfx, Float_t pfy, Float_t pfz, Float_t px, Float_t py, Float_t pz){
	Float_t thetapq;
    TVector3 Vpi(px,py,pz);
    TVector3 Vvirt(pfx,pfy,pfz);
    thetapq = Vvirt.Angle(Vpi)*180./(TMath::Pi());
    
    return thetapq;
}

Float_t PhiPQ(Float_t pfx, Float_t pfy, Float_t pfz, Float_t px, Float_t py, Float_t pz){
	Float_t phipq;
	TVector3 Vpi(px,py,pz);
	TVector3 Vvirt(pfx,pfy,pfz);
	Double_t phi_z = TMath::Pi()-Vvirt.Phi();
    Vvirt.RotateZ(phi_z);
    Vpi.RotateZ(phi_z);
    TVector3 Vhelp(0.,0.,1.);
    Double_t phi_y = Vvirt.Angle(Vhelp);
    Vvirt.RotateY(phi_y);
    Vpi.RotateY(phi_y);
    phipq=Vpi.Phi() * 180./(TMath::Pi());	
    
    return phipq;
}

Float_t PmaxCM(Float_t W){
	Float_t pmaxcm;
	
	pmaxcm = TMath::Sqrt(TMath::Power(W * W - kMntr * kMntr + kMassPi * kMassPi, 2) - 4. * kMassPi * kMassPi * W * W) / 2. / W;
	
	return pmaxcm;
}

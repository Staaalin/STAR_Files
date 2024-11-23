#include "v0dst.h"
#include <bitset>
//ROOT header files
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TChain.h"
// C/C++ heder
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <cmath>
using namespace std;

TChain* ChainThem(const char* filelist, const char* treename, int nlist = 0, int block = 100);

const int kCentBin = 9;
const int kDeltamombin = 0;
const int kyBin = 3;//signal rap
//const int kyBin = 9;//pair rap
const int kktBin = 4;
const int totalv0 = 200;
#define kZBin  58
#define kcentBin 9
#define kEvent 6
//#define totalv0 50 
TLorentzVector Dosmearplus(TLorentzVector Four_mom, int index_E);
TLorentzVector Dosmearminus(TLorentzVector Four_mom, int index_E);
double getQ_LCMS(TLorentzVector Four_mom1, TLorentzVector Four_mom2);
TVector3 getQosl_LCMS(TLorentzVector Four_mom1, TLorentzVector Four_mom2);
TLorentzVector boostfourmom_to_AuAuCMS(TLorentzVector Four_mom,double energy);
TLorentzVector getnewfourmom(TLorentzVector Four_mom,double delta_mom);
double getphistar(TLorentzVector Four_mom1, TLorentzVector Four_mom2, int q1, int q2,double Bz, double tpcR);
Int_t getZbin(Float_t Zvert , Int_t index_E);
Int_t getCentBin(Int_t nrefmult , Int_t index_E);
Int_t getsinglerapbin(Float_t rap ,Int_t index_E);
Int_t getpairrapbin(Float_t rap ,Int_t index_E);
double getSL(Int_t padRow1To24Track1 ,Int_t padRow25To45Track1 ,ULong64_t IpadRow1 ,Int_t nhits1 , Int_t padRow1To24Track2 , Int_t padRow25To45Track2 ,ULong64_t IpadRow2 ,Int_t nhits2 ,Int_t index_E);
double getshiftplus(double p, int index_E);
double getshiftminus(double p, int index_E);
struct Lm_mixed{
	Int_t id;

	Double_t Px;
	Double_t Py;
	Double_t Pz;

	Double_t mass;
	//Int_t ZBin;
	//Int_t CentBin;
	//unsigned long padrow1;
	//unsigned long padrow2;
	Int_t padrow1;
	Int_t padrow2;
	ULong64_t ipadrow;
	Int_t nhits;
	Int_t totalv0;

};
Lm_mixed pionplus[kZBin][kcentBin][kEvent][totalv0],pionminus[kZBin][kcentBin][kEvent][totalv0];

int main(int argc, char** argv){
	TString InputFileList(argv[1]);
	TString OutputDir(argv[2]);
	Int_t cutID = stoi(string(argv[3]));//0 defualt, 1 2 3.......sys_err
	Int_t index_E = stoi(string(argv[4]));//30-->3.0GeV 32-->3.2Gev...... 77-->7.7GeV
	Float_t pdgmass_pionplus = 0.13957039;
	Float_t pdgmass_pionminus = 0.13957039;
	int qbin=50;
	double qmin=-0.25;
	double qmax=0.25;
	int smear_index=1;//1->do smear 0->do not smear
	if(cutID!=0)smear_index=0;
	int index_1DCF=0;
	int index_3DCF=1;
	int index_Hphistar=0;
	//momshift scan only use for rapindex==0 currently
	int index_1Dmomshift=0;
	int index_3Dmomshift=0;
	double deltamom=0.002;
	int rapindex=-99;//0-->single ybin 1-->pair ybin;
	if(kyBin==9) rapindex=1;//0-->single ybin 1-->pair ybin;
	if(kyBin==3) rapindex=0;//0-->single ybin 1-->pair ybin;
	if(rapindex<0){
		cout<<"Erorr: rapindex can not < 0"<<endl;
		return 0;
	}
	//some cuts
	double rapcutmax = 100;
	double rapcutmin = -100;
	if(rapindex==0){
		rapcutmax = 0;
		rapcutmin = -1;
		if(index_E==77){
			rapcutmax = 0.5;
			rapcutmin = -0.5;
		}
	}

	double ptcutmax = 1.5;
	double ptcutmin = 0.15;

	double ktcutmax=10.6; 
	double ktcutmin=-0.15; 
	double ktcut1=-0.15; 
	double ktcut2=10.25; 
	double ktcut3=10.35; 
	double ktcut4=10.60; 
	double ktcut5=11.0; 

	if(kktBin==4){
		ktcutmax=0.6; 
		ktcutmin=0.15; 
		ktcut1=0.15; 
		ktcut2=0.25; 
		ktcut3=0.35; 
		ktcut4=0.45; 
		ktcut5=0.60; 
	}
	int nhitsfitmax=15;
	if(cutID==1)nhitsfitmax=17;
	if(cutID==2)nhitsfitmax=13;

	//form TPC_R 0->0.6m 13->1.9m   i<nRmin||i>nRmax
	int nRmin=6;
	int nRmax=10;
	//double dphistarcut=0.0;
	//double detacut2=0.0;
	double dphistarcut=0.06;
	double detacut2=0.04;
	if(cutID==3){
		dphistarcut=0.065;//sys
		detacut2=0.045;//sys
	}
	if(cutID==4){
		dphistarcut=0.055;//sys
		detacut2=0.035;//sys
	}

	double slcutmin=-0.5;
	double slcutmax=0.6;
	if(cutID==5)slcutmax=0.8;//sys
	if(cutID==6)slcutmax=0.4;//sys

	double dcamax=3.0;
	if(cutID==7)dcamax=2.0;//sys
	if(cutID==8)dcamax=2.5;//sys

	//double pcut=0.0;
	double pcut=0.55;
	double pcutmin=0.15;
	double pcutmax=1.5;
	double m2max=0.08;
	double m2min=-0.05;
	double etacutmin=-20;
	double etacutmax=0;
	if(index_E==77){
		etacutmin=-1.5;
		etacutmax=1.5;
	}
	double ycm=-999;
	if(index_E==30)ycm=1.045;
	if(index_E==32)ycm=1.139;
	if(index_E==35)ycm=1.254;
	if(index_E==39)ycm=1.375;
	if(index_E==45)ycm=1.52;
	if(index_E==52)ycm=1.683;
	if(index_E==77)ycm=0;

	int find_pair_cut_index=0;//0->don't need to find pair cut, 1->to find sl cut, 2->to find dphidetacut
	int qoutqside_index=1;//add qout*qside cut to rm merging effect 0_->don't add, 1->add
	int qoutqside_index2=0;//add qout*qside cut to rm merging effect 0_->don't add, 1->add(2bad octant)

	TFile * fcoul;
	fcoul= new TFile("/star/u/qyq/book1_like_sign_pion_coulomb.root");//from CAT,r=5fm
	if(cutID==9) fcoul = new TFile("/star/u/qyq/cl/book1_like_sign_pion_coul_R=7.0.root");//from CAT,r=7fm
	if(cutID==10) fcoul = new TFile("/star/u/qyq/cl/book1_like_sign_pion_coul_R=3.0.root");//from CAT,r=3fm
	TH1D *hcoul;
	if(cutID!=9&&cutID!=10)hcoul= (TH1D*)fcoul->Get("hCk");
	if(cutID==9) hcoul = (TH1D*)fcoul->Get("hckR_7.0");
	if(cutID==10) hcoul = (TH1D*)fcoul->Get("hckR_3.0");
	Int_t runID = 3;//not use 
	//add some things for smear
	double lamplus[4][2][4],outplus[4][2][4],sideplus[4][2][4],longplus[4][2][4],ol2plus[4][2][4];//[centbin][ybin][ktbin]
	double lamminus[4][2][4],outminus[4][2][4],sideminus[4][2][4],longminus[4][2][4],ol2minus[4][2][4];
	TString ename;
	if(index_E==30)ename=Form("3p0");
	if(index_E==32)ename=Form("3p2");
	if(index_E==35)ename=Form("3p5");
	if(index_E==39)ename=Form("3p9");
	if(index_E==45)ename=Form("4p5");
	if(index_E==52)ename=Form("5p2");
	if(index_E==77)ename=Form("7p7");
	ifstream plusout(Form("/star/u/qyq/data01/rawHBTpionPar/%spluspluspout2rap_momshiftscan_dmom0.txt",ename.Data())); 
	ifstream minusout(Form("/star/u/qyq/data01/rawHBTpionPar/%sminusminuspout2rap_momshiftscan_dmom0.txt",ename.Data())); 
	//ifstream plusout("/star/u/qyq/data01/rawHBTpionPar/3p2pluspluspout2rap_momshiftscan_dmom0.txt"); 
	//ifstream minusout("/star/u/qyq/data01/rawHBTpionPar/3p2minusminuspout2rap_momshiftscan_dmom0.txt"); 
	for(int i=0;i<4;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<4;k++){
				plusout>>lamplus[i][j][k]>>outplus[i][j][k]>>sideplus[i][j][k]>>longplus[i][j][k]>>ol2plus[i][j][k];
				minusout>>lamminus[i][j][k]>>outminus[i][j][k]>>sideminus[i][j][k]>>longminus[i][j][k]>>ol2minus[i][j][k];
			}
		}
	}

	cout<<"Run18 3GeV FXT now"<<endl;

	TChain * t = ChainThem(InputFileList.Data(),"V0PicoDst",1,1);
	if(!t){ cout<<"ERROR: no files are added to the chain!"<<endl; return 0; }
	v0dst v0dst(t);
	TFile ohm(OutputDir,"recreate");
	cout<<"to define  QA histograms "<<endl;

	TH1F *hrefmult = new TH1F("hrefmult","reference multiplicity",800,0,800);
	TH1F *hgrefmult = new TH1F("hgrefmult","reference multiplicity",800,0,800);
	TH1F *hSelectNRefMultCorr  = new TH1F( "SelectRefMultCorr", "Corrected Reference Multiplicity of selected events;refmult;counts", 800, 0.0,800.0 ) ;
	TH2F *hrefmultvstofmatch_before_cut = new TH2F ("refmultvstofmatch_before_cut","refmult vs tofmatch;refmult;tofmatch",150,0,500,100,0,200);
	TH2F *hrefmultvstofmatch = new TH2F ("refmultvstofmatch","refmult vs tofmatch;refmult;tofmatch",150,0,500,100,0,200);
	TH1F *hvz = new TH1F("hvz","vertex z position;vz(cm); ",200,190,210);
	TH2F *hvxvy = new TH2F("hvxvy","vx vs vy;vx(cm);vy(cm)",100,-5,5,100,-5,5);
	TH1F *hevtZbin =  new TH1F("hevtZbin","",10,-1,9);
	TH1F *hevtcentbin =  new TH1F("hevtcentbin","",11,-1,10);

	cout<<"to define Track histograms "<<endl;
	//pionminus

	TH1F *hmass2pionminus = new TH1F("hmass2pionminus",";m^{2}(GeV/c^{2});",100,-0.2,0.2);
	TH1F *hdcaminus = new TH1F("hdcaminus","",50,0,5); 
	TH1F *hnsigmaminus =  new TH1F("hnsigmaminus","",100,-3,3);
	TH2F *hTofPIDminus = new TH2F("htofpidminus",";p(GeV/c);m^{2}(GeV/c^{2})",300,0.,3,300,-0.5,2);
	TH1F *hdedxminus =  new TH1F("hdedxminus","",100,0,10);
	TH2F *hpdedxminus =  new TH2F("hpdedxminus","",100,0,2,100,0,10);
	TH2F *hptdedxminus =  new TH2F("hptdedxminus","",100,0,2,100,0,10);
	TH2F *hnsigma_m2minus = new TH2F("hnsigma_m2minus","",80,-4,4,300,-1,2);
	TH3F *hnsigma_P2_3Dminus = new TH3F("hnsigma_P2_3Dminus","nsigma vs cent vs p; nsigma; centnumber;pt",300,-10,10,10,-1,9,600,0,2);
	TH3F *hmass2_P_3Dminus = new TH3F("hmass2_P_3Dminus","mass2 vs cent vs p; mass2; centnumber;p",300,-1,2.5,10,-1,9,600,0,4);
	TH2F *hy_pt_pionminus = new TH2F("hy_pt_pionminus","",400,-2.5,2.5,1000,0,2);
	TH2F *hy_pt_pionminus_same = new TH2F("hy_pt_pionminus_same","",400,-2.5,2.5,1000,0,2);
	TH2F *hy_pt_pionminus_mix = new TH2F("hy_pt_pionminus_mix","",400,-2.5,2.5,1000,0,2);
	TH1F *hpxminus = new TH1F("hpxminus","",100,-2,2);
	TH1F *hpyminus = new TH1F("hpyminus","",100,-2,2);
	TH1F *hpzminus = new TH1F("hpzminus","",100,-2,2);
	TH1F *hpxminus_same = new TH1F("hpxminus_same","",100,-2,2);
	TH1F *hpyminus_same = new TH1F("hpyminus_same","",100,-2,2);
	TH1F *hpzminus_same = new TH1F("hpzminus_same","",100,-2,2);
	TH1F *hpxminus_mix = new TH1F("hpxminus_mix","",100,-2,2);
	TH1F *hpyminus_mix = new TH1F("hpyminus_mix","",100,-2,2);
	TH1F *hpzminus_mix = new TH1F("hpzminus_mix","",100,-2,2);
	TH1F *hphiminus = new TH1F("hphiminus","",300,-7,7);
	TH1F *hphiminus_same = new TH1F("hphiminus_same","",300,-7,7);
	TH1F *hphiminus_mix = new TH1F("hphiminus_mix","",300,-7,7);
	TH1F *hetaminus = new TH1F("hetaminus","",100,-2,2);
	TH1F *hetaminus_same = new TH1F("hetaminus_same","",100,-2,2);
	TH1F *hetaminus_mix = new TH1F("hetaminus_mix","",100,-2,2);
	//pionplus	
	TH1F *hmass2pionplus = new TH1F("hmass2pionplus",";m^{2}(GeV/c^{2});",100,-0.2,0.2);
	TH1F *hdcaplus = new TH1F("hdcaplus","",50,0,5); 
	TH1F *hnsigmaplus =  new TH1F("hnsigmaplus",";p(GeV/c);m^{2}(GeV/c^{2})",100,-3,3);
	TH2F *hTofPIDplus = new TH2F("htofpidplus","",300,0.,3,300,-0.5,2);
	TH1F *hdedxplus =  new TH1F("hdedxplus","",100,0,10);
	TH2F *hpdedxplus =  new TH2F("hpdedxplus","",100,0,2,100,0,10);
	TH2F *hptdedxplus =  new TH2F("hptdedxplus","",100,0,2,100,0,10);
	TH2F *hnsigma_m2plus = new TH2F("hnsigma_m2plus","",80,-4,4,300,-1,2);
	TH3F *hnsigma_P2_3Dplus = new TH3F("hnsigma_P2_3Dplus","nsigma vs cent vs p; nsigma; centnumber;pt",300,-10,10,10,-1,9,600,0,2);
	TH3F *hmass2_P_3Dplus = new TH3F("hmass2_P_3Dplus","mass2 vs cent vs p; mass2; centnumber;p",300,-1,2.5,10,-1,9,600,0,4);
	TH2F *hy_pt_pionplus = new TH2F("hy_pt_pionplus","",400,-2.5,2.5,1000,0,2);
	TH2F *hy_pt_pionplus_same = new TH2F("hy_pt_pionplus_same","",400,-2.5,2.5,1000,0,2);
	TH2F *hy_pt_pionplus_mix = new TH2F("hy_pt_pionplus_mix","",400,-2.5,2.5,1000,0,2);
	TH1F *hpxplus = new TH1F("hpxplus","",100,-2,2);
	TH1F *hpyplus = new TH1F("hpyplus","",100,-2,2);
	TH1F *hpzplus = new TH1F("hpzplus","",100,-2,2);
	TH1F *hpxplus_same = new TH1F("hpxplus_same","",100,-2,2);
	TH1F *hpyplus_same = new TH1F("hpyplus_same","",100,-2,2);
	TH1F *hpzplus_same = new TH1F("hpzplus_same","",100,-2,2);
	TH1F *hpxplus_mix = new TH1F("hpxplus_mix","",100,-2,2);
	TH1F *hpyplus_mix = new TH1F("hpyplus_mix","",100,-2,2);
	TH1F *hpzplus_mix = new TH1F("hpzplus_mix","",100,-2,2);
	TH1F *hphiplus = new TH1F("hphiplus","",300,-7,7);
	TH1F *hphiplus_same = new TH1F("hphiplus_same","",300,-7,7);
	TH1F *hphiplus_mix = new TH1F("hphiplus_mix","",300,-7,7);
	TH1F *hetaplus = new TH1F("hetaplus","",100,-2,2);
	TH1F *hetaplus_same = new TH1F("hetaplus_same","",100,-2,2);
	TH1F *hetaplus_mix = new TH1F("hetaplus_mix","",100,-2,2);
	TH1F * NPionminus[kCentBin+1];
	TH1F * NPionplus[kCentBin+1];

	cout<<"to define Pair histograms "<<endl;
	TH2F *hqinv_sl_plusplus_same = new TH2F("hqinv_sl_plusplus_same","same event;q_{inv} (GeV);SL",150,0,0.8,80,-0.6,1);
	TH2F *hqinv_sl_minusminus_same = new TH2F("hqinv_sl_minusminus_same","same event;q_{inv} (GeV);SL",150,0,0.8,80,-0.6,1);

	TH2F *hqinv_sl_plusplus_mix = new TH2F("hqinv_sl_plusplus_mix","mixed event;q_{inv} (GeV);SL",150,0,0.8,80,-0.6,1);
	TH2F *hqinv_sl_minusminus_mix = new TH2F("hqinv_sl_minusminus_mix","mixed event;q_{inv} (GeV);SL",150,0,0.8,80,-0.6,1);
	TH2F *hpairy_kt_pionplusplus = new TH2F("hpairy_kt_pionplusplus","",400,-2.5,2.5,1000,0,3);
	TH2F *hpairy_kt_pionminusminus = new TH2F("hpairy_kt_pionminusminus","",400,-2.5,2.5,1000,0,3);
	TH2F *hpairy_kt_pionplusminus = new TH2F("hpairy_kt_pionplusminus","",400,-2.5,2.5,1000,0,3);
	TProfile *haverage_kt_plusminus[kCentBin][kyBin];
	TProfile *haverage_kt_plusplus[kCentBin][kyBin];
	TProfile *haverage_kt_minusminus[kCentBin][kyBin];
	for(int i=0;i<kCentBin;i++){
		for(int j=0;j<kyBin;j++){
			TString hnameplusminus=Form("haverage_kt_plusminus_cent%d_y%d",i,j);
			TString hnameplus=Form("haverage_kt_plusplus_cent%d_y%d",i,j);
			TString hnameminus=Form("haverage_kt_minusminus_cent%d_y%d",i,j);
			haverage_kt_plusminus[i][j] = new TProfile(hnameplusminus,"",kktBin,0,kktBin);
			haverage_kt_plusplus[i][j] = new TProfile(hnameplus,"",kktBin,0,kktBin);
			haverage_kt_minusminus[i][j] = new TProfile(hnameminus,"",kktBin,0,kktBin);
		}
	}
	TH1F *hdphiplus = new TH1F("hdphiplus","",100,-0.5,0.5);
	TH1F *hdetaplus = new TH1F("hdetaplus","",100,-0.5,0.5);
	TH1F *hdmomplus = new TH1F("hdmomplus","",100,-0.5,0.5);
	TH1F *hdphiminus = new TH1F("hdphiminus","",100,-0.5,0.5);
	TH1F *hdetaminus = new TH1F("hdetaminus","",100,-0.5,0.5);
	TH1F *hdmomminus = new TH1F("hdmomminus","",100,-0.5,0.5);
	TH1F *hkt_plusminus[kCentBin][kyBin];
	TH1F *hkt_plusplus[kCentBin][kyBin];
	TH1F *hkt_minusminus[kCentBin][kyBin];
	for(int i=0;i<kCentBin;i++){
		for(int j=0;j<kyBin;j++){
			TString hnameplusminus=Form("hkt_plusminus_cent%d_y%d",i,j);
			TString hnameplus=Form("hkt_plusplus_cent%d_y%d",i,j);
			TString hnameminus=Form("hkt_minusminus_cent%d_y%d",i,j);
			hkt_plusminus[i][j] = new TH1F(hnameplusminus,"",100,0,2);
			hkt_plusplus[i][j] = new TH1F(hnameplus,"",100,0,2);
			hkt_minusminus[i][j] = new TH1F(hnameminus,"",100,0,2);
		}
	}
	const int nTpcR=14;
	double TpcR[nTpcR]={0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};
	TH2F * hdphistardeta_plusplus[kyBin][kktBin][nTpcR];
	TH2F * hdphistardeta_plusplus_mix[kyBin][kktBin][nTpcR];
	TH2F * hdphistardeta_minusminus[kyBin][kktBin][nTpcR];
	TH2F * hdphistardeta_minusminus_mix[kyBin][kktBin][nTpcR];
	if(index_Hphistar==1){
		cout<<"to define dphistardeta histograms Loop"<<endl;
		for(int iybin=0;iybin<kyBin;iybin++){
			for(int i=0;i<kktBin;i++){
				for(int j=0;j<nTpcR;j++){
					TString hnameplusplus=Form("hdphistardeta_plusplusy%dkt%dR%.2f",iybin,i,TpcR[j]);
					TString hnameplusplus_mix=Form("hdphistardeta_plusplus_mixy%dkt%dR%.2f",iybin,i,TpcR[j]);
					TString hnameminusminus=Form("hdphistardeta_minusminusy%dkt%dR%.2f",iybin,i,TpcR[j]);
					TString hnameminusminus_mix=Form("hdphistardeta_minusminus_mixy%dkt%dR%.2f",iybin,i,TpcR[j]);
					hdphistardeta_plusplus[iybin][i][j] = new TH2F(hnameplusplus,"",150,-0.2,0.2,150,-0.2,0.2);
					hdphistardeta_plusplus_mix[iybin][i][j] = new TH2F(hnameplusplus_mix,"",150,-0.2,0.2,150,-0.2,0.2);
					hdphistardeta_minusminus[iybin][i][j] = new TH2F(hnameminusminus,"",150,-0.2,0.2,150,-0.2,0.2);
					hdphistardeta_minusminus_mix[iybin][i][j] = new TH2F(hnameminusminus_mix,"",150,-0.2,0.2,150,-0.2,0.2);
				}
			}
		}

	}
	cout<<"to define 1DCF histograms "<<endl;
	TH3F * hQ_plus_minus[kCentBin+1];
	TH3F * hQ_plus_plus[kCentBin+1];
	TH3F * hQ_minus_minus[kCentBin+1];
	TH3F * hLevyQ_plus_plus[kCentBin+1];
	TH3F * hLevyQ_minus_minus[kCentBin+1];

	TH3F * hQ_plus_minus_mix[kCentBin+1];
	TH3F * hQ_plus_plus_mix[kCentBin+1];
	TH3F * hQ_minus_minus_mix[kCentBin+1];
	TH3F * hLevyQ_plus_plus_mix[kCentBin+1];
	TH3F * hLevyQ_minus_minus_mix[kCentBin+1];

	if(index_1DCF==1){
		cout<<"to define 1DCF histograms Loop"<<endl;
		for(Int_t i=0;i<kCentBin;i++)
		{ 
			char *hName1="hQ";
			char hName[200];
			char hName2[200];
			char hName3[200];
			char hName4[200];

			char hName2_mix[200];
			char hName3_mix[200];
			char hName4_mix[200];
			char *hTitle0="Centrality";    
			char hTitle[200];
			sprintf(hTitle,"%s%d",hTitle0,i);
			sprintf(hName2_mix,"%s_Rebinplusminus_mix%d",hName1,i);
			sprintf(hName3_mix,"%s_Rebinplusplus_mix%d",hName1,i);
			sprintf(hName4_mix,"%s_Rebinminusminus_mix%d",hName1,i);

			sprintf(hName2,"%s_Rebinplusminus%d",hName1,i);
			sprintf(hName3,"%s_Rebinplusplus%d",hName1,i);
			sprintf(hName4,"%s_Rebinminusminus%d",hName1,i);
			hQ_plus_minus[i] =  new TH3F(hName2,hTitle,800,0,1,15,0,15,15,0,15);             
			hQ_plus_plus[i] =  new TH3F(hName3,hTitle,800,0,1,15,0,15,15,0,15);             
			hQ_minus_minus[i] =  new TH3F(hName4,hTitle,800,0,1,15,0,15,15,0,15);             

			hQ_plus_minus_mix[i] =  new TH3F(hName2_mix,hTitle,800,0,1,15,0,15,15,0,15);             
			hQ_plus_plus_mix[i] =  new TH3F(hName3_mix,hTitle,800,0,1,15,0,15,15,0,15);             
			hQ_minus_minus_mix[i] =  new TH3F(hName4_mix,hTitle,800,0,1,15,0,15,15,0,15);             

		} 
		for(Int_t i=0;i<kCentBin;i++)
		{ 
			char *hName1="hQLevy";
			char hName3[200];
			char hName4[200];

			char hName3_mix[200];
			char hName4_mix[200];
			char *hTitle0="Centrality";    
			char hTitle[200];
			sprintf(hTitle,"%s%d",hTitle0,i);
			sprintf(hName3_mix,"%s_Rebinplusplus_mix%d",hName1,i);
			sprintf(hName4_mix,"%s_Rebinminusminus_mix%d",hName1,i);

			sprintf(hName3,"%s_Rebinplusplus%d",hName1,i);
			sprintf(hName4,"%s_Rebinminusminus%d",hName1,i);
			hLevyQ_plus_plus[i] =  new TH3F(hName3,hTitle,800,0,1,15,0,15,15,0,15);             
			hLevyQ_minus_minus[i] =  new TH3F(hName4,hTitle,800,0,1,15,0,15,15,0,15);             

			hLevyQ_plus_plus_mix[i] =  new TH3F(hName3_mix,hTitle,800,0,1,15,0,15,15,0,15);             
			hLevyQ_minus_minus_mix[i] =  new TH3F(hName4_mix,hTitle,800,0,1,15,0,15,15,0,15);             
		}
	}
	cout<<"to define 3DCF histograms"<<endl;
	TH3F *hq_same_plusplus[kCentBin][kyBin][kktBin];
	TH3F *hq_mix_plusplus[kCentBin][kyBin][kktBin];
	TH3F *hq_same_minusminus[kCentBin][kyBin][kktBin];
	TH3F *hq_mix_minusminus[kCentBin][kyBin][kktBin];
	TProfile3D *hcoul_plusplus_mix[kCentBin][kyBin][kktBin];
	TProfile3D *hcoul_minusminus_mix[kCentBin][kyBin][kktBin];
	if(index_3DCF==1){
		cout<<"to define 3DCF histograms Loop"<<endl;
		for(int i=0;i<kCentBin;i++){
			for(int j=0;j<kyBin;j++){
				for(int k=0;k<kktBin;k++){
					TString hname21=Form("hqSame_plusplus_cent%dy%dkt%d",i,j,k);
					TString hname22=Form("hqMix_plusplus_cent%dy%dkt%d",i,j,k);
					TString title21=Form("%s;qout;qside;qlong",hname21.Data());
					TString title22=Form("%s;qout;qside;qlong",hname22.Data());
					hq_same_plusplus[i][j][k] = new TH3F(hname21,title21,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_mix_plusplus[i][j][k] = new TH3F(hname22,title22,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					TString hname31=Form("hqSame_minusminus_cent%dy%dkt%d",i,j,k);
					TString hname32=Form("hqMix_minusminus_cent%dy%dkt%d",i,j,k);
					TString title31=Form("%s;qout;qside;qlong",hname31.Data());
					TString title32=Form("%s;qout;qside;qlong",hname32.Data());
					hq_same_minusminus[i][j][k] = new TH3F(hname31,title31,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_mix_minusminus[i][j][k] = new TH3F(hname32,title32,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
				}
			}
		}
		for(int i=0;i<kCentBin;i++){
			for(int j=0;j<kyBin;j++){
				for(int k=0;k<kktBin;k++){
					TString hname12=Form("hcoul_plusplus_mix_cent%dy%dkt%d",i,j,k);
					hcoul_plusplus_mix[i][j][k] = new TProfile3D(hname12,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					TString hname22=Form("hcoul_minusminus_mix_cent%dy%dkt%d",i,j,k);
					hcoul_minusminus_mix[i][j][k] = new TProfile3D(hname22,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
				}
			}
		}
	}

	cout<<"to define 1DmomshiftCF histograms "<<endl;
	TH3F * hdelta_mom_cms_Q_plus_minus[kCentBin+1][kDeltamombin];
	TH3F * hdelta_mom_cms_Q_plus_plus[kCentBin+1][kDeltamombin];
	TH3F * hdelta_mom_cms_Q_minus_minus[kCentBin+1][kDeltamombin];
	TH3F * hdelta_mom_cms_Q_plus_minus_mix[kCentBin+1][kDeltamombin];
	TH3F * hdelta_mom_cms_Q_plus_plus_mix[kCentBin+1][kDeltamombin];
	TH3F * hdelta_mom_cms_Q_minus_minus_mix[kCentBin+1][kDeltamombin];
	if(index_1Dmomshift==1){
		cout<<"to define 1DmomshiftCF histograms Loop"<<endl;
		for(Int_t i=0;i<kCentBin;i++)
		{ 
			for(Int_t j=0;j<kDeltamombin;j++){
				char *hName1="hdelta_mom_cms_Q";
				char hName[200];
				char hName2[200];
				char hName3[200];
				char hName4[200];

				char hName2_mix[200];
				char hName3_mix[200];
				char hName4_mix[200];
				char *hTitle0="Centrality";    
				char hTitle[200];
				sprintf(hTitle,"%s%d",hTitle0,i);
				sprintf(hName2_mix,"%s_Rebinplusminus_mix%d%d",hName1,i,j);
				sprintf(hName3_mix,"%s_Rebinplusplus_mix%d%d",hName1,i,j);
				sprintf(hName4_mix,"%s_Rebinminusminus_mix%d%d",hName1,i,j);

				sprintf(hName2,"%s_Rebinplusminus%d%d",hName1,i,j);
				sprintf(hName3,"%s_Rebinplusplus%d%d",hName1,i,j);
				sprintf(hName4,"%s_Rebinminusminus%d%d",hName1,i,j);
				hdelta_mom_cms_Q_plus_minus[i][j] =  new TH3F(hName2,hTitle,800,0,1,kktBin,0,kktBin,kyBin+1,0,kyBin+1);             
				hdelta_mom_cms_Q_plus_plus[i][j] =  new TH3F(hName3,hTitle,800,0,1,kktBin,0,kktBin,kyBin+1,0,kyBin+1);             
				hdelta_mom_cms_Q_minus_minus[i][j] =  new TH3F(hName4,hTitle,800,0,1,kktBin,0,kktBin,kyBin+1,0,kyBin+1);             

				hdelta_mom_cms_Q_plus_minus_mix[i][j] =  new TH3F(hName2_mix,hTitle,800,0,1,kktBin,0,kktBin,kyBin+1,0,kyBin+1);             
				hdelta_mom_cms_Q_plus_plus_mix[i][j] =  new TH3F(hName3_mix,hTitle,800,0,1,kktBin,0,kktBin,kyBin+1,0,kyBin+1);             
				hdelta_mom_cms_Q_minus_minus_mix[i][j] =  new TH3F(hName4_mix,hTitle,800,0,1,kktBin,0,kktBin,kyBin+1,0,kyBin+1);             
			}
		}
	}
	cout<<"to define 3DmomshiftCF histograms"<<endl;
	TH3F *hdelta_mom_cms_q_same_plusplus[kCentBin][kyBin][kktBin][kDeltamombin];
	TH3F *hdelta_mom_cms_q_mix_plusplus[kCentBin][kyBin][kktBin][kDeltamombin];
	TH3F *hdelta_mom_cms_q_same_minusminus[kCentBin][kyBin][kktBin][kDeltamombin];
	TH3F *hdelta_mom_cms_q_mix_minusminus[kCentBin][kyBin][kktBin][kDeltamombin];
	TProfile3D *hdelta_mom_cms_coul_plusplus_mix[kCentBin][kyBin][kktBin][kDeltamombin];
	TProfile3D *hdelta_mom_cms_coul_minusminus_mix[kCentBin][kyBin][kktBin][kDeltamombin];
	if(index_3Dmomshift==1){
		cout<<"to define 3DmomshiftCF histograms Loop"<<endl;
		for(int i=0;i<kCentBin;i++){
			for(int j=0;j<kyBin;j++){
				for(int k=0;k<kktBin;k++){
					for(int m=0;m<kDeltamombin;m++){
						TString hname21=Form("hdelta_mom_cms_qSame_plusplus_cent%dy%dkt%ddmom%d",i,j,k,m);
						TString hname22=Form("hdelta_mom_cms_qMix_plusplus_cent%dy%dkt%ddmom%d",i,j,k,m);
						TString title21=Form("%s;qout;qside;qlong",hname21.Data());
						TString title22=Form("%s;qout;qside;qlong",hname22.Data());
						hdelta_mom_cms_q_same_plusplus[i][j][k][m] = new TH3F(hname21,title21,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
						hdelta_mom_cms_q_mix_plusplus[i][j][k][m] = new TH3F(hname22,title22,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
						TString hname31=Form("hdelta_mom_cms_qSame_minusminus_cent%dy%dkt%ddmom%d",i,j,k,m);
						TString hname32=Form("hdelta_mom_cms_qMix_minusminus_cent%dy%dkt%ddmom%d",i,j,k,m);
						TString title31=Form("%s;qout;qside;qlong",hname31.Data());
						TString title32=Form("%s;qout;qside;qlong",hname32.Data());
						hdelta_mom_cms_q_same_minusminus[i][j][k][m] = new TH3F(hname31,title31,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
						hdelta_mom_cms_q_mix_minusminus[i][j][k][m] = new TH3F(hname32,title32,qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					}
				}
			}
		}
		for(int i=0;i<kCentBin;i++){
			for(int j=0;j<kyBin;j++){
				for(int k=0;k<kktBin;k++){
					for(int m=0;m<kDeltamombin;m++){
						TString hname12=Form("hdelta_mom_cms_coul_plusplus_mix_cent%dy%dkt%ddmom%d",i,j,k,m);
						hdelta_mom_cms_coul_plusplus_mix[i][j][k][m] = new TProfile3D(hname12,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
						TString hname22=Form("hdelta_mom_cms_coul_minusminus_mix_cent%dy%dkt%ddmom%d",i,j,k,m);
						hdelta_mom_cms_coul_minusminus_mix[i][j][k][m] = new TProfile3D(hname22,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					}
				}
			}
		}
	}
	//add some hist for smear
	cout<<"to define smear histograms"<<endl;
	TH3F *hq_A_ideal_plusplus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_B_ideal_plusplus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_A_smear_plusplus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_B_smear_plusplus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_A_ideal_minusminus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_B_ideal_minusminus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_A_smear_minusminus[4][kyBin][kktBin];//centbin has been rebined
	TH3F *hq_B_smear_minusminus[4][kyBin][kktBin];//centbin has been rebined
	if(smear_index==1){//smear open
		cout<<"to define smear histograms Loop"<<endl;
		for(int i=0;i<4;i++){
			for(int j=0;j<kyBin;j++){
				for(int k=0;k<kktBin;k++){
					TString hAidealplusplusname=Form("hq_A_ideal_pluspluscent%dy%dkt%d",i,j,k);
					TString hBidealplusplusname=Form("hq_B_ideal_pluspluscent%dy%dkt%d",i,j,k);
					TString hAsmearplusplusname=Form("hq_A_smear_pluspluscent%dy%dkt%d",i,j,k);
					TString hBsmearplusplusname=Form("hq_B_smear_pluspluscent%dy%dkt%d",i,j,k);
					TString hAidealminusminusname=Form("hq_A_ideal_minusminuscent%dy%dkt%d",i,j,k);
					TString hBidealminusminusname=Form("hq_B_ideal_minusminuscent%dy%dkt%d",i,j,k);
					TString hAsmearminusminusname=Form("hq_A_smear_minusminuscent%dy%dkt%d",i,j,k);
					TString hBsmearminusminusname=Form("hq_B_smear_minusminuscent%dy%dkt%d",i,j,k);
					hq_A_ideal_plusplus[i][j][k] = new TH3F(hAidealplusplusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_B_ideal_plusplus[i][j][k] = new TH3F(hBidealplusplusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_A_smear_plusplus[i][j][k] = new TH3F(hAsmearplusplusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_B_smear_plusplus[i][j][k] = new TH3F(hBsmearplusplusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_A_ideal_minusminus[i][j][k] = new TH3F(hAidealminusminusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_B_ideal_minusminus[i][j][k] = new TH3F(hBidealminusminusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_A_smear_minusminus[i][j][k] = new TH3F(hAsmearminusminusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
					hq_B_smear_minusminus[i][j][k] = new TH3F(hBsmearminusminusname,"",qbin,qmin,qmax,qbin,qmin,qmax,qbin,qmin,qmax);
				}
			}
		}
	}
	for(Int_t i=0;i<kCentBin;i++)
	{
		char centname[9][100]={"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
		char hName[100];
		char *hName0="hmCent";
		char *hName1="hQ";
		char hName2[200];
		char *hTitle0="Centrality:";    
		char hTitle[200];
		char hTitle2[200];
		sprintf(hTitle,"%s%s;Npionplus;Counts",hTitle0,centname[i]);
		sprintf(hTitle2,"%s%s;Npionminus;Counts",hTitle0,centname[i]);
		sprintf(hName,"%sNPionminus%d",hName0,i);
		NPionminus[i] = new TH1F(hName,hTitle2,40,0,40);
		sprintf(hName,"%sNPionplus%d",hName0,i);
		NPionplus[i] = new TH1F(hName,hTitle,40,0,40);
	}

	Long64_t nentries = t->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	Int_t zvert = -1;
	Int_t evt[kZBin][kCentBin] ={0};
	Int_t centbin=-1;
	cout<<"nentries "<<nentries<<endl;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = v0dst.LoadTree(jentry);
		if (ientry < 0) break;
		nb = t->GetEntry(jentry);   nbytes += nb;
		if(jentry%1000==0)cout<<jentry<<" "<<nentries<<endl;
		//cout<<"evt= "<<jentry<<" "<<nentries<<endl;

		if(index_E==30){
			if(v0dst.runnumber == 19151029 || v0dst.runnumber ==19151045 || v0dst.runnumber ==19152001 || v0dst.runnumber ==19152078 || v0dst.runnumber ==19153023 || v0dst.runnumber ==19153032 || v0dst.runnumber ==19153065 || v0dst.runnumber ==19154012 || v0dst.runnumber ==19154013 || v0dst.runnumber ==19154014 || v0dst.runnumber ==19154015 || v0dst.runnumber ==19154016 || v0dst.runnumber ==19154017 || v0dst.runnumber ==19154018 || v0dst.runnumber ==19154019 || v0dst.runnumber ==19154020 || v0dst.runnumber ==19154021 || v0dst.runnumber ==19154022 || v0dst.runnumber ==19154023 || v0dst.runnumber ==19154024 || v0dst.runnumber ==19154026 || v0dst.runnumber ==19154046 || v0dst.runnumber ==19154051 || v0dst.runnumber ==19154056)continue;
		}
		if(index_E==35){
			if(v0dst.runnumber ==20355020 || v0dst.runnumber ==20355021 || v0dst.runnumber ==21044023 || v0dst.runnumber ==21045024 || v0dst.runnumber ==21045025 || v0dst.runnumber ==21043027 || v0dst.runnumber ==21044035 || v0dst.runnumber ==21045004)continue;//official
		}
		if(index_E==39){
			if(v0dst.runnumber == 21035011 || v0dst.runnumber ==21036012)continue;//official
		}
		if(index_E==45){
			if(v0dst.runnumber ==21032001)continue;//official
		}
		if(index_E==52){
			if(v0dst.runnumber ==21034002||v0dst.runnumber ==21034007)continue;//official
		}
		//if(v0dst.runnumber == 20180005 || v0dst.runnumber ==20180006 || v0dst.runnumber ==20180019 || v0dst.runnumber ==20180025|| v0dst.runnumber ==20181016|| v0dst.runnumber ==20182034|| v0dst.runnumber ==20183001|| v0dst.runnumber ==20183013|| v0dst.runnumber ==20183014|| v0dst.runnumber ==20183019)continue;
		const UShort_t refmult = v0dst.nrefmult;
		const Double_t vz = v0dst.primvertexZ;

		///event cuts
		hrefmult->Fill(v0dst.nrefmult);         
		hgrefmult->Fill(v0dst.grefmult);  
		hrefmultvstofmatch_before_cut ->Fill(v0dst.nrefmult,v0dst.ntofmatch);
		//pile up correction
		double tofmatch=v0dst.ntofmatch;
		double mult=v0dst.nrefmult;
		if(index_E==32){
			double c[5]={-13.59,1.515,0.02816,-1.195E-4,-9.639E-7};
			double b[5]={19.48,5.428,-0.007,-2.428E-4,1.197E-7};
			if(mult<(c[0]+c[1]*tofmatch+c[2]*pow(tofmatch,2)+c[3]*pow(tofmatch,3)+c[4]*pow(tofmatch,4)))continue;
			if(mult>(b[0]+b[1]*tofmatch+b[2]*pow(tofmatch,2)+b[3]*pow(tofmatch,3)+b[4]*pow(tofmatch,4)))continue;
		}
		if(index_E==35){
			//Erik Loyd
			double a0 = -0.3646, a1 = 3.294, a2 = 0.03156, a3 = -4.758e-4, a4 = 9.955e-7;
			double b0 = 23.28, b1 = 5.247, b2 = 0.04037, b3 = -1.206e-3, b4 = 5.792e-06;
			double c0 =  -14.82, c1 = 1.583,c2 = 0.02684, c3 =4.605e-5,   c4 = -2.410e-06;
			double paraPileup[15]={a0, a1, a2, a3, a4,b0, b1, b2, b3, b4,c0, c1, c2, c3, c4,};
			double refmultcutmode = paraPileup[0] + paraPileup[1]*(tofmatch) + paraPileup[2]*pow(tofmatch,2) + paraPileup[3]*pow(tofmatch,3) + paraPileup[4]*pow(tofmatch,4);
			double refmultcutmax  = paraPileup[5] + paraPileup[6]*(tofmatch) + paraPileup[7]*pow(tofmatch,2) + paraPileup[8]*pow(tofmatch,3) + paraPileup[9]*pow(tofmatch,4);
			double refmultcutmin  = paraPileup[10] + paraPileup[11]*(tofmatch) + paraPileup[12]*pow(tofmatch,2) + paraPileup[13]*pow(tofmatch,3) + paraPileup[14]*pow(tofmatch,4);
			if( mult > refmultcutmax || mult < refmultcutmin )continue;
			if( tofmatch>=100&&mult>325 )continue;
		}
		if(index_E==39){
			//Erik Loyd run20
			double a0 = -1.573, a1 = 3.576, a2 = 0.02297, a3 = -2.734e-4, a4 = -2.380e-7;
			double b0 = 29.74, b1 = 4.421, b2 = 0.09139, b3 = -1.977e-3, b4 = 9.435e-06;
			double c0 =  -20.53, c1 = 2.557,c2 = -0.02094, c3 = 8.943e-4, c4 = -6.879e-06;
			double paraPileup[15]={a0, a1, a2, a3, a4,b0, b1, b2, b3, b4,c0, c1, c2, c3, c4,};
			double refmultcutmode = paraPileup[0] + paraPileup[1]*(tofmatch) + paraPileup[2]*pow(tofmatch,2) + paraPileup[3]*pow(tofmatch,3) + paraPileup[4]*pow(tofmatch,4);
			double refmultcutmax  = paraPileup[5] + paraPileup[6]*(tofmatch) + paraPileup[7]*pow(tofmatch,2) + paraPileup[8]*pow(tofmatch,3) + paraPileup[9]*pow(tofmatch,4);
			double refmultcutmin  = paraPileup[10] + paraPileup[11]*(tofmatch) + paraPileup[12]*pow(tofmatch,2) + paraPileup[13]*pow(tofmatch,3) + paraPileup[14]*pow(tofmatch,4);
			if( mult > refmultcutmax || mult < refmultcutmin )continue;
			if( tofmatch>=86&&mult>344 )continue;
		}
		if(index_E==45){
			//Erik Loyd
			double a0 = -2.031, a1 = 3.657, a2 = 0.02087, a3 = -1.771e-4, a4 = -8.434e-7;
			double b0 = 35.02, b1 = 3.586, b2 = 0.1368, b3 = -2.578e-3, b4 = 1.189e-05;
			double c0 =  -24.84, c1 = 3.289,c2 = -0.05722, c3 =1.491e-3,   c4 = -9.678e-06;
			double paraPileup[15]={a0, a1, a2, a3, a4,b0, b1, b2, b3, b4,c0, c1, c2, c3, c4,};
			double refmultcutmode = paraPileup[0] + paraPileup[1]*(tofmatch) + paraPileup[2]*pow(tofmatch,2) + paraPileup[3]*pow(tofmatch,3) + paraPileup[4]*pow(tofmatch,4);
			double refmultcutmax  = paraPileup[5] + paraPileup[6]*(tofmatch) + paraPileup[7]*pow(tofmatch,2) + paraPileup[8]*pow(tofmatch,3) + paraPileup[9]*pow(tofmatch,4);
			double refmultcutmin  = paraPileup[10] + paraPileup[11]*(tofmatch) + paraPileup[12]*pow(tofmatch,2) + paraPileup[13]*pow(tofmatch,3) + paraPileup[14]*pow(tofmatch,4);
			if( mult > refmultcutmax || mult < refmultcutmin )continue;
		}
		if(index_E==52){
			double a0=-4.94829,a1=4.90393,a2=0.000928248,a3=-0.000284087,a4=1.88985e-06;
			double b0=18.6707,b1=6.92307,b2=-0.0293523,b3=0.000412261,b4=-4.74922e-06;
			double c0=-14.4436,c1=-0.047413,c2=0.100793,c3=-0.00121203,c4=5.59521e-06;
			double paraPileup[15]={a0, a1, a2, a3, a4,b0, b1, b2, b3, b4,c0, c1, c2, c3, c4,};
			double refmultcutmode = paraPileup[0] + paraPileup[1]*(tofmatch) + paraPileup[2]*pow(tofmatch,2) + paraPileup[3]*pow(tofmatch,3) + paraPileup[4]*pow(tofmatch,4);
			double refmultcutmax  = paraPileup[5] + paraPileup[6]*(tofmatch) + paraPileup[7]*pow(tofmatch,2) + paraPileup[8]*pow(tofmatch,3) + paraPileup[9]*pow(tofmatch,4);
			double refmultcutmin  = paraPileup[10] + paraPileup[11]*(tofmatch) + paraPileup[12]*pow(tofmatch,2) + paraPileup[13]*pow(tofmatch,3) + paraPileup[14]*pow(tofmatch,4);
			if( mult > refmultcutmax || mult < refmultcutmin )continue;
		}
		hSelectNRefMultCorr->Fill(v0dst.nrefmult) ;
		hrefmultvstofmatch ->Fill(v0dst.nrefmult,v0dst.ntofmatch);
		Double_t vr = sqrt(v0dst.primvertexX*v0dst.primvertexX+v0dst.primvertexY*v0dst.primvertexY);
		hvxvy->Fill(v0dst.primvertexX,v0dst.primvertexY);
		hvz->Fill(v0dst.primvertexZ); 
		///find zbin
		zvert = getZbin(v0dst.primvertexZ ,index_E);
		if(zvert<0)continue;
		hevtZbin->Fill(zvert); 
		//find centbin
		centbin = getCentBin(v0dst.nrefmult , index_E);        
		if(index_E!=52&&index_E!=45&&index_E!=30)centbin = v0dst.cent;
		hevtcentbin->Fill(centbin); 
		if(centbin<0)continue;
		// assume delta_P_side=delta_P_long
		double half_sidefinaldeltamom3p0[4]={0.0051 ,0.0049 ,0.0037 ,0.0016};//3.0GeV Vinh
		double half_sidefinaldeltamom3p2[4]={0.0041 ,0.0037 ,0.0031 ,0.0015};//3.2GeV Vinh
		double half_sidefinaldeltamom3p5[4]={0.0036 ,0.0031 ,0.0021 ,0.0026};//3.5GeV Vinh
		double half_sidefinaldeltamom3p9[4]={0.0034 ,0.0025 ,0.0017 ,0.0011};//3.9GeV Vinh
		double half_sidefinaldeltamom4p5[4]={0.0030 ,0.0023 ,0.0011 ,0.0000};//4.5GeV Vinh
		double half_sidefinaldeltamom5p2[4]={0.0022 ,0.0017 ,0.0007 ,-0.0005};//5.2GeV Vinh
		double half_sidefinaldeltamom7p7[4]={0.0013 ,0.0010 ,0.0003 ,-0.0011};//7.7GeV Vinh
		if(kDeltamombin==1){
			if(centbin==8||centbin==7){//0-10%
				if(index_E==30)deltamom=half_sidefinaldeltamom3p0[0];
				if(index_E==32)deltamom=half_sidefinaldeltamom3p2[0];
				if(index_E==35)deltamom=half_sidefinaldeltamom3p5[0];
				if(index_E==39)deltamom=half_sidefinaldeltamom3p9[0];
				if(index_E==45)deltamom=half_sidefinaldeltamom4p5[0];
				if(index_E==52)deltamom=half_sidefinaldeltamom5p2[0];
				if(index_E==77)deltamom=half_sidefinaldeltamom7p7[0];
			}
			if(centbin==5||centbin==6){//10-30%
				if(index_E==30)deltamom=half_sidefinaldeltamom3p0[1];
				if(index_E==32)deltamom=half_sidefinaldeltamom3p2[1];
				if(index_E==35)deltamom=half_sidefinaldeltamom3p5[1];
				if(index_E==39)deltamom=half_sidefinaldeltamom3p9[1];
				if(index_E==45)deltamom=half_sidefinaldeltamom4p5[1];
				if(index_E==52)deltamom=half_sidefinaldeltamom5p2[1];
				if(index_E==77)deltamom=half_sidefinaldeltamom7p7[1];
			}
			if(centbin==3||centbin==4){//30-50%
				if(index_E==30)deltamom=half_sidefinaldeltamom3p0[2];
				if(index_E==32)deltamom=half_sidefinaldeltamom3p2[2];
				if(index_E==35)deltamom=half_sidefinaldeltamom3p5[2];
				if(index_E==39)deltamom=half_sidefinaldeltamom3p9[2];
				if(index_E==45)deltamom=half_sidefinaldeltamom4p5[2];
				if(index_E==52)deltamom=half_sidefinaldeltamom5p2[2];
				if(index_E==77)deltamom=half_sidefinaldeltamom7p7[2];
			}
			if(centbin==0||centbin==1||centbin==2){//50-80%
				if(index_E==30)deltamom=half_sidefinaldeltamom3p0[3];
				if(index_E==32)deltamom=half_sidefinaldeltamom3p2[3];
				if(index_E==35)deltamom=half_sidefinaldeltamom3p5[3];
				if(index_E==39)deltamom=half_sidefinaldeltamom3p9[3];
				if(index_E==45)deltamom=half_sidefinaldeltamom4p5[3];
				if(index_E==52)deltamom=half_sidefinaldeltamom5p2[3];
				if(index_E==77)deltamom=half_sidefinaldeltamom7p7[3];
			}
		}

		evt[zvert][centbin]=evt[zvert][centbin]+1;
		int NMIX=5;
		int myrandom = gRandom->Integer(NMIX)+1;
		Int_t counterpionminus =0;

		for(int i = 0; i<v0dst.nPionminus;i++)
		{
			double nsigmapion=v0dst.nsigmapion_Pionminus[i];
			double nsigmaproton=v0dst.nsigmaproton_Pionminus[i];
			double nsigmakaon=v0dst.nsigmakaon_Pionminus[i];
			double nsigmaelectron=v0dst.nsigmaelectron_Pionminus[i];
			double p_px= v0dst.px_Pionminus[i];
			double p_py= v0dst.py_Pionminus[i];
			double p_pz= v0dst.pz_Pionminus[i];
			double pt_p= sqrt(p_px*p_px+p_py*p_py);
			double p= sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz);
			double m2= v0dst.mass2pion_Pionminus[i];
			double beta=v0dst.Betapion_Pionminus[i];
			if(index_E==35||index_E==39||index_E==45||index_E==52){
				double em2= v0dst.emass2pion_Pionminus[i];
				if(em2>-99)m2=em2;
				double ebeta=v0dst.eBetapion_Pionminus[i];
				if(ebeta>-99)beta=ebeta;
			}
			double beta_expected=p/sqrt(pdgmass_pionminus*pdgmass_pionminus+p*p);
			double shift=getshiftminus(p,index_E);
			TLorentzVector P_Pionminus;
			P_Pionminus.SetXYZM(p_px,p_py,p_pz,pdgmass_pionminus);
			double Rap =P_Pionminus.Rapidity();
			Rap =-(Rap+ycm);
			if(index_E==77)Rap=-Rap;
			double eta= P_Pionminus.Eta();
			if(eta>etacutmax||eta<etacutmin)continue;
			if(v0dst.nhitsFit_Pionminus[i]<=nhitsfitmax)continue;
			hnsigmaminus->Fill(v0dst.nsigmapion_Pionminus[i]);
			hdedxminus->Fill(v0dst.dedx_Pionminus[i]);
			hdcaminus->Fill(v0dst.dca_Pionminus[i]);
			if(v0dst.dca_Pionminus[i]>dcamax)continue;
			if(p<pcut){
				hnsigma_P2_3Dminus->Fill(v0dst.nsigmapion_Pionminus[i],centbin,sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
			}
			if(p>pcut){
				if(m2>m2min&&m2<m2max&&fabs(1.0/beta-1/beta_expected)<0.015){
					hnsigma_P2_3Dminus->Fill(v0dst.nsigmapion_Pionminus[i],centbin,sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
				}
			}
			if(p<pcutmin)continue;
			if(p>pcutmax)continue;
			if(p<pcut){
				if(fabs(nsigmapion-shift)>2.0)continue;
				if(fabs(nsigmaproton)<2.0)continue;
				if(fabs(nsigmakaon)<2.0)continue;
				if(fabs(nsigmaelectron)<2.0)continue;
			}
			if(p>pcut){
				if(fabs(nsigmapion-shift)>3.0)continue;
				hmass2_P_3Dminus->Fill(m2,centbin,sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
				if(fabs(1.0/beta-1/beta_expected)>0.015)continue;
				if(m2 > m2max || m2 < m2min)continue;
			}
			hy_pt_pionminus->Fill(Rap,pt_p);  
			if(pt_p>ptcutmax||pt_p<ptcutmin)continue;
			if(Rap<rapcutmin||Rap>rapcutmax)continue; 

			hTofPIDminus->Fill(p, m2);
			hpdedxminus->Fill(p,v0dst.dedx_Pionminus[i]);                   
			hptdedxminus->Fill(pt_p,v0dst.dedx_Pionminus[i]);
			hmass2pionminus->Fill(m2);
			hnsigma_m2minus->Fill(v0dst.nsigmapion_Pionminus[i] , m2);
			hpxminus->Fill(p_px);
			hpyminus->Fill(p_py);
			hpzminus->Fill(p_pz);
			hphiminus->Fill(P_Pionminus.Phi());
			hetaminus->Fill(P_Pionminus.Eta());
			pionminus[zvert][centbin][0][counterpionminus].id = v0dst.id_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].Px = v0dst.px_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].Py = v0dst.py_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].Pz = v0dst.pz_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].mass = pdgmass_pionminus;
			pionminus[zvert][centbin][0][counterpionminus].padrow1 = v0dst.PadRow1_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].padrow2 = v0dst.PadRow2_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].ipadrow = v0dst.IPadRow_Pionminus[i];
			pionminus[zvert][centbin][0][counterpionminus].nhits = v0dst.nhitsFit_Pionminus[i];

			if(evt[zvert][centbin]>NMIX){

				pionminus[zvert][centbin][myrandom][counterpionminus].id = v0dst.id_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].Px = v0dst.px_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].Py = v0dst.py_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].Pz = v0dst.pz_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].mass = pdgmass_pionminus;
				pionminus[zvert][centbin][myrandom][counterpionminus].padrow1 = v0dst.PadRow1_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].padrow2 = v0dst.PadRow2_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].ipadrow = v0dst.IPadRow_Pionminus[i];
				pionminus[zvert][centbin][myrandom][counterpionminus].nhits = v0dst.nhitsFit_Pionminus[i];
			}
			if(evt[zvert][centbin]<=NMIX){
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].id = v0dst.id_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].Px = v0dst.px_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].Py = v0dst.py_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].Pz = v0dst.pz_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].mass = pdgmass_pionminus;
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].padrow1 = v0dst.PadRow1_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].padrow2 = v0dst.PadRow2_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].ipadrow = v0dst.IPadRow_Pionminus[i];
				pionminus[zvert][centbin][evt[zvert][centbin]][counterpionminus].nhits = v0dst.nhitsFit_Pionminus[i];
			}
			counterpionminus++;
		}

		NPionminus[centbin]->Fill(counterpionminus);

		pionminus[zvert][centbin][0][0].totalv0 = counterpionminus;   
		if(evt[zvert][centbin]>NMIX){
			pionminus[zvert][centbin][myrandom][0].totalv0 = counterpionminus;   
		}
		if(evt[zvert][centbin]<=NMIX){
			pionminus[zvert][centbin][evt[zvert][centbin]][0].totalv0 = counterpionminus;   
		}
		///Select pionplus
		Int_t counterpionplus =0;
		for(int i = 0; i<v0dst.nPionplus;i++)
		{
			double nsigmapion=v0dst.nsigmapion_Pionplus[i];
			double nsigmaproton=v0dst.nsigmaproton_Pionplus[i];
			double nsigmakaon=v0dst.nsigmakaon_Pionplus[i];
			double nsigmaelectron=v0dst.nsigmaelectron_Pionplus[i];
			double p_px= v0dst.px_Pionplus[i];
			double p_py= v0dst.py_Pionplus[i];
			double p_pz= v0dst.pz_Pionplus[i];
			double p= sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz);
			double pt_p= sqrt(p_px*p_px+p_py*p_py);
			double m2= v0dst.mass2pion_Pionplus[i];
			double beta=v0dst.Betapion_Pionplus[i];
			if(index_E==35||index_E==39||index_E==45||index_E==52){
				double em2= v0dst.emass2pion_Pionplus[i];
				if(em2>-99)m2=em2;
				double ebeta=v0dst.eBetapion_Pionplus[i];
				if(ebeta>-99)beta=ebeta;
			}
			double beta_expected=p/sqrt(pdgmass_pionplus*pdgmass_pionplus+p*p);
			double shift=getshiftplus(p,index_E);
			TLorentzVector P_Pionplus;
			P_Pionplus.SetXYZM(p_px,p_py,p_pz,pdgmass_pionplus);
			double Rap =P_Pionplus.Rapidity();
			Rap =-(Rap+ycm);
			if(index_E==77)Rap=-Rap;
			double eta= P_Pionplus.Eta();
			if(eta>etacutmax||eta<etacutmin)continue;
			if(v0dst.nhitsFit_Pionplus[i]<=nhitsfitmax)continue;
			hnsigmaplus->Fill(v0dst.nsigmapion_Pionplus[i]);
			hdedxplus->Fill(v0dst.dedx_Pionplus[i]);
			hdcaplus->Fill(v0dst.dca_Pionplus[i]);
			if(v0dst.dca_Pionplus[i]>dcamax)continue;
			if(p<pcut){
				hnsigma_P2_3Dplus->Fill(v0dst.nsigmapion_Pionplus[i],centbin,sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
			}
			if(p>pcut){
				if(m2>m2min&&m2<m2max&&fabs(1.0/beta-1/beta_expected)<0.015){
					hnsigma_P2_3Dplus->Fill(v0dst.nsigmapion_Pionplus[i],centbin,sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
				}
			}
			if(p<pcutmin)continue;
			if(p>pcutmax)continue;
			if(p<pcut){
				if(fabs(nsigmapion-shift)>2)continue;
				if(fabs(nsigmaproton)<2.0)continue;
				if(fabs(nsigmakaon)<2.0)continue;
				if(fabs(nsigmaelectron)<2.0)continue;
			}
			if(p>pcut){
				if(fabs(nsigmapion-shift)>3.0)continue;
				hmass2_P_3Dplus->Fill(m2,centbin,sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
				if(fabs(1.0/beta-1/beta_expected)>0.015)continue;
				if(m2 > m2max || m2 < m2min)continue;
			}
			hy_pt_pionplus->Fill(Rap,pt_p);  
			if(pt_p>ptcutmax||pt_p<ptcutmin)continue;
			if(Rap<rapcutmin||Rap>rapcutmax)continue; 
			hTofPIDplus->Fill(p,m2);
			hpdedxplus->Fill(p,v0dst.dedx_Pionplus[i]);                   
			hptdedxplus->Fill(pt_p,v0dst.dedx_Pionplus[i]);
			hmass2pionplus->Fill(m2);
			hnsigma_m2plus->Fill(v0dst.nsigmapion_Pionplus[i] , m2);
			hpxplus->Fill(p_px);
			hpyplus->Fill(p_py);
			hpzplus->Fill(p_pz);
			hphiplus->Fill(P_Pionplus.Phi());
			hetaplus->Fill(P_Pionplus.Eta());

			pionplus[zvert][centbin][0][counterpionplus].id = v0dst.id_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].Px = v0dst.px_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].Py = v0dst.py_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].Pz = v0dst.pz_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].mass = pdgmass_pionplus;
			pionplus[zvert][centbin][0][counterpionplus].padrow1 = v0dst.PadRow1_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].padrow2 = v0dst.PadRow2_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].ipadrow = v0dst.IPadRow_Pionplus[i];
			pionplus[zvert][centbin][0][counterpionplus].nhits = v0dst.nhitsFit_Pionplus[i];
			if(evt[zvert][centbin]>NMIX){
				pionplus[zvert][centbin][myrandom][counterpionplus].id = v0dst.id_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].Px = v0dst.px_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].Py = v0dst.py_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].Pz = v0dst.pz_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].mass = pdgmass_pionplus;
				pionplus[zvert][centbin][myrandom][counterpionplus].padrow1 = v0dst.PadRow1_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].padrow2 = v0dst.PadRow2_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].ipadrow = v0dst.IPadRow_Pionplus[i];
				pionplus[zvert][centbin][myrandom][counterpionplus].nhits = v0dst.nhitsFit_Pionplus[i];
			}
			if(evt[zvert][centbin]<=NMIX){
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].id = v0dst.id_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].Px = v0dst.px_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].Py = v0dst.py_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].Pz = v0dst.pz_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].mass = pdgmass_pionplus;
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].padrow1 = v0dst.PadRow1_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].padrow2 = v0dst.PadRow2_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].ipadrow = v0dst.IPadRow_Pionplus[i];
				pionplus[zvert][centbin][evt[zvert][centbin]][counterpionplus].nhits = v0dst.nhitsFit_Pionplus[i];
			}
			counterpionplus++;
		}

		NPionplus[centbin]->Fill(counterpionplus);
		pionplus[zvert][centbin][0][0].totalv0 = counterpionplus;   
		if(evt[zvert][centbin]>NMIX){
			pionplus[zvert][centbin][myrandom][0].totalv0 = counterpionplus;   
		}
		if(evt[zvert][centbin]<=NMIX){
			pionplus[zvert][centbin][evt[zvert][centbin]][0].totalv0 = counterpionplus;   
		}

		//star to do same event
		Int_t N1=0,N2=0;//N1 is pionplus N2 is pionminus

		N1 = pionplus[zvert][centbin][0][0].totalv0;
		N2 = pionminus[zvert][centbin][0][0].totalv0;

		//pion+ pion-
		for(Int_t iPLUS = 0; iPLUS<N1;iPLUS++)
		{
			for(Int_t iMINUS=0;iMINUS<N2;iMINUS++)
			{
				//make sure pionplus and pionminus are coming from same event												
				if(pionplus[zvert][centbin][0][iPLUS].id==pionminus[zvert][centbin][0][iMINUS].id)continue;

				TVector3 v1,v1t,v2,v2t;
				TLorentzVector LA1,LA2,Qvect;
				Double_t kstar=-999,qout=-999,qside=-999,qlong=-999,Q1=-999;
				double rand=gRandom->Uniform(0,1);
				if(rand<0.5){
					v1.SetXYZ(pionminus[zvert][centbin][0][iMINUS].Px,pionminus[zvert][centbin][0][iMINUS].Py,pionminus[zvert][centbin][0][iMINUS].Pz);
					v2.SetXYZ(pionplus[zvert][centbin][0][iPLUS].Px,pionplus[zvert][centbin][0][iPLUS].Py,pionplus[zvert][centbin][0][iPLUS].Pz);
					LA1.SetVectM(v1,pionminus[zvert][centbin][0][iMINUS].mass);
					LA2.SetVectM(v2,pionplus[zvert][centbin][0][iPLUS].mass);
				}

				if(rand>0.5){
					v2.SetXYZ(pionminus[zvert][centbin][0][iMINUS].Px,pionminus[zvert][centbin][0][iMINUS].Py,pionminus[zvert][centbin][0][iMINUS].Pz);
					v1.SetXYZ(pionplus[zvert][centbin][0][iPLUS].Px,pionplus[zvert][centbin][0][iPLUS].Py,pionplus[zvert][centbin][0][iPLUS].Pz);
					LA2.SetVectM(v2,pionminus[zvert][centbin][0][iMINUS].mass);
					LA1.SetVectM(v1,pionplus[zvert][centbin][0][iPLUS].mass);
				}


				double rap1=-(LA1.Rapidity()+ycm),rap2=-(LA2.Rapidity()+ycm),pairrap=-((LA1+LA2).Rapidity()+ycm);
				if(index_E>=77){
					rap1=-rap1;
					rap2=-rap2;
					pairrap=-pairrap;
				}
				double kt = 0.5*sqrt((v1.x()+v2.x())*(v1.x()+v2.x())+(v1.y()+v2.y())*(v1.y()+v2.y()));
				int ybin1,ybin2,ybinpair;
				ybin1=getsinglerapbin(rap1,index_E);
				ybin2=getsinglerapbin(rap2,index_E);
				ybinpair=getpairrapbin(pairrap,index_E);
				hpairy_kt_pionplusminus ->Fill(pairrap,kt);
				if(rapindex==1&&ybinpair<0)continue;
				if(rapindex==0&&ybin1<0)continue;
				if(rapindex==0&&ybin2<0)continue;

				double ktbin=-1;
				if(kt>ktcutmax||kt<ktcutmin)continue;
				if(kt>ktcut1&&kt<ktcut2)ktbin=0.5;
				if(kt>ktcut2&&kt<ktcut3)ktbin=1.5;
				if(kt>ktcut3&&kt<ktcut4)ktbin=2.5;
				if(kt>ktcut4&&kt<ktcut5)ktbin=3.5;
				int Centrality=centbin;
				int KTBin=ktbin-0.5;

				Qvect = (LA1-LA2);
				double qinv2=Qvect.Mag2();
				double qinv=sqrt(fabs(qinv2));
				kstar= qinv/2.0;
				if(rapindex==1){
					if(index_1DCF==1)hQ_plus_minus[centbin]->Fill(kstar,ktbin,ybinpair+0.5);//last bin -->> allrap range
					if(qinv<0.2)haverage_kt_plusminus[centbin][ybinpair]->Fill(ktbin,kt);
					if(qinv<0.2)hkt_plusminus[centbin][ybinpair]->Fill(kt);
				}
				if(rapindex==0){
					if(index_1DCF==1)hQ_plus_minus[centbin]->Fill(kstar,ktbin,14.5);//last bin -->> allrap range
					if(qinv<0.2)haverage_kt_plusminus[centbin][2]->Fill(ktbin,kt);//last bin -->> allrap range
					if(qinv<0.2)hkt_plusminus[centbin][2]->Fill(kt);//last bin -->> allrap range
					if(ybin1==ybin2){
						if(index_1DCF==1)hQ_plus_minus[centbin]->Fill(kstar,ktbin,ybin1+0.5);
						if(qinv<0.2)haverage_kt_plusminus[centbin][ybin1]->Fill(ktbin,kt);
						if(qinv<0.2)hkt_plusminus[centbin][ybin1]->Fill(kt);
						for(int i=0;i<kDeltamombin;i++){
							TLorentzVector newLA1_cms,newLA2_cms;
							if(rand<0.5){
								newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),(i+1)*deltamom);
								newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),-(i+1)*deltamom);
								if(index_E>=77){
									newLA1_cms=getnewfourmom(LA1,(i+1)*deltamom);
									newLA2_cms=getnewfourmom(LA2,-(i+1)*deltamom);
								}
								double newkstarcms=0.5*sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
								if(index_1Dmomshift==1)hdelta_mom_cms_Q_plus_minus[centbin][i]->Fill(newkstarcms,ktbin,ybin1);

							}
							if(rand>0.5){
								newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),-(i+1)*deltamom);
								newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),(i+1)*deltamom);
								if(index_E>=77){
									newLA1_cms=getnewfourmom(LA1,-(i+1)*deltamom);
									newLA2_cms=getnewfourmom(LA2,(i+1)*deltamom);
								}
								double newkstarcms=0.5*sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
								if(index_1Dmomshift==1)hdelta_mom_cms_Q_plus_minus[centbin][i]->Fill(newkstarcms,ktbin,ybin1);

							}
						}
					}
				}

			}
		}

		////pion+ pion+
		for(Int_t iPLUS1 = 0; iPLUS1<N1;iPLUS1++)
		{
			TVector3 v22;
			TLorentzVector LA22;
			v22.SetXYZ(pionplus[zvert][centbin][0][iPLUS1].Px,pionplus[zvert][centbin][0][iPLUS1].Py,pionplus[zvert][centbin][0][iPLUS1].Pz);
			LA22.SetVectM(v22,pionplus[zvert][centbin][0][iPLUS1].mass);
			double rap22=-(LA22.Rapidity()+ycm);
			hy_pt_pionplus_same->Fill(rap22,v22.Pt());  
			hpxplus_same->Fill(v22.Px());
			hpyplus_same->Fill(v22.Py());
			hpzplus_same->Fill(v22.Pz());
			hphiplus_same->Fill(v22.Phi());
			hetaplus_same->Fill(v22.Eta());
			for(Int_t iPLUS2=iPLUS1+1;iPLUS2<N1;iPLUS2++)
			{
				if(iPLUS1==iPLUS2)continue;
				//make sure pionplus and pionplus are coming from same event												
				if(pionplus[zvert][centbin][0][iPLUS1].id==pionplus[zvert][centbin][0][iPLUS2].id)continue;

				TVector3 v1,v1t,v2,v2t;
				TLorentzVector LA1,LA2,Qvect;
				Double_t kstar=-999,qout=-999,qside=-999,qlong=-999,Q1=-999;
				double rand=gRandom->Uniform(0,1);
				if(rand<0.5){
					v1.SetXYZ(pionplus[zvert][centbin][0][iPLUS2].Px,pionplus[zvert][centbin][0][iPLUS2].Py,pionplus[zvert][centbin][0][iPLUS2].Pz);
					v2.SetXYZ(pionplus[zvert][centbin][0][iPLUS1].Px,pionplus[zvert][centbin][0][iPLUS1].Py,pionplus[zvert][centbin][0][iPLUS1].Pz);
					LA1.SetVectM(v1,pionplus[zvert][centbin][0][iPLUS2].mass);
					LA2.SetVectM(v2,pionplus[zvert][centbin][0][iPLUS1].mass);
				}
				if(rand>0.5){
					v2.SetXYZ(pionplus[zvert][centbin][0][iPLUS2].Px,pionplus[zvert][centbin][0][iPLUS2].Py,pionplus[zvert][centbin][0][iPLUS2].Pz);
					v1.SetXYZ(pionplus[zvert][centbin][0][iPLUS1].Px,pionplus[zvert][centbin][0][iPLUS1].Py,pionplus[zvert][centbin][0][iPLUS1].Pz);
					LA2.SetVectM(v2,pionplus[zvert][centbin][0][iPLUS2].mass);
					LA1.SetVectM(v1,pionplus[zvert][centbin][0][iPLUS1].mass);
				}



				double rap1=-(LA1.Rapidity()+ycm),rap2=-(LA2.Rapidity()+ycm),pairrap=-((LA1+LA2).Rapidity()+ycm);
				if(index_E>=77){
					rap1=-rap1;
					rap2=-rap2;
					pairrap=-pairrap;
				}
				double kt = 0.5*sqrt((v1.x()+v2.x())*(v1.x()+v2.x())+(v1.y()+v2.y())*(v1.y()+v2.y()));
				int ybin1,ybin2,ybinpair;
				ybin1=getsinglerapbin(rap1,index_E);
				ybin2=getsinglerapbin(rap2,index_E);
				ybinpair=getpairrapbin(pairrap,index_E);
				hpairy_kt_pionplusplus ->Fill(pairrap,kt);
				if(rapindex==1&&ybinpair<0)continue;
				if(rapindex==0&&ybin1<0)continue;
				if(rapindex==0&&ybin2<0)continue;
				Qvect = (LA1-LA2);
				double qinv2=Qvect.Mag2();
				double qinv=sqrt(fabs(qinv2));
				kstar= qinv/2.0;
				int NHITS1=pionplus[zvert][centbin][0][iPLUS1].nhits;
				int NHITS2=pionplus[zvert][centbin][0][iPLUS2].nhits;
				double sl= getSL(pionplus[zvert][centbin][0][iPLUS1].padrow1,pionplus[zvert][centbin][0][iPLUS1].padrow2,pionplus[zvert][centbin][0][iPLUS1].ipadrow,pionplus[zvert][centbin][0][iPLUS1].nhits,pionplus[zvert][centbin][0][iPLUS2].padrow1,pionplus[zvert][centbin][0][iPLUS2].padrow2,pionplus[zvert][centbin][0][iPLUS2].ipadrow,pionplus[zvert][centbin][0][iPLUS2].nhits,index_E);
				hqinv_sl_plusplus_same->Fill(qinv,sl); 
				if(sl<=slcutmin||sl>=slcutmax)continue;
				float dphi = v1.Phi()-v2.Phi();               
				dphi = atan2(sin(dphi),cos(dphi)) ;
				float deta = v1.Eta()-v2.Eta();
				double ktbin=-1;
				if(kt>ktcutmax||kt<ktcutmin)continue;
				if(kt>ktcut1&&kt<ktcut2)ktbin=0.5;
				if(kt>ktcut2&&kt<ktcut3)ktbin=1.5;
				if(kt>ktcut3&&kt<ktcut4)ktbin=2.5;
				if(kt>ktcut4&&kt<ktcut5)ktbin=3.5;
				int Centrality=centbin;
				int KTBin=ktbin-0.5;

				int phistar_index=0;
				for(int i=0;i<nTpcR;i++){
					if(i<nRmin||i>nRmax)continue;
					double dphistar=getphistar(LA1,LA2,1,1,-0.5,TpcR[i]);
					if(fabs(dphistar)<dphistarcut){
						phistar_index=1;
						break;
					}
					if(rapindex==1&&index_E==52&&ybinpair>=8&&fabs(dphistar)<(dphistarcut+0.03)){
						phistar_index=1;
						break;
					}
				}
				if(phistar_index==1&&fabs(deta)<detacut2)continue;
				if(index_Hphistar==1){
					for(int i=0;i<nTpcR;i++){
						double dphistar=getphistar(LA1,LA2,1,1,-0.5,TpcR[i]);
						if(rapindex==0){
							if(qinv<0.1)hdphistardeta_plusplus[2][KTBin][i]->Fill(dphistar,deta);
							if(ybin1==ybin2&&qinv<0.1)hdphistardeta_plusplus[ybin1][KTBin][i]->Fill(dphistar,deta);
						}
						if(rapindex==1){
							if(qinv<0.1)hdphistardeta_plusplus[ybinpair][KTBin][i]->Fill(dphistar,deta);
						}
					}
				}
				if(index_1DCF==1){
					if(find_pair_cut_index==1){
						if(sl>-0.5&&sl<0.9)hQ_plus_plus[centbin]->Fill(qinv,ktbin,1.5);
						if(sl>-0.5&&sl<0.8)hQ_plus_plus[centbin]->Fill(qinv,ktbin,2.5);
						if(sl>-0.5&&sl<0.7)hQ_plus_plus[centbin]->Fill(qinv,ktbin,3.5);
						if(sl>-0.5&&sl<0.6)hQ_plus_plus[centbin]->Fill(qinv,ktbin,4.5);
						if(sl>-0.5&&sl<0.5)hQ_plus_plus[centbin]->Fill(qinv,ktbin,5.5);
						if(sl>-0.5&&sl<0.4)hQ_plus_plus[centbin]->Fill(qinv,ktbin,6.5);
						if(sl>-0.5&&sl<0.3)hQ_plus_plus[centbin]->Fill(qinv,ktbin,7.5);
						if(sl>-0.5&&sl<0.2)hQ_plus_plus[centbin]->Fill(qinv,ktbin,8.5);
						if(sl>-0.5&&sl<0.1)hQ_plus_plus[centbin]->Fill(qinv,ktbin,9.5);
						if(sl>-0.5&&sl<0.0)hQ_plus_plus[centbin]->Fill(qinv,ktbin,10.5);
					}
				}
				TVector3 Qosl=getQosl_LCMS(LA1,LA2);
				qout=Qosl.X();
				qside=Qosl.Y();
				qlong=Qosl.Z();
				//for plusplus qout*qside>0
				//for minusminus qout*qside<0
				if(qoutqside_index==1&&(qout*qside)>0)continue;
				if(qoutqside_index2==1){
					if(qout<0&&qside<0&&qlong>0)continue;
					if(qout>0&&qside>0&&qlong<0)continue;
				}
				double Q_LCMS=getQ_LCMS(LA1,LA2);
				if(rapindex==1){
					if(index_1DCF==1)hQ_plus_plus[centbin]->Fill(qinv,ktbin,ybinpair+0.5);//last bin -->> allrap range
					if(index_1DCF==1)hLevyQ_plus_plus[centbin]->Fill(Q_LCMS,ktbin,ybinpair+0.5);
					if(qinv<0.2)haverage_kt_plusplus[centbin][ybinpair]->Fill(ktbin,kt);
					if(qinv<0.2)hkt_plusplus[centbin][ybinpair]->Fill(kt);
					if(index_3DCF==1)hq_same_plusplus[centbin][ybinpair][KTBin]->Fill(qout,qside,qlong);
				}
				if(rapindex==0){
					if(index_1DCF==1)hQ_plus_plus[centbin]->Fill(qinv,ktbin,14.5);//last bin -->> allrap range
					if(index_1DCF==1)hLevyQ_plus_plus[centbin]->Fill(Q_LCMS,ktbin,14.5);
					if(qinv<0.2)haverage_kt_plusplus[centbin][2]->Fill(ktbin,kt);//last bin -->> allrap range
					if(qinv<0.2)hkt_plusplus[centbin][2]->Fill(kt);//last bin -->> allrap range
					if(index_3DCF==1)hq_same_plusplus[centbin][2][KTBin]->Fill(qout,qside,qlong);
					if(ybin1==ybin2){
						if(index_1DCF==1)hQ_plus_plus[centbin]->Fill(qinv,ktbin,ybin1+0.5);
						if(index_1DCF==1)hLevyQ_plus_plus[centbin]->Fill(Q_LCMS,ktbin,ybin1+0.5);
						if(qinv<0.2)haverage_kt_plusplus[centbin][ybin1]->Fill(ktbin,kt);
						if(qinv<0.2)hkt_plusplus[centbin][ybin1]->Fill(kt);
						if(index_3DCF==1)hq_same_plusplus[centbin][ybin1][KTBin]->Fill(qout,qside,qlong);
						for(int i=0;i<kDeltamombin;i++){
							TLorentzVector newLA1_cms,newLA2_cms;
							newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),-(i+1)*deltamom);
							newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),-(i+1)*deltamom);
							if(index_E>=77){
								newLA1_cms=getnewfourmom(LA1,-(i+1)*deltamom);
								newLA2_cms=getnewfourmom(LA2,-(i+1)*deltamom);
							}
							double newqinvcms=sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
							if(index_1Dmomshift==1)hdelta_mom_cms_Q_plus_plus[centbin][i]->Fill(newqinvcms,ktbin,ybin1);
							TVector3 newQoslcms=getQosl_LCMS(newLA1_cms,newLA2_cms);
							double newqoutcms=newQoslcms.X();
							double newqsidecms=newQoslcms.Y();
							double newqlongcms=newQoslcms.Z();
							if(index_3Dmomshift==1)hdelta_mom_cms_q_same_plusplus[centbin][ybin1][KTBin][i]->Fill(newqoutcms,newqsidecms,newqlongcms);
						}
						if(qinv<0.1)hdphiplus->Fill(dphi);
						if(qinv<0.1)hdetaplus->Fill(deta);
						if(qinv<0.1)hdmomplus->Fill(v1.Mag()-v2.Mag());
					}
				}
			}
		}

		//pion- pion-
		for(Int_t iMINUS1 = 0; iMINUS1<N2;iMINUS1++)
		{
			TVector3 v22;
			TLorentzVector LA22;
			v22.SetXYZ(pionminus[zvert][centbin][0][iMINUS1].Px,pionminus[zvert][centbin][0][iMINUS1].Py,pionminus[zvert][centbin][0][iMINUS1].Pz);
			LA22.SetVectM(v22,pionminus[zvert][centbin][0][iMINUS1].mass);
			double rap22=-(LA22.Rapidity()+ycm);
			hy_pt_pionminus_same->Fill(rap22,v22.Pt());  
			hpxminus_same->Fill(v22.Px());
			hpyminus_same->Fill(v22.Py());
			hpzminus_same->Fill(v22.Pz());
			hphiminus_same->Fill(v22.Phi());
			hetaminus_same->Fill(v22.Eta());
			for(Int_t iMINUS2=iMINUS1+1;iMINUS2<N2;iMINUS2++)
			{
				if(iMINUS1==iMINUS2)continue;
				//make sure pionminus and pionminus are coming from same event												
				if(pionminus[zvert][centbin][0][iMINUS1].id==pionminus[zvert][centbin][0][iMINUS2].id)continue;

				TVector3 v1,v1t,v2,v2t;
				TLorentzVector LA1,LA2,Qvect;
				Double_t kstar=-999,qout=-999,qside=-999,qlong=-999,Q1=-999;
				double rand=gRandom->Uniform(0,1);
				if(rand<0.5){
					v1.SetXYZ(pionminus[zvert][centbin][0][iMINUS2].Px,pionminus[zvert][centbin][0][iMINUS2].Py,pionminus[zvert][centbin][0][iMINUS2].Pz);
					v2.SetXYZ(pionminus[zvert][centbin][0][iMINUS1].Px,pionminus[zvert][centbin][0][iMINUS1].Py,pionminus[zvert][centbin][0][iMINUS1].Pz);
					LA1.SetVectM(v1,pionminus[zvert][centbin][0][iMINUS2].mass);
					LA2.SetVectM(v2,pionminus[zvert][centbin][0][iMINUS1].mass);
				}
				if(rand>0.5){
					v2.SetXYZ(pionminus[zvert][centbin][0][iMINUS2].Px,pionminus[zvert][centbin][0][iMINUS2].Py,pionminus[zvert][centbin][0][iMINUS2].Pz);
					v1.SetXYZ(pionminus[zvert][centbin][0][iMINUS1].Px,pionminus[zvert][centbin][0][iMINUS1].Py,pionminus[zvert][centbin][0][iMINUS1].Pz);
					LA2.SetVectM(v2,pionminus[zvert][centbin][0][iMINUS2].mass);
					LA1.SetVectM(v1,pionminus[zvert][centbin][0][iMINUS1].mass);
				}


				double rap1=-(LA1.Rapidity()+ycm),rap2=-(LA2.Rapidity()+ycm),pairrap=-((LA1+LA2).Rapidity()+ycm);
				if(index_E>=77){
					rap1=-rap1;
					rap2=-rap2;
					pairrap=-pairrap;
				}
				double kt = 0.5*sqrt((v1.x()+v2.x())*(v1.x()+v2.x())+(v1.y()+v2.y())*(v1.y()+v2.y()));
				int ybin1,ybin2,ybinpair;
				ybin1=getsinglerapbin(rap1,index_E);
				ybin2=getsinglerapbin(rap2,index_E);
				ybinpair=getpairrapbin(pairrap,index_E);
				hpairy_kt_pionminusminus ->Fill(pairrap,kt);
				if(rapindex==1&&ybinpair<0)continue;
				if(rapindex==0&&ybin1<0)continue;
				if(rapindex==0&&ybin2<0)continue;

				Qvect = (LA1-LA2);
				double qinv2=Qvect.Mag2();
				double qinv=sqrt(fabs(qinv2));
				kstar= qinv/2.0;
				int NHITS1=pionminus[zvert][centbin][0][iMINUS1].nhits;
				int NHITS2=pionminus[zvert][centbin][0][iMINUS2].nhits;
				double sl= getSL(pionminus[zvert][centbin][0][iMINUS1].padrow1,pionminus[zvert][centbin][0][iMINUS1].padrow2,pionminus[zvert][centbin][0][iMINUS1].ipadrow,pionminus[zvert][centbin][0][iMINUS1].nhits,pionminus[zvert][centbin][0][iMINUS2].padrow1,pionminus[zvert][centbin][0][iMINUS2].padrow2,pionminus[zvert][centbin][0][iMINUS2].ipadrow,pionminus[zvert][centbin][0][iMINUS2].nhits,index_E);
				hqinv_sl_minusminus_same->Fill(qinv,sl); 
				if(sl<=slcutmin||sl>=slcutmax)continue;
				float dphi = v1.Phi()-v2.Phi();               
				dphi = atan2(sin(dphi),cos(dphi)) ;
				float deta = v1.Eta()-v2.Eta();
				double ktbin=-1;
				if(kt>ktcutmax||kt<ktcutmin)continue;
				if(kt>ktcut1&&kt<ktcut2)ktbin=0.5;
				if(kt>ktcut2&&kt<ktcut3)ktbin=1.5;
				if(kt>ktcut3&&kt<ktcut4)ktbin=2.5;
				if(kt>ktcut4&&kt<ktcut5)ktbin=3.5;
				int Centrality=centbin;
				int KTBin=ktbin-0.5;
				int phistar_index=0;
				for(int i=0;i<nTpcR;i++){
					if(i<nRmin||i>nRmax)continue;
					double dphistar=getphistar(LA1,LA2,-1,-1,-0.5,TpcR[i]);
					if(fabs(dphistar)<dphistarcut){
						phistar_index=1;
						break;
					}
					if(rapindex==1&&index_E==52&&ybinpair>=8&&fabs(dphistar)<(dphistarcut+0.03)){
						phistar_index=1;
						break;
					}
				}
				if(phistar_index==1&&fabs(deta)<detacut2)continue;
				if(index_Hphistar==1){
					for(int i=0;i<nTpcR;i++){
						double dphistar=getphistar(LA1,LA2,-1,-1,-0.5,TpcR[i]);
						if(rapindex==0){
							if(qinv<0.1)hdphistardeta_minusminus[2][KTBin][i]->Fill(dphistar,deta);
							if(ybin1==ybin2&&qinv<0.1)hdphistardeta_minusminus[ybin1][KTBin][i]->Fill(dphistar,deta);
						}
						if(rapindex==1){
							if(qinv<0.1)hdphistardeta_minusminus[ybinpair][KTBin][i]->Fill(dphistar,deta);
						}
					}
				}
				if(index_1DCF==1){
					if(find_pair_cut_index==1){
						if(sl>-0.5&&sl<0.9)hQ_minus_minus[centbin]->Fill(qinv,ktbin,1.5);
						if(sl>-0.5&&sl<0.8)hQ_minus_minus[centbin]->Fill(qinv,ktbin,2.5);
						if(sl>-0.5&&sl<0.7)hQ_minus_minus[centbin]->Fill(qinv,ktbin,3.5);
						if(sl>-0.5&&sl<0.6)hQ_minus_minus[centbin]->Fill(qinv,ktbin,4.5);
						if(sl>-0.5&&sl<0.5)hQ_minus_minus[centbin]->Fill(qinv,ktbin,5.5);
						if(sl>-0.5&&sl<0.4)hQ_minus_minus[centbin]->Fill(qinv,ktbin,6.5);
						if(sl>-0.5&&sl<0.3)hQ_minus_minus[centbin]->Fill(qinv,ktbin,7.5);
						if(sl>-0.5&&sl<0.2)hQ_minus_minus[centbin]->Fill(qinv,ktbin,8.5);
						if(sl>-0.5&&sl<0.1)hQ_minus_minus[centbin]->Fill(qinv,ktbin,9.5);
						if(sl>-0.5&&sl<0.0)hQ_minus_minus[centbin]->Fill(qinv,ktbin,10.5);
					}
				}
				TVector3 Qosl=getQosl_LCMS(LA1,LA2);
				qout=Qosl.X();
				qside=Qosl.Y();
				qlong=Qosl.Z();
				//for plusplus qout*qside>0
				//for minusminus qout*qside<0
				//if(qout*qside<0)hQ_minus_minus[centbin]->Fill(qinv,ktbin,ptbin1);
				if(qoutqside_index==1&&(qout*qside)<0)continue;
				if(qoutqside_index2==1){
					if(qout<0&&qside>0&&qlong>0)continue;
					if(qout>0&&qside<0&&qlong<0)continue;
				}
				double Q_LCMS=getQ_LCMS(LA1,LA2);
				if(rapindex==1){
					if(index_1DCF==1)hQ_minus_minus[centbin]->Fill(qinv,ktbin,ybinpair+0.5);//last bin -->> allrap range
					if(index_1DCF==1)hLevyQ_minus_minus[centbin]->Fill(Q_LCMS,ktbin,ybinpair+0.5);
					if(qinv<0.2)haverage_kt_minusminus[centbin][ybinpair]->Fill(ktbin,kt);
					if(qinv<0.2)hkt_minusminus[centbin][ybinpair]->Fill(kt);
					if(index_3DCF==1)hq_same_minusminus[centbin][ybinpair][KTBin]->Fill(qout,qside,qlong);
				}
				if(rapindex==0){
					if(index_1DCF==1)hQ_minus_minus[centbin]->Fill(qinv,ktbin,14.5);//last bin -->> allrap range
					if(index_1DCF==1)hLevyQ_minus_minus[centbin]->Fill(Q_LCMS,ktbin,14.5);
					if(qinv<0.2)haverage_kt_minusminus[centbin][2]->Fill(ktbin,kt);//last bin -->> allrap range
					if(qinv<0.2)hkt_minusminus[centbin][2]->Fill(kt);//last bin -->> allrap range
					if(index_3DCF==1)hq_same_minusminus[centbin][2][KTBin]->Fill(qout,qside,qlong);
					if(ybin1==ybin2){
						if(index_1DCF==1)hQ_minus_minus[centbin]->Fill(qinv,ktbin,ybin1+0.5);
						if(index_1DCF==1)hLevyQ_minus_minus[centbin]->Fill(Q_LCMS,ktbin,ybin1+0.5);
						if(qinv<0.2)haverage_kt_minusminus[centbin][ybin1]->Fill(ktbin,kt);
						if(qinv<0.2)hkt_minusminus[centbin][ybin1]->Fill(kt);
						if(index_3DCF==1)hq_same_minusminus[centbin][ybin1][KTBin]->Fill(qout,qside,qlong);
						for(int i=0;i<kDeltamombin;i++){
							TLorentzVector newLA1_cms,newLA2_cms;
							newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),(i+1)*deltamom);
							newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),(i+1)*deltamom);
							if(index_E>=77){
								newLA1_cms=getnewfourmom(LA1,(i+1)*deltamom);
								newLA2_cms=getnewfourmom(LA2,(i+1)*deltamom);
							}
							double newqinvcms=sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
							if(index_1Dmomshift==1)hdelta_mom_cms_Q_minus_minus[centbin][i]->Fill(newqinvcms,ktbin,ybin1);
							TVector3 newQoslcms=getQosl_LCMS(newLA1_cms,newLA2_cms);
							double newqoutcms=newQoslcms.X();
							double newqsidecms=newQoslcms.Y();
							double newqlongcms=newQoslcms.Z();
							if(index_3Dmomshift==1)hdelta_mom_cms_q_same_minusminus[centbin][ybin1][KTBin][i]->Fill(newqoutcms,newqsidecms,newqlongcms);
						}
						if(qinv<0.1)hdphiminus->Fill(dphi);
						if(qinv<0.1)hdetaminus->Fill(deta);
						if(qinv<0.1)hdmomminus->Fill(v1.Mag()-v2.Mag());
					}
				}
			}
		}


		//star to do mix
		if(evt[zvert][centbin]>=NMIX)
		{
			Int_t Centrality1 = -1;
			for(Int_t iev = 0; iev< 1;iev++)//current event
			{
				for(Int_t jev = 1; jev< (NMIX+1);jev++)//different event
				{
					if(evt[zvert][centbin]==NMIX){
						if(jev==NMIX)continue;//first time to do mix,at this time [0]=[NMIX](they are same event,so we should reject it),make sure here is mixed event 
					}
					if(evt[zvert][centbin]>NMIX){
						if(jev==myrandom)continue;//now [0]=[myrandom](they are same event, so we should reject it)make sure here is mixed envent
					}
					Centrality1 = centbin;

					if(Centrality1<0)continue;
					int N1evt1 = pionplus[zvert][centbin][iev][0].totalv0;
					int N1evt2 = pionplus[zvert][centbin][jev][0].totalv0;
					int N2evt1 = pionminus[zvert][centbin][iev][0].totalv0;
					int N2evt2 = pionminus[zvert][centbin][jev][0].totalv0;

					//mix pion+ pion-
					for(Int_t iPLUS=0;iPLUS<N1evt1;iPLUS++)
					{
						for(Int_t iMINUS=0;iMINUS<N2evt2;iMINUS++)
						{
							TVector3 v1,v1t,v2,v2t;
							TLorentzVector LA1,LA2,Qvect;
							Double_t kstar=-999,qout=-999,qside=-999,qlong=-999,Q1=-999;
							double rand=gRandom->Uniform(0,1);
							if(rand<0.5){
								v1.SetXYZ(pionminus[zvert][centbin][jev][iMINUS].Px,pionminus[zvert][centbin][jev][iMINUS].Py,pionminus[zvert][centbin][jev][iMINUS].Pz);
								v2.SetXYZ(pionplus[zvert][centbin][iev][iPLUS].Px,pionplus[zvert][centbin][iev][iPLUS].Py,pionplus[zvert][centbin][iev][iPLUS].Pz);
								LA1.SetVectM(v1,pionminus[zvert][centbin][jev][iMINUS].mass);
								LA2.SetVectM(v2,pionplus[zvert][centbin][iev][iPLUS].mass);
							}
							if(rand>0.5){
								v2.SetXYZ(pionminus[zvert][centbin][jev][iMINUS].Px,pionminus[zvert][centbin][jev][iMINUS].Py,pionminus[zvert][centbin][jev][iMINUS].Pz);
								v1.SetXYZ(pionplus[zvert][centbin][iev][iPLUS].Px,pionplus[zvert][centbin][iev][iPLUS].Py,pionplus[zvert][centbin][iev][iPLUS].Pz);
								LA2.SetVectM(v2,pionminus[zvert][centbin][jev][iMINUS].mass);
								LA1.SetVectM(v1,pionplus[zvert][centbin][iev][iPLUS].mass);
							}


							double rap1=-(LA1.Rapidity()+ycm),rap2=-(LA2.Rapidity()+ycm),pairrap=-((LA1+LA2).Rapidity()+ycm);
							if(index_E>=77){
								rap1=-rap1;
								rap2=-rap2;
								pairrap=-pairrap;
							}
							int ybin1,ybin2,ybinpair;
							ybin1=getsinglerapbin(rap1,index_E);
							ybin2=getsinglerapbin(rap2,index_E);
							ybinpair=getpairrapbin(pairrap,index_E);
							if(rapindex==1&&ybinpair<0)continue;
							if(rapindex==0&&ybin1<0)continue;
							if(rapindex==0&&ybin2<0)continue;
							double kt = 0.5*sqrt((v1.x()+v2.x())*(v1.x()+v2.x())+(v1.y()+v2.y())*(v1.y()+v2.y()));
							double ktbin=-1;
							if(kt>ktcutmax||kt<ktcutmin)continue;
							if(kt>ktcut1&&kt<ktcut2)ktbin=0.5;
							if(kt>ktcut2&&kt<ktcut3)ktbin=1.5;
							if(kt>ktcut3&&kt<ktcut4)ktbin=2.5;
							if(kt>ktcut4&&kt<ktcut5)ktbin=3.5;
							int KTBin=ktbin-0.5;
							if(iPLUS==0){
								hy_pt_pionminus_mix->Fill(rap1,v1.Pt());  
								hpxminus_mix->Fill(v1.Px());
								hpyminus_mix->Fill(v1.Py());
								hpzminus_mix->Fill(v1.Pz());
								hphiminus_mix->Fill(v1.Phi());
								hetaminus_mix->Fill(v1.Eta());
							}
							if(iMINUS==0){
								hy_pt_pionplus_mix->Fill(rap2,v2.Pt());  
								hpxplus_mix->Fill(v2.Px());
								hpyplus_mix->Fill(v2.Py());
								hpzplus_mix->Fill(v2.Pz());
								hphiplus_mix->Fill(v2.Phi());
								hetaplus_mix->Fill(v2.Eta());
							}

							Qvect = (LA1-LA2);
							double qinv2=Qvect.Mag2();
							double qinv=sqrt(fabs(qinv2));
							kstar= qinv/2.0;
							if(rapindex==1){
								if(index_1DCF==1)hQ_plus_minus_mix[centbin]->Fill(kstar,ktbin,ybinpair+0.5);//last bin -->> allrap range
							}
							if(rapindex==0){
								if(index_1DCF==1)hQ_plus_minus_mix[centbin]->Fill(kstar,ktbin,14.5);//last bin -->> allrap range
								if(ybin1==ybin2){
									if(index_1DCF==1)hQ_plus_minus_mix[centbin]->Fill(kstar,ktbin,ybin1+0.5);
									for(int i=0;i<kDeltamombin;i++){
										TLorentzVector newLA1_cms,newLA2_cms;
										if(rand<0.5){
											newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),(i+1)*deltamom);
											newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),-(i+1)*deltamom);
											if(index_E>=77){
												newLA1_cms=getnewfourmom(LA1,(i+1)*deltamom);
												newLA2_cms=getnewfourmom(LA2,-(i+1)*deltamom);
											}
											double newkstarcms=0.5*sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
											if(index_1Dmomshift==1)hdelta_mom_cms_Q_plus_minus_mix[centbin][i]->Fill(newkstarcms,ktbin,ybin1);

										}
										if(rand>0.5){
											newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),-(i+1)*deltamom);
											newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),(i+1)*deltamom);
											if(index_E>=77){
												newLA1_cms=getnewfourmom(LA1,-(i+1)*deltamom);
												newLA2_cms=getnewfourmom(LA2,(i+1)*deltamom);
											}
											double newkstarcms=0.5*sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
											if(index_1Dmomshift==1)hdelta_mom_cms_Q_plus_minus_mix[centbin][i]->Fill(newkstarcms,ktbin,ybin1);

										}
									}
								}
							}
						}
					}
					//mix pion+ pion+
					for(Int_t iPLUS1=0;iPLUS1<N1evt1;iPLUS1++)
					{
						for(Int_t iPLUS2=0;iPLUS2<N1evt2;iPLUS2++)
						{
							TVector3 v1,v1t,v2,v2t;
							TLorentzVector LA1,LA2,Qvect;
							Double_t kstar=-999,qout=-999,qside=-999,qlong=-999,Q1=-999;
							double rand=gRandom->Uniform(0,1);
							if(rand<0.5){
								v1.SetXYZ(pionplus[zvert][centbin][jev][iPLUS2].Px,pionplus[zvert][centbin][jev][iPLUS2].Py,pionplus[zvert][centbin][jev][iPLUS2].Pz);
								v2.SetXYZ(pionplus[zvert][centbin][iev][iPLUS1].Px,pionplus[zvert][centbin][iev][iPLUS1].Py,pionplus[zvert][centbin][iev][iPLUS1].Pz);
								LA1.SetVectM(v1,pionplus[zvert][centbin][jev][iPLUS2].mass);
								LA2.SetVectM(v2,pionplus[zvert][centbin][iev][iPLUS1].mass);
							}
							if(rand>0.5){
								v2.SetXYZ(pionplus[zvert][centbin][jev][iPLUS2].Px,pionplus[zvert][centbin][jev][iPLUS2].Py,pionplus[zvert][centbin][jev][iPLUS2].Pz);
								v1.SetXYZ(pionplus[zvert][centbin][iev][iPLUS1].Px,pionplus[zvert][centbin][iev][iPLUS1].Py,pionplus[zvert][centbin][iev][iPLUS1].Pz);
								LA2.SetVectM(v2,pionplus[zvert][centbin][jev][iPLUS2].mass);
								LA1.SetVectM(v1,pionplus[zvert][centbin][iev][iPLUS1].mass);
							}

							double rap1=-(LA1.Rapidity()+ycm),rap2=-(LA2.Rapidity()+ycm),pairrap=-((LA1+LA2).Rapidity()+ycm);
							if(index_E>=77){
								rap1=-rap1;
								rap2=-rap2;
								pairrap=-pairrap;
							}
							int ybin1,ybin2,ybinpair;
							ybin1=getsinglerapbin(rap1,index_E);
							ybin2=getsinglerapbin(rap2,index_E);
							ybinpair=getpairrapbin(pairrap,index_E);
							if(rapindex==1&&ybinpair<0)continue;
							if(rapindex==0&&ybin1<0)continue;
							if(rapindex==0&&ybin2<0)continue;

							Qvect = (LA1-LA2);
							double qinv2=Qvect.Mag2();
							double qinv=sqrt(fabs(qinv2));
							kstar= qinv/2.0;
							int NHITS1=pionplus[zvert][centbin][iev][iPLUS1].nhits;
							int NHITS2=pionplus[zvert][centbin][jev][iPLUS2].nhits;
							double sl= getSL(pionplus[zvert][centbin][iev][iPLUS1].padrow1,pionplus[zvert][centbin][iev][iPLUS1].padrow2,pionplus[zvert][centbin][iev][iPLUS1].ipadrow,pionplus[zvert][centbin][iev][iPLUS1].nhits,pionplus[zvert][centbin][jev][iPLUS2].padrow1,pionplus[zvert][centbin][jev][iPLUS2].padrow2,pionplus[zvert][centbin][jev][iPLUS2].ipadrow,pionplus[zvert][centbin][jev][iPLUS2].nhits,index_E);
							hqinv_sl_plusplus_mix->Fill(qinv,sl); 
							if(sl<=slcutmin||sl>=slcutmax)continue;
							float dphi = v1.Phi()-v2.Phi();               
							dphi = atan2(sin(dphi),cos(dphi)) ;
							float deta = v1.Eta()-v2.Eta();
							double kt = 0.5*sqrt((v1.x()+v2.x())*(v1.x()+v2.x())+(v1.y()+v2.y())*(v1.y()+v2.y()));
							double ktbin=-1;
							if(kt>ktcutmax||kt<ktcutmin)continue;
							if(kt>ktcut1&&kt<ktcut2)ktbin=0.5;
							if(kt>ktcut2&&kt<ktcut3)ktbin=1.5;
							if(kt>ktcut3&&kt<ktcut4)ktbin=2.5;
							if(kt>ktcut4&&kt<ktcut5)ktbin=3.5;
							int KTBin=ktbin-0.5;
							int phistar_index=0;
							for(int i=0;i<nTpcR;i++){
								if(i<nRmin||i>nRmax)continue;
								double dphistar=getphistar(LA1,LA2,1,1,-0.5,TpcR[i]);
								if(fabs(dphistar)<dphistarcut){
									phistar_index=1;
									break;
								}
								if(rapindex==1&&index_E==52&&ybinpair>=8&&fabs(dphistar)<(dphistarcut+0.03)){
									phistar_index=1;
									break;
								}
							}
							if(phistar_index==1&&fabs(deta)<detacut2)continue;
							if(index_Hphistar==1){
								for(int i=0;i<nTpcR;i++){
									double dphistar=getphistar(LA1,LA2,1,1,-0.5,TpcR[i]);
									if(rapindex==0){
										if(qinv<0.1)hdphistardeta_plusplus_mix[2][KTBin][i]->Fill(dphistar,deta);
										if(ybin1==ybin2&&qinv<0.1)hdphistardeta_plusplus_mix[ybin1][KTBin][i]->Fill(dphistar,deta);
									}
									if(rapindex==1){
										if(qinv<0.1)hdphistardeta_plusplus_mix[ybinpair][KTBin][i]->Fill(dphistar,deta);
									}
								}
							}

							if(index_1DCF==1){
								if(find_pair_cut_index==1){
									if(sl>-0.5&&sl<0.9)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,1.5);
									if(sl>-0.5&&sl<0.8)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,2.5);
									if(sl>-0.5&&sl<0.7)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,3.5);
									if(sl>-0.5&&sl<0.6)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,4.5);
									if(sl>-0.5&&sl<0.5)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,5.5);
									if(sl>-0.5&&sl<0.4)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,6.5);
									if(sl>-0.5&&sl<0.3)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,7.5);
									if(sl>-0.5&&sl<0.2)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,8.5);
									if(sl>-0.5&&sl<0.1)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,9.5);
									if(sl>-0.5&&sl<0.0)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,10.5);
								}
							}

							TVector3 Qosl=getQosl_LCMS(LA1,LA2);
							qout=Qosl.X();
							qside=Qosl.Y();
							qlong=Qosl.Z();
							//for plusplus qout*qside>0
							//for minusminus qout*qside<0
							//if(qout*qside>0)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,ptbin1);
							if(qoutqside_index==1&&(qout*qside)>0)continue;
							if(qoutqside_index2==1){
								if(qout<0&&qside<0&&qlong>0)continue;
								if(qout>0&&qside>0&&qlong<0)continue;
							}
							double weight=hcoul->GetBinContent(hcoul->FindBin(qinv));
							if(qinv>0.6)weight=1;
							double Q_LCMS=getQ_LCMS(LA1,LA2);
							if(rapindex==1){
								if(index_1DCF==1)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,ybinpair+0.5);//last bin -->> allrap range
								if(index_1DCF==1)hLevyQ_plus_plus_mix[centbin]->Fill(Q_LCMS,ktbin,ybinpair+0.5);
								if(index_3DCF==1)hq_mix_plusplus[centbin][ybinpair][KTBin]->Fill(qout,qside,qlong);
								if(index_3DCF==1)hcoul_plusplus_mix[centbin][ybinpair][KTBin]->Fill(qout,qside,qlong,weight);
							}
							if(rapindex==0){
								if(index_1DCF==1)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,14.5);//last bin -->> allrap range
								if(index_1DCF==1)hLevyQ_plus_plus_mix[centbin]->Fill(Q_LCMS,ktbin,14.5);
								if(index_3DCF==1)hq_mix_plusplus[centbin][2][KTBin]->Fill(qout,qside,qlong);
								if(index_3DCF==1)hcoul_plusplus_mix[centbin][2][KTBin]->Fill(qout,qside,qlong,weight);
								if(ybin1==ybin2){
									if(index_1DCF==1)hQ_plus_plus_mix[centbin]->Fill(qinv,ktbin,ybin1+0.5);
									if(index_1DCF==1)hLevyQ_plus_plus_mix[centbin]->Fill(Q_LCMS,ktbin,ybin1+0.5);
									if(index_3DCF==1)hq_mix_plusplus[centbin][ybin1][KTBin]->Fill(qout,qside,qlong);
									if(index_3DCF==1)hcoul_plusplus_mix[centbin][ybin1][KTBin]->Fill(qout,qside,qlong,weight);
									for(int i=0;i<kDeltamombin;i++){
										TLorentzVector newLA1_cms,newLA2_cms;
										newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),-(i+1)*deltamom);
										newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),-(i+1)*deltamom);
										if(index_E>=77){
											newLA1_cms=getnewfourmom(LA1,-(i+1)*deltamom);
											newLA2_cms=getnewfourmom(LA2,-(i+1)*deltamom);
										}
										double newqinvcms=sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
										if(index_1Dmomshift==1)hdelta_mom_cms_Q_plus_plus_mix[centbin][i]->Fill(newqinvcms,ktbin,ybin1);
										TVector3 newQoslcms=getQosl_LCMS(newLA1_cms,newLA2_cms);
										double newqoutcms=newQoslcms.X();
										double newqsidecms=newQoslcms.Y();
										double newqlongcms=newQoslcms.Z();
										double newweightcms=hcoul->GetBinContent(hcoul->FindBin(newqinvcms));
										if(newqinvcms>0.6)newweightcms=1;
										if(index_3Dmomshift==1)hdelta_mom_cms_q_mix_plusplus[centbin][ybin1][KTBin][i]->Fill(newqoutcms,newqsidecms,newqlongcms);
										if(index_3Dmomshift==1)hdelta_mom_cms_coul_plusplus_mix[centbin][ybin1][KTBin][i]->Fill(newqoutcms,newqsidecms,newqlongcms,newweightcms);
									}
								}
							}
							if(smear_index==1){
								double smearweight;//for smear
								int finalcentbin=-1;
								if(Centrality1==8||Centrality1==7)finalcentbin=0;
								if(Centrality1==6||Centrality1==5)finalcentbin=1;
								if(Centrality1==4||Centrality1==3)finalcentbin=2;
								if(Centrality1==2||Centrality1==1||Centrality1==0)finalcentbin=3;
								if(finalcentbin<0)continue;
								smearweight=1-lamplus[finalcentbin][ybin1][KTBin]+lamplus[finalcentbin][ybin1][KTBin]*weight*(1+exp((-qout*qout*outplus[finalcentbin][ybin1][KTBin]*outplus[finalcentbin][ybin1][KTBin]-qside*qside*sideplus[finalcentbin][ybin1][KTBin]*sideplus[finalcentbin][ybin1][KTBin]-qlong*qlong*longplus[finalcentbin][ybin1][KTBin]*longplus[finalcentbin][ybin1][KTBin]-2*qout*qlong*ol2plus[finalcentbin][ybin1][KTBin])/0.038937929230));
								hq_A_ideal_plusplus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong,smearweight);
								hq_B_ideal_plusplus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong);
								LA1=Dosmearplus(LA1,index_E);
								LA2=Dosmearplus(LA2,index_E);
								TVector3 Qsmear=getQosl_LCMS(LA1,LA2);
								qout=Qsmear.X();
								qside=Qsmear.Y();
								qlong=Qsmear.Z();
								hq_A_smear_plusplus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong,smearweight);
								hq_B_smear_plusplus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong);
							}
						}
					}


					//mix pion- pion-
					for(Int_t iMINUS1=0; iMINUS1<N2evt1;iMINUS1++)
					{
						for(Int_t iMINUS2=0;iMINUS2<N2evt2;iMINUS2++)
						{
							TVector3 v1,v1t,v2,v2t;
							TLorentzVector LA1,LA2,Qvect;
							Double_t kstar=-999,qout=-999,qside=-999,qlong=-999,Q1=-999;
							double rand=gRandom->Uniform(0,1);
							if(rand<0.5){
								v1.SetXYZ(pionminus[zvert][centbin][jev][iMINUS2].Px,pionminus[zvert][centbin][jev][iMINUS2].Py,pionminus[zvert][centbin][jev][iMINUS2].Pz);
								v2.SetXYZ(pionminus[zvert][centbin][iev][iMINUS1].Px,pionminus[zvert][centbin][iev][iMINUS1].Py,pionminus[zvert][centbin][iev][iMINUS1].Pz);
								LA1.SetVectM(v1,pionminus[zvert][centbin][jev][iMINUS2].mass);
								LA2.SetVectM(v2,pionminus[zvert][centbin][iev][iMINUS1].mass);
							}
							if(rand>0.5){
								v2.SetXYZ(pionminus[zvert][centbin][jev][iMINUS2].Px,pionminus[zvert][centbin][jev][iMINUS2].Py,pionminus[zvert][centbin][jev][iMINUS2].Pz);
								v1.SetXYZ(pionminus[zvert][centbin][iev][iMINUS1].Px,pionminus[zvert][centbin][iev][iMINUS1].Py,pionminus[zvert][centbin][iev][iMINUS1].Pz);
								LA2.SetVectM(v2,pionminus[zvert][centbin][jev][iMINUS2].mass);
								LA1.SetVectM(v1,pionminus[zvert][centbin][iev][iMINUS1].mass);
							}

							double rap1=-(LA1.Rapidity()+ycm),rap2=-(LA2.Rapidity()+ycm),pairrap=-((LA1+LA2).Rapidity()+ycm);
							if(index_E>=77){
								rap1=-rap1;
								rap2=-rap2;
								pairrap=-pairrap;
							}
							int ybin1,ybin2,ybinpair;
							ybin1=getsinglerapbin(rap1,index_E);
							ybin2=getsinglerapbin(rap2,index_E);
							ybinpair=getpairrapbin(pairrap,index_E);
							if(rapindex==1&&ybinpair<0)continue;
							if(rapindex==0&&ybin1<0)continue;
							if(rapindex==0&&ybin2<0)continue;

							Qvect = (LA1-LA2);
							double qinv2=Qvect.Mag2();
							double qinv=sqrt(fabs(qinv2));
							kstar= qinv/2.0;
							int NHITS1=pionminus[zvert][centbin][iev][iMINUS1].nhits;
							int NHITS2=pionminus[zvert][centbin][jev][iMINUS2].nhits;
							double sl= getSL(pionminus[zvert][centbin][iev][iMINUS1].padrow1,pionminus[zvert][centbin][iev][iMINUS1].padrow2,pionminus[zvert][centbin][iev][iMINUS1].ipadrow,pionminus[zvert][centbin][iev][iMINUS1].nhits,pionminus[zvert][centbin][jev][iMINUS2].padrow1,pionminus[zvert][centbin][jev][iMINUS2].padrow2,pionminus[zvert][centbin][jev][iMINUS2].ipadrow,pionminus[zvert][centbin][jev][iMINUS2].nhits,index_E);
							hqinv_sl_minusminus_mix->Fill(qinv,sl); 
							if(sl<=slcutmin||sl>=slcutmax)continue;
							float dphi = v1.Phi()-v2.Phi();               
							dphi = atan2(sin(dphi),cos(dphi)) ;
							float deta = v1.Eta()-v2.Eta();
							double kt = 0.5*sqrt((v1.x()+v2.x())*(v1.x()+v2.x())+(v1.y()+v2.y())*(v1.y()+v2.y()));
							double ktbin=-1;
							if(kt>ktcutmax||kt<ktcutmin)continue;
							if(kt>ktcut1&&kt<ktcut2)ktbin=0.5;
							if(kt>ktcut2&&kt<ktcut3)ktbin=1.5;
							if(kt>ktcut3&&kt<ktcut4)ktbin=2.5;
							if(kt>ktcut4&&kt<ktcut5)ktbin=3.5;
							int KTBin=ktbin-0.5;
							int phistar_index=0;
							for(int i=0;i<nTpcR;i++){
								if(i<nRmin||i>nRmax)continue;
								double dphistar=getphistar(LA1,LA2,-1,-1,-0.5,TpcR[i]);
								if(fabs(dphistar)<dphistarcut){
									phistar_index=1;
									break;
								}
								if(rapindex==1&&index_E==52&&ybinpair>=8&&fabs(dphistar)<(dphistarcut+0.03)){
									phistar_index=1;
									break;
								}
							}
							if(phistar_index==1&&fabs(deta)<detacut2)continue;
							if(index_Hphistar==1){
								for(int i=0;i<nTpcR;i++){
									double dphistar=getphistar(LA1,LA2,-1,-1,-0.5,TpcR[i]);
									if(rapindex==0){
										if(qinv<0.1)hdphistardeta_minusminus_mix[2][KTBin][i]->Fill(dphistar,deta);
										if(ybin1==ybin2&&qinv<0.1)hdphistardeta_minusminus_mix[ybin1][KTBin][i]->Fill(dphistar,deta);
									}
									if(rapindex==1){
										if(qinv<0.1)hdphistardeta_minusminus_mix[ybinpair][KTBin][i]->Fill(dphistar,deta);
									}

								}
							}

							if(index_1DCF==1){
								if(find_pair_cut_index==1){
									if(sl>-0.5&&sl<0.9)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,1.5);
									if(sl>-0.5&&sl<0.8)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,2.5);
									if(sl>-0.5&&sl<0.7)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,3.5);
									if(sl>-0.5&&sl<0.6)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,4.5);
									if(sl>-0.5&&sl<0.5)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,5.5);
									if(sl>-0.5&&sl<0.4)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,6.5);
									if(sl>-0.5&&sl<0.3)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,7.5);
									if(sl>-0.5&&sl<0.2)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,8.5);
									if(sl>-0.5&&sl<0.1)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,9.5);
									if(sl>-0.5&&sl<0.0)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,10.5);
								}
							}
							TVector3 Qosl=getQosl_LCMS(LA1,LA2);
							qout=Qosl.X();
							qside=Qosl.Y();
							qlong=Qosl.Z();
							//for plusplus qout*qside>0
							//for minusminus qout*qside<0
							//if(qout*qside<0)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,ptbin1);
							if(qoutqside_index==1&&(qout*qside)<0)continue;
							if(qoutqside_index2==1){
								if(qout<0&&qside>0&&qlong>0)continue;
								if(qout>0&&qside<0&&qlong<0)continue;
							}
							double weight=hcoul->GetBinContent(hcoul->FindBin(qinv));
							if(qinv>0.6)weight=1;
							double Q_LCMS=getQ_LCMS(LA1,LA2);
							if(rapindex==1){
								if(index_1DCF==1)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,ybinpair+0.5);//last bin -->> allrap range
								if(index_1DCF==1)hLevyQ_minus_minus_mix[centbin]->Fill(Q_LCMS,ktbin,ybinpair+0.5);
								if(index_3DCF==1)hq_mix_minusminus[centbin][ybinpair][KTBin]->Fill(qout,qside,qlong);
								if(index_3DCF==1)hcoul_minusminus_mix[centbin][ybinpair][KTBin]->Fill(qout,qside,qlong,weight);
							}
							if(rapindex==0){
								if(index_1DCF==1)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,14.5);//last bin -->> allrap range
								if(index_1DCF==1)hLevyQ_minus_minus_mix[centbin]->Fill(Q_LCMS,ktbin,14.5);
								if(index_3DCF==1)hq_mix_minusminus[centbin][2][KTBin]->Fill(qout,qside,qlong);
								if(index_3DCF==1)hcoul_minusminus_mix[centbin][2][KTBin]->Fill(qout,qside,qlong,weight);
								if(ybin1==ybin2){
									if(index_1DCF==1)hQ_minus_minus_mix[centbin]->Fill(qinv,ktbin,ybin1+0.5);
									if(index_1DCF==1)hLevyQ_minus_minus_mix[centbin]->Fill(Q_LCMS,ktbin,ybin1+0.5);
									if(index_3DCF==1)hq_mix_minusminus[centbin][ybin1][KTBin]->Fill(qout,qside,qlong);
									if(index_3DCF==1)hcoul_minusminus_mix[centbin][ybin1][KTBin]->Fill(qout,qside,qlong,weight);
									for(int i=0;i<kDeltamombin;i++){
										TLorentzVector newLA1_cms,newLA2_cms;
										newLA1_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA1,index_E/20.0),(i+1)*deltamom);
										newLA2_cms=getnewfourmom(boostfourmom_to_AuAuCMS(LA2,index_E/20.0),(i+1)*deltamom);
										if(index_E>=77){
											newLA1_cms=getnewfourmom(LA1,(i+1)*deltamom);
											newLA2_cms=getnewfourmom(LA2,(i+1)*deltamom);
										}
										double newqinvcms=sqrt(fabs((newLA1_cms-newLA2_cms).Mag2()));
										if(index_1Dmomshift==1)hdelta_mom_cms_Q_minus_minus_mix[centbin][i]->Fill(newqinvcms,ktbin,ybin1);
										TVector3 newQoslcms=getQosl_LCMS(newLA1_cms,newLA2_cms);
										double newqoutcms=newQoslcms.X();
										double newqsidecms=newQoslcms.Y();
										double newqlongcms=newQoslcms.Z();
										double newweightcms=hcoul->GetBinContent(hcoul->FindBin(newqinvcms));
										if(newqinvcms>0.6)newweightcms=1;
										if(index_3Dmomshift==1)hdelta_mom_cms_q_mix_minusminus[centbin][ybin1][KTBin][i]->Fill(newqoutcms,newqsidecms,newqlongcms);
										if(index_3Dmomshift==1)hdelta_mom_cms_coul_minusminus_mix[centbin][ybin1][KTBin][i]->Fill(newqoutcms,newqsidecms,newqlongcms,newweightcms);
									}
								}
							}
							if(smear_index==1){
								double smearweight;//for smear
								int finalcentbin=-1;
								if(Centrality1==8||Centrality1==7)finalcentbin=0;
								if(Centrality1==6||Centrality1==5)finalcentbin=1;
								if(Centrality1==4||Centrality1==3)finalcentbin=2;
								if(Centrality1==2||Centrality1==1||Centrality1==0)finalcentbin=3;
								if(finalcentbin<0)continue;
								smearweight=1-lamminus[finalcentbin][ybin1][KTBin]+lamminus[finalcentbin][ybin1][KTBin]*weight*(1+exp((-qout*qout*outminus[finalcentbin][ybin1][KTBin]*outminus[finalcentbin][ybin1][KTBin]-qside*qside*sideminus[finalcentbin][ybin1][KTBin]*sideminus[finalcentbin][ybin1][KTBin]-qlong*qlong*longminus[finalcentbin][ybin1][KTBin]*longminus[finalcentbin][ybin1][KTBin]-2*qout*qlong*ol2minus[finalcentbin][ybin1][KTBin])/0.038937929230));
								hq_A_ideal_minusminus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong,smearweight);
								hq_B_ideal_minusminus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong);
								LA1=Dosmearminus(LA1,index_E);
								LA2=Dosmearminus(LA2,index_E);
								TVector3 Qsmear=getQosl_LCMS(LA1,LA2);
								qout=Qsmear.X();
								qside=Qsmear.Y();
								qlong=Qsmear.Z();
								hq_A_smear_minusminus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong,smearweight);
								hq_B_smear_minusminus[finalcentbin][ybin1][KTBin]->Fill(qout,qside,qlong);
							}
						}
					}
				}
			}
		}//mixend
	}
	ohm.Write();	//save all booked histogram
	ohm.Close();
	delete t;
	return 0;
}

Int_t getZbin(Float_t Zvert , Int_t index_E )
{
	Float_t zbin[kZBin+1]={0};
	for(int i=0;i<kZBin+1;i++){
		zbin[i]=-145+i*5;
	}
	Int_t z;
	if(index_E>=77){
		if(Zvert < zbin[0] || Zvert >= zbin[kZBin])return -1;
		for(z =0;z<kZBin;z++){
			if(zbin[z]<=Zvert && zbin[z+1]>Zvert)return z;
		}
	}
	if(index_E<77){
		if(Zvert<198||Zvert>202)return -1;
		return 1;
	}
}

//void rmpileup(Double_t nrefmult , Double_t tofmatch){
//	double c3p0[5]={0};
//	double c3p2[5]={-13.59,1.515,0.02816,-1.195E-4,-9.639E-7};
//	double c3p5[5]={-13.59,1.515,0.02816,-1.195E-4,-9.639E-7};
//	double c3p9[5]={-13.59,1.515,0.02816,-1.195E-4,-9.639E-7};
//	double c4p5[5]={-13.59,1.515,0.02816,-1.195E-4,-9.639E-7};
//	double c5p2[5]={-13.59,1.515,0.02816,-1.195E-4,-9.639E-7};
//	double c7p7[5]={0};
//	double b3p0[5]={0};
//	double b3p2[5]={19.48,5.428,-0.007,-2.428E-4,1.197E-7};
//	double b3p5[5]={19.48,5.428,-0.007,-2.428E-4,1.197E-7};
//	double b3p9[5]={19.48,5.428,-0.007,-2.428E-4,1.197E-7};
//	double b4p5[5]={19.48,5.428,-0.007,-2.428E-4,1.197E-7};
//	double b5p2[5]={19.48,5.428,-0.007,-2.428E-4,1.197E-7};
//	double b7p7[5]={0};
//	double c[5]={0};
//	double b[5]={0};
//	for(int i=0;i<5;i++){
//		if(index_E==30){
//			b[i]=b3p0[i];
//			c[i]=c3p0[i];
//		}
//		if(index_E==32){
//			b[i]=b3p2[i];
//			c[i]=c3p2[i];
//		}
//		if(index_E==35){
//			b[i]=b3p5[i];
//			c[i]=c3p5[i];
//		}
//		if(index_E==39){
//			b[i]=b3p9[i];
//			c[i]=c3p9[i];
//		}
//		if(index_E==45){
//			b[i]=b4p5[i];
//			c[i]=c4p5[i];
//		}
//		if(index_E==52){
//			b[i]=b5p2[i];
//			c[i]=c5p2[i];
//		}
//		if(index_E==77){
//			b[i]=b7p7[i];
//			c[i]=c7p7[i];
//		}
//	}
//	if(nrefmult<(c[0]+c[1]*tofmatch+c[2]*pow(tofmatch,2)+c[3]*pow(tofmatch,3)+c[4]*pow(tofmatch,4)))continue;
//	if(nrefmult>(b[0]+b[1]*tofmatch+b[2]*pow(tofmatch,2)+b[3]*pow(tofmatch,3)+b[4]*pow(tofmatch,4)))continue;
//}
Int_t getCentBin(Int_t nrefmult , Int_t index_E){
	Float_t cent3p0[kCentBin+1]={5,9,16,26,41,60,86,119,142,195};
	Float_t cent3p2[kCentBin+1]={5,11,20,33,53,81,118,166,197,287};
	Float_t cent3p5[kCentBin+1]={6,12,21,37,59,89,128,181,216,999};
	Float_t cent3p9[kCentBin+1]={6,13,23,40,64,97,141,198,236,999};
	Float_t cent4p5[kCentBin+1]={7,14,27,45,72,108,154,216,257,999};
	Float_t cent5p2[kCentBin+1]={7,15,28,47,76,114,166,233,277,999};//Takahito Todoroki
	Float_t cent7p7[kCentBin+1]={13,17,22,27,35,46,60,80,95,260};
	Float_t centbd[kCentBin+1]={0};
	for(int i=0;i<kCentBin+1;i++){
		if(index_E==30)centbd[i]=cent3p0[i];
		if(index_E==32)centbd[i]=cent3p2[i];
		if(index_E==35)centbd[i]=cent3p5[i];
		if(index_E==39)centbd[i]=cent3p9[i];
		if(index_E==45)centbd[i]=cent4p5[i];
		if(index_E==52)centbd[i]=cent5p2[i];
		if(index_E==77)centbd[i]=cent7p7[i];
	}

	Int_t s;

	if(nrefmult<centbd[0] || nrefmult>centbd[kCentBin])return -1;
	if(nrefmult==centbd[kCentBin])return 8;

	for(s=0;s<kCentBin;s++){
		if(centbd[s]<=nrefmult && centbd[s+1]>nrefmult)
		{                                              
			return s;
		}   
	}   
}   
Int_t getsinglerapbin(Float_t rap ,Int_t index_E)
{
	Float_t ybinFXT[3]={-1.0,-0.5,0};
	Float_t ybinCOL[3]={-0.5, 0,0.5};
	Float_t ybin[3]={0};
	if(index_E>=77){
		for(int i=0;i<3;i++){
			ybin[i]=ybinCOL[i];
		}
	}
	if(index_E<77){
		for(int i=0;i<3;i++){
			ybin[i]=ybinFXT[i];
		}
	}
	if(rap < ybin[0] || rap > ybin[2])return -1;
	Int_t ibin;
	for(ibin =0;ibin<2;ibin++){
		if(ybin[ibin]<rap && ybin[ibin+1]>rap)return ibin;
	}
}

Int_t getpairrapbin(Float_t rap ,Int_t index_E)
{
	Float_t ybin3p0[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin3p2[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin3p5[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin3p9[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin4p5[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin5p2[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin7p7[10]={-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2}; 
	Float_t ybin[10]={0};
	for(int i=0;i<10;i++){
		if(index_E==30)ybin[i]=ybin3p0[i];
		if(index_E==32)ybin[i]=ybin3p2[i];
		if(index_E==35)ybin[i]=ybin3p5[i];
		if(index_E==39)ybin[i]=ybin3p9[i];
		if(index_E==45)ybin[i]=ybin4p5[i];
		if(index_E==52)ybin[i]=ybin5p2[i];
		if(index_E==77)ybin[i]=ybin7p7[i];
	}
	if(rap < ybin[0] || rap > ybin[9])return -1;
	Int_t ibin;
	for(ibin =0;ibin<9;ibin++){
		if(ybin[ibin]<rap && ybin[ibin+1]>rap)return ibin;
	}
}
double getshiftplus(double p, int index_E) {
	double pbd[19]={0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,2.0};
	double shift3p0[18]={0.5};
	double shift3p2[18]={-0.1,0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.2};
	double shift3p5[18]={-0.2,-0.1,-0.0,-0.0,-0.0,0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
	double shift3p9[18]={-0.1,-0.0,-0.0,-0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
	double shift4p5[18]={-0.1,0.0,-0.0,-0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
	double shift5p2[18]={-0.1,0.0,0.0,0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
	double shift7p7[18]={-1.0,-0.7,-0.3,-0.1,0.0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,0.0};
	double shift[18]={0};
	for(int i=0;i<18;i++){
		if(index_E==30)shift[i]=shift3p0[i];
		if(index_E==32)shift[i]=shift3p2[i];
		if(index_E==35)shift[i]=shift3p5[i];
		if(index_E==39)shift[i]=shift3p9[i];
		if(index_E==45)shift[i]=shift4p5[i];
		if(index_E==52)shift[i]=shift5p2[i];
		if(index_E==77)shift[i]=shift7p7[i];
	}
	if(p<pbd[0])return shift[0];
	if(p>pbd[18])return shift[17];
	for(int i=0;i<18;i++){
		if(pbd[i]<=p && pbd[i+1]>p)
		{                                              
			return shift[i];
		}   
	}
}
double getshiftminus(double p, int index_E) {
	double pbd[19]={0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,2.0};
	double shift3p0[18]={0.5};
	double shift3p2[18]={-0.2,-0.0,-0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3};
	double shift3p5[18]={-0.3,-0.1,-0.1,-0.1,-0.0,-0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2};
	double shift3p9[18]={-0.2,-0.1,-0.1,-0.0,-0.0,0.0,0.0,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.3};
	double shift4p5[18]={-0.2,-0.0,-0.0,-0.0,0.0,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.4,0.4,0.4};
	double shift5p2[18]={-0.2,-0.0,-0.0,-0.0,0.0,0.0,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.4,0.4};
	double shift7p7[18]={-1.0,-0.7,-0.3,-0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
	double shift[18]={0};
	for(int i=0;i<18;i++){
		if(index_E==30)shift[i]=shift3p0[i];
		if(index_E==32)shift[i]=shift3p2[i];
		if(index_E==35)shift[i]=shift3p5[i];
		if(index_E==39)shift[i]=shift3p9[i];
		if(index_E==45)shift[i]=shift4p5[i];
		if(index_E==52)shift[i]=shift5p2[i];
		if(index_E==77)shift[i]=shift7p7[i];
	}
	if(p<pbd[0])return shift[0];
	if(p>pbd[18])return shift[17];
	for(int i=0;i<18;i++){
		if(pbd[i]<=p && pbd[i+1]>p)
		{                                              
			return shift[i];
		}   
	}
}
double getQ_LCMS(TLorentzVector Four_mom1, TLorentzVector Four_mom2){
	double x1 = Four_mom1.X();  double y1 = Four_mom1.Y();  double z1 = Four_mom1.Z();  double t1 = Four_mom1.T();
	double x2 = Four_mom2.X();  double y2 = Four_mom2.Y();  double z2 = Four_mom2.Z();  double t2 = Four_mom2.T();

	double dx = x1-x2;
	double dy = y1-y2;
	double dz2_LCMS = 4*(pow(z1*t2-z2*t1,2))/(pow(t1+t2,2)-pow(z1+z2,2));
	double Q_LCMS=sqrt(dx*dx+dy*dy+dz2_LCMS);
	return Q_LCMS;

}
TLorentzVector boostfourmom_to_AuAuCMS(TLorentzVector Four_mom,double energy){
	double mass=0.938272081;
	TLorentzVector p_CM(0,0,sqrt(energy*energy-mass*mass),energy);
	TVector3 beta = 1*p_CM.BoostVector();
	Four_mom.Boost(beta);
	TLorentzVector  newFour_mom(Four_mom.X(),Four_mom.Y(),Four_mom.Z(),Four_mom.T());
	return newFour_mom;

}
TLorentzVector getnewfourmom(TLorentzVector Four_mom,double delta_mom){
	TLorentzVector newFour_mom;
	double p0,p1,px1,py1,pz1;
	p0=Four_mom.P();
	p1=p0+delta_mom;
	double Theta=Four_mom.Theta();
	double Phi=Four_mom.Phi();
	pz1=p1*cos(Theta);
	px1=p1*sin(Theta)*cos(Phi);
	py1=p1*sin(Theta)*sin(Phi);
	TVector3 newThree_mom;
	newThree_mom.SetXYZ(px1,py1,pz1);
	double pionmass=0.13957039;
	newFour_mom.SetVectM(newThree_mom,pionmass);
	return newFour_mom;
}
TVector3 getQosl_LCMS(TLorentzVector Four_mom1, TLorentzVector Four_mom2){
	double x1 = Four_mom1.X();  double y1 = Four_mom1.Y();  double z1 = Four_mom1.Z();  double t1 = Four_mom1.T();
	double x2 = Four_mom2.X();  double y2 = Four_mom2.Y();  double z2 = Four_mom2.Z();  double t2 = Four_mom2.T();

	double dx = x1-x2;
	double xt = x1+x2;

	double dy = y1-y2;
	double yt = y1+y2;

	double dz = z1-z2;
	double zt = z1+z2;

	double dt = t1-t2;
	double tt = t1+t2;

	double k1 = sqrt(xt*xt+yt*yt);
	double k2 = (dx*xt+dy*yt);
	double beta = zt/tt;
	double gamma = 1.0/sqrt(1.0 - beta*beta);
	double Qout = k2/k1;
	double Qside = 2.0*(x2*y1-x1*y2)/k1;
	double Qlong = gamma*(dz - beta*dt);
	TVector3 q3;
	q3.SetXYZ(Qout,Qside,Qlong);
	return q3;
}
TLorentzVector Dosmearplus(TLorentzVector Four_mom,int index_E){

	double a3p0[4]={0.00128838,0.00527254,-0.00150014,0.0043203};//par of dpt/pt
	double b3p0[4]={0.00205018,0.0010024,-0.000396667,-8.51193e-05};//par of deta
	double c3p0[4]={0.00198175,0.00099719,-0.000338764,-0.000732126};//par of dphi

	double a3p2[4]={0.00165891,0.00649343,-0.00199302,0.00336917};//par of dpt/pt
	double b3p2[4]={0.00230002,0.00189034,-0.000703714,-0.000669533};//par of deta
	double c3p2[4]={0.00221697,0.00142921,-0.000460777,-0.00128699};//par of dphi

	double a3p5[4]={0.00165891,0.00649343,-0.00199302,0.00336917};//par of dpt/pt
	double b3p5[4]={0.00230002,0.00189034,-0.000703714,-0.000669533};//par of deta
	double c3p5[4]={0.00221697,0.00142921,-0.000460777,-0.00128699};//par of dphi

	double a3p9[4]={0.00165891,0.00649343,-0.00199302,0.00336917};//par of dpt/pt
	double b3p9[4]={0.00230002,0.00189034,-0.000703714,-0.000669533};//par of deta
	double c3p9[4]={0.00221697,0.00142921,-0.000460777,-0.00128699};//par of dphi

	double a4p5[4]={0.00165891,0.00649343,-0.00199302,0.00336917};//par of dpt/pt
	double b4p5[4]={0.00230002,0.00189034,-0.000703714,-0.000669533};//par of deta
	double c4p5[4]={0.00221697,0.00142921,-0.000460777,-0.00128699};//par of dphi

	double a5p2[4]={0.00165891,0.00649343,-0.00199302,0.00336917};//par of dpt/pt
	double b5p2[4]={0.00230002,0.00189034,-0.000703714,-0.000669533};//par of deta
	double c5p2[4]={0.00221697,0.00142921,-0.000460777,-0.00128699};//par of dphi

	double a7p7[4]={0.00165891,0.00649343,-0.00199302,0.00336917};//par of dpt/pt
	double b7p7[4]={0.00230002,0.00189034,-0.000703714,-0.000669533};//par of deta
	double c7p7[4]={0.00221697,0.00142921,-0.000460777,-0.00128699};//par of dphi

	//double a7p7[4]={0.00177075,0.00787503,-0.0024217,0.00112562};//par of dpt/pt
	//double b7p7[4]={0.0017918,0.00499288,-0.002184,-0.00289791};//par of deta
	//double c7p7[4]={0.00147329,0.00144741,-0.000472098,-0.00122561};//par of dphi

	double a[4]={0};//par of dpt/pt
	double b[4]={0};//par of deta
	double c[4]={0};//par of dphi
	for(int i=0;i<4;i++){
		if(index_E==30){
			a[i]=a3p0[i];
			b[i]=b3p0[i];
			c[i]=c3p0[i];
		}
		if(index_E==32){
			a[i]=a3p2[i];
			b[i]=b3p2[i];
			c[i]=c3p2[i];
		}
		if(index_E==35){
			a[i]=a3p5[i];
			b[i]=b3p5[i];
			c[i]=c3p5[i];
		}
		if(index_E==39){
			a[i]=a3p9[i];
			b[i]=b3p9[i];
			c[i]=c3p9[i];
		}
		if(index_E==45){
			a[i]=a4p5[i];
			b[i]=b4p5[i];
			c[i]=c4p5[i];
		}
		if(index_E==52){
			a[i]=a5p2[i];
			b[i]=b5p2[i];
			c[i]=c5p2[i];
		}
		if(index_E==77){
			a[i]=a7p7[i];
			b[i]=b7p7[i];
			c[i]=c7p7[i];
		}
	}
	double px=Four_mom.Px();
	double py=Four_mom.Py();
	double pz=Four_mom.Pz();
	double pt=Four_mom.Perp();
	double eta=Four_mom.Eta();
	double phi=Four_mom.Phi();
	double E=Four_mom.E();
	double M=sqrt(E*E-px*px-py*py-pz*pz);
	double nsigma=1;
	double sigma_dptoverpt=nsigma*(a[0]*1.0/pt+a[1]*pt+a[2]*pt*pt+a[3]);
	double sigma_deta=nsigma*(b[0]*1.0/pt+b[1]*pt+b[2]*pt*pt+b[3]);
	double sigma_dphi=nsigma*(c[0]*1.0/pt+c[1]*pt+c[2]*pt*pt+c[3]);
	double dpt=gRandom->Gaus(0,sigma_dptoverpt)*pt;
	double deta=gRandom->Gaus(0,sigma_deta);
	double dphi=gRandom->Gaus(0,sigma_dphi);
	double pt_s=pt+dpt;
	double eta_s=eta+deta;
	double phi_s=phi+dphi;
	TLorentzVector Four_mom_s;
	Four_mom_s.SetPtEtaPhiM(pt_s,eta_s,phi_s,M);
	return Four_mom_s;
}
TLorentzVector Dosmearminus(TLorentzVector Four_mom, int index_E){

	double a3p0[4]={0.00128381,0.00497691,-0.0013227,0.00438203};//par of dpt/pt
	double b3p0[4]={0.0020486,0.00109525,-0.000439418,-0.000139845};//par of deta
	double c3p0[4]={0.00194986,0.000812095,-0.000271659,-0.000583314};//par of dphi

	double a3p2[4]={0.00172727,0.00676877,-0.00203881,0.00304546};//par of dpt/pt
	double b3p2[4]={0.00234304,0.00215651,-0.000816452,-0.000897605};//par of deta
	double c3p2[4]={0.00223963,0.00157047,-0.000514462,-0.00139405};//par of dphi

	double a3p5[4]={0.00172727,0.00676877,-0.00203881,0.00304546};//par of dpt/pt
	double b3p5[4]={0.00234304,0.00215651,-0.000816452,-0.000897605};//par of deta
	double c3p5[4]={0.00223963,0.00157047,-0.000514462,-0.00139405};//par of dphi

	double a3p9[4]={0.00172727,0.00676877,-0.00203881,0.00304546};//par of dpt/pt
	double b3p9[4]={0.00234304,0.00215651,-0.000816452,-0.000897605};//par of deta
	double c3p9[4]={0.00223963,0.00157047,-0.000514462,-0.00139405};//par of dphi

	double a4p5[4]={0.00172727,0.00676877,-0.00203881,0.00304546};//par of dpt/pt
	double b4p5[4]={0.00234304,0.00215651,-0.000816452,-0.000897605};//par of deta
	double c4p5[4]={0.00223963,0.00157047,-0.000514462,-0.00139405};//par of dphi

	double a5p2[4]={0.00172727,0.00676877,-0.00203881,0.00304546};//par of dpt/pt
	double b5p2[4]={0.00234304,0.00215651,-0.000816452,-0.000897605};//par of deta
	double c5p2[4]={0.00223963,0.00157047,-0.000514462,-0.00139405};//par of dphi

	double a7p7[4]={0.00172727,0.00676877,-0.00203881,0.00304546};//par of dpt/pt
	double b7p7[4]={0.00234304,0.00215651,-0.000816452,-0.000897605};//par of deta
	double c7p7[4]={0.00223963,0.00157047,-0.000514462,-0.00139405};//par of dphi

	//double a7p7[4]={0.00175894,0.00781781,-0.00240634,0.0011798};//par of dpt/pt
	//double b7p7[4]={0.00148796,0.00187602,-0.000686221,-0.000971511};//par of deta
	//double c7p7[4]={0.00148742,0.00148931,-0.00048614,-0.00126682};//par of dphi

	double a[4]={0};//par of dpt/pt
	double b[4]={0};//par of deta
	double c[4]={0};//par of dphi
	for(int i=0;i<4;i++){
		if(index_E==30){
			a[i]=a3p0[i];
			b[i]=b3p0[i];
			c[i]=c3p0[i];
		}
		if(index_E==32){
			a[i]=a3p2[i];
			b[i]=b3p2[i];
			c[i]=c3p2[i];
		}
		if(index_E==35){
			a[i]=a3p5[i];
			b[i]=b3p5[i];
			c[i]=c3p5[i];
		}
		if(index_E==39){
			a[i]=a3p9[i];
			b[i]=b3p9[i];
			c[i]=c3p9[i];
		}
		if(index_E==45){
			a[i]=a4p5[i];
			b[i]=b4p5[i];
			c[i]=c4p5[i];
		}
		if(index_E==52){
			a[i]=a5p2[i];
			b[i]=b5p2[i];
			c[i]=c5p2[i];
		}
		if(index_E==77){
			a[i]=a7p7[i];
			b[i]=b7p7[i];
			c[i]=c7p7[i];
		}
	}
	double px=Four_mom.Px();
	double py=Four_mom.Py();
	double pz=Four_mom.Pz();
	double pt=Four_mom.Perp();
	double eta=Four_mom.Eta();
	double phi=Four_mom.Phi();
	double E=Four_mom.E();
	double M=sqrt(E*E-px*px-py*py-pz*pz);
	double nsigma=1;
	double sigma_dptoverpt=nsigma*(a[0]*1.0/pt+a[1]*pt+a[2]*pt*pt+a[3]);
	double sigma_deta=nsigma*(b[0]*1.0/pt+b[1]*pt+b[2]*pt*pt+b[3]);
	double sigma_dphi=nsigma*(c[0]*1.0/pt+c[1]*pt+c[2]*pt*pt+c[3]);
	double dpt=gRandom->Gaus(0,sigma_dptoverpt)*pt;
	double deta=gRandom->Gaus(0,sigma_deta);
	double dphi=gRandom->Gaus(0,sigma_dphi);
	double pt_s=pt+dpt;
	double eta_s=eta+deta;
	double phi_s=phi+dphi;
	TLorentzVector Four_mom_s;
	Four_mom_s.SetPtEtaPhiM(pt_s,eta_s,phi_s,M);
	return Four_mom_s;
}
double getphistar(TLorentzVector Four_mom1, TLorentzVector Four_mom2, int q1, int q2,double Bz, double tpcR){
	double deltaphistar = Four_mom1.Phi()-Four_mom2.Phi() + TMath::ASin(-0.15*(q1)*Bz*tpcR/Four_mom1.Perp())-TMath::ASin(-0.15*(q2)*Bz*tpcR/Four_mom2.Perp());
	deltaphistar = atan2(sin(deltaphistar),cos(deltaphistar));
	return deltaphistar;
}
double getSL(Int_t padRow1To24Track1 ,Int_t padRow25To45Track1 ,ULong64_t IpadRow1 ,Int_t nhits1 , Int_t padRow1To24Track2 , Int_t padRow25To45Track2 ,ULong64_t IpadRow2 ,Int_t nhits2, Int_t index_E) {
	// AND logic
	unsigned long bothPads1To24 = padRow1To24Track1 & padRow1To24Track2;
	unsigned long bothPads25To45 = padRow25To45Track1 & padRow25To45Track2;
	ULong64_t     bothIPads = IpadRow1 & IpadRow2;
	// XOR logic
	unsigned long onePad1To24 = padRow1To24Track1 ^ padRow1To24Track2;
	unsigned long onePad25To45 = padRow25To45Track1 ^ padRow25To45Track2;
	ULong64_t     oneIPads = IpadRow1 ^ IpadRow2;
	unsigned long bitI;
	int ibits;
	int Quality = 0;
	double normQual = 0.0;
	int MaxQuality = nhits1+nhits2;
	for (ibits=8;ibits<=31;ibits++) {
		bitI = 0;
		bitI |= 1UL<<(ibits);
		if ( onePad1To24 & bitI ) {
			Quality++;
			continue;
		}
		else{
			if ( bothPads1To24 & bitI ) Quality--;
		}
	}
	for (ibits=0;ibits<=20;ibits++) {
		bitI = 0;
		bitI |= 1UL<<(ibits);
		if ( onePad25To45 & bitI ) {
			Quality++;
			continue;
		}
		else{
			if ( bothPads25To45 & bitI ) Quality--;
		}
	}
	if(index_E==30){
		normQual = (double)Quality/( (double) MaxQuality );
		return normQual;
	}
	for (ibits=0;ibits<=40;ibits++) {
		bitI = 0;
		bitI |= 1UL<<(ibits);
		if ( oneIPads & bitI ) {
			Quality++;
			continue;
		}
		else{
			if ( bothIPads & bitI ) Quality--;
		}
	}
	normQual = (double)Quality/( (double) MaxQuality );
	return normQual;
}        

TChain* ChainThem(const char* filelist,const char* treename, int nlist, int block){
	TChain *globChain = new TChain(treename);

	cout << ">>> Load Chain from file: " << filelist << endl;
	ifstream fList(filelist);

	//  globChain->Add(filelist,0);
	int Ncount = 0;
	int Nfiles = 0;
	char lineFromFile[255];

	if (!fList)
	{
		cout << "!!! Can't open file " << filelist << endl;
		return nullptr;
	}
	while(fList.getline(lineFromFile, 250))
	{
		Ncount++;
		// if(Ncount<=nlist*block)continue;
		//  if(Ncount>(nlist+1)*block)break;
		TFile tempf(lineFromFile);
		if(!tempf.IsZombie() && tempf.GetNkeys())
		{
			globChain->Add(lineFromFile,0);
			cout <<Nfiles<< ">> File '" << lineFromFile << "' has been loaded" << endl;
			Nfiles ++;
		}
		//     else
		//cout << ">> Can't load file '" << lineFromFile << "'" << endl;
	}
	if(Nfiles == 0){
		delete globChain;
		return NULL;
	}

	cout << ">> Total number of entries: " << globChain->GetEntriesFast() << endl;
	fList.close();

	return globChain;
}




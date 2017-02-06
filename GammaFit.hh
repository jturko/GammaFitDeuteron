

#ifndef GAMMAFIT_H
#define GAMMAFIT_H
#endif

#include "TCanvas.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TFile.h"
#include "TColor.h"
#include "TApplication.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom3.h"

#include "TMath.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

class GammaFit 
{

public:
    GammaFit(int cal=3, std::string source="24Na");
    ~GammaFit();

    void PrintParameters() {
        std::cout << "from GammaFit class... " << std::endl;
        std::cout << "      A = " << fParameters[0] << std::endl;
        std::cout << "      B = " << fParameters[1] << std::endl;
        std::cout << "      C = " << fParameters[2] << std::endl;
        std::cout << "   gain = " << fParameters[3] << " keVee/channel" << std::endl;
        std::cout << " offset = " << fParameters[4] << " keVee " << std::endl;
    }

    void Sort(double * par);
    void Sort(double A=0.1, double B=0.05, double C=1e-4, double gain=0.5, double offset=0.);

    void SetParameters(double * par);

    double LightOutput(double E, double * par) {
        return ( par[0]*E-par[1]*(1.0-TMath::Exp(-par[2]*TMath::Power(E,par[3]))) );
        //return ( fLightOffset/1000. + par[0]*E-par[1]*(1.0-TMath::Exp(-par[2]*TMath::Power(E,par[3]))) );
    }
    double Resolution(double E, double * par) {
        return (E*TMath::Sqrt(TMath::Power(par[0],2)+TMath::Power(par[1],2)/E+TMath::Power(par[2]/E,2)))/(2.*TMath::Sqrt(2.*TMath::Log(2.))) ;
    }

    void ApplyCutoffLow(double cutoff, std::string str) {
        if(str == "sim")      for(int i=0; i<fSimHist->FindBin(cutoff); i++) fSimHist->SetBinContent(i,0.);
        else if(str == "exp") for(int i=0; i<fExpHist->FindBin(cutoff); i++) fExpHist->SetBinContent(i,0.);
        else                  std::cout << "--->  Error with ApplyCutoffLow!" << std::endl;
    }
    void ApplyCutoffHigh(double cutoff, std::string str) {
        if(str == "sim")      for(int i=fSimHist->FindBin(cutoff); i<fSimHist->GetNbinsX(); i++) fSimHist->SetBinContent(i,0.);
        else if(str == "exp") for(int i=fExpHist->FindBin(cutoff); i<fExpHist->GetNbinsX(); i++) fExpHist->SetBinContent(i,0.);
        else                  std::cout << "--->  Error with ApplyCutoffHigh!" << std::endl;
    }
    
    double DoChi2() {
        TH1F * h_e = (TH1F*)fExpHist->Clone();
        TH1F * h_s = (TH1F*)fSimHist->Clone();
        for(int i=h_e->FindBin(fCutoffHigh); i<h_e->GetNbinsX(); i++) h_e->SetBinContent(i,0.);
        for(int i=h_s->FindBin(fCutoffHigh); i<h_s->GetNbinsX(); i++) h_s->SetBinContent(i,0.);
        for(int i=0; i<h_e->FindBin(fCutoffLow); i++) h_e->SetBinContent(i,0.);
        for(int i=0; i<h_s->FindBin(fCutoffLow); i++) h_s->SetBinContent(i,0.);
        fChi2 = h_e->Chi2Test(h_s,"CHI2/NDF");
        delete h_e;
        delete h_s;

        if(fChi2 < 0.0001) fChi2 = 1e5;
        return fChi2;
    }

    void Draw() {
        if(fExpHist) fExpHist->SetLineColor(kBlack);
        if(fSimHist) fSimHist->SetLineColor(kRed);
        else { std::cout << "no hists sorted yet!" << std::endl; return; }
        fExpHist->GetXaxis()->SetRangeUser(fCutoffLow-30,fCutoffHigh+150);
        double ymax = 2*fExpHist->GetBinContent(fExpHist->GetMaximumBin());
        fExpHist->GetYaxis()->SetRangeUser(0.1,ymax);

        fExpHist->Draw();
        fSimHist->Draw("same");   
    }

    void SetCutoffHigh(double cut) { fCutoffHigh = cut; }
    void SetCutoffLow(double cut) { fCutoffLow = cut; }
    double GetCutoffHigh() { return fCutoffHigh; }    
    double GetCutoffLow() { return fCutoffLow; }    

    TH1F * GetSimHist() { return fSimHist; }
    TH1F * GetExpHist() { return fExpHist; }
    
    double GetSimEntries() { return fSimEntries; }
    double GetExpEntries() { return fExpEntries; }

    std::string GetSource() { return fSource; }

    double fCutoffHigh;
    double fCutoffLow;
    
    TFile * fSimFile;
    TFile * fExpFile;

    TH1F * fSimHist;
    TH1F * fExpHist;

    TTree * fSimTree;
    TTree * fExpTree;

    TBranch * fEdepBranch;
    TBranch * fEkinBranch;
    TBranch * fPtypeBranch;
    std::vector<double> * fEdepVector;
    std::vector<double> * fEkinVector;
    std::vector<int> * fPtypeVector;

    TBranch * fExpValueBranch;
    UInt_t fExpValue[6];
    std::vector<UInt_t> fExpValueVector;

    double fProtonCoeff[4];
    double fElectronCoeff[4];
    double fDeuteronCoeff[4];
    double fCarbonCoeff[4];
    double fAlphaCoeff[4];
    double fBeCoeff[4];
    double fBCoeff[4];
    double fSmearingCoeff[3];
    double fParameters[5];
    double fGain;
    double fOffset;

    int fCal;
    std::string fSource;

    int fSimEntries;
    int fExpEntries;
    TRandom3 fRandom;

    double fChi2;
    
    int fBinNum;
    double fBinHigh;
    double fBinLow;

    int fSimMaxEntry;


};

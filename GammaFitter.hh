
#ifndef GAMMAFITTER_H
#define GAMMAFITTER_H

#include "TCanvas.h"
#include "TH1F.h"
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

#include "GammaFit.cc"
#include "vec.hh"

class GammaFitter
{

public:
    GammaFitter();
    ~GammaFitter();

    GammaFitter(int cal, std::string s1);
    GammaFitter(int cal, std::string s1, std::string s2);
    GammaFitter(int cal, std::string s1, std::string s2, std::string s3);
    GammaFitter(int cal, std::string s1, std::string s2, std::string s3, std::string s4);

    void Draw();
    void Run(double resolution, double gain, double offset);
    void Run() {
        Run(fParameters[0],fParameters[1],fParameters[2]);
    }
    
    bool Check(int i) { if(i<=-1||i>=GetNumberOfGammaFits()) return false; else return true; }
    
    GammaFit * GetGammaFit(int i) { if(Check(i)) return &fGammaFitVector[i]; else{ std::cout << "out of bounds!" << std::endl; return NULL;} }
    
    int GetNumberOfGammaFits() { return fGammaFitVector.size(); }
 
    void SetNextGammaFit(GammaFit gfit) { fGammaFitVector.push_back(gfit); fSourceVector.push_back(gfit.GetSource()); } 
    void SetNextGammaFit(int cal, std::string source) {
        GammaFit * gfit = new GammaFit(cal,source);
        SetNextGammaFit(*gfit);
    }
   
    void DrawGammaFit(int i) { 
        if(!Check(i)) { std::cout << "out of bounds!" << std::endl; return; }
        if(fCanvas) delete fCanvas; 
        fGammaFitVector.at(i).Draw(); 
    }
    
    void SortRun(int num) { 
        if(Check(num)) { 
            fGammaFitVector.at(num).Sort(fParameters); 
        }
        else{ std::cout << "out of bounds!" << std::endl;} 
    }
    
    void SortAllRuns() { 
        //PrintParameters();
        for(int num=0; num<GetNumberOfGammaFits(); num++) SortRun(num); 
        DoChi2();
    }
    
    void Print() { for(int num=0; num<GetNumberOfGammaFits(); num++) std::cout << "Source = " << fSourceVector.at(num) << std::endl; }

    void SetParameters(double resolution, double gain, double offset) {
        fParameters[0] = resolution;
        fParameters[1] = gain;
        fParameters[2] = offset;
        for(int i=0; i<GetNumberOfGammaFits(); i++) fGammaFitVector.at(i).SetParameters(fParameters);
    }
    void SetParameters(double * par) { // expects a par array w/ 3 elements
        for(int i=0; i<3; i++) fParameters[i] = par[i];
        for(int i=0; i<GetNumberOfGammaFits(); i++) fGammaFitVector.at(i).SetParameters(fParameters);
    }
    
    void PrintParameters() { 
        std::cout << "from GammaFitter class... " << std::endl;    
        std::cout << " resolution = " << 100.*fParameters[0] << " %" << std::endl;
        std::cout << "       gain = " << fParameters[1] << " keVee/channel" << std::endl;
        std::cout << "     offset = " << fParameters[2] << " keVee" << std::endl;
    }

    void InitializeParameters();

    //void NelderMead(double A=0.1, double B=0.05, double C=1e-4, double gain=0.6, double offset=0, int itermax = 50);
    void NelderMead();
    void NelderMead(double resolution, double gain, double offset, int itermax);
    
    double DoChi2() { 
        fSum = 0.;
        fSum2 = 0.;
        for(int i=0;i<GetNumberOfGammaFits();i++) fSum += fGammaFitVector.at(i).DoChi2();
        for(int i=0;i<GetNumberOfGammaFits();i++) fSum2 += fGammaFitVector.at(i).DoChi2() * fGammaFitVector.at(i).DoChi2();
        fSum /= double(GetNumberOfGammaFits());
        fSum2 /= double(GetNumberOfGammaFits());    
    
        return fSum; // CHANGE THIS TO SET WHAT WE WILL MINIMIZE [ chi2 or (chi2)^2 ]
    }
        
    double nm_val(double * par) {
        SetParameters(par);
        SortAllRuns();
        return DoChi2();
    }
    double nm_val(vec v) {
        SetParameters(v.par_array());
        SortAllRuns();
        return DoChi2();
    }
    
    void DrawToFile(std::string name);

    std::vector<GammaFit> fGammaFitVector;   
    std::vector<std::string> fSourceVector;

    double fParameters[3];   
 
    TCanvas * fCanvas;
    
    double fSum;
    double fSum2;

    int fIterMax;

};

#endif

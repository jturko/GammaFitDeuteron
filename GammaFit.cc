



#include "GammaFit.hh"

GammaFit::GammaFit(int cal, std::string source) :
    fCal(cal),
    fSource(source),
    fSimFile(NULL),
    fExpFile(NULL),
    fSimHist(NULL),
    fExpHist(NULL),
    fSimTree(NULL),
    fEdepBranch(NULL),
    fEkinBranch(NULL),
    fPtypeBranch(NULL),
    fEdepVector(NULL),
    fEkinVector(NULL),
    fPtypeVector(NULL),
    fEvtTimeBranch(NULL),
    fEvtTime(0)
{
    fBinNum  = 1000;
    fBinHigh = 5000;
    fBinLow  = 0; 
    fSimMaxEntry = 250000;

    //std::string path = "~/data/testcan_source_raw/cal" + std::to_string(fCal) + "/" + fSource + ".root";
    std::string path = "../../../../testcan_source_raw/cal" + std::to_string(fCal) + "/" + fSource + ".root";
    fExpFile = TFile::Open(path.c_str());
    fExpTree = (TTree*)(fExpFile->Get("tree"));
    fExpTree->SetBranchAddress("value",&fExpValue,&fExpValueBranch);
    fExpEntries = fExpTree->GetEntries();    
    for(int i=0; i<GetExpEntries(); i++) 
    {
        fExpTree->GetEntry(i);
        UInt_t val = fExpValue[0];
        fExpValueVector.push_back(val);
        val = fExpValue[2];
        fExpPSDVector.push_back(val);
    }

    std::string name = "../G4_RAW/" + fSource + ".root";
    fSimFile = TFile::Open(name.c_str());     
    fSimTree = (TTree*)(fSimFile->Get("ntuple/ntuple")); 
    fSimTree->SetBranchAddress("eDepVector",&fEdepVector);
    fSimTree->SetBranchAddress("eKinVector",&fEkinVector);
    fSimTree->SetBranchAddress("particleTypeVector",&fPtypeVector);
    //fSimTree->SetBranchAddress("eDepVector",&fEdepVector,&fEdepBranch);
    //fSimTree->SetBranchAddress("eKinVector",&fEkinVector,&fEkinBranch);
    //fSimTree->SetBranchAddress("particleTypeVector",&fPtypeVector,&fPtypeBranch);
    //fSimTree->SetBranchAddress("time",&fEvtTime,&fEvtTimeBranch);
    fSimTree->SetBranchAddress("posz",&fEvtTime);
    fSimEntries = fSimTree->GetEntries();
    
    fProtonCoeff[0] = 0.74; fProtonCoeff[1] = 3.2; fProtonCoeff[2] = 0.20; fProtonCoeff[3] = 0.97;
    fDeuteronCoeff[0] = 0.75; fDeuteronCoeff[1] = 2.80; fDeuteronCoeff[2] = 0.25; fDeuteronCoeff[3] = 0.93;
    fCarbonCoeff[0] = 0.05; fCarbonCoeff[1] = 0.0; fCarbonCoeff[2] = 0.0;fCarbonCoeff[3] = 0.0;
    fAlphaCoeff[0] = 0.14; fAlphaCoeff[1] = 0.59; fAlphaCoeff[2] = 0.065; fAlphaCoeff[3] = 1.01;
    fBeCoeff[0] = 0.0821; fBeCoeff[1] = 0.0; fBeCoeff[2] = 0.0; fBeCoeff[3] = 0.0;
    fBCoeff[0] = 0.0375; fBCoeff[1] = 0.0; fBCoeff[2] = 0.0; fBCoeff[3] = 0.0;
    fElectronCoeff[0] = 1; fElectronCoeff[1] = 0; fElectronCoeff[2] = 0; fElectronCoeff[3] = 0;    
    
    //
    fSmearingCoeff[0] = 0.1; fSmearingCoeff[1] = 0.05; fSmearingCoeff[2] = 1e-4;
    fGain = 0.5;
    fOffset = 0;
    
    fParameters[0] = fSmearingCoeff[0];
    fParameters[1] = fSmearingCoeff[1]; 
    fParameters[2] = fSmearingCoeff[2];
    fParameters[3] = fGain;
    fParameters[4] = fOffset;
    
    SetParameters(fParameters);

    fCutoffHigh = 5000;
    if(fSource == "24Na") { fCutoffLow = 500; fCutoffHigh = 3000; }
    else if(fSource == "137Cs") { fCutoffLow = 300; fCutoffHigh = 700; }
    else if(fSource == "241Am"||fSource == "241Am_1"||fSource=="241Am_2"||fSource=="241Am_3") { fCutoffLow = 30; fCutoffHigh = 100; fBinNum = 10000; }
    else fCutoffLow = 0;
    
    if(fSource == "241Am" && fCal == 0) fCutoffLow = 50;
    
    std::cout << "Source: " << fSource << " ; cal: " << fCal << std::endl;
}

GammaFit::~GammaFit() {}

void GammaFit::SetParameters(double * par)
{
    fSmearingCoeff[0] = par[0];
    fSmearingCoeff[1] = par[1];
    fSmearingCoeff[2] = par[2];
    fGain             = par[3];
    fOffset           = par[4];
 
    fParameters[0] = par[0];
    fParameters[1] = par[1];
    fParameters[2] = par[2];
    fParameters[3] = par[3];
    fParameters[4] = par[4];
}

void GammaFit::Sort(double A, double B, double C, double gain, double offset)
{
    double par[5];
    par[0] = A;
    par[1] = B;
    par[2] = C;
    par[3] = gain;
    par[4] = offset;
    
    Sort(par);
}

void GammaFit::Sort(double * par)
{
    fRandom.SetSeed(1);
    gErrorIgnoreLevel = kError;    
    
    SetParameters(par);
    //PrintParameters();

    if(fSimHist) { delete fSimHist; fSimHist = NULL; }
    fSimHist = new TH1F("fSimHist","fSimHist",fBinNum,fBinLow,fBinHigh); 
    int nHits = 0;
    double light = 0.;
    double centroidEkin = 0.;    
    double centroidEres = 0.;    
    
    int counter = 0;

    for(int i=0; i<fSimEntries; i++)
    //for(int i=0; i<fSimMaxEntry; i++)
    {
        counter++;
        if( counter%50000==0 ) std::cout << "sorting evt " << counter << "/" << fSimEntries << "; " << double(counter)/double(fSimEntries)*100 << "% complete \r"  << std::flush; 
     
        //fEdepBranch->GetEntry(i);   
        //fEkinBranch->GetEntry(i);   
        //fPtypeBranch->GetEntry(i);   
        //fEvtTimeBranch->GetEntry(i);
        fSimTree->GetEntry(i);
        nHits = fEdepVector->size();
        light = 0.;
        for(int j=0; j<nHits; j++)
        {
            if(fPtypeVector->at(j) == 2 || fPtypeVector->at(j) == 3) {
                //centroidEkin = fEkinVector->at(j);
                //centroidEres = fEkinVector->at(j)-fEdepVector->at(j);
                centroidEkin = LightOutput(fEkinVector->at(j), fElectronCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fElectronCoeff);
            }
            else if(fPtypeVector->at(j) == 4) {
                centroidEkin = LightOutput(fEkinVector->at(j), fProtonCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fProtonCoeff);
            }
            else if(fPtypeVector->at(j) == 6) {
                centroidEkin = LightOutput(fEkinVector->at(j), fDeuteronCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fDeuteronCoeff);
            }
            else if(fPtypeVector->at(j) == 7) {
                centroidEkin = LightOutput(fEkinVector->at(j), fCarbonCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fCarbonCoeff);
            }
            else if(fPtypeVector->at(j) == 8) {
                centroidEkin = LightOutput(fEkinVector->at(j), fAlphaCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fAlphaCoeff);
            }
            else if(fPtypeVector->at(j) == 9) {
                centroidEkin = LightOutput(fEkinVector->at(j), fBeCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fBeCoeff);
            }
            else if(fPtypeVector->at(j) == 10) {
                centroidEkin = LightOutput(fEkinVector->at(j), fBCoeff );
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fBCoeff);
            }
            else { 
                centroidEkin = 0.; 
                centroidEres = 0.; 
            } 

            if(centroidEkin>0.){
                light += 1000.*fRandom.Gaus(centroidEkin, Resolution(centroidEkin,fSmearingCoeff));
            }
            if(centroidEres>0.){
                light -= 1000.*fRandom.Gaus(centroidEres, Resolution(centroidEres,fSmearingCoeff));
            } 
        }//end scatters loop       
        
        //std::cout << "fEvtTime = " << fEvtTime << " vector size = " << fEkinVector->size() << std::endl;    
        if(light>0. && fEvtTime < 1e10) {
            fSimHist->Fill(light);
        }
    }//end event loop


    // gain and offset for experimental histogram
    if(fExpHist) { delete fExpHist; fExpHist = NULL; }
    fExpHist = new TH1F("fExpHist","fExpHist",fBinNum,fBinLow,fBinHigh);
    
    double val = 0;
    for(int i=0; i<fExpEntries; i++) 
    {
        //fExpTree->GetEntry(i);
        //val = double(fExpValue[1])*fGain + fOffset;
        //if(fExpValue[1]>10 && fExpValue[1]<15000) fExpHist->Fill(val);
        
        //val = double(fExpValueVector.at(i))*fGain + fOffset;
        //if(fExpValueVector.at(i)>10 && fExpValueVector.at(i)<15000) fExpHist->Fill(val);
        val = fExpValueVector.at(i)*fGain + fOffset;
        val = fRandom.Gaus(val,0.2*fExpHist->GetBinWidth(0));
        
        if(fExpValueVector.at(i)>10 && fExpValueVector.at(i)<15000) fExpHist->Fill(val);
        //if(fExpValueVector.at(i)>10 && fExpValueVector.at(i)<15000 && fExpPSDVector.at(i)>6000) fExpHist->Fill(val);
    }
    fExpHist->Scale(10000/fExpHist->Integral());
    
    fExpHist->SetLineColor(kBlack);    
    fSimHist->SetLineColor(kRed);    

    //ApplyCutoffLow(fCutoffLow,"sim");    
    //ApplyCutoffLow(fCutoffLow,"exp");    
    
    fSimHist->Scale(fExpHist->Integral(fExpHist->FindBin(fCutoffLow),fExpHist->FindBin(fCutoffHigh),"width")/fSimHist->Integral(fSimHist->FindBin(fCutoffLow),fSimHist->FindBin(fCutoffHigh),"width"));
    
    fSimHist->SetStats(false);
    fExpHist->SetStats(false);
    
    std::cout << "sorting... done!                                                   " << std::endl;
    
    std::string title = fSource + " cal " + std::to_string(fCal) + " ; Chi2 = " + std::to_string(DoChi2());
    fExpHist->SetTitle(title.c_str());

}


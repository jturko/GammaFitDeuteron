
#include "GammaFitter.hh"

GammaFitter::GammaFitter() 
{ 
    fCanvas = NULL;
    SetNextGammaFit(3,"24Na");
    InitializeParameters(); 
} 

GammaFitter::~GammaFitter() {}

GammaFitter::GammaFitter(int cal, std::string source)
{
    fCanvas = NULL;
    SetNextGammaFit(cal,source);
    InitializeParameters();
}

GammaFitter::GammaFitter(int cal, std::string source1, std::string source2)
{
    fCanvas = NULL;
    SetNextGammaFit(cal,source1);
    SetNextGammaFit(cal,source2);
    InitializeParameters();
}

GammaFitter::GammaFitter(int cal, std::string source1, std::string source2, std::string source3)
{
    fCanvas = NULL;
    SetNextGammaFit(cal,source1);
    SetNextGammaFit(cal,source2);
    SetNextGammaFit(cal,source3);
    InitializeParameters();
}

GammaFitter::GammaFitter(int cal, std::string source1, std::string source2, std::string source3, std::string source4)
{
    fCanvas = NULL;
    SetNextGammaFit(cal,source1);
    SetNextGammaFit(cal,source2);
    SetNextGammaFit(cal,source3);
    SetNextGammaFit(cal,source4);
    InitializeParameters();
}

void GammaFitter::InitializeParameters()
{
    if(GetNumberOfGammaFits() > 0) {
        fParameters[0] = fGammaFitVector.at(0).fResolution;
        fParameters[1] = fGammaFitVector.at(0).fGain;
        fParameters[2] = fGammaFitVector.at(0).fOffset;
        std::cout << "parameters initialized to (" << fParameters[0] << "," << fParameters[1] << "," << fParameters[2] << ")" << std::endl;
    }    
    
    fSum = 0;
    fSum2 = 0;

    fIterMax = 50;
}

void GammaFitter::Draw()
{
    if(!fCanvas) {
        fCanvas = new TCanvas();
        if(GetNumberOfGammaFits() == 1) fCanvas->Divide(1);
        else if(GetNumberOfGammaFits() == 2) fCanvas->Divide(2);
        else if(GetNumberOfGammaFits() == 3) fCanvas->Divide(3);
        else if(GetNumberOfGammaFits() == 4) fCanvas->Divide(2,2);
        else { std::cout << "more/less GammaFits that I know how to draw!" << std::endl; return; }
    }
    for(int i=0; i<GetNumberOfGammaFits(); i++) {
        fCanvas->cd(i+1);
        gPad->Clear();
        fGammaFitVector.at(i).Draw();
        gPad->Update();
    }
}

//void GammaFitter::Run(double A, double B, double C, double gain, double offset)
void GammaFitter::Run(double resolution, double gain, double offset)
{
    SetParameters(resolution,gain,offset);
    PrintParameters();
    if(!fCanvas) {
        fCanvas = new TCanvas();
        if(GetNumberOfGammaFits() == 1) fCanvas->Divide(1);
        else if(GetNumberOfGammaFits() == 2) fCanvas->Divide(2);
        else if(GetNumberOfGammaFits() == 3) fCanvas->Divide(3);
        else if(GetNumberOfGammaFits() == 4) fCanvas->Divide(2,2);
        else { std::cout << "more/less GammaFits that I know how to draw!" << std::endl; return; }
    }
    for(int i=0; i<GetNumberOfGammaFits(); i++) {
        SortRun(i);
        fCanvas->cd(i+1);
        gPad->Clear();
        fGammaFitVector.at(i).Draw();
        gPad->Update();
    }
    fSum = 0.;
    fSum2 = 0.;
    for(int i=0;i<GetNumberOfGammaFits();i++) fSum += fGammaFitVector.at(i).DoChi2();
    for(int i=0;i<GetNumberOfGammaFits();i++) fSum2 += fGammaFitVector.at(i).DoChi2() * fGammaFitVector.at(i).DoChi2();
    fSum /= double(GetNumberOfGammaFits());
    fSum2 /= double(GetNumberOfGammaFits());
    std::cout << "sum(chi2)/nfits = " << fSum << " | sum((chi2)^2)/nfits = " << fSum2 << std::endl;
}

void GammaFitter::DrawToFile(std::string input)
{
    std::cout << "drawing all GammaFits to output file \"" << input << "\" ... " << std::flush;
    TCanvas * out = new TCanvas();
    for(int i=0; i<GetNumberOfGammaFits(); i++) 
    {
        std::string name = input;
        fGammaFitVector.at(i).Draw();
        if(i==0) {
            name += "(";
            out->Print(name.c_str(),"pdf");
        }
        else if(i==GetNumberOfGammaFits()-1) {
            name += ")";
            out->Print(name.c_str(),"pdf");
        }
        else out->Print(name.c_str(),"pdf");
    }
    delete out;
    std::cout << " done!" << std::endl;

}

//void GammaFitter::NelderMead(double A, double B, double C, double gain, double offset, int itermax)
void GammaFitter::NelderMead(double resolution, double gain, double offset, int itermax)
{
    std::cout << "starting Nelder Mead method... " << std::endl;
    
    //double A_inc = 0.2;
    //double B_inc = 0.2;
    //double C_inc = 1e-4;
    //double gain_inc = 0.5;
    //double offset_inc = 100;
    double res_inc = 0.02;
    double gain_inc = 0.01;
    double offset_inc = 1;

    //    ( A   , B    , C     , a1   , a2  , a3  , a4   , carbon)
    vec v1(resolution,gain,offset);
    vec v2(v1); v2.set(0,v1.at(0)+res_inc);
    vec v3(v1); v3.set(1,v1.at(1)+gain_inc);
    vec v4(v1); v4.set(2,v1.at(2)+offset_inc);

    std::vector<vec> nmvec;
    nmvec.push_back(v1);
    nmvec.push_back(v2);
    nmvec.push_back(v3);
    nmvec.push_back(v4);

    std::cout << "calculating chi2's for the initial simplex..." << std::endl;
    std::vector<double> chi2vec;
    for(int i=0; i<int(nmvec.size()); i++) {
        std::cout << "simplex " << i+1 << std::endl;
        SetParameters(nmvec.at(i).par_array());         
        SortAllRuns();
        chi2vec.push_back(DoChi2());
    } 
        
    std::cout << "starting the Nelder-Mead iterations..." << std::endl;
    //////////////////////////////////////////////////////////////////
    for(int iter=1; iter<=itermax; iter++) {

        std::vector<vec> temp_par;
        std::vector<double> temp_chi2;
        temp_par.resize(4);
        temp_chi2.resize(4);

        // reordering...
        double test = 1e100;
        int val = 0;
        for(int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                if(chi2vec.at(j) < test) {
                    test = chi2vec.at(j);
                    temp_chi2.at(i) = test;
                    temp_par.at(i) = nmvec.at(j);
                    val = j;
                }
            }
            chi2vec.at(val) = 1e100;
            test = 1e100;
        }
        nmvec = temp_par;
        chi2vec = temp_chi2;

        std::cout << "printing the reordered variables..." << std::endl;
        for(int i=0; i<4; i++) {
            std::cout << " chi2 = " << chi2vec.at(i);
            std::cout << " pars = ";
            for(int j=0; j<2; j++) std::cout << nmvec.at(i).at(j) << " , "; std::cout << nmvec.at(i).at(2);
            std::cout << std::endl;
        }
        
    
        vec B(nmvec.at(0)); double B_chi2 = chi2vec.at(0);
        vec G(nmvec.at(1)); double G_chi2 = chi2vec.at(1);
        vec W(nmvec.at(3)); double W_chi2 = chi2vec.at(3);
        vec M = B.midpoint(G); double M_chi2 = nm_val(M);
        vec R = M.scalar_multiply(2.); R.subtract(W); double R_chi2 = nm_val(R);
        vec E = R.scalar_multiply(2.); E.subtract(M); double E_chi2 = 0;
        vec C; double C_chi2 = 0;        
        vec S; double S_chi2 = 0;    

        // now with the logical decisions....
        if(R_chi2 < G_chi2) {  // case 1
            if(B_chi2 < R_chi2) {
                W = R; W_chi2 = R_chi2;
            }
            else {
                E_chi2 = nm_val(E);
                if(E_chi2 < B_chi2) {
                    W = E; W_chi2 = E_chi2;   
                }
                else {
                    W = R; W_chi2 = R_chi2;
                }        
            }
        }
        else {  // case 2
            if(R_chi2 < W_chi2) {
                W = R; W_chi2 = R_chi2;
            }
            vec C1 = W.midpoint(M); double C1_chi2 = nm_val(C1);
            vec C2 = M.midpoint(R); double C2_chi2 = nm_val(C2);
            if(C1_chi2 < C2_chi2) { C = C1; C_chi2 = C1_chi2; }
            else                  { C = C2; C_chi2 = C2_chi2; }
        
            if(C_chi2 < W_chi2) {
                W = C; W_chi2 = C_chi2;
            }
            else {
                S = B.midpoint(W); S_chi2 = nm_val(S);
                W = S; W_chi2 = S_chi2;
                G = M; G_chi2 = M_chi2;
            }
        }
        nmvec.at(0) = B; chi2vec.at(0) = B_chi2;
        nmvec.at(1) = G; chi2vec.at(1) = G_chi2;
        nmvec.at(3) = W; chi2vec.at(3) = W_chi2;
    
        std::cout << std::endl << "finished iteration # " << iter << "/" << itermax << std::endl << std::endl;
        
        // end of logical loop
    }
    //////////////////////////////////////////////////////////////////
    Run(nmvec.at(0).at(0), nmvec.at(0).at(1), nmvec.at(0).at(2));
}

void GammaFitter::NelderMead() 
{
    if(fIterMax < 1) fIterMax = 50;
    NelderMead(fParameters[0],fParameters[1],fParameters[2],fIterMax);
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "TSystem.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TH1F.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "../../interface/RooParametricHist.h"
//#include "/afs/cern.ch/user/v/vmilosev/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "RooAddition.h"

enum PROCESS{
  data = 0,
  VBFH = 1,
  ggH = 2,
  QCDZnunu = 3,
  EWKZnunu = 4,
  QCDW = 5,
  EWKW = 6,
  QCDDYll = 7,
  EWKZll = 8
};


int makeWS_percategory(std::string year="2017", std::string cat="MTR"){


    const bool is2017 = year=="2017";
    // As usual, load the combine library to get access to the RooParametricHist
    gSystem->Load("libHiggsAnalysisCombinedLimit.so");

    ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //Hardcoded input parameters
    //To be adjusted to channel !

    std::string lChannel = "VBF";
    std::string lCategory = cat+"_";
    std::string lYear = year+"_";
    std::string lOutFileName = "param_ws_"+lCategory+lYear+lChannel+".root";

    //define the variable that is fitted
    RooRealVar lVarFit(("mjj_"+cat+"_"+year).c_str(),"M_{jj} (GeV)",200,5000);
    std::string lVarLabel = "Mjj";

    // Regions
    // Use same order! SR=region 0.
    const unsigned nR = 5;
    std::string lRegions[5] = {"SR","Wenu","Wmunu","Zee","Zmumu"};
    //processes
    //use same order: data= process 0, signal = process 1, QCD Z+Jets in SR = 2, etc....
    const unsigned nP = 12;
    std::string lProcs[nP] = {"data_obs","VBFHtoInv","GluGluHtoInv",
			      "ZJETS","EWKZNUNU",
			      "WJETS","EWKW",
			      "DY","EWKZll",
			      "TOP","VV","QCD"};

    const unsigned nN = 20;
    std::string lNuis[20] = {"bjet_veto","pileup","tau_veto",
			     "eventVetoVEleIdIso","eventVetoLMuId","eventVetoLMuIso",
			     "eventSelTEleIdIso","eventSelTMuId","eventSelTMuIso",
      			     "eventSelVEleIdIso","eventSelLMuId","eventSelLMuIso",
                             "fnlo_SF_QCD_corr_EWK_proc", "fnlo_SF_EWK_corr",
                             "fnlo_SF_QCD_corr_QCD_proc_muF","fnlo_SF_QCD_corr_QCD_proc_muR","fnlo_SF_QCD_corr_QCD_proc_pdf",
                             "eventSelVEleReco", "eventSelTEleReco", "eventVetoVEleReco"
    };

    const bool corrCat[20] = {1,1,1,
			      1,1,1,
			      1,1,1,
			      1,1,1,
			      0,0,
			      0,0,0,
			      1,1,1
    };
    const bool corrYear[20] = {1,1,0,
			       1,1,1,
			       1,1,1,
			       1,1,1,
			       1,1,
			       1,1,1,
			       1,1,1
    };

    const unsigned nS = 2*nN+1;
    std::string lSysts[41];
    for (unsigned iS(0); iS<nS; ++iS){
      if (iS==0) lSysts[iS] = "";
      else lSysts[iS] = (iS-1)%2==0? lNuis[(iS-1)/2]+"Up" : lNuis[(iS-1)/2]+"Down";
      //std::cout << lSysts[iS] << std::endl;
    }
    const bool isSRsyst[41] = {1,1,1,1,1,
			       1,1,1,1,1,
			       1,1,1,0,0,
			       0,0,0,0,0,
			       0,0,0,0,0,
                               1,1,1,1,
                               1,1,1,1,1,1,
                               0,0,0,0,1,1};

    const bool isCRWsyst[41] = {1,1,1,1,1,
				1,1,0,0,0,
				0,0,0,1,1,
				1,1,1,1,0,
				0,0,0,0,0,
                                1,1,1,1,
                                1,1,1,1,1,1,
                                0,0,1,1,0,0};

    const bool isCRZsyst[41] = {1,1,1,1,1,
                                1,1,0,0,0,
                                0,0,0,1,1,
                                1,1,1,1,1,
                                1,1,1,1,1,
                                1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,0,0};			       

    //possibility to override nuisances with hardcoded values
    double hardCodeNuisance[nR][nS];
    
    for (unsigned iR(0); iR<nR; ++iR){     
      for (unsigned iS(0); iS < nS; ++iS){
	hardCodeNuisance[iR][iS] = -1.;
      }
    }
    //@FIXME correct lepton IDISO for now
    //hardCodeNuisance[0][7] = 0.97;
    //hardCodeNuisance[0][8] = 1.03;
    for (unsigned iR(1); iR<nR; ++iR){     
      //VetoVEleIdIso
      //just for W regions
      //inverted for veto, but acting on SR = denominator -> inverted twice...
      if (iR<3) hardCodeNuisance[iR][7] = 1.03;
      if (iR<3) hardCodeNuisance[iR][8] = 0.97;
      //skip muon regions
      if (iR==2 || iR==4) continue;
      //SelTEleIdIso
      hardCodeNuisance[iR][13] = iR==1? 1.03 : 1.06;
      hardCodeNuisance[iR][14] = iR==1? 0.97 : 0.94;
      //SelVEleIdIso // as counted double for tight already for Zll, just ignore this one...
      hardCodeNuisance[iR][19] = iR==1? 1. : 1.;
      hardCodeNuisance[iR][20] = iR==1? 1. : 1.;
    }
    //input file path and name
    //input file is expected to contain one directory per region with names as in lRegions,
    //and one histogram per process with name as in lProcs with shape of the variable that is fitted.
    std::string lInFileName[nR];
    for (unsigned iR(0); iR<nR; ++iR){     
      lInFileName[iR] = "out_VBF_ana_"+lRegions[iR]+"_"+year+"_v"+cat+"_"+year+"_200109/VBF_shapes.root";
    }

    std::string lInFileName_Sam = "WJetsToLNu.root";
     
    //indices of QCD or EWK Z/W processes in SR/CR in lProcs array. 
    const unsigned nT = 2;
    std::string lType[nT] = {"QCD","EWK"};
    unsigned vproc[nT][5] = {
      {PROCESS::QCDZnunu,PROCESS::QCDW,PROCESS::QCDW,PROCESS::QCDDYll,PROCESS::QCDDYll},
      {PROCESS::EWKZnunu,PROCESS::EWKW,PROCESS::EWKW,PROCESS::EWKZll,PROCESS::EWKZll}
    };


    //values of nuisances
    double WZratioSyst = 1.12;

    double jesWWSyst[nT] = {is2017 ? 1.015 : 1.01, 1.01};
    double jerWWSyst[nT] = {is2017 ? 1.015 : 1.01, 1.01};

    double jesZZSyst = 1.01;
    double jerZZSyst = 1.01;

    double jesWZSyst[nT] = {1/1.02, 1/1.01};
    double jerWZSyst[nT] = {is2017 ? 1/1.025 : 1/1.01,is2017 ? 1/1.015 : 1/1.01};

    //inverted because vetos...acting on numerator W.
    double eleRecoVetoWZ = 1/1.01;
    double eleIdIsoVetoWZ = 1/1.03;
    double muIdVetoWZ = 1/1.005;
    double muIsoVetoWZ = 1/1.001;
    double tauVetoWZ = 1/1.01;

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Output file and workspace

    TFile *fOut = new TFile(lOutFileName.c_str(),"RECREATE");
    RooWorkspace wspace("wspace","wspace");

    RooArgList vars(lVarFit);

    //-- Define the shapes and binning -- read from input files

    TFile *finput[nR];
    TFile *finput_sam;

    TH1F *histos[nR][nP][nS];
    for (unsigned iR(0); iR<nR; ++iR){
      for (unsigned iP(0); iP<nP; ++iP){
	for (unsigned iS(0); iS < nS; ++iS){
	  histos[iR][iP][iS] = 0;
	}
      }
    }
    
    //TH1D *sam_qcd_w_histo_nom;
    //TH1D *sam_qcd_w_histo_muR;
    //TH1D *sam_qcd_w_histo_muF;
    //TH1D *sam_qcd_w_histo_Syst_pdf;
    
    //std::cout<< "Opening Sam's." << std::endl;
    //finput_sam = TFile::Open(lInFileName_Sam.c_str());
    //if (!finput_sam) std::cout<< "Sam's file could not be opened." << std::endl;
    //sam_qcd_w_histo_nom = (TH1D*)gDirectory->Get("mjj"); 
    //if (!sam_qcd_w_histo_nom) std::cout<< "Sam's mjj_nom could not be opened." << std::endl;
    //sam_qcd_w_histo_muR = (TH1D*)gDirectory->Get("mjj_Renorm_Up");          
    //sam_qcd_w_histo_muF = (TH1D*)gDirectory->Get("mjj_Fact_Up");          
    //sam_qcd_w_histo_Syst_pdf = (TH1D*)gDirectory->Get("mjj_PDF_Up");          
      
  
    for (unsigned iR(0); iR<nR; ++iR){
      
      finput[iR] = TFile::Open(lInFileName[iR].c_str());
      finput[iR]->cd(lRegions[iR].c_str());
      
      for (unsigned iP(0); iP<nP; ++iP){
	for (unsigned iS(0); iS < nS; ++iS){
	  if (iS==0 || 
	      (iR==0 && !isSRsyst[iS]) ||
	      (iR>0 && iR<3 && !isCRWsyst[iS]) ||
	      (iR>2 && !isCRZsyst[iS])
	      ) histos[iR][iP][iS] = (TH1F*)gDirectory->Get(lProcs[iP].c_str());
	  else if (iP>0) {
	    histos[iR][iP][iS] = (TH1F*)gDirectory->Get((lProcs[iP]+"_"+lSysts[iS]).c_str());
	  }
	  if (!histos[iR][iP][iS]) {
	    if (iP>0){
	      std::cout<< " Histo not found for region " << lRegions[iR] << " process " << lProcs[iP] << " syst " << lSysts[iS] << std::endl;
	      histos[iR][iP][iS] = (TH1F*)gDirectory->Get(lProcs[iP].c_str());
	    }
	    if (iS==0) return 1;
	    continue;
	  }
	  std::cout << " --- histos " << histos[iR][iP][iS]->GetName() << " " << histos[iR][iP][iS]->GetEntries() << " " << histos[iR][iP][iS]->Integral() << std::endl;
	}
      }
    }
    
    const unsigned nB = histos[0][PROCESS::data][0]->GetNbinsX();//-4;
    double bins[nB+1];
    for (unsigned iB(0); iB<nB; ++iB){
      bins[iB] = histos[0][PROCESS::data][0]->GetXaxis()->GetBinLowEdge(iB+1);//+2);
      //std::cout << " bin " << iB << " Low edge " <<  bins[iB] << std::endl;
    }
    bins[nB] = histos[0][PROCESS::data][0]->GetXaxis()->GetBinLowEdge(nB+1);
    //std::cout << " bin " << nB << " Low edge " <<  bins[nB] << std::endl;

    TH1F dummyHist("dummyHist","Dummy hist for binning",nB,bins);

    std::cout << nB << " - Check Binning: ";
    for (unsigned iB(0); iB<nB+1; ++iB){
      std::cout << bins[iB] << " ";
    }
    std::cout << std::endl;

    RooDataHist *data_hist[nR];
    RooDataHist *bkg_hist[nR][nP][nS];
    for (unsigned iR(0); iR<nR; ++iR){
      
      data_hist[iR] = new RooDataHist(("data_obs_"+lRegions[iR]).c_str(),"Data observed",vars,histos[iR][PROCESS::data][0]);
      wspace.import(*data_hist[iR]);
      
      for (unsigned iP(1); iP<nP; ++iP){
	for (unsigned iS(0); iS < nS; ++iS){
	  std::ostringstream label;
	  label << lProcs[iP] << "_hist_" << lRegions[iR];
	  if (iS>0) {
	    label << "_" << lSysts[iS];
	  }
	  if (histos[iR][iP][iS] && histos[iR][iP][iS]->GetEntries()>0){
	    bkg_hist[iR][iP][iS] = new RooDataHist(label.str().c_str(),"Background",vars,histos[iR][iP][iS]);
	    wspace.import(*bkg_hist[iR][iP][iS]);
	  }
	}
      }
    }

    // In the signal region, W and Z are tied together by their ratio from MC + theory uncertainty.
    // SR/CR remaining non-cancellations from JES uncertainties.
    RooRealVar *wzratio[nT];
    //JES fully correlated across years, JER fully uncorrelated
    RooRealVar *jes = new RooRealVar("CMS_VBF_scale_j", "JES nuisance parameter", 0);
    RooRealVar *jer = new RooRealVar(("CMS_res_j_"+year).c_str(), "JER nuisance parameter", 0);
    //RooRealVar *jesZZ[nT];
    //RooRealVar *jerZZ[nT];
    //RooRealVar *wzratioSyst_jes[nT], *wzratioSyst_jer[nT];
    RooRealVar *wzratioSyst_muR[nT], *wzratioSyst_muF[nT], *wzratioSyst_pdf[nT];
    for (unsigned iT(0); iT<nT; ++iT){     
      wzratio[iT] = new RooRealVar(("wzratio"+lType[iT]).c_str(), (lType[iT]+" W/Z ratio nuisance parameter").c_str(),0);
      std::ostringstream lname;
      //lname.str("");
      //lname << lCategory << lYear;
      //lname << lType[iT];
      //lname << "wzratioSyst_jes";
      //wzratioSyst_jes[iT] = new RooRealVar(lname.str().c_str(), (lType[iT]+" W/Z JES ratio nuisance parameter per bin").c_str(),0);
      //lname.str("");
      //lname << lCategory << lYear;
      //lname << lType[iT];
      //lname << "wzratioSyst_jer";
      //wzratioSyst_jer[iT] = new RooRealVar(lname.str().c_str(), (lType[iT]+" W/Z JER ratio nuisance parameter per bin").c_str(),0);
      lname.str("");
      lname << lCategory;
      lname << lType[iT];
      lname << "wzratioQCDcorrSyst_muR";
      wzratioSyst_muR[iT] = new RooRealVar(lname.str().c_str(), (lType[iT]+" W/Z muR ratio nuisance parameter per bin").c_str(),0);
      lname.str("");
      lname << lCategory;
      lname << lType[iT];
      lname << "wzratioQCDcorrSyst_muF";
      wzratioSyst_muF[iT] = new RooRealVar(lname.str().c_str(), (lType[iT]+" W/Z muF ratio nuisance parameter per bin").c_str(),0);
      lname.str("");
      lname << lCategory;
      lname << lType[iT];
      lname << "wzratioQCDcorrSyst_pdf";
      wzratioSyst_pdf[iT] = new RooRealVar(lname.str().c_str(), (lType[iT]+" W/Z Syst pdf ratio nuisance parameter per bin").c_str(),0);
    
 
      //jesZZ[iT] = new RooRealVar((lCategory+lYear+"jesZZ"+lType[iT]).c_str(), (lType[iT]+" Z/Z JES nuisance parameter").c_str(),0);
      //jerWW[iT] = new RooRealVar((lCategory+lYear+"jerWW"+lType[iT]).c_str(), (lType[iT]+" W/W JER nuisance parameter").c_str(),0);
      //jerZZ[iT] = new RooRealVar((lCategory+lYear+"jerZZ"+lType[iT]).c_str(), (lType[iT]+" Z/Z JER nuisance parameter").c_str(),0);
    }


    RooRealVar *ewkqcdratiostat[nB];
    RooRealVar *wzratiostat[nT][nB];
    RooRealVar *wzratioEWK_on_strong[nT][nB];
    RooRealVar *TFstat[nT][nB][nR-1];
    RooRealVar *TFsysts[nN];
    //same nuisance name for different CR, correlated accross CR or given different systs names already....??
    for (unsigned iN(0); iN < nN; ++iN){
      std::ostringstream lname;
      lname.str("");
      if (!corrCat[iN]) lname << lCategory;
      if (!corrYear[iN]) lname << lYear;
     // changing to match the naming convention used in the datacard
     // lname << "TF_syst_" << lNuis[iN];
      lname << "" << lNuis[iN];
      TFsysts[iN] = new RooRealVar(lname.str().c_str(),"CR/SR ratio syst nuisance parameter",0);
    }


    // Create one parameter per bin representing the yield. (note of course we can have multiple processes like this)
    for (unsigned iB(1); iB<nB+1; ++iB){
      std::cout << " -- Processing bin " << iB << std::endl;
      std::ostringstream lname;

      RooFormulaVar *EWKQCDbin = 0;
      
      for (unsigned iT(0); iT<nT; ++iT){     

	std::cout << " --- Processing type " << lType[iT] << std::endl;

	lname.str("");
        lname << lCategory << lYear;
	lname << lType[iT] << "Z_SR_bin" << iB;
	RooRealVar binParZ(lname.str().c_str(),(lType[iT]+" Z+jets yield in signal region, per bin").c_str(),histos[0][iT==0?PROCESS::QCDZnunu:PROCESS::EWKZnunu][0]->GetBinContent(iB),0,10*histos[0][iT==0?PROCESS::QCDZnunu:PROCESS::EWKZnunu][0]->GetBinContent(iB));

	std::cout << "Importing " << binParZ.GetName() <<std::endl;
	wspace.import(binParZ,RooFit::RecycleConflictNodes());
	std::cout << " ..... done " <<std::endl;

	if (iT==0) {
	  
	  //program link between QCD and EWK yields in SR
	  lname.str("");
          lname << lCategory << lYear;
	  lname << "ewkqcdratio_stat_bin" << iB;
	  ewkqcdratiostat[iB] = new RooRealVar(lname.str().c_str()," EWK/QCD ratio stat nuisance parameter",0);
	  lname.str("");
          lname << lCategory << lYear;
	  lname << "TF_EWKQCDSR_bin" << iB;
	  std::ostringstream lFormula;
	  double ratio = histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB) / histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB); 
	  double ratiostat = 1+sqrt(pow(histos[0][PROCESS::EWKZnunu][0]->GetBinError(iB)/histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB),2)+pow(histos[0][PROCESS::QCDZnunu][0]->GetBinError(iB)/histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB),2));
	  lFormula << ratio;
	  lFormula << "*TMath::Power(" << ratiostat << ",@0)";
	  
	  RooFormulaVar TFEWKQCD(lname.str().c_str(),"Transfer factor EWK/QCD Z",lFormula.str().c_str(),RooArgList(*(ewkqcdratiostat[iB])) );
	  wspace.import(TFEWKQCD,RooFit::RecycleConflictNodes());
	  lname.str("");
          lname << lCategory << lYear;
	  lname << "EWKQCD_SR_bin" << iB;
	  EWKQCDbin = new RooFormulaVar(lname.str().c_str(),"EWK Z+jets yield in signal regions from QCD Z yield, per bin","@0*@1",RooArgList(TFEWKQCD,binParZ));
	  wspace.import((*EWKQCDbin),RooFit::RecycleConflictNodes());


      

	} 
	//lname.str("");
	//lname << "EWKQCD_SR_bin" << iB;
	//RooFormulaVar EWKQCDbin = *wspace.function(lname.str().c_str());
	lname.str("");
        lname << lCategory << lYear;
	lname << lType[iT] << "wzratio_stat_bin" << iB;
	wzratiostat[iT][iB] = new RooRealVar(lname.str().c_str(),"W/Z ratio stat nuisance parameter",0);
        lname.str("");
        lname << lCategory;
        lname << lType[iT] << "wzratio_EWK_corr_on_Strong_bin" << iB;        
        wzratioEWK_on_strong[iT][iB] = new RooRealVar(lname.str().c_str(),"W/Z ratio EWK correction on Strong processes nuisance parameter",0);
         

	lname.str("");
        lname << lCategory << lYear;
	lname << lType[iT] << "TF_WZSR_bin" << iB;
	std::ostringstream lFormula;
	lFormula.str("");
        std::cout<< "Sam's histogram bin ."<< iB << std::endl;
        //double WZratioSyst_nom = sam_qcd_w_histo_nom->GetBinContent(iB+1);
        //double WZratioSyst_muF = sam_qcd_w_histo_muF->GetBinContent(iB+1);
        //double WZratioSyst_muR = sam_qcd_w_histo_muR->GetBinContent(iB+1);
        //double WZratioSyst_pdf = sam_qcd_w_histo_Syst_pdf->GetBinContent(iB+1);

        double WZratioSyst_nom = histos[0][PROCESS::QCDW][0]->GetBinContent(iB);
        double WZratioSyst_muF = histos[0][PROCESS::QCDW][29]->GetBinContent(iB);
        double WZratioSyst_muR = histos[0][PROCESS::QCDW][31]->GetBinContent(iB);
        double WZratioSyst_pdf = histos[0][PROCESS::QCDW][33]->GetBinContent(iB);

       
        WZratioSyst_muF = 1.0*WZratioSyst_muF/WZratioSyst_nom;
        WZratioSyst_muR = 1.0*WZratioSyst_muR/WZratioSyst_nom;
        WZratioSyst_pdf = 1.0*WZratioSyst_pdf/WZratioSyst_nom;
        std::cout<< "WZratioSyst_muF :  " << WZratioSyst_muF << std::endl;
        std::cout<< "WZratioSyst_muR :  " << WZratioSyst_muR << std::endl;
        std::cout<< "WZratioSyst_pdf :  " << WZratioSyst_pdf << std::endl;
 
	double ratio = 0;
	double ratiostat = 0;
        double ratio_EWK_corr_on_Strong_proc = (histos[0][PROCESS::QCDW][0]->GetBinContent(iB) / histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB))/(histos[0][PROCESS::QCDW][27]->GetBinContent(iB) / histos[0][PROCESS::QCDZnunu][27]->GetBinContent(iB));
        std::cout<<"Test EWK corr "<<ratio_EWK_corr_on_Strong_proc<<std::endl;


	if (iT==0) {
	  ratio = histos[0][PROCESS::QCDW][0]->GetBinContent(iB) / histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB);
	  ratiostat = 1+sqrt(pow(histos[0][PROCESS::QCDW][0]->GetBinError(iB)/histos[0][PROCESS::QCDW][0]->GetBinContent(iB),2)+pow(histos[0][PROCESS::QCDZnunu][0]->GetBinError(iB)/histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB),2));
	} else{
	  ratio = histos[0][PROCESS::EWKW][0]->GetBinContent(iB) / histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB);
	  ratiostat = 1+sqrt(pow(histos[0][PROCESS::EWKW][0]->GetBinError(iB)/histos[0][PROCESS::EWKW][0]->GetBinContent(iB),2)+pow(histos[0][PROCESS::EWKZnunu][0]->GetBinError(iB)/histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB),2));
	}
	lFormula << ratio;
	if (iT!=0){
         lFormula << "*TMath::Power(" << WZratioSyst_muR << ",@0)*TMath::Power(" << WZratioSyst_muF << ",@1)*TMath::Power(" << WZratioSyst_pdf << ",@2)*TMath::Power(" << ratiostat << ",@3)*TMath::Power(" << ratio_EWK_corr_on_Strong_proc << ",@4)";
	} else {
	 lFormula << "*TMath::Power(" << WZratioSyst_muR << ",@0)*TMath::Power(" << WZratioSyst_muF << ",@1)*TMath::Power(" << WZratioSyst_pdf << ",@2)*TMath::Power(" << ratiostat << ",@3)*TMath::Power(" << ratio_EWK_corr_on_Strong_proc << ",@4)";
	}
	 //add JES/JER
	 lFormula << "*TMath::Power(" << jesWZSyst[iT] << ",@5)";
	 lFormula << "*TMath::Power(" << jerWZSyst[iT] << ",@6)";
	 //add also lepton veto uncertainties
	 //lFormula << "*TMath::Power(" << eleRecoVetoWZ << ",@7)";
	 lFormula << "*TMath::Power(" << eleIdIsoVetoWZ << ",@7)";
	 //lFormula << "*TMath::Power(" << muIdVetoWZ << ",@9)";
	 //lFormula << "*TMath::Power(" << muIsoVetoWZ << ",@10)";
	 lFormula << "*TMath::Power(" << tauVetoWZ << ",@8)";


	std::cout << " ---- Check stat error " << lType[iT] << " WZratio: " << ratiostat << std::endl;
        RooFormulaVar *TFWZ;
	//C-AMM: why the if-else loop below ? Everything is given with [iT] ...
        if (iT != 0){	
	  std::cout << " Ok -> iT=1 "<< lname.str().c_str() << ", " << lFormula.str().c_str() << std::endl;
          //RooArgList variables( *(wzratioSyst_muR[iT]), *(wzratioSyst_muF[iT]), *(wzratioSyst_pdf[iT]),*(wzratiostat[iT][iB]), *(wzratioEWK_on_strong[iT][iB]), *(wzratioSyst_jes[iT]), *(wzratioSyst_jer[iT]));
          RooArgList variables( *(wzratioSyst_muR[iT]), *(wzratioSyst_muF[iT]), *(wzratioSyst_pdf[iT]),*(wzratiostat[iT][iB]), *(wzratioEWK_on_strong[iT][iB]), *(jes), *(jer), *(TFsysts[3]), *(TFsysts[2]));
          variables.Print();
          TFWZ = new RooFormulaVar(lname.str().c_str(),"Transfer factor W/Z",lFormula.str().c_str(),variables );
        
	  wspace.import(*TFWZ,RooFit::RecycleConflictNodes());
        }
        else{
	  std::cout << " Ok -> iT=0 " << lname.str().c_str() << ", " << lFormula.str().c_str() << std::endl;
	  RooArgList variables( *(wzratioSyst_muR[iT]), *(wzratioSyst_muF[iT]), *(wzratioSyst_pdf[iT]),*(wzratiostat[iT][iB]), *(wzratioEWK_on_strong[iT][iB]), *(jes), *(jer), *(TFsysts[3]), *(TFsysts[2]));
	  variables.Print();
          TFWZ = new RooFormulaVar(lname.str().c_str(),"Transfer factor W/Z",lFormula.str().c_str(),variables );
          //std::cout<<RooArgList( *(wzratioSyst_muR[iB]), *(wzratioSyst_muF[iB]), *(wzratioSyst_pdf[iB]),*(wzratiostat[iT][iB]))<<std::endl;
	  TFWZ->Print("V");
          wspace.import(*TFWZ,RooFit::RecycleConflictNodes());
        }
        std::cout<<"Proslo ok"<<std::endl;
        lname.str("");
        lname << lCategory << lYear;
        lname << lType[iT] << "WZ_SR_bin" << iB;
        RooFormulaVar WZbin(lname.str().c_str(),(lType[iT]+" W+jets yield in signal regions from Z yield, per bin").c_str(),"@0*@1",iT==0?RooArgList(*TFWZ,binParZ):RooArgList(*TFWZ,*(EWKQCDbin)));
        wspace.import(WZbin,RooFit::RecycleConflictNodes());

	
	if (iT==0) std::cout << " --- SR QCD Z yield = " << histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB) << " SR QCD W yield = " << histos[0][PROCESS::QCDW][0]->GetBinContent(iB) << " SR signal yield = " << histos[0][PROCESS::VBFH][0]->GetBinContent(iB) << std::endl;
	else std::cout << " --- SR EWK Z yield = " << histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB) << " SR EWK W yield = " << histos[0][PROCESS::EWKW][0]->GetBinContent(iB) << " SR signal yield = " << histos[0][PROCESS::VBFH][0]->GetBinContent(iB) << std::endl;
	
	for (unsigned iR(1); iR<nR; ++iR){
	  
	  std::cout << " ---- Doing control regions: " << lRegions[iR] << std::endl;
	  
	  lname.str("");
          lname << lCategory << lYear;
	  lname << lType[iT] << "TF_" << lRegions[iR] << "_stat_bin" << iB;
	  TFstat[iT][iB][iR-1] = new RooRealVar(lname.str().c_str(),"CR/SR ratio stat nuisance parameter",0);
	  ratio = 0;
	  ratiostat = 0;
	  std::cout << " ----- Check process names for CR/SR ratio: " << lProcs[vproc[iT][iR]] << " with " ;
	  if (iR<3) std::cout << lProcs[vproc[iT][iR]];
	  else std::cout << lProcs[vproc[iT][0]];
	  std::cout << std::endl;
	  if (iR<3) {
	    //WCR / WSR
	    ratio = histos[iR][vproc[iT][iR]][0]->GetBinContent(iB) / histos[0][vproc[iT][iR]][0]->GetBinContent(iB);
	    ratiostat = 1+sqrt(pow(histos[iR][vproc[iT][iR]][0]->GetBinError(iB)/histos[iR][vproc[iT][iR]][0]->GetBinContent(iB),2)+pow(histos[0][vproc[iT][iR]][0]->GetBinError(iB)/histos[0][vproc[iT][iR]][0]->GetBinContent(iB),2));
	  }
	  else {
	    //ZCR / ZSR
	    ratio = histos[iR][vproc[iT][iR]][0]->GetBinContent(iB) / histos[0][vproc[iT][0]][0]->GetBinContent(iB);
	    ratiostat = 1+sqrt(pow(histos[iR][vproc[iT][iR]][0]->GetBinError(iB)/histos[iR][vproc[iT][iR]][0]->GetBinContent(iB),2)+pow(histos[0][vproc[iT][0]][0]->GetBinError(iB)/histos[0][vproc[iT][0]][0]->GetBinContent(iB),2));
	  }
	  lFormula.str("");
	  lFormula << ratio;
	  lFormula << "*TMath::Power(";
	  if (iR<3) lFormula << jesWWSyst[iT];
	  else lFormula << jesZZSyst;
	  lFormula << ",@0)";
	  lFormula << "*TMath::Power(";
	  if (iR<3) lFormula << jerWWSyst[iT];
	  else lFormula << jerZZSyst;
	  lFormula << ",@1)";
	  lFormula << "*TMath::Power(" << ratiostat << ",@2)";

	  RooArgList nuisances;
	  //nuisances.add(iR<3? *(jesWW[iT]) : *(jesZZ[iT]));
	  //nuisances.add(iR<3? *(jerWW[iT]) : *(jerZZ[iT]));
	  nuisances.add(*(jes));
	  nuisances.add(*(jer));
	  nuisances.add(*(TFstat[iT][iB][iR-1]));
	  for (unsigned iN(0); iN < nN; ++iN){
	    // 3 = JES+JER+stat: change if adding syst beforehand !!
	    unsigned iSyst = 3+iN;
	    double ratiovar[2];
	    double ratiosyst[2];
	    //get up and down variations
	    for (unsigned iV(0); iV<2; ++iV){
	      unsigned iS = 2*iN+1+iV;
	      ratiovar[iV] = iR<3 ? histos[iR][vproc[iT][iR]][iS]->GetBinContent(iB) / histos[0][vproc[iT][iR]][iS]->GetBinContent(iB):
		histos[iR][vproc[iT][iR]][iS]->GetBinContent(iB) / histos[0][vproc[iT][0]][iS]->GetBinContent(iB);
	      if (hardCodeNuisance[iR][iS]<0) ratiosyst[iV] = 1+(ratiovar[iV]-ratio)/ratio;
	      else ratiosyst[iV] = hardCodeNuisance[iR][iS];
	      if (ratiosyst[iV] < 0){
		std::cout << " -- ERROR in systematics variations! For process " << lProcs[vproc[iT][iR]] << " region " << lRegions[iR] << " syst " << lSysts[iS] << " bin " << iB << ": ratio = " << ratio << " ratiovar = " << ratiovar[iV] << std::endl;
		return 1;
	      }
	      std::cout << std::setprecision(10) << " ------ bin " << iB << " type " << lType[iT] << " region " << lRegions[iR] << " Check syst " << lSysts[iS] << " " << iV << " " << ratiosyst[iV] << std::endl;
	    }
	    if (ratiovar[0]<ratio){
	      std::cout << " -- INFO: up variations is actually giving smaller ratio...." << lProcs[vproc[iT][iR]] << " region " << lRegions[iR] << " syst " << lNuis[iN] << " bin " << iB << ": ratio = " << ratio << " ratiovar = " << ratiovar[0] << std::endl;
	    }
	    if (ratiovar[1]>ratio){
	      std::cout << " -- INFO: down variations is actually giving larger ratio...." << lProcs[vproc[iT][iR]] << " region " << lRegions[iR] << " syst " << lNuis[iN] << " bin " << iB << ": ratio = " << ratio << " ratiovar = " << ratiovar[1] << std::endl;
	    }
	    //take sign of @i to decide up / down...
	    lFormula << "*( (@" << iSyst << ">=0)*TMath::Power(" << ratiosyst[0] << ",@" << iSyst << ")+(@" << iSyst << "<0)*TMath::Power(" << 1./ratiosyst[1] << ",@" << iSyst << "))";
	    nuisances.add(*(TFsysts[iN]));
	  }


	  lname.str("");
          lname << lCategory << lYear;
	  lname << lType[iT] << "TF_" << lRegions[iR] << "_bin" << iB;
	  RooFormulaVar TF(lname.str().c_str(),"Transfer factor CR/SR",lFormula.str().c_str(),nuisances);
	  wspace.import(TF,RooFit::RecycleConflictNodes());
	  lname.str("");
          lname << lCategory << lYear;
	  lname << lType[iT] << "V_" << lRegions[iR] << "_bin" << iB;
	  RooFormulaVar CRbin(lname.str().c_str(),(lType[iT]+" V+jets yield in control regions, per bin").c_str(),"@0*@1",iR<3?RooArgList(TF,WZbin):iT==0?RooArgList(TF,binParZ):RooArgList(TF,*(EWKQCDbin)));
	  wspace.import(CRbin,RooFit::RecycleConflictNodes());
	  
	  if (iR==3) std::cout << " ---- CR " << lRegions[iR] << " V yield = " << histos[iR][vproc[iT][iR]][0]->GetBinContent(iB) << " data yield " << histos[iR][0][0]->GetBinContent(iB) << std::endl;
	  
	  
	}//loop on regions
      }//loop on types


    }//loop on bins

    std::cout << " - Adding to lists: " << std::endl;
    RooArgList Z_SR_bins[nT];
    RooArgList W_SR_bins[nT];
    RooArgList V_CR_bins[nT][nR-1];
    for (unsigned iT(0); iT<nT; ++iT){     
      for (unsigned iB(0); iB<nB; ++iB){
	std::ostringstream lname;
	if (iT==0){
          lname << lCategory << lYear;
	  lname << lType[iT] << "Z_SR_bin" << iB+1;
	} else {
          lname << lCategory << lYear;
	  lname << "EWKQCD_SR_bin" << iB+1;
	}
	if ( (iT==0 && !wspace.var(lname.str().c_str())) ||
	     (iT==1 && !wspace.function(lname.str().c_str()))
	     ) {
	  std::cout << "Error for " << lType[iT] << " Z bin " << iB << " " << lname.str() << std::endl;
	  return 1;
	}
	if (iT==0) Z_SR_bins[iT].add(*wspace.var(lname.str().c_str()));
	else Z_SR_bins[iT].add(*wspace.function(lname.str().c_str()));
	lname.str("");
        lname << lCategory << lYear;
	lname << lType[iT] << "WZ_SR_bin" << iB+1;
	if (!wspace.function(lname.str().c_str())) {
	  std::cout << "Error for " << lType[iT] << " W bin " << iB << " " << lname.str() << std::endl;
	  return 1;
	}
	W_SR_bins[iT].add(*wspace.function(lname.str().c_str()));
	for (unsigned iR(1); iR<nR; ++iR){
	  lname.str("");
          lname << lCategory << lYear;
	  lname << lType[iT] << "V_" << lRegions[iR] << "_bin" << iB+1;
	  if (!wspace.function(lname.str().c_str())) {
	    std::cout << "Error for  " << lType[iT] << " Z bin " << iB << " region " << lRegions[iR] << " " << lname.str() << std::endl;
	    return 1;
	  }
	  V_CR_bins[iT][iR-1].add(*wspace.function(lname.str().c_str()));
	}
      }//loop on bins
    }//loop on types
    
    std::cout << " - Creating the SR parametric hists" << std::endl;
    
    // Create a RooParametericHist which contains those yields, last argument is just for the binning,
    // can use the data TH1 for that

    for (unsigned iT(0); iT<nT; ++iT){     
      std::cout << " -- Processing type " << lType[iT] << std::endl;

      RooParametricHist p_Z((lType[iT]+"Z_SR").c_str(), (lType[iT]+"Z+jets PDF in signal region").c_str(),lVarFit,Z_SR_bins[iT],dummyHist);
      // Always include a _norm term which should be the sum of the yields (thats how combine likes to play with pdfs)
      RooAddition p_Z_norm((lType[iT]+"Z_SR_norm").c_str(),("Total Number of events from "+lType[iT]+" Z+jets in signal region").c_str(),Z_SR_bins[iT]);
      
      RooParametricHist p_W((lType[iT]+"W_SR").c_str(), (lType[iT]+"W+jets PDF in signal region").c_str(),lVarFit,W_SR_bins[iT],dummyHist);
      RooAddition p_W_norm((lType[iT]+"W_SR_norm").c_str(),("Total Number of events from "+lType[iT]+" W+jets in signal region").c_str(),W_SR_bins[iT]);
    
      std::cout << " -- Importing the parametric hists" << std::endl;
    
      // import the pdfs
      wspace.import(p_Z);
      wspace.import(p_Z_norm,RooFit::RecycleConflictNodes());
      wspace.import(p_W);
      wspace.import(p_W_norm,RooFit::RecycleConflictNodes());
      
      std::cout << " -- Creating and importing the CR parametric hists" << std::endl;
    
      for (unsigned iR(1); iR<nR; ++iR){
	std::ostringstream lname;
        lname << lCategory << lYear;
	lname << lType[iT] << "V_" << lRegions[iR];
	RooParametricHist p_CRV(lname.str().c_str(), "Background PDF in control region",lVarFit,V_CR_bins[iT][iR-1],dummyHist);
	lname << "_norm";
	RooAddition p_CRV_norm(lname.str().c_str(),("Total Number of events from "+lType[iT]+" V+jets background in control region").c_str(),V_CR_bins[iT][iR-1]);
	wspace.import(p_CRV);
	wspace.import(p_CRV_norm,RooFit::RecycleConflictNodes());
      }
    }//loop on types
    
    std::cout << " - Printing workspace." << std::endl;
    wspace.Print();
    
    fOut->cd();
    wspace.Write();
    
    std::cout << " - Workspace " << wspace.GetName() << " written." << std::endl;
    
    // Clean up
    fOut->Close();
    fOut->Delete();
    
    return 0;
    
}

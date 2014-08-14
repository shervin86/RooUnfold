//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfoldExample()
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
  double weight=10;
  cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse response (40, -10.0, 10.0);

  // Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i= 0; i<100000; i++) {
    Double_t xt= gRandom->BreitWigner (0.3, 2.5);
    Double_t x= smear (xt);
    if (x!=cutdummy)
      response.Fill (x, xt,weight);
    else
      response.Miss (xt, weight);
  }

  cout << "==================================== TEST =====================================" << endl;
  TH1D* hTrue= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  // Test with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<10000; i++) {
    Double_t xt= gRandom->Gaus (0.0, 2.0), x= smear (xt);
    hTrue->Fill(xt, weight);
    if (x!=cutdummy) hMeas->Fill(x, weight);
  }

  cout << "==================================== UNFOLD ===================================" << endl;
  RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
//RooUnfoldTUnfold unfold (&response, hMeas);

  TH1D* hReco= (TH1D*) unfold.Hreco();

  unfold.PrintTable (cout, hTrue);
  double yMax= (hMeas->GetMaximum() > hTrue->GetMaximum()) ?  hMeas->GetMaximum() : hTrue->GetMaximum();
  
  hMeas->SetFillColor(kGray);
  hMeas->SetFillStyle(1001);
  hMeas->SetMaximum(yMax*1.2);
  hMeas->Draw("");
  hTrue->SetLineColor(2);
  hTrue->SetLineWidth(2);
  hTrue->Draw("SAME");
  hReco->Draw("same");
  
  
  TLegend *leg = new TLegend(0.7,0.7,1,1);
  leg->SetBorderSize(1);
  leg->AddEntry(hReco, "reco","p");
  leg->AddEntry(hMeas, "meas","lf");
  leg->AddEntry(hTrue, "true", "l");
  leg->Draw();
}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif

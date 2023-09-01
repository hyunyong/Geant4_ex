//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1Analysis.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4INCLGlobals.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }
  auto analysisManager = G4AnalysisManager::Instance();
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  G4Track* tr = step->GetTrack();

  if (tr->GetParticleDefinition()->GetParticleName() == "gamma") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(0, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "e-") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(1, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "e+") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(2, kE/MeV, theta);
  }
 
    //auto kEV = tr->GetVertexKineticEnergy();
    //auto kE = tr->GetDynamicParticle()->GetKineticEnergy();
    //auto theta = tr->GetMomentumDirection().theta();
    //analysisManager->FillH2(0, kE/MeV, theta);
    //auto parentID = tr->GetParentID();
    //auto trID = tr->GetTrackID();
    //G4cout << "stepAct trID: " << trID << ", pName: " << tr->GetParticleDefinition()->GetParticleName() << ", kE: " <<kE << ", kEV: " << kEV <<", parentID: " << parentID << G4endl;

    /*
    if (z < 0) return;
    if (z > 140) return;
    int zIdx = int(z);
    if (kE/GeV < 1) {
      analysisManager->FillH2(zIdx, kE/MeV, theta);}
    */
    /*else {
      zIdx+=14;
      analysisManager->FillH2(zIdx, kE/GeV, theta);}*/
    /*
    if (kE/GeV < 1) {
      analysisManager->FillH2(0, kE/MeV, z);
      analysisManager->FillH2(2, kE/MeV, theta);}
    
    else {
      analysisManager->FillH2(1, kE/GeV, z);
      analysisManager->FillH2(3, kE/GeV, theta);}
   
    analysisManager->FillNtupleDColumn(0, kE/MeV);
    analysisManager->FillNtupleDColumn(1, theta);
    analysisManager->FillNtupleDColumn(2, z);
    analysisManager->AddNtupleRow();
    */
  //}
  /*
  if (tr->GetParticleDefinition()->GetParticleName() == "gamma") {
    analysisManager->FillNtupleDColumn(4, tr->GetKineticEnergy());
    analysisManager->AddNtupleRow();
  }*/
  /*
  if ((tr->GetParticleDefinition()->GetParticleName() == "pi+" or tr->GetParticleDefinition()->GetParticleName() == "pi-") and tr->GetTrackStatus() == 2) {
    G4StepPoint* pStep = step->GetPreStepPoint();
    const G4String& pName = tr->GetParticleDefinition()->GetParticleName();
    //G4cout << pName << G4endl;
    if (pName == "pi+") analysisManager->FillNtupleDColumn(0, pStep->GetKineticEnergy());
    if (pName == "pi-") analysisManager->FillNtupleDColumn(2, pStep->GetKineticEnergy());
    bool muFlag = false;
    bool nuFlag = false;
    const G4TrackVector* dps = step->GetSecondary();
    for (auto dp : *dps) {
      if (dp->GetParticleDefinition()->GetParticleName() == "mu+" or dp->GetParticleDefinition()->GetParticleName() == "mu-") {
        const G4String& mName = dp->GetParticleDefinition()->GetParticleName();
        if (mName == "mu+" and pName == "pi+") analysisManager->FillNtupleDColumn(1, dp->GetKineticEnergy());
        if (mName == "mu-" and pName == "pi-") analysisManager->FillNtupleDColumn(3, dp->GetKineticEnergy());
        muFlag = true;
      }
      //G4cout << dp->GetParticleDefinition()->GetParticleName() << G4endl;
      if  ((dp->GetParticleDefinition()->GetParticleName() == "nu_mu" and pName == "pi+") or (dp->GetParticleDefinition()->GetParticleName() == "anti_nu_mu" and pName == "pi-")) nuFlag = true;
    }
    if (muFlag and nuFlag) analysisManager->AddNtupleRow(); 
  }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


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

#include "B1TrackingAction.hh"
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

B1TrackingAction::B1TrackingAction(B1EventAction* eventAction)
: G4UserTrackingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1TrackingAction::~B1TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1TrackingAction::PreUserTrackingAction(const G4Track* tr)
{

  if (!fScoringVolume) {
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();
  }

  G4LogicalVolume* volume
    = tr->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  if (volume != fScoringVolume) return;
  auto analysisManager = G4AnalysisManager::Instance();

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
  if (tr->GetParticleDefinition()->GetParticleName() == "pi-") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(3, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "pi+") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(4, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "kaon-") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(5, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "kaon+") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(6, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "pi0") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(7, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "eta") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(8, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "mu-") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(9, kE/MeV, theta);
  }
  if (tr->GetParticleDefinition()->GetParticleName() == "mu+") {
    auto kE = tr->GetKineticEnergy();
    auto theta = tr->GetMomentumDirection().theta();
    analysisManager->FillH2(10, kE/MeV, theta);
  }





    //auto kEV = tr->GetVertexKineticEnergy();
    //auto kE = tr->GetDynamicParticle()->GetKineticEnergy();
    //auto parentID = tr->GetParentID();
    //auto trID = tr->GetTrackID();
    //G4cout << "trAct trID: " << trID << ", pName: " << tr->GetParticleDefinition()->GetParticleName() << ", kE: " <<kE << ", kEV: " << kEV <<", parentID: " << parentID << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


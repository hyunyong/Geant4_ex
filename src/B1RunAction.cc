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
// $Id: B1RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1Run.hh"
#include "B1Analysis.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//#include "G4Pow.hh"
#include <math.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction()
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);        

  auto analysisManager = G4AnalysisManager::Instance();
 
  G4cout << "Using " << analysisManager->GetType() << G4endl;
  analysisManager->SetVerboseLevel(1);
  
  //analysisManager->CreateH2("EvsZ1", "photon E (1 MeV to 1000 MeV) vs. Z position (mm)", 999, 1, 1000, 120,0, 1200);
  //analysisManager->CreateH2("EvsZ2", "photon E (1 GeV to 1000 GeV) vs. Z position (mm)", 999, 1, 1000, 120,0, 1200);

  //analysisManager->CreateH3("E1_Theta_Z","photon E (1 Mev to 1000 MeV)  Theta (rad) z (cm)", 27, 1, 1000, 314, 0,3.14, 140, 0, 140, "MeV", "rad", "cm", "", "", "", "log10", "linear", "linear");
  //analysisManager->CreateH3("E2_Theta_Z","photon E (1 Gev to 1000 GeV)  Theta (rad) z (cm)", 27, 1, 1000, 314, 0,3.14, 140, 0, 140, "GeV", "rad", "cm", "", "", "", "log10", "linear", "linear");
  std::vector<double> xAxis; 
  for (double i = 1E0; i < 1E10; i *= pow(10,0.1)) {xAxis.push_back(i);} 
  std::vector<double> yAxis;
  for (double i = 1E-6; i < 1E0; i *= pow(10,0.1)) {yAxis.push_back(i);}
  analysisManager->CreateH2("gamma","photon E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("e","e E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("ep","e+ E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("pim","pi- E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("pip","pi+ E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("km","k- E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("kp","k+ E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("pi0","pi0 E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("eta","eta E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("mum","mu- E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  analysisManager->CreateH2("mup","mu+ E (1 Mev to 10 TeV)  Theta (rad) z", xAxis, yAxis, "MeV","radian");
  

  /*
  for (int i=0;i<140;i++) {
  char chIndex[20];
  sprintf(chIndex,"E1vsTheta_Z_%03d",i);
  G4String tmp = chIndex;
  analysisManager->CreateH2(tmp, "photon E (1 Mev to 1000 MeV) vs. Theta (rad) z index", 100, 0, 1000, 314, 0, 3.14);
  }*/
  /*for (int i=0;i<140;i++) {
  char chIndex[20];
  sprintf(chIndex,"E2vsTheta_Z_%02d",i);
  G4String tmp = chIndex;
  analysisManager->CreateH2(tmp, "photon E (1 GeV to 1000 GeV) vs. Theta (rad) z index", 999, 1, 1000, 314, 0, 3.14);
  }*/

  /*
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetActivation(true);
  analysisManager->CreateNtuple("gamma", "gamma"); 
  analysisManager->CreateNtupleDColumn("kE"); 
  analysisManager->CreateNtupleDColumn("theta"); 
  analysisManager->CreateNtupleDColumn("z"); 
  analysisManager->FinishNtuple(); 
  */


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* B1RunAction::GenerateRun()
{
  return new B1Run; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile(); 
 
  const B1Run* b1Run = static_cast<const B1Run*>(run);

  // Compute dose
  //
  G4double edep  = b1Run->GetEdep();
  G4double edep2 = b1Run->GetEdep2();
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Dose in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " +- " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

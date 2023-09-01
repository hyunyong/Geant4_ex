#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

B1DetectorConstruction::~B1DetectorConstruction()
{ }

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  

  G4NistManager* nist = G4NistManager::Instance();

  G4double env_sizeXY = 1*m, env_sizeZ = 10*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Fe");
  G4bool checkOverlaps = true;
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");


  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     

  G4Box* solidEnv1 =    
    new G4Box("Envelope1",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv1 =                         
    new G4LogicalVolume(solidEnv1,            //its solid
                        env_mat,             //its material
                        "Envelope1");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,env_sizeZ/2.0),         //at (0,0,0)
                    logicEnv1,                //its logical volume
                    "Envelope1",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  fScoringVolume = logicEnv1;
  return physWorld;
}

 



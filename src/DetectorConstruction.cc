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
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

using namespace B2;

namespace B2a
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::DetectorConstruction()
  {
    G4NistManager *nistManager = G4NistManager::Instance();
    // Air defined using NIST Manager
    nistManager->FindOrBuildMaterial("G4_AIR");

    nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_Ar");
    nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    


  
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::~DetectorConstruction()
  {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume *DetectorConstruction::Construct()
  {

    G4String trackerChamberSDname = "/straw";
    TrackerSD* aTrackerSD = new TrackerSD(trackerChamberSDname
                                              );
    G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
    // Setting aTrackerSD to all logical volumes with the same name
    // of "Chamber_LV".

    G4double gas_density = 0.001788*g/cm3;
    G4double nGasComp = 2.;

    G4Material *air = G4Material::GetMaterial("G4_AIR");
    G4Material *Pb = G4Material::GetMaterial("G4_Pb");
    G4Material *Ar = G4Material::GetMaterial("G4_Ar");
    G4Material *CO2 = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
    G4Material* straw_gas = new G4Material("StrawGas", gas_density, nGasComp);
    straw_gas->AddMaterial(Ar, 70.*perCent);
	  straw_gas->AddMaterial(CO2, 30.*perCent);


    double worldLength = 1000;
    G4Box *worldS = new G4Box("World",                                            // its name
                              worldLength / 2, worldLength / 2, worldLength / 2); // its size

    G4LogicalVolume *worldLV = new G4LogicalVolume(
        worldS,   // its solid
        air,      // its material
        "World"); // its name

    //  Must place the World Physical volume unrotated at (0,0,0).
    //
    G4VPhysicalVolume *worldPV = new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        worldLV,         // its logical volume
        "World",         // its name
        0,               // its mother  volume
        false,           // no boolean operations
        0,               // copy number
        false);          // checking overlaps







    double TestboxLenth = 110;
    double TestboxWidth = 50;

    G4Box *testboxS = new G4Box("testbox", TestboxLenth, TestboxWidth, TestboxWidth);

    G4LogicalVolume *testboxLV = new G4LogicalVolume(testboxS, Pb, "Testbox");

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), testboxLV, "Testbox", worldLV, false, 0, false);

    const int n = 20;
    double TubeRad = 2.5;
    double TubePlaneWidth = TubeRad * (n - 1);

    G4Tubs *tubeS[n];
    G4LogicalVolume *tubeLV[n];
    for (int i = 0; i<n; i++)
      {
        tubeS[i] = new G4Tubs("tube", 0., TubeRad, TestboxWidth, 0. * deg, 360. * deg);
        tubeLV[i] = new G4LogicalVolume (tubeS[i], straw_gas, "tube");
        new G4PVPlacement (0, G4ThreeVector(-15 - TestboxLenth, 2 * TubeRad * i - TubePlaneWidth, 0), tubeLV[i], "tube", worldLV, false, 0, false);
      };

    G4Tubs *tube1S[n];
    G4LogicalVolume *tube1LV[n];
    for (int i = 0; i<n; i++)
      {
        tube1S[i] = new G4Tubs("tube1", 0., TubeRad, TestboxWidth, 0. * deg, 360. * deg);
        tube1LV[i] = new G4LogicalVolume (tube1S[i], straw_gas, "tube1");
        new G4PVPlacement (0, G4ThreeVector(15 +  TestboxLenth, 2 * TubeRad * i - TubePlaneWidth, 0), tube1LV[i], "tube1", worldLV, false, 0, false);
      };

    SetSensitiveDetector("tube", aTrackerSD, true);
    SetSensitiveDetector("tube1", aTrackerSD, true);
  
  //  G4VisAttributes *boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));

    //testboxLV->SetVisAttributes(boxVisAtt);
    return worldPV;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::ConstructSDandField()
  {
    // Sensitive detectors
  }

    

}
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
/// \file TrackerSD.cc
/// \brief Implementation of the B2::TrackerSD class

#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name)
 : G4VSensitiveDetector(name)
{

  G4cout<<"TrackerSD"<<G4endl;

  f = new TFile("straw.root", "RECREATE");
  h2d = new TH2F("h2d", "histogram", 100, -50, 50, 100, -50, 50);
  h2d1 = new TH2F("h2d1", "histogram of 1-st straw", 100, -50, 50, 100, -50, 50);
  h2d2 = new TH2F("h2d2", "histogram of 2-nd straw", 100, -50, 50, 100, -50, 50);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD()
{

 G4cout<<"Artem"<<G4endl;
 h2d -> Write();
 h2d1 -> Write();
 h2d2 -> Write();
 f->Close();



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory*)
{

 G4double edep = aStep->GetTotalEnergyDeposit();

 
G4ThreeVector pos  = aStep->GetPostStepPoint()->GetPosition();
  G4cout<<pos[0]<<"  "<<pos[1]<<" "<<pos[2]<<G4endl;

  if (abs(pos[0]+40)<5)
  {
    h2d1 -> Fill(pos[1], pos[2]);
  }

  if (abs(pos[0]-40)<5)
  {
    h2d2 -> Fill(pos[1], pos[2]);
  }

  h2d->Fill(pos[1],pos[2]);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


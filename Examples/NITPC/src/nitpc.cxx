#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <TApplication.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>

#include "MediumMagboltz.hh"
#include "ComponentElmer.hh"
#include "Sensor.hh"
#include "GarfieldConstants.hh"
#include "Random.hh"
#include "AvalancheNIMicroscopic.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "ViewDrift.hh"
#include "ViewFEMesh.hh"

using namespace Garfield;


int main(int argc, char *argv[]){
  

  TApplication *app = new TApplication("app", &argc, argv);

  double r_val=0.31;

  MediumMagboltz *gas = new MediumMagboltz();
  gas->SetComposition("sf6", 100.);
  gas->SetTemperature(293.15);
  gas->SetPressure(76.);
  //if (argc==4) 
  gas->EnablePenningTransfer(r_val, 0, "sf6");
  gas->SetMaxElectronEnergy(300.);
  gas->EnableDrift();
  gas->Initialise(true);
  gas->LoadIonMobility("/work/shimada/Garfield++/Mobility/IonMobility_SF6-_SF6.txt");

  TFile* file = new TFile("Tracking.root", "recreate");
  TTree* tree = new TTree("tree","tree");
  int nStep,status;
  std::vector<double> xpos,ypos,zpos,time,changePoint;
  tree->Branch("nStep",&nStep);
  tree->Branch("status",&status);
  tree->Branch("xpos",&xpos);
  tree->Branch("ypos",&ypos);
  tree->Branch("zpos",&zpos);
  tree->Branch("time",&time);
  tree->Branch("changePoint",&changePoint);

  std::string data_dir = "/work/shimada/Garfield++/GEM";
  std::string header = data_dir + "/gemcell/mesh.header";
  std::string element = data_dir + "/gemcell/mesh.elements";
  std::string node = data_dir + "/gemcell/mesh.nodes";
  std::string eps = data_dir + "/gemcell/dielectrics.dat";
  std::string volt = data_dir + "/gemcell/gemcell.result";
  ComponentElmer *elm
    = new ComponentElmer(header.c_str(), element.c_str(),
			 node.c_str(), eps.c_str(),
			 volt.c_str(), "cm");
  elm->EnablePeriodicityX();
  elm->EnableMirrorPeriodicityY();
  elm->SetMedium(0, gas);

  Sensor *upic = new Sensor();
  upic->AddComponent(elm);
  upic->SetArea(-0.03, -0.03, -0.03, 0.03, 0.03, 0.03);
  
  ViewDrift *viewDrift = new ViewDrift();
	viewDrift->SetArea(-0.03,-0.03,-0.03,0.03,0.03,0.03);
  
  //view mesh
  TCanvas* c = new TCanvas();
  ViewFEMesh* vMesh = new ViewFEMesh();
  vMesh->SetCanvas(c);
  vMesh->SetComponent(elm);
  vMesh->SetPlane(0,-1,0,0,0,0);
  vMesh->SetFillMesh(true);
  vMesh->SetViewDrift(viewDrift);
  vMesh->SetArea(-.03, -.03, -.03, .03, .03, .03);
	vMesh->SetColor(0, kAzure);
	vMesh->SetColor(1, kOrange);
	vMesh->SetColor(2, kRed);
	vMesh->SetColor(3, kGreen);
  
  AvalancheNIMicroscopic *aval = new AvalancheNIMicroscopic();
  aval->SetSensor(upic);
  aval->SetElectronTransportCut(1e-20);
  aval->SetNegativeIonMass(146);
  //aval->SetDistanceSteps(1e-5); // cm
  aval->SetDistanceSteps(1e-4); // cm
  aval->EnablePlotting(viewDrift);
  aval->InputDetachCrossSectionData("DetachCrossSection.txt");
  aval->SetDetachModel(1);
  
  
  // Counting
  double random_x = 0;
  double random_y = 0;
  double random_z = 0.01;
  //double random_z = -0.004;
  aval->AvalancheNegativeIon(random_x,random_y,random_z,0,0,0,0,0);
  Int_t nd = aval->GetNumberOfElectronEndpoints();
  for(int i=0; i<nd; i++){
    xpos.clear();
    ypos.clear();
    zpos.clear();
    time.clear();
    changePoint.clear();
    int Stat;
    double x0,y0,z0,t0,e0;
    double x1,y1,z1,t1,e1;
    aval->GetElectronEndpoint(i,
                              x0,y0,z0,t0,e0,
                              x1,y1,z1,t1,e1,
                              Stat);
    unsigned int nstep = aval->GetNumberOfElectronDriftLinePoints(i);
    status = Stat;
    nStep  = nstep;
    for(int step=0;step<nstep;step++){
      double x,y,z,t;
      int change;
      aval->GetElectronDriftLinePoint(x,y,z,t,change,step,i);
      xpos.push_back(x);
      ypos.push_back(y);
      zpos.push_back(z);
      time.push_back(t);
      changePoint.push_back(change);
    }
    tree->Fill();
  }

  tree->Write();
  file->Close();

  std::cout << "... End Negative Ion Gain Simulation" <<std::endl;

  std::cout << "Plotting Start ... " <<std::endl;
  vMesh->EnableAxes();
  vMesh->Plot();

  c->Update();
  c->Print("track.png");


}

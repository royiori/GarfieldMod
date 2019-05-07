//------------------------------------------------------------------//
// Author: Qian LIU
//------------------------------------------------------------------//


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "ComponentComsol.hh"
#include "ViewField.hh"
#include "ViewFEMesh.hh"
#include "ViewGeometry.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewSignal.hh"

#include <time.h>

using namespace Garfield;
using namespace std;

void ReadFiledMap(TString);
void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2);
void PlotField(ComponentComsol*);
void PlotDriftLine(AvalancheMC* driftIon);
void PlotDriftSetup(AvalancheMicroscopic* aval, AvalancheMC* driftIon);
void PlotDrift();
void PlotMesh();
void PlotSignal(Sensor* sensor);

ComponentComsol *fm;
MediumMagboltz  *gas;
ViewDrift *driftView;

// Dimensions of the GEM, in !!!! CM !!!!
const double radus  = 0.005;
const double pitch  = 4*radus;
const double diam = 2*radus;
const double thick = 0.02;
const double metal  = 0.001;
const double induce = 0.2;
const double driftz  = 0.2;

// Parameters for simulation
const int    nSteps = 100;
const double tStart = 0;
const double tStop  = 100;
const double tStep  = (tStop - tStart) /nSteps;
const string label = "readout";


// Dimensions for area selection for DriftView & Sensor etc.
const double xmin = -pitch*3/2;
const double xmax =  pitch*3/2;
const double ymin = -pitch*3/2;
const double ymax =  pitch*3/2;
const double zmin =  -induce-thick/2-metal; //anode
const double zmax =   driftz+thick/2+metal; //cathode

// Control flag
const bool debug         = 1;
const bool plotField     = 0;
const bool plotDriftLine = 0;
const bool plotDrift     = 0;
const bool plotMesh      = 0;
const bool plotSignal    = 0;
const bool plotHistogram = 1;

const double zprim  = 0.1;//thick/2;//thick/2+2*metal+driftz/2;

char Path[100] = "./simdata/thgem/";
TString RootFile = TString(Path)+TString("/result.root");
TString TrackFile= TString(Path)+TString("/track.root");


//------------------------------------------------------------------//
// Main 
//------------------------------------------------------------------//
int main(int argc, char * argv[]) {

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    ReadFiledMap(TString(Path));
    DefineGas("ar", 80, "DME", 20);

    if(plotField) { PlotField(fm); app.Run(kTRUE); }

    //------------------------------------------------------
    // Create the sensor. Initialize
    fm->SetWeightingFieldLabel(0, label);
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
    //sensor->EnableDebugging();

    if(plotSignal) {
        sensor->WriteSensorFile(TrackFile.Data());
        sensor->AddElectrode(fm, label);
        sensor->SetTimeWindow(-1, tStep, 200);
        sensor->ClearSignal();
    }

    //------------------------------------------------------
    // Now avalanche and track electron
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    aval->EnableSignalCalculation();
    //aval->EnableDebugging();

    // Now track ions
    AvalancheMC* driftIon = new AvalancheMC();
    driftIon->SetSensor(sensor);
    driftIon->SetDistanceSteps(2.e-4);
    //drift->EnableDebugging();

    //------------------------------------------------------
    // plot avalanch drift line
    if(plotDriftLine) { PlotDriftLine(driftIon); app.Run(kTRUE); }

    if(plotDrift) { PlotDriftSetup(aval, driftIon); }

    //------------------------------------------------------
    // define root file
    int ne, ni;
    vector<double> e1hit;
    vector<double> e2hit;
    vector<double> i1hit;
    vector<double> i2hit;
  
    TFile *f = TFile::Open(RootFile,"RECREATE");
    TTree *t = new TTree("tree","trees");
    t->Branch("e1hit",&e1hit);
    t->Branch("e2hit",&e2hit);
    t->Branch("i1hit",&i1hit);
    t->Branch("i2hit",&i2hit);
    t->Branch("ne", &ne, "ne/I");
    t->Branch("ni", &ni, "ni/I");


    //------------------------------------------------------
    // define Histograms according to your needs
    TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", 100, 0, 1000);
    TH1F* hIons      = new TH1F("hIons",      "Number of ions",      100, 0, 1000);

    //------------------------------------------------------
    // define parameters
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;

    const int nEvents   = 10;

    //------------------------------------------------------
    // start simulation
    clock_t start, stop;
    start = clock();
 
    for (int i = nEvents; i--;) {
        if (debug && i%100==0) std::cout << "-->Simulating "<< i << " out of " << nEvents << " events.\n";

        double smear = diam /2 *3;
        double x0 = - smear + 2 * RndmUniform() * smear;
        double y0 = 0;
        double z0 = zprim;
        double t0 = 0.;
        double e0 = 0.1;
        
        sensor->Reset();
        aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);    // Calculate an electron avalanche
        aval->GetAvalancheSize(ne, ni);
        hElectrons->Fill(ne);
        hIons     ->Fill(ni);

        e1hit.clear();
        e2hit.clear();
        i1hit.clear();
        i2hit.clear();

        const int np = aval->GetNumberOfElectronEndpoints();
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
                                      xe2, ye2, ze2, te2, e2, status);

            e1hit.push_back(xe1);
            e1hit.push_back(ye1);
            e1hit.push_back(ze1);
            e1hit.push_back(te1);
            e1hit.push_back( e1);

            e2hit.push_back(xe2);
            e2hit.push_back(ye2);
            e2hit.push_back(ze2);
            e2hit.push_back(te2);
            e2hit.push_back( e2);

            driftIon->DriftIon(xe1, ye1, ze1, te1);
            driftIon->GetIonEndpoint(0, xi1, yi1, zi1, ti1,
                                  xi2, yi2, zi2, ti2, status);
 
            i1hit.push_back(xi1);
            i1hit.push_back(yi1);
            i1hit.push_back(zi1);
            i1hit.push_back(ti1);
                      
            i2hit.push_back(xi2);
            i2hit.push_back(yi2);
            i2hit.push_back(zi2);
            i2hit.push_back(ti2);
        }
        t->Fill();
        sensor->Fill();

	if(plotSignal) { PlotSignal(sensor); app.Run(kTRUE); }
    }

    f->Write();
    sensor->Write();

    stop = clock();
    cout<<" Time consuming: "<<(double)(stop - start)/CLOCKS_PER_SEC<<endl;


    //------------------------------------------------------
    // plotting Histograms according to your needs

    if(plotDrift) {PlotDrift();}

    if (plotHistogram) {
        TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
        cH->Divide(2, 2);
        cH->cd(1);
        hElectrons->Draw();
        cH->cd(2);
        hIons->Draw();
        cH->Modified();
        cH->Update();
    }

    app.Run(kTRUE);

}

//------------------------------------------------------------------//
// Read COMSOL output
//------------------------------------------------------------------//
void ReadFiledMap(TString path)
{   
    // Loading the Comsol files
    fm = new ComponentComsol();
    TString mfile = path+TString("/mesh_THGEM.mphtxt");
    TString dfile = path+TString("/dielectrics.dat");
    TString ffile = path+TString("/Field_THGEM.txt");

    fm->Initialise(mfile.Data(), dfile.Data(), ffile.Data());
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
    fm->SetMagneticField(0.,0.,0.);
    fm->PrintRange();
}


//------------------------------------------------------------------//
// define gas compounents
//------------------------------------------------------------------//
void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2)
{
    gas = new MediumMagboltz();
    gas->SetComposition(gas1, rat1, gas2, rat2);
    gas->SetTemperature(293.15);
    gas->SetPressure(760*0.8);
    gas->Initialise();
    gas->DisableDebugging();
    //gas->EnableDebugging();

    gas->SetMaxElectronEnergy(200.); //eV
    gas->SetMaxPhotonEnergy(200.);
    gas->EnableEnergyRangeAdjustment(true);
    gas->EnableAnisotropicScattering();

    gas->Initialise(true);

    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");

    const std::string path = getenv("GARFIELD_HOME");
    gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");

    // Associate the gas with the corresponding field map material.
    const int nMaterials = fm->GetNumberOfMaterials();
    for (int i = 0; i < nMaterials; ++i) {
        const double eps = fm->GetPermittivity(i);         // Get permittivity (episilon)
        if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);    // set gas as medium if eps is arround 1
    }
    fm->PrintMaterials();
}


//------------------------------------------------------------------//
// 画电场，电力线等图
//------------------------------------------------------------------//
void PlotField(ComponentComsol* fm0)
{
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm0);

    fieldView->SetNumberOfSamples1d(2000);           //default is 1000, the density of the plotting grid
    fieldView->SetNumberOfSamples2d(100,100);        //default is 200, 200. 二维图的nX, nY
    fieldView->SetNumberOfContours(50);              //default is 100, 等高线的次数
    //fieldView->SetVoltageRange(-3000., 3000.);       //the min/max of potential in the Component is used as default, but it can also be defined by user.
    //fieldView->SetElectricFieldRange(0,20000);       //default is (0, 100). it must be defined by user for electric field.
    //fieldView->SetWeightingFieldRange(0,10);         //it must be defined by user for weighting field.

    fieldView->SetPlane(0, 0, 0, 0, 0., thick/2);        //(fx,fy,fz,x0,y0,z0), (fx,fy,fz) 指定的是平面的法线方向。
    fieldView->SetArea(xmin, ymin, xmax, ymax);      //指定这个平面所画的区域大小
    TCanvas* cF1 = new TCanvas();
    fieldView->SetCanvas(cF1);
    fieldView->PlotContour("e");                     // v/p/phi/volt/voltage/pot/potential 都是画电势, e/field 画电场, or ex/ey/ez


    TCanvas* cF2 = new TCanvas();
    fieldView->SetCanvas(cF2);
    fieldView->PlotProfile(0, 0, -induce, 0, 0, thick+2*metal+driftz, "v");

    //double x=0, y=0, z=0;
    //double ex, ey, ez, v;
    //Medium* m;
    //int status;
    //for(int iz=0;iz<10;iz++) {
    //  z = -induce + iz * (thick+2*metal+driftz+induce)/10000;
    //  fm0->ElectricField(x, y, z, ex, ey, ez, v, m, status);
    //  //fm0->EnableDebugging();
    //  //fm0->EnableCheckMapIndices();
    //  cout << "Field at ("<<x<<","<<y<<","<<z<< "): "<<ex<<" "<<ey<<" "<<ez<<" "<<v<<" "<<status<< "\n\n\n";
    //  //fm0->DisableDebugging();
    //  //fm0->DisableCheckMapIndices();
    //}

    return;
}

//------------------------------------------------------------------//
// Plot Mesh
//------------------------------------------------------------------//
void PlotMesh() 
{
    ViewFEMesh* meshView = new ViewFEMesh();
    meshView->SetComponent(fm);
    meshView->SetPlane(0., 1., 0., 0., 0., 0.);
    meshView->SetArea(xmin, zmin, zmin, xmax, zmax, zmax);
    meshView->SetFillColor(0,6); //matid=0 is Cu  
    meshView->SetFillColor(1,7); //matid=1 is FR4, 1=black, 2=Red, 3=Green, 4=blue, 5=yellow
    meshView->SetFillMesh(true);
    meshView->SetViewDrift(driftView);
    TCanvas* cM = new TCanvas();
    meshView->SetCanvas(cM);
    meshView->Plot();
}


//------------------------------------------------------------------//
// Plot drift lines
//------------------------------------------------------------------//
void PlotDriftLine(AvalancheMC* driftIon)
{
    driftView = new ViewDrift();
    driftView->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
    driftIon->EnablePlotting(driftView);
    driftIon->DisableDiffusion();

    const int Nbin=50;
    for(int i=0;i<Nbin;i++) {
        double xtmp = xmin+(i+0.5)*(xmax-xmin)/Nbin;
        driftIon->DriftIon(xtmp, 0, zmin,0);
    }

    TCanvas* cD = new TCanvas();
    driftView->SetCanvas(cD);
    driftView->Plot();

    if (plotMesh) PlotMesh();
}


//------------------------------------------------------------------//
// Plot drift lines of avalanche
//------------------------------------------------------------------//
void PlotDriftSetup(AvalancheMicroscopic *aval, AvalancheMC* driftIon)
{
    driftView = new ViewDrift();
    driftView->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);

    // Plot every 10 collisions (in microscopic tracking).
    aval->SetCollisionSteps(10);
    aval->EnablePlotting(driftView);
    driftIon->EnablePlotting(driftView);
}

void PlotDrift()
{
    TCanvas* cD = new TCanvas();
    driftView->SetCanvas(cD);
    driftView->Plot();
        
    if (plotMesh) PlotMesh();
}

//------------------------------------------------------------------//
// Plot signal
//------------------------------------------------------------------//
void PlotSignal(Sensor* sensor)
{
    cout<<"------------------------------------------------"<<endl;
    cout<<"->                                            <-"<<endl;
    cout<<"-> Plotting signal will only simulate 1 event <-"<<endl;
    cout<<"->                                            <-"<<endl;
    cout<<"------------------------------------------------"<<endl;

    ViewSignal* signalView = new ViewSignal();
    signalView->SetSensor(sensor);
    TCanvas* c1 = new TCanvas();
    signalView->SetCanvas(c1);
    signalView->PlotSignal(label); 
    signalView->PlotSignal(label,false,false,true);
    sensor->ClearSignal();  
}

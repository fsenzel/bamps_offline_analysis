//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h> // for getenv()
#include <stdlib.h> // for getenv()
#include <string>
#include <sstream>
#include "configuration.h"
#include "scattering22.h"
#include "particle.h"
// #include "interpolation_gQ.h"


/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<Particle> ns_casc::particles;
std::vector<Particle> ns_casc::particles_init;
std::vector<Particle> ns_casc::particles_atTimeNow;
std::vector<Particle> ns_casc::particles_atTimeNowCopy;
std::vector<Particle> ns_casc::addedParticles;
std::vector<Particle> ns_casc::addedParticlesCopy;


using namespace std;

// global variables which are changed in config.cpp (input parameter)
double Mcharm = 1.5;       // GeV, mass of charm quark
double Mbottom = 4.8;      // GeV, mass of bottom quark
double Kggccb = 1.0;       // K factor for gg -> ccb cross section
double KgQgQ = 1.0;        // K factor for gQ -> gQ cross section
double kappa_gQgQ = 1.0;   // kappa factor for debye mass, usually 0.2 (Peshier,Gossiaux,Aichelin)
bool isCouplingRunningForGQ = false; // true if a running coupling is used for gQ -> gQ 

double dt = 0.1;

bool isIsotropicCrossSec = false;
bool isConstantCrossSec = false;
double constantCrossSec = 10.0; // mb

// interpolation_gQ theIgc;
// interpolation_gQ theIgb;


config::config(const int argc, const char * const argv[])
{
  //--set default values---------------------------------------
  pathdirCascadeData = "data_cascade";
  name = "default";
  runtime = 5.1;
  P0 = 2.9;
  delta_t = 0.0;
  dt_specified = false;
  originalSeed = 0;
  outputLevel = 5;
  numberOfParticlesToAdd = 0;
  minimumPT = 0;
  freezeOutEnergyDensity = 0.6;
  factor_dt = 0.8;
  centralRingRadius = 1.5;
  ringNumber = 10;
  deltaR = 1.0;
  pythiaParticleFile = "-";
  localCluster = true;
  movieOutputBackground = false;
  movieOutputJets = false;
//   pythiaParticleFileCharmTestparticles = 35;
//   KIniCharm = 1.0;
  //these values are used in case no input file is specified 
  //or in case the input file lacks certain statements
  //----------------------------------------------------------- 
//   pathdirData = "data";
  
  switch (argc)
    {
    case 1: 
      cout << "WARNING: Argument specifying input file is missing. Defaults will be used." << endl;
      break;
    case 2:
      inputFile = argv[1];
      readInputfile();
      cout << "Reading input file: " << inputFile << endl;
      break;
    case 3:
      inputFile = argv[1];
      readInputfile();
      cout << "Reading input file: " << inputFile << endl;
      break;
    default:
      cout << "WARNING: Too many arguments passed. Defaults will be used." << endl;
      break;
    }  

  cout << "Cascade data folder: " << pathdirCascadeData << endl;
  
  if(dt_specified)
    cout << "dt = " << dt << endl;
  else
    cout << "Get dt for each step from cascade data." << endl;
    
//   if(isConstantCrossSec)
//     cout << "Constant cross section of " << constantCrossSec << " mb for gQ-> gQ." << endl;
//   else
//   {
//     cout << "pQCD cross section for gQ-> gQ with K = " << KgQgQ;
//     if(isCouplingRunningForGQ)
//       cout << " and running coupling";
//     else
//       cout << " and constant coupling";
//     cout << ", kappa_gQ = " << kappa_gQgQ << endl;
//   }
// 
//   if(isIsotropicCrossSec)
//     cout << "Sample angle isotropic." << endl;
// 
//   theIgc.configure(charm, isCouplingRunningForGQ);
//   theIgb.configure(bottom, isCouplingRunningForGQ);
}


 
int config::getNameLength() const
{
  return name.length();
}



void config::readInputfile()
{  
  extern double Mcharm;       // GeV, mass of charm quark
  extern double Mbottom;       // GeV, mass of charm quark
  extern double Kggccb;       // K factor for gg -> ccb cross section
  extern double KgQgQ;        // K factor for gc -> gc cross section
  extern double kappa_gQgQ;   // kappa factor for debye mass, usually 0.2 (Peshier,Gossiaux,Aichelin)
  extern bool isCouplingRunningForGQ; // true if a running coupling is used for gQ -> gQ 
  extern double dt;
  
  fstream file(inputFile.c_str(), ios::in);
  string content = "";
  
  do
    {
      file >> content;

      if (content == "#name")
      {
        file >> name;
      }
      if (content == "#outputName")
      {
        file >> outputName;
      }
      else if (content == "#time")
      {
        file >> runtime;
      }
      else if (content == "#P0")
      {
        file >> P0;
      }
      else if ( content == "#jetMovie")
      {
        int temp;
        file >> temp;
        if ( temp != 0 )
        {
          movieOutputJets = true;
        }
      }
      else if ( content == "#backgroundMovie")
      {
        int temp;
        file >> temp;
        if ( temp != 0 )
        {
          movieOutputBackground = true;
        }
      }
      else if (content == "#seed")        
      {
        file >> seed;
      }
      else if (content == "#outLevel")
      {
        file >> outputLevel; 
      }
      else if (content == "#localCluster")
      {
        file >> localCluster;
      }
      else if (content == "#dt")
      {
        file >> delta_t;
        dt = delta_t;
        dt_specified = true;
      }
      else if (content == "#pathdirCascadeData")
        file >> pathdirCascadeData;  
      else if (content == "#isotropicCrossSec")
        file >> isIsotropicCrossSec;
      else if (content == "#constantCrossSec")
      {
        file >> constantCrossSec;
        isConstantCrossSec = true;
      }
      else if ( content == "#nAdded" )
      {
        file >> numberOfParticlesToAdd;
      }
      else if ( content == "#factor_dt" )
      {
        file >> factor_dt;
      }
      else if ( content == "#pythiaParticleFile" )
      {
        file >> pythiaParticleFile;
      }
      else if ( content == "#minPT_added" )
      {
        file >> minimumPT;
      }
//       else if (content == "#pythiaParticleFileCharm")
//         file >> pythiaParticleFileCharm;
//       else if (content == "#pythiaParticleFileCharmTestparticles")
//         file >> pythiaParticleFileCharmTestparticles;

//       else if (content == "#Mcharm")
//         file >> Mcharm;
//       else if (content == "#Mbottom")
//         file >> Mbottom;
//       else if (content == "#Kggccb")
//         file >> Kggccb;
//       else if (content == "#Kgcgc")
//         file >> KgQgQ;
//       else if (content == "#KIniCharm")
//         file >> KIniCharm;   

//       else if (content == "#kappa_gQgQ")
//         file >> kappa_gQgQ;  
//       else if (content == "#Kgcgc")
//         file >> kappa_gQgQ;  
//       else if (content == "#asRunGQ")
//         file >> isCouplingRunningForGQ;   
    }
  while(content != "#end");
}



int config::getJobID()
{
  int lastJobDigits;
  
  char * temp;
  
  if(localCluster)
  {
    temp = getenv("JOB_ID");
    if (temp != NULL)
    {
      string jobID = temp;
      stringstream ss(&jobID[jobID.size()-2]);
      ss >> lastJobDigits;
    }
    else
      lastJobDigits = 0;
  }  
  else
  {
    temp = getenv("PBS_JOBID");
    if (temp != NULL)
    {
      string jobID = temp;
      string::size_type pos = jobID.find(".",0);
      jobID.erase(pos);
      stringstream ss(&jobID[jobID.size()-2]);
      ss >> lastJobDigits;
    }
    else
      lastJobDigits = 0;
  }

  cout << "last job digits: " << lastJobDigits << endl;

  return lastJobDigits;
}

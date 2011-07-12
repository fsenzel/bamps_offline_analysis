//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_box/main/branches/quark_branch/src/binning.cpp $
//$LastChangedDate: 2008-02-22 15:28:55 +0100 (Fri, 22 Feb 2008) $
//$LastChangedRevision: 21 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <fstream>
#include <iostream>

#include "binning.h"

using namespace std;

/*
  default constructor
*/
binning::binning()
{
  if (!theList.empty())
    theList.clear();
  if (!theBins.empty())
    theBins.clear();

  nBins = 100;              //set the default number of bins
  theBins.resize(nBins,0);  //resize vector "theBins" to nBins and fill with 0
  btype = default_binning;  //default binning type
  binnedEvents = 0;         //number of collected events = 0  
  finished = false;  
  setFilename();            //set filename to default
}

/*
  default constructor + filename set
*/
binning::binning(const string name)
{
  if (!theList.empty())
    theList.clear();
  if (!theBins.empty())
    theBins.clear();

  nBins = 100;               //set the default number of bins
  theBins.resize(nBins,0);   //resize vector "theBins" to nBins and fill with 0
  btype = default_binning;   //default binning type
  binnedEvents = 0;          //number of collected events = 0   
  finished = false; 
  setFilename(name);         //set filename
}


/*
  constructor for specifying the binning area and the number of bins
*/
binning::binning(const double min_arg, const double max_arg, const int n)
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize(nBins,0);
  binWidth = (x_max - x_min) / nBins;
  btype = fixed_min_max_number;
  binnedEvents = 0;
  finished = true;  
  setFilename();
}

/*
  constructor for specifying the binning area and the number of bins + filename set
*/
binning::binning(const string name, const double min_arg, const double max_arg, const int n)
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize(nBins,0);
  binWidth = (x_max - x_min) / nBins;
  btype = fixed_min_max_number;
  binnedEvents = 0; 
  finished = true; 
  setFilename(name);
}


/*
  constructor for specifying the number of bins
*/
binning::binning(const int n)
{
  nBins = n;
  theBins.resize(nBins,0);
  btype = fixed_number;
  binnedEvents = 0;
  finished = false;
  setFilename();
}

/*
  constructor for specifying the number of bins + filename set
*/
binning::binning(const string name, const int n)
{
  nBins = n;
  theBins.resize(nBins,0);
  btype = fixed_number;
  binnedEvents = 0;
  finished = false;
  setFilename(name);
}


/*
  constructor for specifying the bin width
*/
binning::binning(const double width)
{
  binWidth = width;
  btype = fixed_width;
  binnedEvents = 0;
  finished = false;
  setFilename();
}

/*
  constructor for specifying the bin width + filename set
*/
binning::binning(const string name, const double width)
{
  binWidth = width;
  btype = fixed_width;
  binnedEvents = 0;
  finished = false;
  setFilename(name);
}


binning::~binning()
{
}


void binning::setMinMaxWidth(const double min_arg, const double max_arg, const double width_arg)
{  
  x_min = min_arg;
  x_max = max_arg;
  binWidth = width_arg;
  nBins = int((x_max-x_min)/binWidth);
  if( (x_min + nBins*binWidth) < x_max )
    nBins++;
  theBins.resize(nBins,0);
  btype = fixed_min_max_number;
  binnedEvents = 0; 
  finished = true; 
  
}


void binning::setMinMaxN(const double min_arg, const double max_arg, const int n)
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize(nBins,0);
  binWidth = (x_max - x_min) / nBins;
  btype = fixed_min_max_number;
  binnedEvents = 0; 
  finished = true; 
}



//set the number of bins
void binning::setNumber(const int n)
{
  nBins = n;
  btype = fixed_number; 
}

//set the bin width
void binning::setWidth(const double w)
{
  binWidth = w;
  btype = fixed_width;
}


//add element x to the binned data
void binning::add(const double x)
{
  //if min and max are known AND data collection hasn't started yet, x can be sorted into the appropriate bin immediately
  if (btype == fixed_min_max_number && theList.size() == 0)   
  {
    int index = int((x - x_min)/binWidth);
    if (index >= nBins)
      index = nBins-1;
    else if (index < 0)
      index = 0;
    ++theBins[index];
  }
  //otherwise all data has to be collected before sorting into bins starts
  else
    theList.push_back(x);   //push x into the temporary list
  
  ++binnedEvents;
}


//sort into bins if necessary and print results to file
void binning::print()
{
  if (!finished)
    //binning still needs to be done
    //in the case where btype==fixed_min_max_number AND theList.size==0 the binning was done directly in binning::add(..)  
  {
    finish();
  } 
   
  fstream file(filename.c_str(),ios::out);
  file << "#bin width: " << binWidth << "   number of binned events: " << binnedEvents << endl;
  for (int i=0; i<nBins; i++)
  {
    file << x_min + (i+0.5)*binWidth << "\t" << theBins[i] << "\t" << theBins[i]/binWidth << endl;
  }   
  file.close();
}


void binning::finish()
{
  if (btype != fixed_min_max_number)
  {
    x_min = x_max = theList[0];
    for(int i=1; i<theList.size(); i++)
    {
      if(theList[i] < x_min)
        x_min = theList[i];
      if(theList[i] > x_max)
        x_max = theList[i];
    }
  }  
  
  if (btype == fixed_width)
  {
    nBins = int((x_max - x_min) / binWidth) + 1;
  }
  else
    binWidth = (x_max - x_min) / nBins; 

  
  int index;
  for(int i=0; i<theList.size(); i++)
  {
    index = int((theList[i] - x_min)/binWidth);
    if (index >= nBins)
      index = nBins-1;
    ++theBins[index];
  }
  
  finished = true;
}



double binning::getBinLabel(const int i)
{
  if (!finished)
  {
    finish();
  } 

  if (i >= 0 && i < nBins)  
    return x_min + (i+0.5)*binWidth;
  else
    return 0;
}


double binning::getBin(const int i)
{
  if (!finished)
  {
    finish();
  } 

  if (i >= 0 && i < nBins)  
    return theBins[i]/binWidth;
  else
    return 0;
}



int binning::getBinRaw(const int i)
{
  if (!finished)
  {
    finish();
  } 

  if (i >= 0 && i < nBins)  
    return theBins[i];
  else
    return 0;
}

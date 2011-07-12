#include <fstream>
#include <iostream>

#include "binning2.h"

using namespace std;




// 2 dimensional binning
binning2d::binning2d(const string name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y)
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  nBins_x = n_x;
  y_min = min_y_arg;
  y_max = max_y_arg;
  nBins_y = n_y;
  
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ ) 
  {
    theBins.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
      theBins[i].push_back ( 0 );
  }
  
  binWidth_x = (x_max - x_min) / nBins_x;
  binWidth_y = (y_max - y_min) / nBins_y;
  btype = fixed_min_max_number;
  binnedEvents = 0; 
  finished = true; 
  setFilename(name);
}


//add element x to the binned data
void binning2d::add(const double x, const double y)
{
  
  int index_x = int((x - x_min)/binWidth_x);
  if (index_x >= nBins_x)
    index_x = nBins_x-1;
  else if (index_x < 0)
    index_x = 0;
  int index_y = int((y - y_min)/binWidth_y);
  if (index_y >= nBins_y)
    index_y = nBins_y-1;
  else if (index_y < 0)
    index_y = 0;
  ++theBins[index_x][index_y];
  
  ++binnedEvents;
}


//sort into bins if necessary and print results to file
void binning2d::print()
{
   
  fstream file(filename.c_str(),ios::out);
  file << "#bin width_x: " << binWidth_x << "   bin width_y: " << binWidth_y <<  "   number of binned events: " << binnedEvents << endl;
  for (int i=0; i<nBins_x; i++)
  {
    for (int j=0; j<nBins_y; j++)
    {
      file.width(15);
      file << x_min + (i+0.5)*binWidth_x;
      file.width(15);
      file << y_min + (j+0.5)*binWidth_y;
      file.width(15);
      file << theBins[i][j];
      file.width(15);
      file << double(theBins[i][j])/binnedEvents/binWidth_x/binWidth_y << endl;
    }
    file << endl;
  }   
  file.close();
}






// do not add 1, but values to bin and average over them
binningValues::binningValues(const string name, const double min_arg, const double max_arg, const int n)
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize(nBins,0.0);
  NumberInBins.resize(nBins,0);
  binWidth = (x_max - x_min) / nBins;
  btype = fixed_min_max_number;
  binnedEvents = 0; 
  finished = true; 
  setFilename(name);
}

//add element x to the binned data
void binningValues::add(const double x, const double v)
{

  int index = int((x - x_min)/binWidth);
  if (index >= nBins)
    index = nBins-1;
  else if (index < 0)
    index = 0;
    
  theBins[index] += v;
  ++NumberInBins[index];

  ++binnedEvents;
}


//sort into bins if necessary and print results to file
void binningValues::print()
{
  fstream file(filename.c_str(),ios::out);
  file << "#bin width: " << binWidth << "   number of binned events: " << binnedEvents << endl;
  for (int i=0; i<nBins; i++)
  {
    file.width(15);
    file << x_min + (i+0.5)*binWidth;
    file.width(15);
    if(NumberInBins[i] != 0)
      file << theBins[i];
    else
      file << 0;
    file.width(15);
    file << NumberInBins[i] << endl;
  }   
  file.close();
}







// 2 dimensional binning
binningValues2d::binningValues2d(const string name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y)
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  nBins_x = n_x;
  y_min = min_y_arg;
  y_max = max_y_arg;
  nBins_y = n_y;
  
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ ) 
  {
    theBins.push_back ( vector<double>() );
    NumberInBins.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0.0 );
      NumberInBins[i].push_back ( 0 );
    }
  }
  
  binWidth_x = (x_max - x_min) / nBins_x;
  binWidth_y = (y_max - y_min) / nBins_y;
  btype = fixed_min_max_number;
  binnedEvents = 0; 
  finished = true; 
  setFilename(name);
}


//add element x to the binned data
void binningValues2d::add(const double x, const double y, const double v)
{
  
  int index_x = int((x - x_min)/binWidth_x);
  if (index_x >= nBins_x)
    index_x = nBins_x-1;
  else if (index_x < 0)
    index_x = 0;
  int index_y = int((y - y_min)/binWidth_y);
  if (index_y >= nBins_y)
    index_y = nBins_y-1;
  else if (index_y < 0)
    index_y = 0;
    
  theBins[index_x][index_y] += v;
  ++NumberInBins[index_x][index_y];
  
  ++binnedEvents;
}


//sort into bins if necessary and print results to file
void binningValues2d::print()
{
   
  fstream file(filename.c_str(),ios::out);
  file << "#bin width_x: " << binWidth_x << "   bin width_y: " << binWidth_y <<  "   number of binned events: " << binnedEvents << endl;
  for (int i=0; i<nBins_x; i++)
  {
    for (int j=0; j<nBins_y; j++)
    {
      file.width(15);
      file << x_min + (i+0.5)*binWidth_x;
      file.width(15);
      file << y_min + (j+0.5)*binWidth_y;
      file.width(15);
      if(NumberInBins[i][j] != 0)
        file << theBins[i][j];
      else
        file << 0;
      file.width(15);
      file << NumberInBins[i][j] << endl;
    }
    file << endl;
  }   
  file.close();
}





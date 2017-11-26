// A simple program to generate exponential random numbers and
// store them in a histogram; also optionally writes the individual
// values to a file.

// Glen Cowan
// RHUL Physics
// 2 December 2006

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>

using namespace std;

int main(int argc, char **argv) {

// Set up output file

  ofstream dataFile;
  dataFile.open("expData.txt");
  if ( dataFile.fail() ) {
    cout << "Could not open output file -- exit program" << endl;
    exit(1);
  }

// Create a TRandom object to generate random numbers uniform in ]0,1]
// Use the "Marsenne Twister" algorithm of TRandom3

  int seed = 12345;
  TRandom* ran = new TRandom3(seed);
 
//  Fill with exponential random numbers.

  const double xi1 = 1.0;              // mean value of the exponential
  const double xi2 = 5.0;              // value of xi2 for exponential
  const double alpha = 0.2;            // coeff
  double x;
  int numVal = 0;
  cout << "Enter number of values to generate: ";
  cin >> numVal;

  for (int i = 0; i<numVal; ++i){
    double r1 = ran->Rndm();
    double r2 = ran->Rndm();
    if (r1 < alpha) {
        x = -xi1 * log(r2);
        cout << "x1 : " << x << endl;
    }
    else {
        x = -xi2 * log(r2);
        cout << "x2 : " << x << endl;
    }

    dataFile << x << endl;
  }

// Save all histograms and close up.

  dataFile.close();

  return 0;

}

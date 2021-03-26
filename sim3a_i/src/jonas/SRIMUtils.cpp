#include "jonas/SRIMUtils.h"
#include "jonas/LineShapes.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

TH1I * readRangeFile(string fileName)
{
  vector<double> ranges;
  ifstream infile;

  /*
  Ok, so we try and open the input file. If not successful, we tell the user and return.
  Otherwise we proceed by reading the file line by line, until we find some lines with useful info.
  */
  infile.open(fileName);
  if(!infile.is_open()){
    cout<<"Could not open file\n";
    return 0;
  }

  int row = 0;
  double maxRange = 0;
  while(!infile.eof())
  {
    row++;
    string line;
    getline(infile, line);
    if(row < 18) continue;      //Not interested in file header (18 lines).
    stringstream ss(line);

    double ion, x, y, z;
    ss >> ion >> x >> y >> z;   //x,y,z in Ångstrøm.
    ranges.push_back(x);
    if(x > maxRange) maxRange = x;
  }
  infile.close();

  TH1I *h = new TH1I("h","",100,0,maxRange);
  for(int i=0; i<ranges.size(); i++){
    double range = ranges.at(i);
    h->Fill(range);
  }

  return h;    
}

unique_ptr<TF1> getRangeProfile(string fileName)
{
  TH1I *h = readRangeFile(fileName);
  if(!h) return 0;
  double max = h->GetXaxis()->GetXmax();
  unique_ptr<TF1> f(new TF1("range",SimpleAlpha,0,max,4));

  double integral = h->Integral();
  double center =  h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  double sigma = h->GetRMS();
  double tail = - h->GetSkewness() * 1.5 * sigma;
  
  f->SetParameters(5*integral,center,sigma,tail); //These guesses are a bit ad hoc and might not be entirely robust.
  TFitResultPtr r = h->Fit(f.get(),"Q0R"); //Q0R
  int status = r;
  if(status != 0) cout << "Fit to SRIM range file failed!" << endl;

  return f;
}

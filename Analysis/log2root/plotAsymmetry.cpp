// ================================== header ===================================
// function plotaAsymmetry
// -------------------
// Author: Declan Millar
// ================================== start ====================================
TH1D* plotAsymmetry(ifstream * logStream){

  // declarations.
  std::string asymName, xTitle,  yTitle; // asymribution labels
  // std::string xVar, xUnits, yVar, yUnits 
  std::string xS, yS, sigpS, sigmS, stop; // for reading in data pairs as strings
  int n_DataPairs; // number of x,y pairs  
  double x, y, sigp, sigm; // for conversion of data pairs to integers
  std::vector<double> xV,yV,sigpV,sigmV; // for storing data pairs
  int nBins, nBinEdges;
  double binWidth;
  std::vector<double> binLowEdges;

  // name strings.
  *logStream >> asymName >> yTitle >> xTitle;
  // error when no units in title
  //*logStream >> asymName >> yVar >> yUnits >> xVar >> xUnits;
  //xTitle = xVar + " " + xUnits;
  //yTitle = yVar + " " + yUnits;
  stop = "END"; // when found, stop reading the histograms

  // loop over data pairs.
  n_DataPairs=0;
  if (logStream->is_open()){
    while(1){        
      *logStream >> xS;      
      // if (!logStream->good()) break; // breaks loop when we stop finding data
      if (xS == stop) // breaks the loop if reads stop
      {
        // printf("Stop string: %s\n",xS.c_str());
        break;
      }
      *logStream >> yS;
      *logStream >> sigpS;
      *logStream >> sigmS;
      x = atof(xS.c_str());
      y = atof(yS.c_str());
      sigp = atof(sigpS.c_str());
      sigm = atof(sigmS.c_str());
      n_DataPairs+=1;
      // printf ("Data pairs number: %i\n",n_DataPairs);
      // std::cout << x << ' ' << y << std::endl;
      xV.push_back(x);
      yV.push_back(y);
      sigpV.push_back(sigp);
      sigmV.push_back(sigm);
    }
  }
  else printf("  Stream failed!");
  // printf ("  Read in %i data pairs.\n",n_DataPairs);

  // histogram bin information.
  nBins = xV.size();
  nBinEdges = nBins +1;
  binWidth = xV[1]-xV[0]; // only works for fixed bin width!
  binLowEdges.resize(nBinEdges);

  // Find lower edges of bins.
  for (int i=0; i<nBins; i++){
    binLowEdges[i]=xV[i]-binWidth/2;
    // printf("Bin Low Edge = %f\n",binLowEdges[i]);
  }
  binLowEdges[nBins]=xV[nBins-1]+binWidth/2; // adds end bin edge to vector.

  // replace hyphens with spaces
  for(int i=0;i <xTitle.size();i++)
    if(xTitle[i]=='-')
        xTitle[i]=' ';
  for(int i=0;i <yTitle.size();i++)
    if(yTitle[i]=='-')
        yTitle[i]=' ';


  // print histogram information.
  printf ("  Name:      %s\n",asymName.c_str());
  printf ("  x-Title:   %s\n",xTitle.c_str());
  printf ("  y-Title:   %s\n",yTitle.c_str());
  printf ("  n_bins:    %i\n",nBins);
  printf ("  bin_width: %f\n",binWidth);  

  // create Histogram.
  TH1D *hist = new TH1D(asymName.c_str(), " " , nBins, &(binLowEdges[0]));
  hist -> GetYaxis() -> SetTitle(yTitle.c_str());
  hist -> GetYaxis() -> CenterTitle();
  hist -> GetXaxis() -> SetTitle(xTitle.c_str());
  hist -> GetXaxis() -> CenterTitle();

  // fill Histogram.
  for (int i=0; i<nBins; i++){
    hist->Fill(xV[i],yV[i]);
  }

  return hist;
}
// =================================== end =====================================
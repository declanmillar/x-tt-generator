TH1D* plotDistribution(double luminosity, double efficiency, ifstream *logStream){

  // declarations.
  std::string distName, xTitle,  yTitle; 
  std::string xS, yS, stop; 
  int n_DataPairs; 
  double x, y; 
  std::vector<double> xV, yV; 
  int nBins, nBinEdges;
  double binWidth;
  std::vector<double> binLowEdges;
  stop = "END"; 

  // name strings.
  *logStream >> distName >> yTitle >> xTitle;

  // replace hyphens with spaces
  for (int i = 0; i < xTitle.size(); i++)
    if (xTitle[i] == '-') xTitle[i] = ' ';
  for(int i = 0; i < yTitle.size(); i++)
    if (yTitle[i] == '-') yTitle[i] = ' ';  

  // loop over data pairs.
  n_DataPairs = 0;
  if (logStream->is_open()) {
    while(1) {        
      *logStream >> xS;      
      if (xS == stop) break;
      *logStream >> yS;
      x = atof(xS.c_str());
      y = atof(yS.c_str());
      n_DataPairs += 1;
      xV.push_back(x);
      yV.push_back(y);
    }
  }
  else printf("Stream failed!");

  // histogram bin information.
  nBins = xV.size();
  nBinEdges = nBins + 1;
  binWidth = xV[1]-xV[0]; // only works for fixed bin width!
  binLowEdges.resize(nBinEdges);

  // Find lower edges of bins.
  for (int i=0; i<nBins; i++){
    binLowEdges[i]=xV[i]-binWidth/2;
    // printf("Bin Low Edge = %f\n",binLowEdges[i]);
  }
  binLowEdges[nBins]=xV[nBins-1]+binWidth/2; // adds end bin edge to vector.

  // multiply by the luminosity if specified
  if (luminosity > 0)
  {
    for (int i=0; i<nBins; i++)
    {
      yV[i]=yV[i]*luminosity*efficiency;
      // printf("Bin Low Edge = %f\n",binLowEdges[i]);
    }
    yTitle="Events";
  }

  // print histogram information.
  printf ("  Name:      %s\n",distName.c_str());
  printf ("  x-Title:   %s\n",xTitle.c_str());
  printf ("  y-Title:   %s\n",yTitle.c_str());
  printf ("  n_bins:    %i\n",nBins);
  printf ("  bin_width: %f\n",binWidth);  

  // create Histogram.
  TH1D *hist = new TH1D(distName.c_str(), " " , nBins, &(binLowEdges[0]));
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
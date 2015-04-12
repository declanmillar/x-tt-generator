TH1D* plotDistribution(double luminosity, double efficiency, ifstream *logStream){

  // declarations.
  std::string distName, xTitle,  yTitle; 
  std::string xS, yS, stop, sumw2S; 
  int n_DataPairs; 
  double x, y, sumw2; 
  std::vector<double> xV, yV, sumw2V; 
  int nBins, nBinEdges;
  double binWidth;
  std::vector<double> binLowEdges;
  stop = "END"; 

  // name strings.
  *logStream >> distName >> yTitle >> xTitle;

  // replace hyphens with spaces
  for (int i = 0; i < (int) xTitle.size(); i++)
    if (xTitle[i] == '-') xTitle[i] = ' ';
  for(int i = 0; i < (int) yTitle.size(); i++)
    if (yTitle[i] == '-') yTitle[i] = ' ';  

  // loop over data pairs.
  n_DataPairs = 0;
  if (logStream->is_open()) {
    while(1) {        
      *logStream >> xS;      
      if (xS == stop) break;
      *logStream >> yS;
      if (distName == "cost5") *logStream >> sumw2S;
      x = atof(xS.c_str());
      y = atof(yS.c_str());
      if (distName == "cost5") sumw2 = atof(sumw2S.c_str());
      n_DataPairs += 1;
      xV.push_back(x);
      yV.push_back(y);
      if (distName == "cost5") sumw2V.push_back(sumw2);
    }
  }
  else printf("Stream failed!");

  // sum all elements in sumw2 array
  double sumsumw2 = std::accumulate( sumw2V.begin(), sumw2V.end(), 0 );

  printf("%f\n", sumsumw2);

  // histogram bin information.
  nBins = xV.size();
  nBinEdges = nBins + 1;
  binWidth = xV[1]-xV[0]; // only works for fixed bin width!
  binLowEdges.resize(nBinEdges);

  // Find lower edges of bins.
  for (int i = 0; i < (int) nBins; i++){
    binLowEdges[i]=xV[i]-binWidth/2;
    // printf("Bin Low Edge = %f\n",binLowEdges[i]);
  }
  binLowEdges[nBins]=xV[nBins-1]+binWidth/2; // adds end bin edge to vector.

  // multiply by the luminosity if specified
  if (luminosity > 0)
  {
    for (int i = 0; i < (int) nBins; i++)
    {
      yV[i]=yV[i]*luminosity*efficiency;
      // printf("Bin Low Edge = %f\n",binLowEdges[i]);
    }
    yTitle="Events";
  }

  // print histogram information.
  printf ("  Name:      %s\n", distName.c_str());
  printf ("  x-Title:   %s\n", xTitle.c_str());
  printf ("  y-Title:   %s\n", yTitle.c_str());
  printf ("  n_bins:    %i\n", nBins);
  printf ("  bin_width: %f\n", binWidth);  

  // create Histogram.
  TH1D *hist = new TH1D(distName.c_str(), " ", nBins, &(binLowEdges[0]));
  hist->GetYaxis()->SetTitle(yTitle.c_str());
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitle(xTitle.c_str());
  hist->GetXaxis()->CenterTitle();

  // fill Histogram.
  for (int i = 1; i < nBins+1; i++) {
    hist->Fill(xV[i-1], yV[i-1]);
    if (distName == "cost5") hist->SetBinError(i, sqrt(sumw2V[i-1]));
  }

  return hist;
}
// =================================== end =====================================
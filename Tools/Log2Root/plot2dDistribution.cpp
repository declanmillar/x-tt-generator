// plot2dDistribution
// author: Declan Millar <d.millar@soton.ac.uk>

TH2D* plot2dDistribution(double luminosity, double efficiency, ifstream * logStream){

  // declarations.
  std::string distName, xTitle, yTitle,  zTitle;
  std::string stop;
  std::string         X,  Y,  Z;
  double              x,  y,  z;
  std::vector<double> xv, yv, zv;
  int nLines;
  int nBinsx ;//, nBinEdgesx;
  int nBinsy ;//, nBinEdgesy;
  double xlow, xup, ylow, yup;
  // double binWidthx, binWidthy;
  // std::vector<double> binLowEdgesx;
  // std::vector<double> binLowEdgesy;
  stop = "END"; // when found, stop reading the histograms

  // Precusor info
  *logStream >> distName >> zTitle
             >> xTitle >> nBinsx >> xlow >> xup
             >> yTitle >> nBinsy >> ylow >> yup ;

  // replace hyphens with spaces
  for(int i=0;i <xTitle.size();i++)
    if(xTitle[i]=='-')
        xTitle[i]=' ';
  for(int i=0;i <yTitle.size();i++)
    if(yTitle[i]=='-')
        yTitle[i]=' ';    
  for(int i=0;i <zTitle.size();i++)
    if(zTitle[i]=='-')
        zTitle[i]=' ';  

  // loop over data pairs.
  nLines=0;
  if (logStream->is_open()){
    while(1){        
      *logStream >> X;      
      if (X == stop) break;
      *logStream >> Y >> Z;
      x = atof(X.c_str());
      y = atof(Y.c_str());
      z = atof(Z.c_str());
      nLines+=1;
      std::cout << x << ' ' << y << ' ' << z << "\n";
      xv.push_back(x);
      yv.push_back(y);
      zv.push_back(z);
    }
  }
  else printf("  Stream failed!");

  // histogram bin information (only works for fixed bin width)
  // nBinEdgesx = nBinsx + 1;  
  // binWidthx = xv[1]-xv[0];
  // binLowEdgesx.resize(nBinEdgesx);

  // find lower edges of bins
  // for (int i = 0; i < nBinsx; i++){
  //   binLowEdgesx[i]=xv[i]-binWidthx/2;
  // }
  // binLowEdgesx[nBinsx]=xv[nBinsx-1]+binWidthx/2;

  // nBinEdgesy = nBinsy + 1;
  // binWidthy = yv[1]-yx[0];
  // binLowEdgesy.resize(nBinEdgesy);
  
  // find lower edges of bins
  // for (int i = 0; i < nBinsy; i++){
  //   binLowEdgesy[i]=yv[i]-binWidthy/2;
  // }
  // binLowEdgesy[nBinsy]=yv[nBinsy-1]+binWidthy/2;
  

  // multiply by the luminosity if specified
  if (luminosity > 0)
  {
    for (int i = 0; i < zv.size(); i++)
    {
      zv[i]=zv[i]*luminosity;
    }
    zTitle="Events";
  }

  // print histogram information
  printf ("  Name:      %s\n",distName.c_str());
  printf ("  z-Title:   %s\n",zTitle.c_str());
  printf ("  x-Title:   %s\n",xTitle.c_str());
  printf ("  nbinsx:    %i\n",nBinsx);
  printf ("  xlow:      %f\n",xlow);
  printf ("  xup:       %f\n",xup);
  printf ("  y-Title:   %s\n",yTitle.c_str());  
  printf ("  nbinsy:    %i\n",nBinsy);
  printf ("  ylow:      %f\n",ylow);
  printf ("  yup:       %f\n",yup);

  // create histogram
  TH2D *hist = new TH2D(distName.c_str(), " " , nBinsx, xlow, xup, nBinsy, ylow, yup);
  hist -> GetYaxis() -> SetTitle(zTitle.c_str());
  hist -> GetYaxis() -> CenterTitle();
  hist -> GetXaxis() -> SetTitle(xTitle.c_str());
  hist -> GetXaxis() -> CenterTitle();

  // fill histogram
  for (int i = 0; i < xv.size(); i++){
    hist->Fill(xv[i],yv[i],zv[i]);
  }

  return hist;
}
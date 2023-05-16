#include <iostream>
#include <iomanip>
#include <utility>
#include <fstream>
#include <tuple>
#include "gsl/gsl_integration.h"
#include <assert.h>
#include "Distributions.hpp"

//x_bins, pin_bins, g, d, n, mX
int main(int argc, char **argv)
{

  assert(argc==9);

  int nbins_x = atoi(argv[1]);
  double step_x = 1./double(nbins_x);

  int nbins_pin = atoi(argv[2]);
  double pin_min = atof(argv[3]);
  double pin_max = atof(argv[4]);
  double step_pin = (pin_max - pin_min)/double(nbins_pin);

  double g = atof(argv[5]);
  double d = atof(argv[6]);
  double n = atof(argv[7]);

  int mX = atoi(argv[8]);

  use_tables = 0;

  std::ostringstream fso;
  fso << "../../tables/m" << mX << "x_g_" << int(g) << "_d_" << int(d) << "_n_" << int(n) << ".dat";
  std::ofstream outfile(fso.str().c_str());

  //Write header
  outfile << "pin_min " << pin_min << " pin_max " << pin_max << " nbins_pin " << nbins_pin
          << " nbins_x " << nbins_x << std::endl;

  //Declare gsl integration workspace
  gsl_integration_workspace *wdk = gsl_integration_workspace_alloc(100000);
  gsl_integration_workspace *wkcm = gsl_integration_workspace_alloc(100000);
  gsl_integration_workspace *wx = gsl_integration_workspace_alloc(100000);

  //Skip x=0, since it is always 0
  for (int i=0; i<=nbins_x; i++) {
    double x = double(i)*step_x;
    
    for (int j=0; j<=nbins_pin; j++) {
      double pin = pin_min + double(j)*step_pin;
    
      double result;
      if (mX==1) result = m1x(x, pin, g, d, n, wdk, wkcm, wx);
      else if (mX==2) result = m2x(x, pin, g, d, n, wdk, wkcm, wx);
      else if (mX==3) result = m3x(x, pin, g, d, n, wdk, wkcm, wx);
      else if (mX==4) result = m4x(x, pin, g, d, n, wdk, wkcm, wx);
      else if (mX==5) result = m5x(x, pin, g, d, n, wdk, wkcm, wx);
      else if (mX==6) result = m6x(x, pin, g, d, n, wdk, wkcm, wx);
      else if (mX==7) result = m7x(x, pin, g, d, n, wdk, wkcm, wx);
      else {
        std::cout << "Unexpected value of mX = " << mX << std::endl;
	exit(EXIT_FAILURE);
      }

      std::cout << std::setprecision(12) << "Result = " << result << std::endl;
      outfile << std::setprecision(12) << result << " ";

    }
    outfile << std::endl;
  }

  outfile.close();

  return 0;
      
  //std::cout << " i= " << i << std::endl;
  //std::vector<double> parti = gen_particles(500.,0.,0.,21,0.01,wdk,wkcm,wx);

}

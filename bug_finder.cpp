#include <iostream>
#include <iomanip>
#include <utility>
#include <tuple>
#include "read_tables.hpp"
//#include "read_kcm_tables.hpp"
#include "gsl/gsl_integration.h"
#include <assert.h>
#include "Distributions.hpp"

int main()
{

  gsl_integration_workspace *wdk = gsl_integration_workspace_alloc(100000);
  gsl_integration_workspace *wkcm = gsl_integration_workspace_alloc(100000);
  gsl_integration_workspace *wx = gsl_integration_workspace_alloc(100000);

  double pin_min = 4.;
  double pin_max = 1400.;
  int pinbins = 500.;
  int xbins=100;
  double step_x = 1./double(xbins);
  double step_pin = (pin_max-pin_min)/double(pinbins);
  
  double g = 1.;
  double d = 1.;
  double n= 9.;
  
  int ijob=0;
  for (int i=0; i<=xbins; i++) {
    double x = double(i)*step_x;
    
    for (int j=0; j<=pinbins; j++) {
      double pin = pin_min + double(j)*step_pin;
  
      ijob++;
      if (ijob<26966) continue;

      std::cout << " x= " << x << " pin= " << pin << std::endl;
      double num_result = m6x(x, pin, g, d, n, wdk, wkcm, wx);
      std::cout << " Result= " << num_result << std::endl;

    }

  }
   
  return 0;
}

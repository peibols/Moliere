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

  std::string tables_path="../table_runner/tables/";
  read_tables(tables_path);
  use_tables = 1;
 
  //std::string kcm_tables_path="../table_runner/tables/";
  //read_kcm_tables(kcm_tables_path);
  //use_kcm_tables = 1;
  
  double pin_min = 4.;
  double pin_max = 1400.;
  int pinbins = 500.;
  double step_pin = (pin_max-pin_min)/double(pinbins);
  
  /*
  for (int i=0; i<1000; i++) {
    double pin=pin_min+(pin_max-pin_min)*((double) rand() / (RAND_MAX));
    std::cout << " p_in= " << pin << std::endl;
    std::vector<double> parti = gen_particles(pin,0.,0.,21,0.01,wdk,wkcm,wx);
    std::cout << " Done i= " << i << std::endl;
  }
  */

  double x = 1.;
  double g = -1.;
  double d = 1.;
  double n= 6.;
  
  for (int i=400; i<450; i++) {
    double pin = double(i)*step_pin+step_pin/7.;
    //pin = double(i)*step_pin+pin_min;
    use_tables = 0;
    double num_result = m1x(x, pin, g, d, n, wdk, wkcm, wx); 
    use_tables = 1;
    double tab_result = m1x(x, pin, g, d, n, wdk, wkcm, wx);
    std::cout << std::setprecision(9) << "Num_result= " << num_result << " Tab_result= " << tab_result << std::endl;
  }
   
  return 0;
}

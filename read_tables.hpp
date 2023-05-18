#pragma once

#include <iostream>
#include <utility>
#include <tuple>
#include <fstream>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

double ****mx_table;
//const gsl_interp2d_type *TmX = gsl_interp2d_bilinear; //SO MUCH FASTER THAN BICUBIC
//const gsl_interp2d_type *TmX = gsl_interp2d_bicubic;
//gsl_interp_accel *xacc;
//gsl_interp_accel *yacc;
//gsl_spline2d *spline;

const int SIZE_PIN = 501;
const int SIZE_X = 101;
const int SIZE_KCM = 501;

double pin_vals[SIZE_PIN];
double pin_step, pin_min, pin_max;

double x_vals[SIZE_X];
double step_x;

bool use_tables = 0;

void read_tables(std::string tables_path) {

  //Allocate tables
  int nX=7, nA=68, nP=SIZE_PIN, nY=SIZE_X;
  mx_table = (double ****)malloc(nX * sizeof(double ***));
  if (mx_table==NULL)
  {
    fprintf(stderr, "out of memory\n");
    exit(0);
  }
  for (int i=0; i<nX; i++)
  {
    mx_table[i]=(double ***)malloc(nA * sizeof(double **));
    if (mx_table[i]==NULL)
    {
      fprintf(stderr, "out of memory\n");
      exit(0);
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      mx_table[i][j]=(double **)malloc(nP * sizeof(double *));
      if (mx_table[i][j]==NULL)
      {
        fprintf(stderr, "out of memory\n");
        exit(0);
      }
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      for (int k=0; k<nP; k++)
      {
        mx_table[i][j][k]=(double *)malloc(nY * sizeof(double));
        if (mx_table[i][j][k]==NULL)
        {
          fprintf(stderr, "out of memory\n");
          exit(0);
        }
      }
    }
  }

  //Read
  int nbins_pin, nbins_x;
  for (int iX=0; iX<nX; iX++) {
  //for (int iX=0; iX<1; iX++) { //DEBUG
    int n = 1;
    for (int iA=0; iA<nA; iA++) {
      int res = (iA % 4);
      int g, d;
      if (res == 0) g = -1, d = -1;
      else if (res == 1) g = -1, d = 1;
      else if (res == 2) g = 1, d = -1;
      else if (res == 3) g = 1, d = 1;
      else {
        std::cout << "WTF res= " << res << std::endl;
	exit(1);
      }
      std::ostringstream fsi;
      fsi << tables_path.c_str() << "m" << iX+1 << "x_g_" << int(g) << "_d_" << int(d) << "_n_" << int(n) << ".dat";
      std::ifstream infile(fsi.str().c_str());
      if (infile.fail()) {
        std::cout << "No file = " << fsi.str().c_str() << std::endl;
	exit(1);
      }

      std::cout << "Reading file= " << fsi.str() << std::endl;
      std::string s;
      getline(infile,s); //Header
      if (iX==0 && iA==0) {
        std::istringstream ifh(s);
	std::string scrap;
        ifh >> scrap >> pin_min >> scrap >> pin_max >> scrap >> nbins_pin >> scrap >> nbins_x;
      }

      for (int iY=0; iY<nY; iY++) {
        getline(infile,s);
        std::istringstream ifi(s);
        for (int iP=0; iP<nP; iP++) ifi >> mx_table[iX][iA][iP][iY];
      }
      infile.close();

      //break; //DEBUG
      if (res==3) n++;

    }
    
    //break; //DEBUG
  }
  std::cout << "Finished reading mX tables" << std::endl;

  //spline = gsl_spline2d_alloc(TmX, nP, nY);
  //xacc = gsl_interp_accel_alloc();
  //yacc = gsl_interp_accel_alloc();

  pin_step = (pin_max - pin_min)/double(nbins_pin);
  for (int iP=0; iP<nP; iP++) pin_vals[iP] = pin_min + pin_step*double(iP);
  //for (int iP=0; iP<nP; iP++) std::cout << " pin val iP= " << iP << " = " << pin_vals[iP] << std::endl; 

  step_x = 1./double(nbins_x);
  for (int iY=0; iY<nY; iY++) x_vals[iY] = step_x*double(iY);
  //for (int iY=0; iY<nY; iY++) std::cout << " x_vals iY= " << iY << " = " << x_vals[iY] << std::endl;

  return;

}

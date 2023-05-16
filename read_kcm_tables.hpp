#pragma once

#include <iostream>
#include <utility>
#include <tuple>
#include <fstream>
#include <sstream>

double *****mkcm_table;

int nbins_kcm;

bool use_kcm_tables = 0;

void read_kcm_tables() {

  //Allocate tables
  int nX=7, nA=68, nP=SIZE_PIN, nY=SIZE_X, nK=SIZE_KCM;
  mkcm_table = (double *****)malloc(nX * sizeof(double ****));
  if (mkcm_table==NULL)
  {
    fprintf(stderr, "out of memory\n");
    exit(0);
  }
  for (int i=0; i<nX; i++)
  {
    mkcm_table[i]=(double ****)malloc(nA * sizeof(double ***));
    if (mkcm_table[i]==NULL)
    {
      fprintf(stderr, "out of memory\n");
      exit(0);
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      mkcm_table[i][j]=(double ***)malloc(nP * sizeof(double **));
      if (mkcm_table[i][j]==NULL)
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
        mkcm_table[i][j][k]=(double **)malloc(nY * sizeof(double *));
        if (mkcm_table[i][j][k]==NULL)
        {
          fprintf(stderr, "out of memory\n");
          exit(0);
        }
      }
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      for (int k=0; k<nP; k++)
      {
        for (int l=0; l<nY; l++)
        {
          mkcm_table[i][j][k][l]=(double *)malloc(nK * sizeof(double));
          if (mkcm_table[i][j][k][l]==NULL)
          {
            fprintf(stderr, "out of memory\n");
            exit(0);
          }
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
      fsi << "../table_runner/tables/m" << iX+1 << "kcm_g_" << int(g) << "_d_" << int(d) << "_n_" << int(n) << ".dat";
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
	nbins_kcm = nbins_pin;
      }

      for (int iY=0; iY<nY; iY++) {
        getline(infile,s);
        std::istringstream ifi(s);
        for (int iP=0; iP<nP; iP++) {
	  for (int iK=0; iK<nK; iK++) {
	    ifi >> mkcm_table[iX][iA][iP][iY][iK];
	  }
	}
      }
      infile.close();

      //break; //DEBUG
      if (res==3) n++;

    }
    
    //break; //DEBUG
  }
  std::cout << "Finished reading mkcm tables" << std::endl;

  return;

}

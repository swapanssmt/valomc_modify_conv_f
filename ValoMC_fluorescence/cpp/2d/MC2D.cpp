
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <time.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <utility>
#include <iomanip>




// define USE_OMP to utilize threading with OpenMP
//#define USE_OMP

#ifdef USE_OMP
#include <omp.h>
#endif

#include "Array.hpp"
#include "MC2D.hpp"
#include "../versionstring.h"
#include "../fileio.hpp"

MC2D MC;

int LoadProblem(char *fin);
int SaveProblem(char *fout, int time);
bool Progress(double perc);
bool suppress_progressbar;


void FinalChecks(int, int) {}; // [AL] finalize is intentionally left empty

int main(int argc, char **argv)
{

  printf("                 ValoMC-2D\n");
  char infobuf[5012];
  version_string(infobuf);
  printf("%s",infobuf);
  suppress_progressbar = false;

  // Display help
  if ((argc < 3))
  {
    printf("Use syntax:\n");
    printf(" MC2D inputfile outputfile\n");
    printf(" or \n");
    printf(" MC2D inputfile outputfile -s\n");
    printf(" to suppress progress bar.\n");
    printf("\n");
    printf("Authors: Aki Pulkkinen, Aleksi Leino and Tanja Tarvainen (2018).\n");
    printf("The simulation code is originally written by Aki Pulkkinen.\n");
    printf("\n");
    return (0);
  }

  if(argc == 4) {
    std::string param = std::string(argv[3]);
    if(param.compare("-s") == 0) suppress_progressbar = true;
  }


  MC.seed = (unsigned long)time(NULL);

  char *fin = argv[1];
  char *fout = argv[2];

  //  printf("Loading problem %s\n", fin);

  if (LoadProblem(fin))
  {
    printf("Error while loading!\n");
    return (1);
  }

  printf("Initializing MC2D\n");

  // [AL]
  try
  {
    MC.ErrorChecks();
    MC.Init();
  }
  catch (mcerror e)
  {
    printf("MC Initialization failed!: Reason %s\n", errorstring(e));
    return 1;
  }

  printf("Computing...\n");
  int tstart = (unsigned long)time(NULL);
  MC.MonteCarlo(Progress, FinalChecks);
  int tend = (unsigned long)time(NULL);

  printf("Saving problem\n");

  if (SaveProblem(fout, tend - tstart))
  {
    printf("Error while saving!\n");
    return (1);
  }

  if (MC.loss)
  {
    printf(" %ld photons lost during computation!\n", MC.loss);
  }
  printf(" %ld fluorescence photons generated during computation!\n", MC.N_F_Photons);

  //#ifdef USE_OMP
  //  double tend = omp_get_wtime();

  printf("Computation took %i seconds\n", tend - tstart);
  //#end

  return (0);
}


int LoadProblem_TXT(char *fin)
{
  /*
    File Structure for text input:

      Ne Nb Nr Nphoton NBin2Dtheta
      f phase0 Qyield_f Tau_f
      nx ny nz [AL]
      H
      BH
      r
      mua_ex_sol mua_ex_f mua_em_sol mus_ex mus_em g n
      BCType
      [BCn]
      [BCLNormal] -- if provided, BCn has to be provided as well
      [BCLightDirectionType] // [AL]
      [BCIntensity] // [AL]
      [GaussianSigma] //[AL]
*/

  int_fast64_t ii;
  int_fast64_t Ne, Nb, Nr, sd1,sd2;
  int fsr;

  FILE *fp = fopen(fin, "r");
  if (fp == NULL)
    return (1);
   

  fsr = fscanf(fp, "%li %li %li %li %li\n", &Ne, &Nb, &Nr, &MC.Nphoton, &MC.NBin2Dtheta);
  fsr = fscanf(fp, "%lf %lf %lf %lf %li %li\n", &MC.f, &MC.Qyield_f, &MC.Tau_f, &MC.phase0, &sd1, &sd2); // [AL]

  char line[5012]; // skip a line
  char *tmpbuf = fgets(line, 5011, fp);

  if(sd2) {
     MC.seed = sd1;
  } else {
     MC.seed = (unsigned long) time(NULL);
  }


  printf("Constants:\n");
  printf("  %10s   (%e)\n", "f", MC.f);
  printf("  %10s   (%e)\n", "Qyield_f", MC.Qyield_f);
  printf("  %10s   (%e)\n", "Tau_f", MC.Tau_f);
  printf("  %10s   (%e)\n", "phase0", MC.phase0);
  printf("  %10s   (%li)\n", "Ne", Ne);
  printf("  %10s   (%li)\n", "Nb", Nb);
  printf("  %10s   (%li)\n", "Nr", Nr);
  printf("  %10s   (%li)\n", "Nphoton", MC.Nphoton);
  printf("  %10s   (%li)\n", "NBin2Dtheta", MC.NBin2Dtheta);
  printf("  %10s   (%li)\n", "seed", MC.seed);
  printf("Arrays:\n");

  // make negative phase0 positive by adding a multiple of 2*pi
  if (MC.phase0 < 0)
  {
    MC.phase0 += 2 * M_PI * ceil(-MC.phase0 / (2 * M_PI));
  }
  // read the arrays
  readAndResize(fp, Ne, 3, true, &MC.H, "H");
  readAndResize(fp, Nb, 2, true, &MC.BH, "BH");
  readAndResize(fp, Nr, 2, true, &MC.r, "r");
  readAndResize(fp, Ne, 1, true, &MC.mua_ex_sol, &MC.mua_ex_f, &MC.mua_em_sol, &MC.mus_ex, &MC.mus_em, &MC.g, &MC.n, "mua_ex_sol", "mua_ex_f", "mua_em_sol", "mus_ex", "mus_em", "g", "n");
  readAndResize(fp, Nb, 1, true, &MC.BCType, "BCType");
  readAndResize(fp, Nb, 1, false, &MC.BCn, "BCn");
  readAndResize(fp, Nb, 2, false, &MC.BCLNormal, "BCLightDirection");
  readAndResize(fp, Nb, 1, false, &MC.BCLightDirectionType, "BCLightDirectionType");
  readAndResize(fp, Nb, 1, false, &MC.BCIntensity, "BCIntensity");
  readAndResize(fp, Nb, 1, false, &MC.GaussianSigma, "GaussianSigma");

  fclose(fp);

  return (0);
}

int LoadProblem(char *fin)
{
  return (LoadProblem_TXT(fin));
}

int SaveProblem_TXT(char *fout, int time)
{
  int ii;

  FILE *fp = fopen(fout, "w");
  if (fp == NULL)
    return (1);

  fprintf(fp, "%i 0\n", time);
  for (ii = 0; ii < MC.ER.Nx; ii++)
    fprintf(fp, "%e %e\n", MC.ER[ii], MC.EI[ii]);
  for (ii = 0; ii < MC.EBR.Nx; ii++)
    fprintf(fp, "%e %e\n", MC.EBR[ii], MC.EBI[ii]);
  //******************MODIFY************************
  for (ii = 0; ii < MC.R_ER.N; ii++)
    fprintf(fp, "%e %e\n", MC.R_ER[ii], MC.R_EI[ii]);
  for (ii = 0; ii < MC.R_EBR.N; ii++)
    fprintf(fp, "%e %e\n", MC.R_EBR[ii], MC.R_EBI[ii]);
  //*********************************************
  //*************** modify fluoroscence *********************
  for (ii = 0; ii < MC.F_ER.Nx; ii++)
    fprintf(fp, "%e %e\n", MC.F_ER[ii], MC.F_EI[ii]);
  for (ii = 0; ii < MC.F_EBR.Nx; ii++)
    fprintf(fp, "%e %e\n", MC.F_EBR[ii], MC.F_EBI[ii]);
  //******************MODIFY************************
  for (ii = 0; ii < MC.F_R_ER.N; ii++)
    fprintf(fp, "%e %e\n", MC.F_R_ER[ii], MC.F_R_EI[ii]);
  for (ii = 0; ii < MC.F_R_EBR.N; ii++)
    fprintf(fp, "%e %e\n", MC.F_R_EBR[ii], MC.F_R_EBI[ii]);
  //**********************************************************
  fclose(fp);

  return (0);
}

int SaveProblem(char *fin, int time)
{
  return (SaveProblem_TXT(fin, time));
}

bool Progress(double perc)
{
  if(!suppress_progressbar) {
     printf("  %f %%\r", perc);
     fflush(stdout);
  }
  return true;
}

/*=========================================================================
                                                                                
Copyright (c) 2011, Los Alamos National Security, LLC

All rights reserved.

Copyright 2011. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

#include "OutGnuplot.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

OutGnuplot::OutGnuplot(PlasmaData* data) : PlasmaOut(data)
{
}

OutGnuplot::~OutGnuplot()
{
}

///////////////////////////////////////////////////////////////////////
//
// Output gnuplot PNG files for a time step
//
///////////////////////////////////////////////////////////////////////

void OutGnuplot::outputPlots(
                        int tIndx,
                        MomentSolution* curMoment,
                        ParticleSolution* curParticle,
                        double* time_vec,
                        double** cc_rms,
                        double* E_tot,
                        double* E_diff,
                        double* mv_diff)
{
  ostringstream tstr;
  double xrange[2], yrange[2];

  // Problem description
  if (this->data->prob_type == LANDAU_DAMPING)
    tstr << " (LANDAU";
  else if (this->data->prob_type == TWO_STREAM_INSTABILITY)
    tstr << " (TSI";
  else
    tstr << " (IASW";
  tstr << " cc=" << this->data->cc_flag
          << " order=" << this->data->temp_order << ")";

  // Particles
  outputParticles(tIndx, curParticle, tstr.str().c_str());

  // Density
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    outputSystemStat(tIndx, sIndx, 
       "out_r_", "Density", tstr.str().c_str(), "X", "rho", 
       this->data->p_color[sIndx].c_str(), this->data->p_name[sIndx].c_str(),
       this->data->nx, this->data->xpos_node, curParticle->getR(sIndx));

  // Fluid momentum
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    outputSystemStat(tIndx, sIndx, 
       "out_ru_", "Fluid momentum", tstr.str().c_str(), "X", "ru", 
       this->data->p_color[sIndx].c_str(), this->data->p_name[sIndx].c_str(),
       this->data->nfx, this->data->xpos_face, curParticle->getRU(sIndx));

  // Stress tensor
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    outputSystemStat(tIndx, sIndx, 
       "out_s2_", "Stress tensor", tstr.str().c_str(), "X", "s2", 
       this->data->p_color[sIndx].c_str(), this->data->p_name[sIndx].c_str(),
       this->data->nx, this->data->xpos_node, curParticle->getS2(sIndx));

  // Electric field
  outputSystemStat(tIndx, 0,
       "out_E_", "Electric field", tstr.str().c_str(), "X", "E", 
       "blue", "TOTAL",
       this->data->nfx, this->data->xpos_face, curMoment->getE());

  // Charge conservation over time
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    outputTimeStat(tIndx, sIndx,
       "out_ccrms_", "Charge conservation", tstr.str().c_str(),
       "time", "cc_rms", 
       this->data->p_color[sIndx].c_str(), this->data->p_name[sIndx].c_str(),
       time_vec, cc_rms[sIndx]);

  // Total energy over time
  outputTimeStat(tIndx, 0,
       "out_Etot_", "E field energy vs time", tstr.str().c_str(), 
       "time", "E_tot", "blue", "E_tot", time_vec, E_tot);

  // Delta energy over time
  outputTimeStat(tIndx, 0,
       "out_Ediff_", "Delta E field vs time", tstr.str().c_str(), 
       "time", "E_diff", "blue", "E_diff", time_vec, E_diff);

  // Total energy over time
  outputTimeStat(tIndx, 0,
       "out_mvdiff_", "Delta momentum vs time", tstr.str().c_str(), 
       "time", "mv_diff", "blue", "mv_diff", time_vec, mv_diff);
}

///////////////////////////////////////////////////////////////////////
//
// Output particles
//
///////////////////////////////////////////////////////////////////////

void OutGnuplot::outputParticles(
                        int tIndx,
                        ParticleSolution* curParticle,
                        const char* title2)
{
  ostringstream fstr, ostr, tstr, cstr;

  // Range on gnuplot
  double xrange[2], yrange[2];
  xrange[0] = 0.0;  xrange[1] = this->data->lx;
  yrange[0] = -0.25; yrange[1] = 0.25;

  // Output all species in one gnuplot
  ostr << "output/out_particles_" << tIndx << ".png";
  tstr << "Phase space" << title2;

  // Script file to write to
  ofstream script("plasma.script");

  script << "set title \"" << tstr.str().c_str() << "\"" << endl;
  script << "set xlabel \"X\"" << endl;
  script << "set ylabel \"V\"" << endl;
  script << "set xrange [" << xrange[0] << ":" << xrange[1] << "]" << endl;
/*
  script << "set yrange [" << yrange[0] << ":" << yrange[1] << "]" << endl;
*/

  // Plot data to get output EPS file
  script << "set term png" << endl;
  script << "set output \"" << ostr.str().c_str() << "\"" << endl;
  script << "plot \\" << endl;

  // Write particles for each species into files
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    fstr.seekp(0); fstr << "output/out_particles_" << sIndx << std::ends;
    ofstream pStream(fstr.str().c_str(), ios::out);

    double* x = curParticle->getX(sIndx);
    double* v = curParticle->getV(sIndx);
    int numberOfParticles = curParticle->getNumberOfParticles(sIndx);

    for (int pIndx = 0; pIndx < numberOfParticles; pIndx++)
      pStream << x[pIndx] << "   " << v[pIndx] << endl;
    pStream.close();

    script << "   \"" << fstr.str().c_str() << "\"";
    script << " pt 0 title \"" << this->data->p_name[sIndx] << "\"";
    if (sIndx < (this->data->p_size - 1))
      script << ", \\";   
    script << endl; 
  }
  script.close();

  // Call gnuplot to execute script and process the output svg file
  cstr << "gnuplot plasma.script";
  system(cstr.str().c_str());
}

///////////////////////////////////////////////////////////////////////
//
// Output a statistic over the length of the system
//
///////////////////////////////////////////////////////////////////////

void OutGnuplot::outputSystemStat(
                        int tIndx,
                        int sIndx,
                        const char* fileBase,
                        const char* title1,
                        const char* title2,
                        const char* xTitle,
                        const char* yTitle,
                        const char* lineColor,
                        const char* lineName,
                        int numData,
                        double* xData,
                        double* yData)
{
  ostringstream fstr, ostr, tstr, lstr;
  tstr << title1 << title2 << std::ends;

  // Range on gnuplot
  double xrange[2], yrange[2];
  xrange[0] = 0.0;  xrange[1] = this->data->lx;
  yrange[0] = yrange[1] = 0.0;

  // Write data for time and species to file
  fstr.seekp(0); fstr << "output/" << fileBase << sIndx << std::ends;
  ofstream dStream(fstr.str().c_str(), ios::out);
  for (int indx = 0; indx < numData; indx++)
    dStream << xData[indx] << "   " << yData[indx] << endl;
  dStream.close();

  // Output file name and line style
  ostr << "output/" << fileBase << sIndx << "_" << tIndx << ".png";
  lstr << " with lines lw 2 lc rgb \"" << lineColor
       << "\" title \"" << lineName << "\"";

  // Plot the data
  runGnuplot(ostr.str().c_str(), tstr.str().c_str(), 
             xTitle, xrange, yTitle, yrange,
             fstr.str().c_str(), lstr.str().c_str());
}

///////////////////////////////////////////////////////////////////////
//
// Output a statistic over time
//
///////////////////////////////////////////////////////////////////////

void OutGnuplot::outputTimeStat(
                        int tIndx,
                        int sIndx,
                        const char* fileBase,
                        const char* title1,
                        const char* title2,
                        const char* xTitle,
                        const char* yTitle,
                        const char* lineColor,
                        const char* lineName,
                        double* time_vec,
                        double* yData)
{
  ostringstream fstr, ostr, tstr, lstr;
  tstr << title1 << title2 << std::ends;

  // Range on gnuplot
  double xrange[2], yrange[2];
  xrange[0] = 0.0;  xrange[1] = this->data->tmax;
  yrange[0] = yrange[1] = 0.0;

  // Write data for time and species to file
  fstr.seekp(0); fstr << "output/" << fileBase << sIndx << std::ends;
  ofstream dStream(fstr.str().c_str(), ios::out);
  for (int indx = 0; indx < tIndx; indx++)
    dStream << time_vec[indx] << "   " << yData[indx] << endl;
  dStream.close();

  // Output file name and line style
  ostr << "output/" << fileBase << sIndx << "_" << tIndx << ".png";
  lstr << " with lines lw 2 lc rgb \"" << lineColor
       << "\" title \"" << lineName << "\"";

  // Plot the data
  runGnuplot(ostr.str().c_str(), tstr.str().c_str(), 
             xTitle, xrange, yTitle, yrange,
             fstr.str().c_str(), lstr.str().c_str());
}

///////////////////////////////////////////////////////////////////////
//
// Write the script file which will drive gnuplot and run gnuplot
//
///////////////////////////////////////////////////////////////////////

void OutGnuplot::runGnuplot(
                        const char* outFile,
                        const char* title,
                        const char* xlabel,
                        double* xrange,
                        const char* ylabel,
                        double* yrange,
                        const char* dataFile,
                        const char* lineStyle)
{
  // Script file to write to
  ofstream script("plasma.script");

  script << "set title \"" << title << "\"" << endl;
  script << "set xlabel \"" << xlabel << "\"" << endl;
  script << "set ylabel \"" << ylabel << "\"" << endl;
  script << "set xrange [" << xrange[0] << ":" << xrange[1] << "]" << endl;
  if (yrange[0] != yrange[1])
    script << "set yrange [" << yrange[0] << ":" << yrange[1] << "]" << endl;

  // Plot data to get output EPS file
  script << "set term png" << endl;
  script << "set output \"" << outFile << "\"" << endl;
  script << "plot \"" << dataFile << "\"" << lineStyle << endl;
  script.close();

  // Call gnuplot to execute script and process the output svg file
  ostringstream ostr;
  ostr << "gnuplot plasma.script";
  system(ostr.str().c_str());
}

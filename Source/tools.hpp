//
//  tools.hpp
//  NewProject - ConsoleApp
//
//  Created by Livio Saracino on 25/03/2019.
//

#ifndef tools_hpp
#define tools_hpp

#include <stdio.h>
#include "../JuceLibraryCode/JuceHeader.h"
#include <math.h>
#include <list>
#include <array>
// libindigo and armadillo have been manually copied in the build folder !!
#include "indigo/signal/dynamics.hpp"
#define ARMA_DONT_USE_WRAPPER //cazzo e' ??
#include "armadillo"


#endif /* tools_hpp */


struct envelope_type {
    arma::vec envelope;
    arma::vec time;
};


//open file wav
void open(std::vector<double> &_inSignal, File &file, double &_fs, double &_inSignalLen);

//compute rms or peak envelope
envelope_type Envelope(arma::vec &x, const double &fs, const double &rms_length, const double &overlap, const String type);
//reduce the dynamic range discarding value lower than min and return a type ready for GMM_AIC usage
arma::mat ReduceDynamicRange(envelope_type &envelope, const double min_value );

//microdynamics
std::vector< std::array<int,2> > microdynamics(envelope_type &envelope_rms_fast, envelope_type &envelope_rms_slow, double rms_fast_length, double rms_slow_length, double overlap_rms_fast, int hystereisi_high_threshold, int hystereisi_low_threshold, double &start);

//compute gmm model
arma::gmm_diag GMM_AIC(arma::mat &rmsdb, int N_components);

//pdf and cdf
arma::vec pdf( double mean, double variance, arma::vec &x);
arma::vec cdf(arma::vec &pdf);
//percentile from gmm model
double getPercentile(arma::gmm_diag model, double percentile, double sample_space_start, double sample_space_end, double sample_space_points);
//getModelPdf
// ..
//getGaussianPercentile
double getGausPercentile(double percentile, double mean, double variance, arma::vec &x);

//linear regression
arma::vec linearFit(arma::vec x, arma::vec y);
//Schroder
double Schroder_BackwardIntegration(arma::vec &envelope, double hop_size, double tx);


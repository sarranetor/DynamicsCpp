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


/**
    Computes signal envelope rms or peak.
 
    @param
    x the signal
    fs the sampling frequency
    rms_length the width of the window in ms
    overlap % of window overlap
    type has to be "peak" or "rms"
    @return the envelope vector
**/
envelope_type Envelope(arma::vec &x, const double &fs, const double &rms_length, const double &overlap, const String type);

//reduce the dynamic range discarding value lower than min and return a type ready for GMM_AIC usage
arma::mat ReduceDynamicRange(envelope_type &envelope, const double min_value );

/**
    thanks to two rms envelope with different window length finds the points in which they cross.
    it computes the difference btw them (fast-slow), and thanks to an hysteresis comparator finds the zero crossings.
 
    @param info about the two envelopes and the hysteresis ramge.
    start index of rms_fast that corresponds in time to the first rms_slow element
    @return vector of array. each array contains the zero crossing index plus a fleg that identify the type of crossing (ascendent or descendent).
**/
std::vector< std::array<int,2> > microdynamics(envelope_type &envelope_rms_fast, envelope_type &envelope_rms_slow, double rms_fast_length, double rms_slow_length, double overlap_rms_fast, int hystereisi_high_threshold, int hystereisi_low_threshold, double &start);

/**
    computes the gaussian mixture model of a matrix/vector.
 
    @param
    signal: each coloumn is an element to fit the model.
    N_components: max number of gaussins.
    @return the output model is chosen with the AIC.
**/
arma::gmm_diag GMM_AIC(arma::mat &signal, int N_components);

/**
    computes percentile from gmm model.
 **/
double getPercentile(arma::gmm_diag model, double percentile, double sample_space_start, double sample_space_end, double sample_space_points);

/**
    computes the probability density function of a gaussian curve.
**/
arma::vec pdf( double mean, double variance, arma::vec &x);

/**
    computes the cumulative distribution function from a pdf.
    note: the input vector mast respect the definition of pdf -> area under the curve = 1.
**/
arma::vec cdf(arma::vec &pdf);

/**
    computes percentile from a gaussian curve.
**/
double getGausPercentile(double percentile, double mean, double variance, arma::vec &x);

/**
    linear regression using Minimum mean square error [MMSE] estimator. y = b1*x + b0.
 
    @return [b1,b0]
**/
arma::vec linearFit(arma::vec x, arma::vec y);

/**
    computes a certain decay time using the Schroder_BackwardIntegration technique.
 
    @param
    signal or envelope of which we want to compute a certain decay. (impulse responce)
    hop_size btw one point and other in ms.
    tx Db ..
    @return Time that the signal needed to reduce its energy by tx Db.
 **/
double Schroder_BackwardIntegration(arma::vec &envelope, double hop_size, double tx);


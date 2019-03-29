//
//  algo_impl.hpp
//  NewProject - ConsoleApp
//
//  Created by Livio Saracino on 29/03/2019.
//

#ifndef algo_impl_hpp
#define algo_impl_hpp

#include <stdio.h>
#include "../JuceLibraryCode/JuceHeader.h"
#include <math.h>
#include <list>
// libindigo and armadillo have been manually copied in the build folder !!
#include "indigo/signal/dynamics.hpp"
#define ARMA_DONT_USE_WRAPPER //cazzo e' ??
#include "armadillo"

#include "tools.hpp"

#endif /* algo_impl_hpp */

// gate and compressor drums smart preset
std::map<std::string, double> snare_kick_dynamics(arma::vec &inSignal, double &_fs);
//gate threshold computation
double GateThreshold(arma::gmm_diag &model);

// compressor melodic smart preset
std::map<std::string, double> voice_dynamics(arma::vec &inSignal, double &_fs);

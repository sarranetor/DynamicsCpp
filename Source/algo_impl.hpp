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

// RATIO
#define COMP_RATIO_VOICE 1.0/2.0
#define COMP_RATIO_KICK_SNARE 1.0/3.0

// DB RANGE
#define MIN_DB_VALUE -100
#define MAX_DB_VALUE 0
#define NUMBER_DB_VALUES 1000

#define Gate_th_Percentile 0.98
#define MaxLevel_Percentile 0.998
#define MaxAfterComp_Percentile 0.99

// ms RANGE
#define MIN_VALUE 0
#define MAX_ATTACK_VALUE 100
#define MAX_RELEASE_VALUE 200
#define MS_SENSIBILITY 0.2

#define Ta_Percentile 0.5
#define Tr_Percentile 0.5

#endif /* algo_impl_hpp */


// gate and compressor drums smart preset
std::map<std::string, double> snare_kick_dynamics(arma::vec &inSignal, double &_fs);
// gate threshold computation
double GateThreshold(arma::gmm_diag &model);

// compressor melodic smart preset
std::map<std::string, double> voice_dynamics(arma::vec &inSignal, double &_fs);

// other instruments
// ..

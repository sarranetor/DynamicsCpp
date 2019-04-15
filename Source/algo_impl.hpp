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


/**
 gate and compressor drums smart preset
 
 @return a map with the following keys.
 gate_threshold, gate_hold, gate_release, compressor_threshold, compressor_ratio, compressor_attack, compressor_release, compressor_gain.
 their value is {-1} when the algorithm wasn't able to compute the value for some reason.
 **/
template<typename T>
std::map<std::string, double> snare_kick_dynamics(T* data, int size, double _fs);


/**
 compressor melodic smart preset
 
 @return a map with the following keys.
 comp_threshold, comp_ratio, comp_makeupgain, comp_attack, comp_release.
 their value is {-1} when the algorithm wasn't able to compute the value for some reason.
 **/
template<typename T>
std::map<std::string, double> voice_dynamics(T* data, int size, double _fs);


// other instruments
// ..

/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic startup code for a JUCE application.

  ==============================================================================
*/

#include "../JuceLibraryCode/JuceHeader.h"
#include <math.h>
#include <list>
// libindigo and armadillo have been manually copied in the build folder !!
#include "indigo/signal/dynamics.hpp"
#define ARMA_DONT_USE_WRAPPER //cazzo e' ??
#include "armadillo"

// my stuff
#include "algo_impl.hpp"


// function declaration ========================================================
// ..

int main (int argc, char* argv[])
{

    //==================================================
    // AUDIO FILE LOAD
    //==================================================
    
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/MoosSnareTop1.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/CommitmentSnareTop1.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/KillerQueen1.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Voice/segments/FunnyValentines1.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Voice/segments/DunningKruger1_10s.wav");
    
    std::string track = "neckdeep_10s.wav";
    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Voice/segments/" + track);
    
//    std::string track = "traffiker2.wav";
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/" + track);
    
    std::vector<double> _inSignal;
    double _fs;
    double _inSignalLen;
    // the incoming signal is stored in a col arma vec
    arma::vec inSignal;
    
    open(_inSignal, file, _fs, _inSignalLen);

    //copy _inSignal to inSignal with iterators
    inSignal.zeros(_inSignalLen);
    arma::vec::iterator ptrarma;
    std::vector<double>::iterator ptr;
    for(ptrarma = inSignal.begin(), ptr = _inSignal.begin(); ptrarma < inSignal.end() && ptr <_inSignal.end(); ptrarma++, ptr++)
        *ptrarma = *ptr;
    
    
    //=================================================
    // SNARE / KICK GATE AND COMPRESSOR SMART PRESET
    //=================================================
//    std::map<std::string, double> m = snare_kick_dynamics(inSignal, _fs);
//    //print map
//    for (auto it = m.cbegin(); it != m.cend(); ++it) {
//        std::cout << "{" << (*it).first << ": " << (*it).second << "}\n";
//    }
    
    //=================================================
    // VOICE COMPRESSOR SMART PRESET
    //=================================================
    std::map<std::string, double> m = voice_dynamics(inSignal, _fs);
    //print map
    for (auto it = m.cbegin(); it != m.cend(); ++it) {
        std::cout << "{" << (*it).first << ": " << (*it).second << "}\n";
    }

    return 0;
}


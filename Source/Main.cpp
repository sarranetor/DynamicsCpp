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
#include "indigo/modules/module_sndfile.h"
#define ARMA_DONT_USE_WRAPPER //cazzo e' ??
#include "armadillo"

// my stuff
#include "algo_impl.hpp"


// function  ========================================================
/**
 open file wav
 **/
void open(std::vector<double> &_inSignal, File &file, double &_fs, double &_inSignalLen);

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
    
//    std::string track = "lawsuits1.wav";
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Voice/segments/" + track);
//
    std::string track = "allflat.wav";
    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/" + track);
    
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

//================================================
// JUCE OPEN A FILE          [audio sample buffer appears to support float only]
//================================================
void open(std::vector<double> &_inSignal, File &file, double &_fs, double &_inSignalLen)
{
    // deposit buffers: can be preallocated if we know the size
    AudioBuffer<float> fileBuffer;
    AudioFormatManager formatManager;
    
    //set the formatManager to deal with the various audio formats
    formatManager.registerBasicFormats();
    
    std::unique_ptr<AudioFormatReader> reader(formatManager.createReaderFor (file));
    
    _fs = reader->sampleRate;
    _inSignalLen = reader->lengthInSamples;
    //give size and fill the vector (still a dynamic operation)
    _inSignal.resize(_inSignalLen);
    std::fill(_inSignal.begin(), _inSignal.begin()+_inSignalLen, 0);
    
    // set size of AudioBuffer juce class and copy the audio there
    fileBuffer.setSize (reader->numChannels, (int) reader->lengthInSamples);
    reader->read (&fileBuffer, 0, _inSignalLen, 0, true, true);
    
    //copy the audioBuffer content in std::vec
    for (int n=0; n<_inSignalLen; ++n)
        _inSignal[n] = fileBuffer.getSample(0, n);
    
    std::cout << std::endl;
    std::cout << "fs = " << _fs << std::endl;
    std::cout << "# samples = " << _inSignalLen << std::endl;
    std::cout << "Upload -> done" << std::endl;
    
    reader.reset ();
}

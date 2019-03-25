/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic startup code for a JUCE application.

  ==============================================================================
*/

#include "../JuceLibraryCode/JuceHeader.h"
#include <math.h>
// libindigo and armadillo have been manually copied in the build folder !!
#include "indigo/signal/dynamics.hpp"
#define ARMA_DONT_USE_WRAPPER //cazzo e' ??
#include "armadillo"


//==============================================================================

struct envelope_type {
    arma::vec envelope;
    arma::vec time;
};

// function declaration ========================================================
//open file wav
void open(std::vector<double> &_inSignal, File &file, double &_fs, double &_inSignalLen);
//compute rms or peak envelope
envelope_type Envelope(arma::vec &x,double &fs, double &rms_length, double &overlap, String type);
//compute gmm model
arma::gmm_diag GMM_AIC(arma::mat &rmsdb, int N_components);
//pdf and cdf
arma::vec pdf( double mean, double variance, arma::vec &x);
arma::vec cdf(arma::vec &pdf);

//gate threshold computation
double GateThreshold(arma::gmm_diag model);


int main (int argc, char* argv[])
{
//    File file("/Users/Livio.Saracino/Developer/DunningKruger1_10sComp.wav");
    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/MoosSnareTop1.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Voice/segments/FunnyValentines1.wav");
    
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
    
    double rms_length = 10;
    double overlap = 0.8;
    envelope_type envelopeRms, envelopePeak;
    envelopeRms = Envelope(inSignal, _fs, rms_length, overlap, "rms");
    envelopePeak = Envelope(inSignal, _fs, rms_length, overlap, "peak");
    
//    envelopePeak.envelope.print();
//    envelopePeak.time.print();
    
    
    //data as arma::mat. each col is an element. //deleate the element below -100 db. or sub setting. silence!
    
    // counting number of silence elements
    int N_silence_elem_rms = 0;
    int N_silence_elem_peak = 0;
    arma::vec::iterator ptr_vec;
    
    for (ptr_vec = envelopeRms.envelope.begin() ; ptr_vec != envelopeRms.envelope.end() ; ptr_vec++)
        if ( *ptr_vec < -100 )
            N_silence_elem_rms++;
    
    for (ptr_vec = envelopePeak.envelope.begin() ; ptr_vec != envelopePeak.envelope.end() ; ptr_vec++)
        if ( *ptr_vec < -100 )
            N_silence_elem_peak++;
    
    // vector for gmm modelling RMS --------
    arma::mat rmsdb(1, envelopeRms.envelope.size() - N_silence_elem_rms );
    for (arma::mat::iterator ptr_mat = rmsdb.begin(), ptr_vec = envelopeRms.envelope.begin(); ptr_mat != rmsdb.end() && ptr_vec != envelopeRms.envelope.end(); ptr_vec++ ){
        if ( *ptr_vec > -100 ){
            *ptr_mat = *ptr_vec;
            ptr_mat++;
        }
    }
    
    
    // vector for gmm modelling Peak --------
    arma::mat peak_db(1, envelopePeak.envelope.size() - N_silence_elem_peak );
    for (arma::mat::iterator ptr_mat = peak_db.begin(), ptr_vec = envelopePeak.envelope.begin(); ptr_mat != peak_db.end() && ptr_vec != envelopePeak.envelope.end(); ptr_vec++ ){
        if ( *ptr_vec > -100 ){
            *ptr_mat = *ptr_vec;
            ptr_mat++;
        }
    }
    
//    peak_db.print();
    

//    arma::mat rmsdb{envelopeRms.envelope.t()};
//    envelopeRms.envelope.print();
//    rmsdb.t().print();
    
    std::cout << "\nsize " << rmsdb.size() << std::endl;
//    rmsdb.print();
    
    
    arma::gmm_diag test; // modelDep deposit variable
    test.learn(rmsdb, 5, arma::eucl_dist, arma::random_subset, 20, 10, 1e-6, true);
    std::cout << "\nmeans" << std::endl;
    test.means.print();
    std::cout << "conv" << std::endl;
    test.dcovs.print();
    std::cout << "wheigth" << std::endl;
    test.hefts.print();
//    model.avg_log_p(envelope.envelope);
    double overall_likelihood = test.avg_log_p(rmsdb);
    printf("\nLikehood:%lf", overall_likelihood);
    
    test.learn(peak_db, 7, arma::eucl_dist, arma::random_subset, 20, 10, 1e-6, true);
    std::cout << "\nmeans" << std::endl;
    test.means.print();
    std::cout << "conv" << std::endl;
    test.dcovs.print();
    std::cout << "wheigth" << std::endl;
    test.hefts.print();
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Rms{GMM_AIC(rmsdb, 7)};
    std::cout << "\nmeans" << std::endl;
    model_Rms.means.print();
    std::cout << "conv" << std::endl;
    model_Rms.dcovs.print();
    std::cout << "wheigth" << std::endl;
    model_Rms.hefts.print();
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Peak{GMM_AIC(peak_db, 7)};
    std::cout << "\nmeans" << std::endl;
    model_Peak.means.print();
    std::cout << "conv" << std::endl;
    model_Peak.dcovs.print();
    std::cout << "wheigth" << std::endl;
    model_Peak.hefts.print();
    
    
    //=========================
    //check if threshold works RMS
    //=========================
    // function needed: generate pdf, generate cdf from pdf, and calculate percentile.
        // arma::normpdf  arma::normcdf
    
    // select second laudest gaussian
        // arma::find  arma::indexmax
    
    arma::uvec indexMinorThanLoudest = arma::find( model_Peak.means < model_Peak.means.max() );
    arma::uvec indexSecondLoudest = arma::find( model_Peak.means == model_Peak.means.elem(indexMinorThanLoudest).max() );
    // [mean, cov, wheight]
    arma::vec SecondLaudestGaussian_Params{model_Peak.means[indexSecondLoudest[0]], model_Peak.dcovs[indexSecondLoudest[0]], model_Peak.hefts[indexSecondLoudest[0]]};
    
//    SecondLaudestGaussian_Params.print();
    arma::vec dbAx{arma::linspace<arma::vec>(-100, 0, 1000)}; // 1000 elements, step 0.1 db
    
//  pdf single gaussian calculation
    arma::vec pdf_secondlaudest{pdf( model_Peak.means[indexSecondLoudest[0]] , model_Peak.dcovs[indexSecondLoudest[0]], dbAx)};
//    arma::vec pdf_secondlaudest{pdf( -24.38309067 , 12.42126215, dbAx)};
    pdf_secondlaudest.print();
//    pdf_secondlaudest = pdf_secondlaudest * SecondLaudestGaussian_Params[2];
    
    // cdf
    arma::vec cdf_secondlaudest{ cdf(pdf_secondlaudest) };
    cdf_secondlaudest.print();
    
    // percentile computation
    arma::uvec indexes = arma::find( cdf_secondlaudest < 0.98 );
    double gate_thershold = dbAx[ indexes.max() ];
    
    //=========================
    //check if threshold works PEAK
    //=========================
    

    
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


//================================================
// GMM + AIC AND RETURN THE CHOSEN MODEL
//================================================
arma::gmm_diag GMM_AIC(arma::mat &signal, int N_components)
{
    arma::gmm_diag model, modelDep; // modelDep deposit variable
    double overall_likelihood, oldAIC = 1e10, AIC = 0, AICtest = 0;
    
    for ( int i=1; i<=N_components; i++){
        
        modelDep.learn(signal, i, arma::eucl_dist, arma::random_spread, 20, 10, 1e-6, false);
        // per sample average log likehood
        overall_likelihood = modelDep.avg_log_p(signal);
        
        arma::rowvec eachsample_likelihood = modelDep.log_p(signal);
        double likehood = arma::sum(eachsample_likelihood);
        
        // log likehood = per sample average log likehood * N_samples
        // AIC = -2*log likehood + 2*N_modelParam
        
        AICtest = -2 * likehood + 2*(modelDep.means.size() + modelDep.dcovs.size() + i - 1);
        AIC = -2 * overall_likelihood * signal.size() + 2*(modelDep.means.size() + modelDep.dcovs.size() + i - 1);
        
        printf("\n iter:%d AIC:%lf", i, AIC);
        printf("\n iter:%d AICtest:%lf", i, AICtest);
        
        if ( AIC < oldAIC ){
            oldAIC = AIC;
            model = modelDep;
        }
    }
    
    return model;
}


//================================================
// RMS OR PEAK ENVELOPE COMPUTATION
//================================================
envelope_type Envelope(arma::vec &x,double &fs, double &rms_length, double &overlap, String type)
{
    double x_length = std::floor(x.size());
    double rms_samples =  std::floor((rms_length * (fs/1000)));
    double hop_size =  std::floor(((1-overlap) * rms_samples));
    double n_loop =  std::floor(((x_length - rms_samples)/hop_size));
    
    // rms envelope output
    envelope_type envelope;
    envelope.envelope.zeros(n_loop);
    envelope.time.zeros(n_loop);
    
    double startWindow, endWindow;
    double EnvelopeValue, time;
    arma::vec::iterator ptr_env, ptr_time;
    //inizialise iterators
    ptr_env = envelope.envelope.begin();
    ptr_time = envelope.time.begin();
    
    //conpute envelope ------------------------------------
    for(int n=0; n<n_loop; n++){
        //window extremes
        startWindow = n * hop_size;
        endWindow = (startWindow+rms_samples)<x_length ? (startWindow+rms_samples) : x_length;
        
        //compute rmsValue or peak and corresponding time of the computed value
        if(type=="peak"){
            EnvelopeValue = arma::abs(x(arma::span(startWindow,endWindow))).max();
        }
        else if(type=="rms"){
//            x(arma::span(startWindow,endWindow)).print();
            
            EnvelopeValue = arma::mean( arma::square( x(arma::span(startWindow,endWindow)) ));
            EnvelopeValue = sqrt(EnvelopeValue);
        }
        
        time = (startWindow + endWindow) / fs / 2.0;
        
        //save the values
        *(ptr_env+n) = EnvelopeValue; // db transformation // if n is not valid ??
        *(ptr_time+n) = time;
        
        //        printf(" length: %d", endWindow - startWindow );
        //        printf(" mean: %f", rmsValue );
        //        printf(" time: %f", time );
    }
    
    // db transformation
    envelope.envelope = 20* arma::log10(envelope.envelope + 0.0000001);
    
    
    return envelope;
}

//================================================
// PDF OF GAISSIAN
//================================================
arma::vec pdf( double mean, double variance, arma::vec &x)
{
    arma::vec pdf(x.size());
    
    for (int i=0; i<x.size(); i++){
//        double  coeff = (1/(std::sqrt(2 * arma::datum::pi * variance)));
//        double indip = std::pow((x(i) - mean) , 2);
        pdf(i) = (1/(std::sqrt(2 * arma::datum::pi * variance))) * std::exp( - std::pow((x(i) - mean) , 2) / (2 * variance) );
    }
    
    return pdf;
}

//================================================
// CDF FROM PDF  [do not like this implementation much]
//================================================
arma::vec cdf(arma::vec &pdf)
{
    //this should be an integral, but not sure if for a normal distrubution has a closed form.
    
    arma::vec cdf{arma::cumsum(pdf)};
    cdf = cdf / cdf.max();
    
    return cdf;
}


//================================================
// CUSTOM ALGO ANALYSIS FUNCTIONS
//================================================

double GateThreshold(arma::gmm_diag model)
{
    
}

//
//  tools.cpp
//  NewProject - ConsoleApp
//
//  Created by Livio Saracino on 25/03/2019.
//

#include "tools.hpp"


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
        //        printf("\n iter:%d AICtest:%lf", i, AICtest);
        
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
// Schroder_BackwardIntegration
//================================================
double Schroder_BackwardIntegration(arma::vec &envelope, double tx)
{
    //take the signal only after the max
    double max = envelope.max();
    
    
    // ^2 and backward integration
    
    //Db trasformation and normalization
    
    //fitting and Tx calculation
    
    return 0;
}

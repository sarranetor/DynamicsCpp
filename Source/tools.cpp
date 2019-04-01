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
// MICRODYNAMICS // transient event list [[index in rms fast, type event], ..]  type event: ascendent[1]/descendent[-1]
//================================================
std::list< std::vector<int> > microdynamics(envelope_type &envelope_rms_fast, envelope_type &envelope_rms_slow, double rms_fast_length, double rms_slow_length, double overlap_rms_fast, int hystereisi_high_threshold, int hystereisi_low_threshold, double &start)
{
    std::list< std::vector<int> > transient_list;
    // first sample of rms_slow is later that rms_fas cause bigger window
    start = std::round( (rms_slow_length/2 - rms_fast_length/2) / (rms_fast_length * (1 - overlap_rms_fast)) );
    
    // rms_fast - rms_slow
    arma::vec DiffRms{ envelope_rms_fast.envelope(arma::span(start, envelope_rms_slow.envelope.size()-1+start)) - envelope_rms_slow.envelope };
    
    //flag variables
    int dep = 0; // 0: below threshold         1: above threshold
    //set flag variable
    if ( DiffRms[0] < 0 )
        dep = -1;
    // transient event [index in rms fast, type event]  type event: ascendent[1]/descendent[-1]
    std::vector<int> transient(2);
    
    
    //================================================ iteration among DiffRms values
    int i=0;
    for ( arma::vec::iterator ptr_DiffRms=DiffRms.begin(); ptr_DiffRms<DiffRms.end(); i++, ptr_DiffRms++ ) //TOO CONVOLUTED!
    {
        if ( *ptr_DiffRms < hystereisi_low_threshold ) //descending
        {
            if (dep == 1)
            {
                transient[0] = i - 1 ; // - 1 cause i want to catch the transient the sample before happening  // CAREFULL FOR BUGS !!
                transient[1] = -1;
                transient_list.push_back(transient);
            }
            //end processing
            dep = -1;
        }
        else if ( *ptr_DiffRms > hystereisi_high_threshold ) // ascending
        {
            if (dep == -1)
            {
                transient[0] = i - 1; // - 1 cause i want to catch the transient the sample before happening
                transient[1] = 1;
                transient_list.push_back(transient);
            }
            //end processing
            dep = 1;
        }
    }
    
    
    return transient_list;
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
// PERCENTILE FROM GMM MODEL -- sample space range should contain the pdf (pdf > 1e-10 in sample space range)
//===============================================
double getPercentile(arma::gmm_diag model, double percentile, double sample_space_start, double sample_space_end, double sample_space_points)
{
    double sample_space_percentile_value;
    
    //====================== sample space considered
    arma::vec x{arma::linspace<arma::vec>( sample_space_start, sample_space_end, sample_space_points)};
    
    // ====================== pdf computation
    arma::vec pdf_model(x.size());
    double mean, variance, weight;
    // iterate though the sample space
    for (int i=0; i<x.size(); i++)
    {
        // initialize pdf(i)
        pdf_model(i) = 0;
        //iterate through gaussians in model
        for ( int j = 0; j < model.n_gaus(); j++)
        {
            mean = model.means[j];
            variance = model.dcovs[j];
            weight = model.hefts[j];
            //cumulative sum for each component in the model
            pdf_model(i) = pdf_model(i) + weight * (1/(std::sqrt(2 * arma::datum::pi * variance))) * std::exp( - std::pow((x(i) - mean) , 2) / (2 * variance) );
        }
    }
    
    // ====================== cdf computation
    arma::vec cdf_model{ cdf(pdf_model) };
    
    // ====================== percentile computation
    arma::uvec indexes = arma::find( cdf_model < percentile );  //value for which p(X<sample_space_percentile_value) < percentile
    sample_space_percentile_value = x[ indexes.max() ];
    
    
    return sample_space_percentile_value;
}


//================================================
// Schroder_BackwardIntegration
// NOTE: some modification are needed if this has to be a general function
// NOTE: index schereder integration sligthly different that python
// time decay computation: NOT T60
//================================================
double Schroder_BackwardIntegration(arma::vec &envelope, double hop_size, double tx)
{
    double Tx;
    arma::vec Pfit(2);
    
    //take the signal only after the max
    double max = envelope.max();
    arma::uvec index_max = arma::find(envelope==max);
    arma::vec envelope_after_max { envelope( arma::span( index_max[0], envelope.size() - 1 )) };
    
    //================================================ CONTROL
    if ( envelope_after_max.min() < - (5 + tx) )
    {
        // Db inverse-trasformation, ^2 and backward integration
        envelope_after_max = arma::exp10( envelope_after_max / 20);
        envelope_after_max = arma::pow(envelope_after_max, 2);
        
            //armadillo 7 doesn't have reverse function and doesn't have reverse iterator
        arma::vec schroder_integration(envelope_after_max.size());
        arma::vec::iterator schroder_iter = schroder_integration.end() - 1; //iterator pointing to the last element
        arma::vec::iterator envelope_iter = envelope_after_max.end() - 1; //iterator pointing to the last element
        double integration_dep = 0;
        
            // != .begin() - 1      cause if just != .begin() won't consider the first element
        for ( ; schroder_iter != schroder_integration.begin() - 1 && envelope_iter != envelope_after_max.begin() - 1; schroder_iter--, envelope_iter--)
        {
            integration_dep = *envelope_iter + integration_dep;
            *schroder_iter = integration_dep;
        }
        //    schroder_integration.print();
        
        //normalization and Db trasformation
        schroder_integration = schroder_integration / schroder_integration.max();
        schroder_integration = 20* arma::log10(schroder_integration);
        // time vector needed for poly fitting
        arma::vec time { arma::linspace<arma::vec>(0, (schroder_integration.size() - 1) * hop_size, schroder_integration.size())};
        
        //fitting and Tx calculation
        arma::uvec index_start_fitting = arma::find( schroder_integration < - 5 );
        arma::uvec index_end_fitting = arma::find( schroder_integration < - (5 + tx) ); // is gonna be empty if there aren't value below - (5 + tx) => impossible to fit
        
//        schroder_integration.print();
        
        //================================================ CONTROL
        if ( index_end_fitting.size() != 0 ) // => possible to fit
        {
            int start_index = index_start_fitting[0];
            int end_index = index_end_fitting[0];
            // linear fitting
            Pfit = linearFit( time(arma::span(start_index,end_index)), schroder_integration(arma::span(start_index,end_index)) );
            
            // time decay computation: NOT T60
            Tx = - tx / Pfit[0];
        }
        else
        {
             Tx = NULL;
        }
        
    }
    else
    {
        Tx = NULL;
    }

    return Tx;
}


//================================================ Linear Fit
arma::vec linearFit(arma::vec x, arma::vec y)  //VALUE NOT PASSED BY REFERENCE
{
    arma::vec b(2);
    
    // X = [x ,1]
    arma::mat X(x.size(), 1, arma::fill::ones ); //X with just a col of ones
    X.insert_cols(0, x); //add x col
    // oore-penrose pseudoinverse computation
    arma::mat PseudoInv = arma::pinv(X); // PseudoInv = ( X.t() * X )^-1 * X.t()
    
    //normal equation solution b = ( X.t() * X )^-1 * X.t() * y
    b = PseudoInv * y;
    
    return b;
}


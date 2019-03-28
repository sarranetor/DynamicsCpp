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
#include "tools.hpp"

//==============================================================================

// function declaration ========================================================

//gate threshold computation
double GateThreshold(arma::gmm_diag model);


int main (int argc, char* argv[])
{

    //==================================================
    // AUDIO FILE LOAD
    //==================================================
    
//    File file("/Users/Livio.Saracino/Developer/DunningKruger1_10sComp.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/MoosSnareTop1.wav");
    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/CommitmentSnareTop1.wav");
//    File file("/Users/Livio.Saracino/Developer/Compressor research/tracks/Snare/segments/KillerQueen1.wav");
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
    
    
    //==================================================
    // ENVELOPE COMPUTATION
    //==================================================
    // envelop settings
    double rms_length = 10;
    double overlap = 0.8;
    // rms slow settings
    double rms_length_slow = 200;
    double step = rms_length * (1 - overlap); //hop size
    double overlap_slow = ((rms_length * (1 - overlap)) - rms_length_slow) / (- rms_length_slow);
    // envelope struct declaration
    envelope_type envelopeRms, envelopePeak, envelopeRms_slow;
    
    envelopeRms = Envelope(inSignal, _fs, rms_length, overlap, "rms"); //rms fast
    envelopeRms_slow = Envelope(inSignal, _fs, rms_length_slow, overlap_slow, "rms"); //rms slow
    
//    envelopeRms_slow.envelope.print();
    
    envelopePeak = Envelope(inSignal, _fs, rms_length, overlap, "peak");

    
    
    //==================================================
    // GMM + AIC [peak and rms]
    //==================================================
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
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Rms{GMM_AIC(rmsdb, 7)};
//    std::cout << "\nmeans" << std::endl;
//    model_Rms.means.print();
//    std::cout << "conv" << std::endl;
//    model_Rms.dcovs.print();
//    std::cout << "wheigth" << std::endl;
//    model_Rms.hefts.print();
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Peak{GMM_AIC(peak_db, 7)};
//    std::cout << "\nmeans" << std::endl;
//    model_Peak.means.print();
//    std::cout << "conv" << std::endl;
//    model_Peak.dcovs.print();
//    std::cout << "wheigth" << std::endl;
//    model_Peak.hefts.print();
    
    
    //==================================================
    // peak gate threshld
    //==================================================

    double gate_threshold{ GateThreshold(model_Peak) };
    std::cout << "\ngate_threshold: " << gate_threshold << std::endl;
    
    //==================================================
    // microdynamic test
    //==================================================
    
    int hystereisi_high_threshold = 1;
    int hystereisi_low_threshold = -1;
    double start = std::round( (rms_length_slow/2 - rms_length/2) / (rms_length * (1 - overlap)) );
    
    arma::vec DiffRms{ envelopeRms.envelope(arma::span(start, envelopeRms_slow.envelope.size()-1+start)) - envelopeRms_slow.envelope };
    
//    DiffRms.print();
    
    //flag variables
    int dep = 0; // 0: below threshold         1: above threshold
    //set flag variable
    if ( DiffRms[0] < 0 )
        dep = -1;
    // transient event [index in rms fast, type event]  type event: ascendent[1]/descendent[-1]
    std::vector<int> transient(2);
    // transient events list
    std::list< std::vector<int> > transient_list;
    // iteration among DiffRms values
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
    
    //==================================================
    // gate timing computations [Tr and Th] -- using PEAK ENVELOPE
    //==================================================
    
    int offset = 1; // ..
    
    //create transient list positive[ascendent] and transient list negative[descendent]
    std::list< std::vector<int> > transient_list_ascendent;
    std::list< std::vector<int> > transient_list_descendent;
    
    for ( std::vector<int> l : transient_list)
    {
        if ( l[1] == 1 )
            transient_list_ascendent.push_back(l);
        else if ( l[1] == -1 )
            transient_list_descendent.push_back(l);
    }
    
    //event discriminator btw main event and spill event -- INDEXes CORRESPOND TO RMS FAST OR PEAK ENVELOPE
    std::list< std::vector<int> > main_event_list; // [ [index_event_start, index_event_end] , ... ]
    std::list< std::vector<int> > spill_event_list; // [ [index_event_start, index_event_end] , ... ]
    std::list< std::vector<int> > event_list; // [ [index_event_start, index_event_end, label] , ... ] // label==1 [main]     label==0 [spill]
    
    // transient_list_ascendent.size should be even. iter transient_list_ascendent till the second to last element.
    std::list<std::vector<int>>::iterator iter;
    int begin_event_index, end_event_index;
    std::vector<int> peak_event_indexes(2);
    std::vector<int> peak_event_indexes_withlabel(3);
    
    for ( iter = transient_list_ascendent.begin(); iter != std::prev(transient_list_ascendent.end()); iter++)
    {
        begin_event_index = start + (*iter)[0] - offset; //start + cause index are referred to rms slow
        iter++; //cause it is a biderectional iterator
        end_event_index = start + (*iter)[0] - offset;
        iter--;
        // ..
        peak_event_indexes[0] = begin_event_index;
        peak_event_indexes[1] = end_event_index;
        // ..
        peak_event_indexes_withlabel[0] = begin_event_index;
        peak_event_indexes_withlabel[1] = end_event_index;
        
        // if event max > gate threshold (=> main)
        if ( envelopePeak.envelope(arma::span( begin_event_index, end_event_index )).max() > gate_threshold )
        {
            main_event_list.push_back(peak_event_indexes);
            
            peak_event_indexes_withlabel[2] = 1;
            event_list.push_back(peak_event_indexes_withlabel);
        }
        // if event max < gate threshold (=> spill)
        else if ( envelopePeak.envelope(arma::span( begin_event_index, end_event_index )).max() < gate_threshold )
        {
            spill_event_list.push_back(peak_event_indexes);
            
            peak_event_indexes_withlabel[2] = 0;
            event_list.push_back(peak_event_indexes_withlabel);
        }
    }
    
    //================================================== Ioi_main_spill computation
    
    arma::vec Ioi_main_spill;
    int counter_row = 0;
    double time_start_main, time_start_spill;
    double Ioi_main_spill_mean;
    
    for ( iter = event_list.begin(); iter != std::prev(event_list.end()); iter++ )
    {
        // if current event is main and next is spill
//        std::cout << (*iter)[2] << std::endl;
//        std::cout << (*std::next(iter))[2] << std::endl;
        
        if ( (*iter)[2] == 1 && (*std::next(iter))[2] == 0 )
        {
            //time start current event: main
            time_start_main = envelopePeak.time[ (*iter)[0] ];
            //time start next event: spill
            time_start_spill = envelopePeak.time[ (*std::next(iter))[0] ];
            
            // append the data to the arma::vec
            arma::rowvec row{time_start_spill - time_start_main}; //create a new row to append
            Ioi_main_spill.insert_rows(counter_row, row); //append
            counter_row++; //increment pointer to row
        }
    }
    
//    Ioi_main_spill.print();
    Ioi_main_spill_mean = arma::mean(Ioi_main_spill);
    
    
    //================================================== schroder
    
    double Tx;
    arma::vec envelop_one_event;
    arma::vec Thold_vector;
    
//    int s = main_event_list.;
//    int e = (main_event_list[0])[1];
    
    counter_row = 0;
    for ( iter = main_event_list.begin(); iter != main_event_list.end(); iter++)
    {
        envelop_one_event = envelopePeak.envelope( arma::span( (*iter)[0], (*iter)[1] ) );
//        envelop_one_event.print();
        Tx = Schroder_BackwardIntegration( envelop_one_event , step, 30); //if it fails Tx == NULL [0]  // Tx in [ms]
        
        if ( Tx != 0)
        {
            arma::rowvec row{ Tx }; //create a new row to append
            Thold_vector.insert_rows(counter_row, row); //append
            counter_row++; //increment pointer to row
        }
        
    }
    
    //================================================== Th and Tr
    double Th = arma::mean(Thold_vector);
    double Tr = Ioi_main_spill_mean * 1e3 - Th;
    
    std::cout << "Thold: " << Th << std::endl;
    std::cout << "Trelease: " << Tr << std::endl;
    
    return 0;
}




//================================================
// CUSTOM ALGO ANALYSIS FUNCTIONS
//================================================

double GateThreshold(arma::gmm_diag model)
{
    //find index for mean,cov, and weigth for second loudest gaussian
    arma::uvec indexMinorThanLoudest = arma::find( model.means < model.means.max() );
    arma::uvec indexSecondLoudest = arma::find( model.means == model.means.elem(indexMinorThanLoudest).max() );
    
    arma::vec SecondLaudestGaussian_Params{model.means[indexSecondLoudest[0]], model.dcovs[indexSecondLoudest[0]], model.hefts[indexSecondLoudest[0]]};

    // value that the process could reach in db
    arma::vec dbAx{arma::linspace<arma::vec>(-100, 0, 1000)}; // 1000 elements, step 0.1 db
    
    //  pdf single gaussian calculation
    arma::vec pdf_secondlaudest{pdf( model.means[indexSecondLoudest[0]] , model.dcovs[indexSecondLoudest[0]], dbAx)};

    // cdf
    arma::vec cdf_secondlaudest{ cdf(pdf_secondlaudest) };

    // percentile computation
    arma::uvec indexes = arma::find( cdf_secondlaudest < 0.98 );  //value for which p(X<gate_thershold)<0.98
    double gate_thershold = dbAx[ indexes.max() ];
    
    return gate_thershold;
}

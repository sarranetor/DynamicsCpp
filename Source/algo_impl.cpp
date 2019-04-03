//
//  algo_impl.cpp
//  NewProject - ConsoleApp
//
//  Created by Livio Saracino on 29/03/2019.
//

#include "algo_impl.hpp"
//#include "tools.hpp"

#define COMP_RATIO 1.0/3.0

//====================================================================================================
// SNARE / KICK GATE AND COMPRESSOR SMART PRESET
//====================================================================================================

double GateThreshold(arma::gmm_diag &model)
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


std::map<std::string, double> snare_kick_dynamics(arma::vec &inSignal, double &_fs)
{
    std::map<std::string, double> m;
    
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
    //    std::cout << "\ngate_threshold: " << gate_threshold << std::endl;
    
    //==================================================
    // microdynamics
    //==================================================
    
    double start;
    std::list< std::vector<int> > transient_list;
    
    transient_list = microdynamics(envelopeRms, envelopeRms_slow, rms_length, rms_length_slow, overlap, 1, -1, start);
    
    //==================================================
    // gate timing computations [Tr and Th] -- using PEAK ENVELOPE
    //==================================================
    
    //================================================== Event labelling
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
    
    double Tx_hold, Tx_attack;
    arma::vec envelop_one_event;
    arma::vec Thold_vector;
    arma::vec Tattack_vector;
    
    //    int s = main_event_list.;
    //    int e = (main_event_list[0])[1];
    
    int counter_row_hold = 0;
    int counter_row_attack = 0;
    
    for ( iter = main_event_list.begin(); iter != main_event_list.end(); iter++)
    {
        envelop_one_event = envelopePeak.envelope( arma::span( (*iter)[0], (*iter)[1] ) );
        
        //=============== computation of Thold gate == Trelease compressor
        Tx_hold = Schroder_BackwardIntegration( envelop_one_event , step, 30); //if it fails Tx == NULL [0]  // Tx in [ms]
        if ( Tx_hold > 0.0)
        {
            arma::rowvec row{ Tx_hold }; //create a new row to append
            Thold_vector.insert_rows(counter_row_hold, row); //append
            counter_row_hold++; //increment pointer to row
        }
        
        //=============== computation of Tattack compressor
        Tx_attack = Schroder_BackwardIntegration( envelop_one_event , step, 10); //if it fails Tx == NULL [0]  // Tx in [ms]
        if ( Tx_attack > 0.0)
        {
            arma::rowvec row{ Tx_attack }; //create a new row to append
            Tattack_vector.insert_rows(counter_row_attack, row); //append
            counter_row_attack++; //increment pointer to row
        }
    }
    
    //================================================== Th and Tr
    double Th_gate = arma::mean(Thold_vector);
    double Tr_gate = Ioi_main_spill_mean * 1e3 - Th_gate;
    
    //    std::cout << "Thold: " << Th_gate << std::endl;
    //    std::cout << "Trelease: " << Tr_gate << std::endl;
    
    
    //==================================================
    // COMPRESSION SETTINGS
    //==================================================
    
    double Tr_compressor = Th_gate;
    double Ta_compressor = arma::mean(Tattack_vector);
    
    double MaxLevel_peak = getPercentile(model_Peak, 0.998, -100, 0, 1000);
    
    double threshold = MaxLevel_peak - 30;
    double ratio = 1.0/3.0;
    
    double MaxAfterComp = threshold + (MaxLevel_peak-threshold)*ratio;
    double MakeUPgain = MaxLevel_peak - MaxAfterComp;
    
    if ( MakeUPgain > std::abs(MaxLevel_peak) )
    {
        MakeUPgain = std::abs(MaxLevel_peak);
    }
    
    //fill map
    m.insert(std::pair<std::string, double>("gate_threshold", gate_threshold));
    m.insert(std::pair<std::string, double>("gate_hold", Th_gate));
    m.insert(std::pair<std::string, double>("gate_release", Tr_gate));
    
    m.insert(std::pair<std::string, double>("compressor_threshold", threshold));
    m.insert(std::pair<std::string, double>("compressor_attack", Ta_compressor));
    m.insert(std::pair<std::string, double>("compressor_release", Tr_compressor));
    m.insert(std::pair<std::string, double>("compressor_gain", MakeUPgain));
    
    return m;
}




//====================================================================================================
// VOICE COMPRESSOR SMART PRESET
//====================================================================================================

std::map<std::string, double> voice_dynamics(arma::vec &inSignal, double &_fs)
{
    std::map<std::string, double> m;
    
    double release_Dbrange_transient = 10;
    double attack_Dbrange_transient = 7;
    double MaxTimeDecay = 0.08;
    double MaxTimeAttack = 0.03;
    
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
//        std::cout << "\nmeans" << std::endl;
//        model_Rms.means.print();
//        std::cout << "conv" << std::endl;
//        model_Rms.dcovs.print();
//        std::cout << "wheigth" << std::endl;
//        model_Rms.hefts.print();
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Peak{GMM_AIC(peak_db, 7)};
    //    std::cout << "\nmeans" << std::endl;
    //    model_Peak.means.print();
    //    std::cout << "conv" << std::endl;
    //    model_Peak.dcovs.print();
    //    std::cout << "wheigth" << std::endl;
    //    model_Peak.hefts.print();
    
    
    //==================================================
    // compressor threshold
    //==================================================
    
    double ratio = COMP_RATIO;
    double MaxLevel_rms = getPercentile(model_Rms, 0.998, -100, 0, 1000);
    double MaxLevel_peak = getPercentile(model_Peak, 0.998, -100, 0, 1000);
    
    /*
     the voice db range is considered to be 30 db thus [MaxLevel_rms, MaxLevel_rms - 30].
     this range is divided in 4 areas each one of 30/4 = 7.5 db, and these areas are labbeled as:
     very loud, loud, soft, very soft.
     */
    
    // find first gaussian in with loud label
    arma::uvec index_not_veryloud = arma::find( model_Rms.means < MaxLevel_rms - 30/4 );
    arma::uvec index_loud_gaussian = arma::find( model_Rms.means == model_Rms.means.elem(index_not_veryloud).max() );
    
    // value that the process could reach in db
    arma::vec x{arma::linspace<arma::vec>(-100, 0, 1000)}; // 1000 elements, step 0.1 db
    //  pdf single gaussian calculation
    arma::vec pdf_loud_gaussian{pdf( model_Rms.means[index_loud_gaussian[0]] , model_Rms.dcovs[index_loud_gaussian[0]], x)};
    // cdf
    arma::vec cdf_loud_gaussian{ cdf(pdf_loud_gaussian) };
    // percentile computation
    arma::uvec indexes = arma::find( cdf_loud_gaussian < 0.99 );  //value for which p(X<gate_thershold)<0.98
    
    double MaxAfterComp = x[ indexes.max() ];
    double threshold = (MaxAfterComp - MaxLevel_rms * ratio) / (1 - ratio);
    double MakeUPgain = MaxLevel_rms - MaxAfterComp;
    
    if ( MakeUPgain > std::abs(MaxLevel_peak) )
    {
        MakeUPgain = std::abs(MaxLevel_peak);
    }
    
    //==================================================
    // compressor timing computations [Tr and Ta] -- using rms envelope
    //==================================================
    
    //================================================== microdynamics
    double start;
    std::list< std::vector<int> > transient_list;
    
    transient_list = microdynamics(envelopeRms, envelopeRms_slow, rms_length, rms_length_slow, overlap, 1, -1, start);
    
    //================================================== transient evaluation
    arma::vec Trelease_vector , Tattack_vector;
    int counter_row_release = 0; //each time the condition are respected a row is added to Trelease_vector
    int counter_row_attack = 0; //each time the condition are respected a row is added to Tattack_vector
    double end_index = start + envelopeRms_slow.envelope.size(); //last elem to consider in envelopeRms [after this rms slow  doesn't exist]
    
    //=========== attack calculation variable
    double attack_time; // store computed attack before to add in Tattack_vector
    double start_event_rms_value; // rms value where the attack start
    double Dbrange_transient_attack;
    
    //=========== release decay calculation variable
    double decay_time; // store computed decay before to add in Trelease_vector
    double end_event_rms_value; // rms value where the decay start
    double Dbrange_transient_release;
    
    for ( std::vector<int> l : transient_list)
    {
        if ( l[1] == 1 ) // == 1 [ascendent]
        {
            start_event_rms_value = envelopeRms.envelope[ start + l[0] ];
            Dbrange_transient_attack = start_event_rms_value - attack_Dbrange_transient;
            
            // index where the signal goes to Dbrange_transient_release. first index.
            arma::uvec index = arma::find(  envelopeRms.envelope(arma::span(0, start + l[0])) < Dbrange_transient_attack );
            
            attack_time = envelopeRms.time[start + l[0]] - envelopeRms.time[ *(index.end()-1) ];
            
            // the attack time is considered only if is less that a max
            if ( attack_time <= MaxTimeAttack )
            {
                arma::rowvec row { attack_time };
                Tattack_vector.insert_rows(counter_row_attack, row);
                counter_row_attack++; // increase row counter
            }
            
        }
        else if ( l[1] == -1 ) // ==-1 [descendent]
        {
            end_event_rms_value = envelopeRms.envelope[ start + l[0] ];
            Dbrange_transient_release = end_event_rms_value - release_Dbrange_transient;
            
            // index where the signal goes to Dbrange_transient_release. first index.
            arma::uvec index = arma::find(  envelopeRms.envelope(arma::span(start + l[0], end_index)) < Dbrange_transient_release );
            
            // start +, cause the transient index refers to the Diff btw fast and slow. but the envelopeRms.time[start] = envelopeRms_slow.time[0]
            decay_time = envelopeRms.time[ start + l[0] + index[0]] - envelopeRms.time[start + l[0]];
            
            // the decay time is considered only if is less that a max
            if ( decay_time <= MaxTimeDecay)
            {
                arma::rowvec row { decay_time };
                Trelease_vector.insert_rows(counter_row_release, row);
                counter_row_release++; // increase row counter
            }
        }
    }
    
//    Trelease_vector.print();
//    Tattack_vector.print();
    
    //==================================================
    // Tattack_vector and Trelease_vector statistics and chosing criteria
    //==================================================
    
    //control when Trelease_vector.size() very low
    // ..
    
    arma::mat Trelease_vector_ms { Trelease_vector * 1000 }; // *1000 for [ms]
    Trelease_vector_ms = Trelease_vector_ms.t();
    
    arma::gmm_diag model_Trelease{ GMM_AIC(Trelease_vector_ms, 4) };
//    std::cout << "\nmeans" << std::endl;
//    model_Trelease.means.print();
//    std::cout << "conv" << std::endl;
//    model_Trelease.dcovs.print();
//    std::cout << "weigth" << std::endl;
//    model_Trelease.hefts.print();
    
    arma::mat Tattack_vector_ms { Tattack_vector * 1000 }; // *1000 for [ms]
    Tattack_vector_ms = Tattack_vector_ms.t();
    
    arma::gmm_diag model_Tattack{ GMM_AIC(Tattack_vector_ms, 4) };
//    std::cout << "\nmeans" << std::endl;
//    model_Tattack.means.print();
//    std::cout << "conv" << std::endl;
//    model_Tattack.dcovs.print();
//    std::cout << "weigth" << std::endl;
//    model_Tattack.hefts.print();
    
    // selection criteria
    // ..
    
    double Tr = getPercentile(model_Trelease, 0.5, 0, 200, 1000);
    double Ta = getPercentile(model_Tattack, 0.5, 0, 100, 500);
    // ..
    Tr = Tr * 2; // cause i want to measure t20 and i approximate t20 = t10 * 2
    
    //fill map
    m.insert(std::pair<std::string, double>("comp_threshold", threshold));
    m.insert(std::pair<std::string, double>("comp_ratio", ratio));
    m.insert(std::pair<std::string, double>("comp_makeupgain", MakeUPgain));
    m.insert(std::pair<std::string, double>("comp_attack", Ta));
    m.insert(std::pair<std::string, double>("comp_release", Tr));
    
    return m;
}

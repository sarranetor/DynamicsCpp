//
//  algo_impl.cpp
//  NewProject - ConsoleApp
//
//  Created by Livio Saracino on 29/03/2019.
//

#include "algo_impl.hpp"
//#include "tools.hpp"


// DB RANGE
constexpr int min_db_value = -100;
constexpr int max_db_value = 0;
constexpr int n_db_pointsInRange = 1000;
// envelop settings
constexpr double rms_length = 10; // ms
constexpr double overlap = 0.8;
// rms slow settings
constexpr double rms_length_slow = 200; // ms
constexpr double step = rms_length * (1 - overlap); //hop size [ms]
constexpr double overlap_slow = ((rms_length * (1 - overlap)) - rms_length_slow) / (- rms_length_slow);

constexpr double hysteresis_dbRange= 1;


//====================================================================================================
// SNARE / KICK GATE AND COMPRESSOR SMART PRESET
//====================================================================================================

// foundamental settings
constexpr double ratio_snare_kick = 1.0/3.0;
constexpr double db_decay_holdAndRelease = 30;
constexpr double db_decay_attack = 10;
constexpr double threshold_offset = 30;

/**
    takes the second loudest gaussian of the model and computes a certain percentile and return it.
    @return the gare threshold setting in db
 **/
double _GateThreshold(arma::gmm_diag &model);

/**
    knowing the gate threshold in db end the position of the transient in the signal, it computes and stores the events start and end indexes.
 
    @param
    main_event_vector is filled with start and end indexes of the events above gate threshold.
    spill_event_vector is filled with start and end indexes of the events below gate threshold.
    event_vector is filled with start and end indexes of each events plus a label  1==main  0==spill.
 **/
void _EventLabelling(std::vector< std::array<int,2> > &transient_vector, double gate_threshold, double start, envelope_type &envelopePeak, std::vector< std::array<int,2> > &main_event_vector, std::vector< std::array<int,2> > &spill_event_vector, std::vector< std::array<int,3> > &event_vector);

/**
    computes the inter-onsets-interval btw a main event and the first spill event that follows.
    @return value in seconds
 **/
double _Ioi_main_spill_computation(std::vector<std::array<int,3>> &event_vector, envelope_type &envelopePeak);

/**
    computes the decay time of each event using the schreoder backward integration.
    @param tx Db ..
    @return a vector filled with the decay time [ms] of each event. if this calulation was possible.
 **/
arma::vec _decayEvents(std::vector< std::array<int,2> > &main_event_vector, envelope_type envelopePeak, double tx);


std::map<std::string, double> snare_kick_dynamics(arma::vec &inSignal, double &_fs)
{
    std::map<std::string, double> m;
    
    //================================================== ENVELOPE COMPUTATION
    envelope_type envelopeRms, envelopePeak, envelopeRms_slow;
    envelopeRms = Envelope(inSignal, _fs, rms_length, overlap, "rms"); //rms fast
    envelopeRms_slow = Envelope(inSignal, _fs, rms_length_slow, overlap_slow, "rms"); //rms slow
    envelopePeak = Envelope(inSignal, _fs, rms_length, overlap, "peak");
    
    //================================================== GMM + AIC [peak and rms]
    //data as arma::mat. each col is an element. //deleate the element below -100 db. or sub setting. silence!
    arma::mat rmsdb = ReduceDynamicRange(envelopeRms, min_db_value);
    arma::mat peak_db = ReduceDynamicRange(envelopePeak, min_db_value);
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Rms{ GMM_AIC(rmsdb, 7) };
    arma::gmm_diag model_Peak{ GMM_AIC(peak_db, 7) };
    
    //================================================== peak gate threshld
    double gate_threshold{ _GateThreshold(model_Peak) };

    //================================================== microdynamics
    double start;
    std::vector< std::array<int,2> > transient_vector;
    transient_vector = microdynamics(envelopeRms, envelopeRms_slow, rms_length, rms_length_slow, overlap, hysteresis_dbRange, -hysteresis_dbRange, start);
    
    //================================================== gate timing computations [Tr and Th] -- using PEAK ENVELOPE

    //======Event labelling
    //event discriminator btw main event and spill event -- INDEXes CORRESPOND TO RMS FAST OR PEAK ENVELOPE
    std::vector< std::array<int,2> > main_event_vector; // [ [index_event_start, index_event_end] , ... ]
    std::vector< std::array<int,2> > spill_event_vector; // [ [index_event_start, index_event_end] , ... ]
    std::vector< std::array<int,3> > event_vector; // [ [index_event_start, index_event_end, label] , ... ] // label==1 [main]     label==0 [spill]
    
    _EventLabelling(transient_vector, gate_threshold, start, envelopePeak, main_event_vector, spill_event_vector, event_vector);
    
    //================================================== Ioi_main_spill computation
    double Ioi_main_spill_mean = _Ioi_main_spill_computation(event_vector, envelopePeak);
    
    //================================================== schroder integration for decay computation
    arma::vec Thold_vector = _decayEvents(main_event_vector, envelopePeak,  db_decay_holdAndRelease);
    arma::vec Tattack_vector = _decayEvents(main_event_vector, envelopePeak,  db_decay_attack);
    
    //================================================== Th and Tr gate
    double Th_gate = arma::mean(Thold_vector);
    double Tr_gate = Ioi_main_spill_mean * 1e3 - Th_gate; // * 1e3 cause Ioi_main_spill_mean is in seconds
    
    //================================================== COMPRESSION SETTINGS
    double Tr_compressor = Th_gate;
    double Ta_compressor = arma::mean(Tattack_vector);
    
    double MaxLevel_peak = getPercentile(model_Peak, 0.998, min_db_value, max_db_value, n_db_pointsInRange);
    
    double threshold = MaxLevel_peak - threshold_offset;
    double ratio = ratio_snare_kick;
    
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
    m.insert(std::pair<std::string, double>("compressor_ratio", ratio));
    m.insert(std::pair<std::string, double>("compressor_attack", Ta_compressor));
    m.insert(std::pair<std::string, double>("compressor_release", Tr_compressor));
    m.insert(std::pair<std::string, double>("compressor_gain", MakeUPgain));
    
    return m;
}


double _GateThreshold(arma::gmm_diag &model)
{
    //find index for mean,cov, and weigth for second loudest gaussian
    arma::uvec indexMinorThanLoudest = arma::find( model.means < model.means.max() );
    arma::uvec indexSecondLoudest = arma::find( model.means == model.means.elem(indexMinorThanLoudest).max() );
    
    arma::vec SecondLaudestGaussian_Params{model.means[indexSecondLoudest[0]], model.dcovs[indexSecondLoudest[0]], model.hefts[indexSecondLoudest[0]]};
    
    // value that the process could reach in db
    arma::vec dbAx{arma::linspace<arma::vec>(min_db_value, max_db_value, n_db_pointsInRange)}; // 1000 elements, step 0.1 db
    
    double gate_thershold = getGausPercentile(0.98, model.means[indexSecondLoudest[0]], model.dcovs[indexSecondLoudest[0]], dbAx);
    
    return gate_thershold;
}


void _EventLabelling(std::vector< std::array<int,2> > &transient_vector, double gate_threshold, double start, envelope_type &envelopePeak, std::vector< std::array<int,2> > &main_event_vector, std::vector< std::array<int,2> > &spill_event_vector, std::vector< std::array<int,3> > &event_vector)
{
    constexpr int MainLabel = 1, SpillLabel = 0;
    int offset = 1; // how many samples before the transient detection i consider the event startiong
    
    //create transient list positive[ascendent] and transient list negative[descendent]
    std::vector< std::array<int,2> > transient_vector_ascendent;
    std::vector< std::array<int,2> > transient_vector_descendent;
    
    for ( std::array<int,2> l : transient_vector)
    {
        if ( l[1] == 1 )
            transient_vector_ascendent.push_back(l);
        else if ( l[1] == -1 )
            transient_vector_descendent.push_back(l);
    }

    
    // transient_list_ascendent.size should be even. iter transient_list_ascendent till the second to last element.
    int begin_event_index, end_event_index;
    std::array<int,2> peak_event_indexes;
    std::array<int,3> peak_event_indexes_withlabel;
    
    for ( std::vector< std::array<int,2> >::iterator iter = transient_vector_ascendent.begin(); iter != std::prev(transient_vector_ascendent.end()); iter++)
    {
        begin_event_index = start + (*iter)[0] - offset; //start + cause index are referred to rms slow
        end_event_index = start + ( *(iter+1) )[0] - offset;
        // array filling before being pushed in vector
        peak_event_indexes[0] = begin_event_index;
        peak_event_indexes[1] = end_event_index;
        peak_event_indexes_withlabel[0] = begin_event_index;
        peak_event_indexes_withlabel[1] = end_event_index;
        
        // if event max > gate threshold (=> main)
        if ( envelopePeak.envelope(arma::span( begin_event_index, end_event_index )).max() > gate_threshold )
        {
            main_event_vector.push_back(peak_event_indexes);
            peak_event_indexes_withlabel[2] = MainLabel;
            event_vector.push_back(peak_event_indexes_withlabel);
        }
        // if event max < gate threshold (=> spill)
        else if ( envelopePeak.envelope(arma::span( begin_event_index, end_event_index )).max() < gate_threshold )
        {
            spill_event_vector.push_back(peak_event_indexes);
            peak_event_indexes_withlabel[2] = SpillLabel;
            event_vector.push_back(peak_event_indexes_withlabel);
        }
    }
    
}


arma::vec _decayEvents(std::vector< std::array<int,2> > &event_vector, envelope_type envelopePeak, double tx)
{
    int counter_row = 0;
    double Tx;
    arma::vec envelop_one_event;
    arma::vec Tx_vector;
    
    for ( std::vector<std::array<int,2>>::iterator iter = event_vector.begin(); iter != event_vector.end(); iter++)
    {
        envelop_one_event = envelopePeak.envelope( arma::span( (*iter)[0], (*iter)[1] ) );
        
        //=============== computation of Thold gate == Trelease compressor
        Tx = Schroder_BackwardIntegration( envelop_one_event , step, tx); //if it fails Tx == NULL [0]  // Tx in [ms]
        if ( Tx > 0.0)
        {
            arma::rowvec row{ Tx }; //create a new row to append
            Tx_vector.insert_rows(counter_row, row); //append
            counter_row++; //increment pointer to row
        }
    }
    
    return Tx_vector;
}


double _Ioi_main_spill_computation(std::vector<std::array<int,3>> &event_vector, envelope_type &envelopePeak)
{
    arma::vec Ioi_main_spill;
    int counter_row = 0;
    double time_start_main, time_start_spill;
    double Ioi_main_spill_mean;
    
    for ( std::vector<std::array<int,3>>::iterator iter = event_vector.begin(); iter != std::prev(event_vector.end()); iter++ )
    {
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
    return Ioi_main_spill_mean = arma::mean(Ioi_main_spill);
}



//====================================================================================================
// VOICE COMPRESSOR SMART PRESET
//====================================================================================================

// foundamental settings
constexpr double release_Dbrange_transient = 10;
constexpr double attack_Dbrange_transient = 7;
constexpr double MaxTimeDecay = 0.08;
constexpr double MaxTimeAttack = 0.03;
constexpr double ratio_voice = 1.0/2.0;

/**
    fills Trelease_vector and Tattack_vector with the time in ms that the transient needed to move of release_Dbrange_transient db and attack_Dbrange_transient db.
 **/
void _transientEvaluation (std::vector< std::array<int,2> > &transient_vector, envelope_type &envelopeRms, envelope_type &envelopeRms_slow, arma::vec &Trelease_vector, arma::vec &Tattack_vector, double start);


std::map<std::string, double> voice_dynamics(arma::vec &inSignal, double &_fs)
{
    std::map<std::string, double> m;
    
    //================================================== ENVELOPE COMPUTATION
    envelope_type envelopeRms, envelopePeak, envelopeRms_slow;
    envelopeRms = Envelope(inSignal, _fs, rms_length, overlap, "rms"); //rms fast
//    envelopeRms.envelope.print();
    envelopeRms_slow = Envelope(inSignal, _fs, rms_length_slow, overlap_slow, "rms"); //rms slow
    envelopePeak = Envelope(inSignal, _fs, rms_length, overlap, "peak");
    
    //================================================== GMM + AIC [peak and rms]
    //data as arma::mat. each col is an element. //deleate the element below -100 db. or sub setting. silence!
    arma::mat rmsdb = ReduceDynamicRange(envelopeRms, min_db_value);
    arma::mat peak_db = ReduceDynamicRange(envelopePeak, min_db_value);
    
    // GMM + AIC -----------------------------------
    arma::gmm_diag model_Rms{GMM_AIC(rmsdb, 7)};
    arma::gmm_diag model_Peak{GMM_AIC(peak_db, 7)};
    
    //================================================== compressor threshold and makeupgain
    double ratio = ratio_voice;
    double MaxLevel_rms = getPercentile(model_Rms, 0.998, min_db_value, max_db_value, n_db_pointsInRange);
    double MaxLevel_peak = getPercentile(model_Peak, 0.998, min_db_value, max_db_value, n_db_pointsInRange);
    
    /*
     the voice db range is considered to be 30 db thus [MaxLevel_rms, MaxLevel_rms - 30].
     this range is divided in 4 areas each one of 30/4 = 7.5 db, and these areas are labbeled as:
     very loud, loud, soft, very soft.
     */
    
    // find first gaussian in with loud label
    arma::uvec index_not_veryloud = arma::find( model_Rms.means < MaxLevel_rms - 30.0/4.0 );
    arma::uvec index_loud_gaussian = arma::find( model_Rms.means == model_Rms.means.elem(index_not_veryloud).max() );
    
    // value that the process could reach in db
    arma::vec x{arma::linspace<arma::vec>(min_db_value, max_db_value, n_db_pointsInRange)}; // 1000 elements, step 0.1 db
    
    //value for which p(X<gate_thershold)<0.98
    double MaxAfterComp = getGausPercentile(0.99, model_Rms.means[index_loud_gaussian[0]], model_Rms.dcovs[index_loud_gaussian[0]], x);
    
    double threshold = (MaxAfterComp - MaxLevel_rms * ratio) / (1 - ratio);
    double MakeUPgain = MaxLevel_rms - MaxAfterComp;
    
    if ( MakeUPgain > std::abs(MaxLevel_peak) )
    {
        MakeUPgain = std::abs(MaxLevel_peak);
    }
    
    //================================================== compressor timing computations [Tr and Ta] -- using rms envelope

    //================================================== microdynamics
    double start;
    std::vector< std::array<int,2> > transient_vector;
    
    transient_vector = microdynamics(envelopeRms, envelopeRms_slow, rms_length, rms_length_slow, overlap, hysteresis_dbRange, -hysteresis_dbRange, start);
    
    //================================================== transient evaluation
    arma::vec Trelease_vector , Tattack_vector;
    _transientEvaluation(transient_vector, envelopeRms, envelopeRms_slow, Trelease_vector, Tattack_vector, start);
 
    //================================================== Tattack_vector and Trelease_vector statistics and chosing criteria
    double Tr = 0, Ta = 0;
    
    //control when Trelease_vector.size()/Tattack_vector.size() very low => not possible to compute GMM
    if ( Trelease_vector.size() > 1 )
    {
        arma::mat Trelease_vector_ms { Trelease_vector * 1000 }; // *1000 for [ms] conversion
        Trelease_vector_ms = Trelease_vector_ms.t();
        arma::gmm_diag model_Trelease{ GMM_AIC(Trelease_vector_ms, 4) };
        // selection criteria
        Tr = getPercentile(model_Trelease, 0.5, 0, 200, 1000);
        Tr = Tr * 2; // cause i want to measure t20 and i approximate t20 = t10 * 2 in order to have more data
    }
    else
    {
        Tr = arma::mean(Trelease_vector);
    }
    
    //control when Trelease_vector.size()/Tattack_vector.size() very low => not possible to compute GMM
    if ( Tattack_vector.size() > 1 )
    {
        arma::mat Tattack_vector_ms { Tattack_vector * 1000 }; // *1000 for [ms] conversion
        Tattack_vector_ms = Tattack_vector_ms.t();
        arma::gmm_diag model_Tattack{ GMM_AIC(Tattack_vector_ms, 4) };
        // selection criteria
        Ta = getPercentile(model_Tattack, 0.5, 0, 100, 500);
    }
    else
    {
        Ta = arma::mean(Tattack_vector);
    }
    
    //fill map
    m.insert(std::pair<std::string, double>("comp_threshold", threshold));
    m.insert(std::pair<std::string, double>("comp_ratio", ratio));
    m.insert(std::pair<std::string, double>("comp_makeupgain", MakeUPgain));
    m.insert(std::pair<std::string, double>("comp_attack", Ta));
    m.insert(std::pair<std::string, double>("comp_release", Tr));
    
    return m;
}


void _transientEvaluation (std::vector< std::array<int,2> > &transient_vector, envelope_type &envelopeRms, envelope_type &envelopeRms_slow, arma::vec &Trelease_vector, arma::vec &Tattack_vector, double start)
{
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
    
    for ( std::array<int,2> l : transient_vector)
    {
        if ( l[1] == 1 ) // == 1 [rising]
        {
            start_event_rms_value = envelopeRms.envelope[ start + l[0] ];
            Dbrange_transient_attack = start_event_rms_value - attack_Dbrange_transient;
            
            // index where the signal goes to Dbrange_transient_release. first index.
            arma::uvec index = arma::find(  envelopeRms.envelope(arma::span(0, start + l[0])) < Dbrange_transient_attack );
            
            if ( index.size() > 0 )
            {
                // start +, cause the transient index refers to the Diff btw fast and slow. but the envelopeRms.time[start] = envelopeRms_slow.time[0]
                attack_time = envelopeRms.time[start + l[0]] - envelopeRms.time[ *(index.end()-1) ];
                // the attack time is considered only if is less that a max
                if ( attack_time <= MaxTimeAttack )
                {
                    arma::rowvec row { attack_time };
                    Tattack_vector.insert_rows(counter_row_attack, row);
                    counter_row_attack++; // increase row counter
                }
            }
        }
        else if ( l[1] == -1 ) // == -1 [descendent]
        {
            end_event_rms_value = envelopeRms.envelope[ start + l[0] ];
            Dbrange_transient_release = end_event_rms_value - release_Dbrange_transient;
            
            // index where the signal goes to Dbrange_transient_release. first index.
            arma::uvec index = arma::find(  envelopeRms.envelope(arma::span(start + l[0], end_index)) < Dbrange_transient_release );
            
            if ( index.size() > 0 )
            {
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
    }
    
    std::cout << "release" << std::endl;
    Trelease_vector.print();
    std::cout << "attack" << std::endl;
    Tattack_vector.print();
}

#include <boost/make_shared.hpp>

#include "euter/objectstore.h"

// Population related
#include "euter/population.h"
#include "euter/population_view.h"
#include "euter/celltypes.h"
#include "pycellparameters/pyparameteraccess.h"
//#include "euter/nativerandomgenerator.h"

// Connector
#include "euter/alltoallconnector.h"
#include "euter/connector.h"
#include "euter/fixednumberpreconnector.h"
#include "euter/fixedprobabilityconnector.h"
#include "euter/fromlistconnector.h"
#include "euter/nativerandomgenerator.h"
#include "euter/onetooneconnector.h"
#include "euter/projection.h"

//logger
#include <log4cxx/basicconfigurator.h>

// marocco
#include "pymarocco/PyMarocco.h"
#include "pymarocco/runtime/Runtime.h"
#include "sthal/Wafer.h"
#include "sthal/HICANN.h"
#include "hal/Coordinate/iter_all.h"

// submit experiments to the wafer
#include "submit.h"

//#include "boost/serialization/serialization.h"
//#include "boost/serialization/ublas.hpp"

using boost::make_shared;
namespace C = HMF::Coordinate;
namespace D = halco::hicann::v2;

/**
 * Setting low lever parameter for the BrainScaleS System
 * 
 * @param wafer pointer to wafer object
 * @param gmax set to 1023
 * @param gmax_div set to 1
 */
void set_stahl_params(boost::shared_ptr<sthal::Wafer> wafer, double gmax, double gmax_div){
    for(auto hicann: wafer->getAllocatedHicannCoordinates()){
        auto fgs = (*wafer)[hicann].floating_gates;
        
        for(auto block: C::iter_all<C::FGBlockOnHICANN>()){
            fgs.setShared(block, HMF::HICANN::shared_parameter::V_gmax0, gmax);
            fgs.setShared(block, HMF::HICANN::shared_parameter::V_gmax1, gmax);
            fgs.setShared(block, HMF::HICANN::shared_parameter::V_gmax2, gmax);
            fgs.setShared(block, HMF::HICANN::shared_parameter::V_gmax3, gmax);
        }
        
        for(auto driver : C::iter_all<D::SynapseDriverOnHICANN>()){
            for(auto row : C::iter_all<C::RowOnSynapseDriver>()){
                (*wafer)[hicann].synapses[driver][row].set_gmax_div(C::left, gmax_div);
            }
        }
        
        for(size_t i = 0; i<fgs.getNoProgrammingPasses(); i++){
            auto cfg = fgs.getFGConfig(C::Enum(i));
            cfg.fg_biasn = 0;
            cfg.fg_bias = 0;
            fgs.setFGConfig(C::Enum(i), cfg);
        }
        for(auto block: C::iter_all<C::FGBlockOnHICANN>()){
            fgs.setShared(block, HMF::HICANN::shared_parameter::V_dllres, 275);
            fgs.setShared(block, HMF::HICANN::shared_parameter::V_ccas, 800);
        }
    }
}

/**
 * Init some of the loggers
 */
void init_logger(){log4cxx::BasicConfigurator::resetConfiguration();
	log4cxx::BasicConfigurator::configure();
    auto logger =  log4cxx::Logger::getRootLogger();
    logger->setLevel(log4cxx::Level::getWarn());
    for (auto logger : {"ESS", "marocco", "calibtic", "stah"}){
        auto logger_inst =  log4cxx::Logger::getLogger(logger);
        logger_inst->setLevel(log4cxx::Level::getWarn());
    }
}

int main(int , const char** )
{
    init_logger();
    
    //Random generator used for random connectors
    boost::shared_ptr<RandomGenerator> rng =
	        boost::make_shared<NativeRandomGenerator>(1234);

    ObjectStore store;
    auto p1 = Population::create(store, 2, CellType::SpikeSourceArray);
    
	//P a = P::getDefault(); Default Parameters
    //P::getNames(); Names of parameters
    
    // Set spike times
    auto& params_2 = reinterpret_cast<TypedCellParameterVector<CellType::SpikeSourceArray>&>(p1->parameters());
     for(auto& i : params_2.parameters()){
        i.spike_times = {30.0,40.0};
    }
    
	auto p2 = Population::create(store, 3, CellType::IF_cond_exp);
    auto& params = reinterpret_cast<TypedCellParameterVector<CellType::IF_cond_exp>&>(p2->parameters());
    
    //Record spike times
    for(auto& i : params.parameters()){
        i.record_spikes = true;
    }
    
	double weight = 5;
	double delay  = 5;
    //AllToAllConnector(true, weight, delay);
    //FixedProbabilityConnector(prob, true, weight, delay);
    auto conn  =boost::make_shared<AllToAllConnector>(true, weight, delay);
    auto proj = Projection::create(
		store, p1, p2,conn, rng, "", "excitatory"); // + Synapse_dynamics , label
    
    // ProjectionMatrix weights = proj->getWeights(); Used to get the set weights
    
    
    auto marocco = pymarocco::PyMarocco::create();
    
    marocco->continue_despite_synapse_loss = true;
    
    // Choose between Hardware, ESS, and None
    marocco->backend = pymarocco::PyMarocco::Backend::Hardware;
    
    marocco->calib_backend = pymarocco::PyMarocco::CalibBackend::Binary;
    marocco->neuron_placement.default_neuron_size(4); // denmems per neuron
    
    // Some low-level defaults we might consider to change
    //marocco->neuron_placement.restrict_rightmost_neuron_blocks(false); 
    //marocco->neuron_placement.minimize_number_of_sending_repeaters(false); 
    //marocco->param_trafo.use_big_capacitors = true;  //default true
    //marocco->input_placement.consider_firing_rate(true);
    //marocco->input_placement.bandwidth_utilization(0.8);
    //marocco->calib_path = "/wang/data/calibration/brainscales/default";
    //marocco->defects.backend = pymarocco::Defects::Backend::XML;
    //marocco->defects.path = "/wang/data/calibration/brainscales/default";
    
    // Choose hicann and wafer
    marocco->manual_placement.on_hicann(p1->id(), 
                                        HMF::Coordinate::HICANNOnWafer(
                                            HMF::Coordinate::Enum(297)));//std::vector<HMF::Coordinate::HICANNOnWafer> const& hicanns
    auto runtime = pymarocco::runtime::Runtime::create(HMF::Coordinate::Wafer(33));
    
    // Save marocco to ObjectStore
    ObjectStore::Settings settings;
    ObjectStore::metadata_map metadata;
    metadata["marocco"] = marocco;
    metadata["marocco_runtime"] = runtime;
    store.setup(settings, metadata); //runtime object
    
    // Run mapping, necessary for changing low-level parameters
    marocco->backend = pymarocco::PyMarocco::Backend::None;
    marocco->skip_mapping = false;
    store.run(100);// ms
    submit(store);
    set_stahl_params(runtime->wafer(), 1023,1);
    
    
    marocco->backend = pymarocco::PyMarocco::Backend::Hardware;
    marocco->skip_mapping = true;

    // Set low-level weights
    auto hw_synapses = runtime->results()->synapse_routing.synapses().find(proj->id());
    for(auto hw_syn : hw_synapses){
        auto syn_hand = hw_syn.hardware_synapse();
        auto hicann = syn_hand->toHICANNOnWafer();
        auto proxy = (*runtime->wafer())[hicann].synapses[*syn_hand];
        proxy.weight = HMF::HICANN::SynapseWeight(15);
    }
    
    // Run the emulation
    submit(store);
    
    // Printing spikes in ms
    for(size_t pop = 0; pop != p2->size(); ++pop) {
        for(auto i : p2->getSpikes(pop)){
		std::cout << i *1e3<< ", ";
        }
        std::cout <<std::endl;
	}
	
	// Get back membrane voltage
	/*
     auto const& voltageTrace = p1->getMembraneVoltageTrace(size_t n)
     for (auto const& tv_pair : voltageTrace) {
         data(pos, 1) = std::get<0>(tv_pair) *1e3;  //time
	     data(pos, 2) = std::get<1>(tv_pair); // mem 
     }
     */
    
    store.reset();
    return 0;
}

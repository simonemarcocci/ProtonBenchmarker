include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lEG -pthread -lm -ldl -rdynamic)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})

art_make( BASENAME_ONLY
          LIBRARY_NAME ubana_ProtonBenchmarker
          LIB_LIBRARIES
            ubana_ProtonBenchmarker_Algos
	    ubana_ProtonBenchmarker_Datatypes
	    larcorealg_Geometry
	    larcorealg_GeoAlgo
	    ubobj_UBXSec
	    ${PANDORASDK}
	    ${PANDORAMONITORING}
	    cetlib cetlib_except
  
          MODULE_LIBRARIES 
	    cetlib cetlib_except
	    larpandora_LArPandoraInterface
	    LArPandoraContent
            ubana_ProtonBenchmarker_Algos
	    ubana_ProtonBenchmarker_Datatypes
	    ubobj_UBXSec
            larcore_Geometry_Geometry_service
	    larcorealg_Geometry
	    larcorealg_GeoAlgo
            larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
	    larsim_MCCheater_BackTrackerService_service
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            lardataobj_RecoBase
            lardata_RecoObjects larreco_Deprecated  
            larreco_RecoAlg
            larpandora_LArPandoraInterface
            SignalShapingServiceMicroBooNE_service
            nusimdata_SimulationBase
            ${LARRECO_LIB}
            ${LARDATA_LIB}
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            art_Persistency_Common canvas
	    #canvas_Persistency_Common 
            art_Persistency_Provenance canvas
	    #canvas_Persistency_Provenance
            art_Utilities canvas
	    #canvas_Utilities
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${ROOT_GEOM}
            ${ROOT_XMLIO}
            ${ROOT_GDML}
            ${ROOT_BASIC_LIB_LIST}
  )   

### install_headers()
install_fhicl()
install_source()

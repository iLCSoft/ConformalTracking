#include "ConformalTrackingV2.h"

#include "ParameterParser.h"
#include "Parameters.h"

using namespace lcio;
using namespace marlin;
using namespace std;

// Static instantiation of the processor
ConformalTrackingV2 aConformalTrackingV2;

/*
 
 Pattern recognition code for the CLIC detector, using conformal mapping and cellular automaton
 
 */

ConformalTrackingV2::ConformalTrackingV2() : ConformalTracking("ConformalTrackingV2") {
  // Processor description
  _description = "ConformalTrackingV2 constructs tracks using a combined conformal mapping and cellular automaton approach.";

  // Input collections and parameters
  registerParameters();

  // Parameters for tracking
  registerProcessorParameter("Steps", "Steps", m_rawSteps, m_rawSteps);
}

void ConformalTrackingV2::parseStepParameters() {
  ParameterParser::parseParameters(_stepParameters, m_rawSteps, m_inputTrackerHitCollections);
}  //parseStepParameters

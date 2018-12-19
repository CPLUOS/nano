// 
// This code is from copying 
//   https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/JetMETCorrections/Modules/src/JetResolution.cc
//   in 4d7b6ef commit
// Copied by Byeonghak Ko
// 


#define STANDALONE


#ifndef STANDALONE
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <CondFormats/DataRecord/interface/JetResolutionRcd.h>
#include <CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h>
#else
//#include "JetResolution.h" // MODIFIED
#include "nano/external/interface/JetResolution.h" // MODIFIED
#endif

namespace JMENano {

    JetResolution::JetResolution(const std::string& filename) {
        m_object = std::make_shared<JetResolutionObject>(filename);
    }

    JetResolution::JetResolution(const JetResolutionObject& object) {
        m_object = std::make_shared<JetResolutionObject>(object);
    }

#ifndef STANDALONE
    const JetResolution JetResolution::get(const edm::EventSetup& setup, const std::string& label) {
        edm::ESHandle<JetResolutionObject> handle;
        setup.get<JetResolutionRcd>().get(label, handle); 

        return *handle;
    }
#endif

    float JetResolution::getResolution(const JetParameters& parameters) const {
        const JetResolutionObject::Record* record = m_object->getRecord(parameters);
        if (! record)
            return 1;

        return m_object->evaluateFormula(*record, parameters);
    }

    JetResolutionScaleFactor::JetResolutionScaleFactor(const std::string& filename) {
        m_object = std::make_shared<JetResolutionObject>(filename);
    }

    JetResolutionScaleFactor::JetResolutionScaleFactor(const JetResolutionObject& object) {
        m_object = std::make_shared<JetResolutionObject>(object);
    }

#ifndef STANDALONE
    const JetResolutionScaleFactor JetResolutionScaleFactor::get(const edm::EventSetup& setup, const std::string& label) {
        edm::ESHandle<JetResolutionObject> handle;
        setup.get<JetResolutionScaleFactorRcd>().get(label, handle); 

        return *handle;
    }
#endif

    float JetResolutionScaleFactor::getScaleFactor(const JetParameters& parameters, Variation variation/* = Variation::NOMINAL*/) const {
        const JetResolutionObject::Record* record = m_object->getRecord(parameters);
        if (! record)
            return 1;

        const std::vector<float>& parameters_values = record->getParametersValues();
        return parameters_values[static_cast<size_t>(variation)];
    }

}

#ifndef STANDALONE
#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(JME::JetResolution);
TYPELOOKUP_DATA_REG(JME::JetResolutionScaleFactor);
#endif

import FWCore.ParameterSet.Config as cms

#### MC PU DISTRIBUTIONS
from SimGeneral.MixingModule.mix_2015_25ns_Startup_PoissonOOTPU_cfi import mix as mix_2015_25ns
from SimGeneral.MixingModule.mix_2015_50ns_Startup_PoissonOOTPU_cfi import mix as mix_2015_50ns
from SimGeneral.MixingModule.mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi import mix as mix_2015_25ns_realistScenario
pu_distribs = { "74X_mcRun2_asymptotic_v2" : mix_2015_25ns.input.nbPileupEvents.probValue,
                "76X_mcRun2_asymptotic_v12" : mix_2015_25ns_realistScenario.input.nbPileupEvents.probValue
                }

#### DATA PU DISTRIBUTIONS
data_pu_distribs = { "Jamboree_golden_JSON" : [5.12e+04,3.66e+05,5.04e+05,4.99e+05,7.5e+05,1.1e+06,2.53e+06,9.84e+06,4.4e+07,1.14e+08,1.94e+08,2.63e+08,2.96e+08,2.74e+08,2.06e+08,1.26e+08,6.38e+07,2.73e+07,1.1e+07,5.2e+06,3.12e+06,1.87e+06,9.35e+05,3.64e+05,1.1e+05,2.64e+04,5.76e+03,1.53e+03,594,278,131,59.8,26,10.8,4.29,1.62,0.587,0.203,0.0669,0.0211,0.00633,0.00182,0.000498,0.00013,3.26e-05,7.77e-06,1.77e-06,3.85e-07,7.99e-08,1.58e-08,3e-09,5.43e-10] }
    
    
pileupProducer = cms.EDProducer("PileupWeightProducer",
                                #hardcodedWeights = cms.untracked.bool(True),
                                pileupInfoTag    = cms.InputTag("slimmedAddPileupInfo"),
                                PileupMC = cms.vdouble(pu_distribs["76X_mcRun2_asymptotic_v12"]),
                                PileupData = cms.vdouble(data_pu_distribs["Jamboree_golden_JSON"]),
                                )

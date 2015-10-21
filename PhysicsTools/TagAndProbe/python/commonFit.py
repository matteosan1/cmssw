import FWCore.ParameterSet.Config as cms

all_pdfs = cms.PSet(
tight_20p0To40p0_0p0To1p5 = cms.vstring(
"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])",
"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])",
"ZGeneratorLineShape::signalPhyPass(mass,\"MCTemplates/templatesID.root\", \"hMass_20.0To40.0_0.0To1.5_Pass\")",
"ZGeneratorLineShape::signalPhyFail(mass,\"MCTemplates/templatesID.root\", \"hMass_20.0To40.0_0.0To1.5_Fail\")",
"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
"FCONV::signalPass(mass, signalPhyPass, signalResPass)",
"FCONV::signalFail(mass, signalPhyFail, signalResFail)",
"efficiency[0.5,0,1]",
"signalFractionInPassing[1.0]"
),

tight_20p0To40p0_1p5To2p5 = cms.vstring(
"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])",
"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])",
"ZGeneratorLineShape::signalPhyPass(mass,\"MCTemplates/templatesID.root\", \"hMass_20.0To40.0_1.5To2.5_Pass\")",
"ZGeneratorLineShape::signalPhyFail(mass,\"MCTemplates/templatesID.root\", \"hMass_20.0To40.0_1.5To2.5_Fail\")",
"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
"FCONV::signalPass(mass, signalPhyPass, signalResPass)",
"FCONV::signalFail(mass, signalPhyFail, signalResFail)",
"efficiency[0.5,0,1]",
"signalFractionInPassing[1.0]"
),

tight_40p0To60p0_0p0To1p5 = cms.vstring(
"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])",
"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])",
"ZGeneratorLineShape::signalPhyPass(mass,\"MCTemplates/templatesID.root\", \"hMass_40.0To60.0_0.0To1.5_Pass\")",
"ZGeneratorLineShape::signalPhyFail(mass,\"MCTemplates/templatesID.root\", \"hMass_40.0To60.0_0.0To1.5_Fail\")",
"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
"FCONV::signalPass(mass, signalPhyPass, signalResPass)",
"FCONV::signalFail(mass, signalPhyFail, signalResFail)",
"efficiency[0.5,0,1]",
"signalFractionInPassing[1.0]"
),

tight_40p0To60p0_1p5To2p5 = cms.vstring(
"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])",
"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])",
"ZGeneratorLineShape::signalPhyPass(mass,\"MCTemplates/templatesID.root\", \"hMass_40.0To60.0_1.5To2.5_Pass\")",
"ZGeneratorLineShape::signalPhyFail(mass,\"MCTemplates/templatesID.root\", \"hMass_40.0To60.0_1.5To2.5_Fail\")",
"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
"FCONV::signalPass(mass, signalPhyPass, signalResPass)",
"FCONV::signalFail(mass, signalPhyFail, signalResFail)",
"efficiency[0.5,0,1]",
"signalFractionInPassing[1.0]"
),

tight_60p0To100p0_0p0To1p5 = cms.vstring(
"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])",
"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])",
"ZGeneratorLineShape::signalPhyPass(mass,\"MCTemplates/templatesID.root\", \"hMass_60.0To100.0_0.0To1.5_Pass\")",
"ZGeneratorLineShape::signalPhyFail(mass,\"MCTemplates/templatesID.root\", \"hMass_60.0To100.0_0.0To1.5_Fail\")",
"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
"FCONV::signalPass(mass, signalPhyPass, signalResPass)",
"FCONV::signalFail(mass, signalPhyFail, signalResFail)",
"efficiency[0.5,0,1]",
"signalFractionInPassing[1.0]"
),

tight_60p0To100p0_1p5To2p5 = cms.vstring(
"RooGaussian::signalResPass(mass, meanP[1.0,-5.000,5.000],sigmaP[0.5,0.001,5.000])",
"RooGaussian::signalResFail(mass, meanF[1.0,-5.000,5.000],sigmaF[0.5,0.001,5.000])",
"ZGeneratorLineShape::signalPhyPass(mass,\"MCTemplates/templatesID.root\", \"hMass_60.0To100.0_1.5To2.5_Pass\")",
"ZGeneratorLineShape::signalPhyFail(mass,\"MCTemplates/templatesID.root\", \"hMass_60.0To100.0_1.5To2.5_Fail\")",
"RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
"FCONV::signalPass(mass, signalPhyPass, signalResPass)",
"FCONV::signalFail(mass, signalPhyFail, signalResFail)",
"efficiency[0.5,0,1]",
"signalFractionInPassing[1.0]"
),

)
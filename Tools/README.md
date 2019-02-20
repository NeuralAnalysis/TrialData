README: TrialData/Tools

This folder contains a set of general MATLAB processing functions for use with TrialData structures. As a convention, every function should take as inputs the trial data structure first, and a parameters struct second. They should return as outputs the processed trial data struct first, and then any other useful outputs can follow. However, there are some exceptions to this rule when it is advantageous.

1) getTDidx will return the trial indices first and the trial data struct second, to allow inline indexing of structures/variables without the need for dummy variables

2) getTDfields simply returns a list of fields

LIST OF AVAILABLE FUNCTIONS

addCorrelatedNoise

addFiringRates

appendTDs

binTD

catTDs

convBasisFunc

dimReduce

dupeAndShift

getCommonUnits

getDifferential

getEnvelope

getMoveOnsetAndPeak

getNorm

getSig

getTDfields

getTDidx

removeBadNeurons

removeBadTrials

reorderTDfields

smoothSignals

stripSpikeSorting

subtractConditionMean

trialAverage

trimTD

zscoreSignals
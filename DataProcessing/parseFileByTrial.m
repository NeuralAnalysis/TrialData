function trial_data = parseFileByTrial(data,inputArgs)
%   data: a CDS file or BDF file
%
% inputArgs:
%   meta: .perturbation, .epoch, .angle_dir, .rotation_angle, .force_magnitude, .force_angle (REQUIRED)
%   trialResults: which reward codes to use ('R','A','F','I')
%   binSize: default 0.01 sec
%   extraTime: [time before, time after] beginning and end of trial (default [0.5 0.3] sec)

if isa(data,'commonDataStructure')
    trial_data = parseFileByTrial_cds(data,inputArgs);
else
    error('BDF not currently supported.');
    trial_data = parseFileByTrial_bdf(data,inputArgs);
end
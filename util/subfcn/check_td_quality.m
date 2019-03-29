function trial_data = check_td_quality(trial_data)

if isempty(trial_data)
    error('No trials in the trial_data struct!');
end

% make  sure it's a struct
if ~isstruct(trial_data)
    error('First input must be trial_data struct!');
end

% make sure they all have the same bin size and that it's useful
if isfield(trial_data,'bin_size')
    unique_bin_sizes = unique([trial_data.bin_size]);
    if length(unique_bin_sizes) ~= 1
        error('Each trial must have the same bin size!');
    end
    if isnan(unique_bin_sizes)
        error('Bin size is NaN! This is a problem');
    end
    if unique_bin_sizes <= 0
        error('Bin size must be positive!');
    end
else
    error('No bin sizes field!');
end


function out = processNSx(filename,signal_info)
% NOTE: By default this returns mV. Alternative is to return the blackrock
% digital values. Maybe someday make this more flexible

convert_to_mv = true;
if isfield(signal_info,'params')
    assignParams(who,signal_info.params)
end

NSx = openNSx_td(filename,'read');

% strip whitespace from labels
labels = {NSx.ElectrodesInfo.Label};
for i = 1:length(labels)
    temp = uint16(labels{i});
    labels{i} = labels{i}(temp ~= 0);
end

data = double(NSx.Data)';
if convert_to_mv
    data = data/4e3;
end

out.t = 0:1/NSx.MetaTags.SamplingFreq:NSx.MetaTags.DataDurationSec-1/NSx.MetaTags.SamplingFreq;
out.data = data;
out.labels = labels;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = processNSx(filename,params)
if nargin > 1
    warning('no parameters supported for processNSx now...')
end

NSx = openNSx_td(filename,'read');

% strip whitespace from labels
labels = {NSx.ElectrodesInfo.Label};
for i = 1:length(labels)
    temp = uint16(labels{i});
    labels{i} = labels{i}(temp ~= 0);
end

out.duration = NSx.MetaTags.DataDurationSec;
out.samprate = NSx.MetaTags.SamplingFreq;
out.data = NSx.Data';
out.labels = labels;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
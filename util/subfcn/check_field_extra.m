function field_extra = check_field_extra(field_extra,signals)

if nargin == 1
    signals = {''};
end

if isempty(field_extra) % default to input signal names
    field_extra = repmat({''},1,size(signals,1));
else
    if ~iscell(field_extra),  field_extra = {field_extra};  end
    if length(field_extra) > 1 && length(field_extra) ~= size(signals,1)
        error('If you specify field_extra, you must have one entry or one entry for every signal.');
    end
    if length(field_extra) == 1 && size(signals,1) > 1
        field_extra = repmat(field_extra,1,size(signals,1));
    end
end
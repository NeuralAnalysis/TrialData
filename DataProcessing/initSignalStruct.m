function signal_info = initSignalStruct(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function serves as a way to initialize the signal struct input to
% make the convertDataToTD processing code simpler. See the code below
% where the struct is instantiated to see a list of parameters that you can
% modify.
%
% INPUTS:
%   varargin : any number of 'name','value' pairs
%
%   Parameters that can be modified:
%     1) filename : Must be specified! Where to find the file.
%     2) routine  : Optional. This contains an inline function handle (@)
%                   to a file that contains the code necessary to process
%                   the data if the default functionality is not sufficient
%                   or desired.
%           Note that filename and routine must only have one entry per
%           struct. Together, these define the unique processing operations
%           taken by the convertDataToTD code.
%     3) name : Must be specified, but can be a cell array of arbitrary
%               length. Each provided name will put one signal in the
%               trial_data struct, in a field corresponding to the name.
%     4) label : which signals to take out of the file. This is most useful
%                when using the built-in processing (e.g. to select analog
%                channels recorded on cerebus) but could have creative uses
%                when interfacing with your custom code.
%     5) routine_params : a struct containing arbitrary parameters for the
%                         routine, if desired.
%     6) operation : an inline function handle (e.g. @diff). Note this is
%                    distinct from routine since it runs AFTER the data has
%                    been loaded and grouped. So for a given name, this
%                    will operate on the continuous data between the
%                    loading routine and any subsequent built-in
%                    processing. Useful for differential EMGs, for example.
%            NOTE: As of now, this function must only take a single matrix
%            as input, where TIME is actually the columns and signals are
%            the rows. This is the transposed version of the trial_data
%            standard, but was done to accomodate some basic matlab
%            functions with simplicity (e.g. diff)
%     7) type : (string) what type of signal we are processing. Options:
%             - 'generic' : any old time-varying signal
%             - 'spikes' : spiking data. Will add a unit_guide [elec, unit]
%                          Doesn't have to be from NEV, you can use this
%                          with your own processing code.
%             - 'lfp' : (NOT IMPLEMENTED) will eventually filter LFP bands
%                       using continuously sampled data
%             - 'emg' : will do a rectifying and filtering (see subfunction
%                       in code. Can pass in parameters in the main code)
%             - 'trigger' : a pulse signal. Code will turn it into an idx_.
%                           Can pass in trigger thresholding parameter
%             - 'event' : bin number that will become an idx_ field
%                           Note for now, it assumes the processing code
%                           will spit out a bin number in the correct time
%                           frame, and won't do any corrections.
%             - 'meta'  : any arbitrary meta parameter
%
%    NOTE: If name is a cell array of length > 1, it will process multiple
%    signals from this file. You have two options for the parameters...
%    1) pass in each parameter as a cell array of the same length of names,
%       where each entry corresponds to the matching name. You should be
%       able to pass in '' or [] if you just want to use defaults for any
%       entry in the cell array.
%    2) fill them in normally or do nothing. This code will replicate them
%       to be of the length of names.
%
% OUTPUTS:
%   signal_info : the struct
%
% Here are a few potentially useful examples of how to initialize:
%
% 1) a way to take split arrays going to one pedestal and separate them
% initSignalStruct( ... % M1 neural data
%     'filename','blackrock_data.nev', ...
%     'name',{'M1','PMd'}, ...
%     'type','spikes', ...
%     'label',{m1_idx,pmd_idx}) % label should be electrode indices here
% 
% 2) a way to do differential EMGs recorded with cerebus
% initSignalStruct( ...
%     'filename','emg_data.ns5', ...
%     'name',{'Biceps','Triceps'}, ...
%     'label',{ {'BicepsElec1','BiceptsElec2'}, {'TricesElec1','TricepsElec2'} }, ...
%     'operation',@diff)
% note that 'label' points to the name in the datafile of the desired
% signals, and 'name' is the field name that will be in the TD
% 
% 3) a way to load signals recorded with some other means
% initSignalStruct( ...
%     'filename','analog_signals.dat', ...
%     'routine',@yourPersonalCodeToProcessData, ...
%     'name',{'pressure','motor_commands'}, ...
%     'type',{'','trigger'}) % for type, '' defaults to 'generic'
% 
% 4) a way to load events that were tracked manually and saved
% initSignalStruct( ...
%     'filename','tracked_events.mat', ...
%     'routine',@addEventsToTD, ...
%     'name',{'FootStrike','FootOff'}, ...
%     'type','event'), ...
%
% Written by Matt Perich. Updated April 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(varargin) == 1 % assume it's a cell array of fields
    varargin = varargin{1};
end
if ~iscell(varargin) || mod(length(varargin),2) ~= 0
    error('Inputs must be in ...''parameter'',''value'',... pairs');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intialize
signal_info = struct( ...
    'filename','', ...
    'routine','', ...
    'routine_params',struct(), ...
    'label',{{}}, ...
    'name','', ...
    'type','generic', ...
    'operation',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add data to signal_info
for i = 1:2:length(varargin)
    if isfield(signal_info,varargin{i})
        signal_info.(varargin{i}) = varargin{i+1};
    else
        error(['Field ' varargin{i} ' not recognized.']);
    end
end

% filename, routine, and routine_params must always be single
if iscell(signal_info.filename)
    if length(signal_info.filename) > 1
        error('Can only use one file per signal_info struct');
    else
        signal_info.filename = signal_info.filename{1};
    end
end
if iscell(signal_info.routine)
    if length(signal_info.routine) > 1
        error('Can only use one routine per signal_info struct');
    else
        signal_info.routine = signal_info.routine{1};
    end
end
if iscell(signal_info.routine_params)
    if length(signal_info.routine_params) > 1
        error('Can only use one routine per signal_info struct');
    else
        signal_info.routine_params = signal_info.routine_params{1};
    end
end

% look at name. If it's a cell  of length > 1, multiple signals are
% requested. In this case, duplicate all of the existing ones
if iscell(signal_info.name)
    num_signals = length(signal_info.name);
else
    num_signals = 1;
end
fn = fieldnames(signal_info);

for i = 1:length(fn)
    % skip filename and routine
    if ~strcmpi(fn{i},'filename') && ~strcmpi(fn{i},'routine') && ~strcmpi(fn{i},'routine_params')
        
        if ~iscell(signal_info.(fn{i}))
            % if it's not a cell, make it a cell of length num_signals
            signal_info.(fn{i}) = repmat({signal_info.(fn{i})},1,num_signals);
            
        else % it's a cell
            if length(signal_info.(fn{i})) == 1
                % if it's a cell of length one, make it length num_signals
                signal_info.(fn{i}) = repmat({signal_info.(fn{i})},1,num_signals);
            elseif strcmpi(fn{i},'label')
                % label is an exception, since it can be a cell anyway
                temp = signal_info.(fn{i});
                if isempty(temp) || ~iscell(temp)
                    % if it's a cell array of cell arrays, assume it's already the
                    % appropriate length for num_signals
                    signal_info.(fn{i}) = repmat({signal_info.(fn{i})},1,num_signals);
                elseif iscell(temp{1}) && length(temp) ~= num_signals
                    error('Label does not match .name field');
                end
            elseif length(signal_info.(fn{i})) > num_signals
                error('Something does not match with the .name field and other parameters.');
            end
        end
        
    end
end

% search through it all and look for empty type or routine_params
temp = signal_info.type;
idx = find(cellfun(@isempty,temp));
for i = idx
    temp{i} = 'generic';
    signal_info.type = temp;
end








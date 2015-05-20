
% modify path
addpath('/Users/renaud/Documents/programming/MATLAB/MCTW/matlab/');

% displayed at start
fprintf('now calculating the efficiency curve...\n\n');

% what to process
WTP         = 'NMDA';
gains       = [0.1:0.1:0.8,0.85:0.05:1.6,2,3,5,9,12];

% list data files
dataFiles	= dir(['output/*' WTP '*']);
nFiles      = length(dataFiles);

% load parameters
parameters

% extra parameters
dt          = 0.05;             %               [ms]
binarization= 3;                %               [ms]
depth       = 40;               % depth         [time bins]
shift       = 0;                % time shift    [ms]

% load files and start the processing
for fls = 12:1:12%nFiles
    
    % load data
    fprintf(['\t Now dealing with file ' dataFiles(fls).name '...\t']);
    load(['output/' dataFiles(fls).name]);
    time    = timeClipped(3:end,:);
    S       = SClipped(3:end,:);
    clear SClipped timeClipped;
    
    % lookup table
    V       = S(:,10);
    
    % read gain factor from the filename
    gain(fls) = str2num(dataFiles(fls).name(14:18));
    if gain(fls) == 1
        gainIndex = fls;
    end
    
    % find APs
    [pks,locs]      = findpeaks(V,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',3/dt);
    if length(pks) == 0
        frequency(fls) = 0;
    else
        frequency(fls) = 1e03*length(pks)/DTA;
    end
    BinOutput = (locs-1)*dt;
    
    % input spike train
    load stimNirenON;
    stim                            = stim(deadTime+1:end);
    BinInput                        = [];
    j                               = 1;
    while j < length(stim)
        if stim(j) > 0
            BinInput = cat(1,BinInput,j);
            j = j+9;
        end
        j = j+1;
    end
    BinInput = (BinInput-1)*dt;
    
    % CTW
    try
        [mi(fls),h(fls),ch(fls),chshuffle(fls)] = mutual_information_ctw(BinInput,BinOutput,0,125000,3,depth,shift);
    catch
        mi(fls) = 0; h(fls) = 0; ch(fls) = 0; chshuffle(fls) = 0;
    end
    
    % clear variables at the end of the loop
    %clear time S pks locs V BinInput BinOutput;
    fprintf(' Done\n');
end
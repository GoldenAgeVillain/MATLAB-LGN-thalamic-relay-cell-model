
% what to process
WTP         = 'NMDA';
gains       = [0.1:0.1:0.8,0.85:0.05:1.6,2,3,5,9,12];

% list data files
dataFiles	= dir(['output/*' WTP '*']);
nFiles      = length(dataFiles);

% load parameters
parameters

% extra parameters
dt          = 0.05e-03;                     % sec
binarization= 3e-03;                        % sec
windowB     = 2040;                         % time steps
windowA     = 4020;                         % time steps
psth        = zeros(windowA+windowB,1);     % generic container
idx         = 1;                            % generic index

% load files and start the processing
for fls = 12:1:12%nFiles
    
    % load data
    fprintf(['\t Now dealing with file ' dataFiles(fls).name '...\t']);
    load(['output/' dataFiles(fls).name]);
    time    = timeClipped(3:end,:);
    S       = SClipped(3:end,:);
    clear SClipped timeClipped;
    
    % lookup table
    m       = S(:,1);
    h       = S(:,2);
    n       = S(:,3);
    m_iA    = S(:,4);
    h_iA    = S(:,5);
    m_iT    = S(:,6);
    h_iT    = S(:,7);
    Ca      = S(:,8);
    O       = S(:,9);
    V       = S(:,10);
    
    % read gain factor from the filename
    gain(fls) = str2num(dataFiles(fls).name(14:18));
    if gain(fls) == 1
        gainIndex = fls;
    end
    
    % find APs
    [pks,locs]      = findpeaks(S(:,end),'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',3e-03/dt);
    if length(pks) == 0
        frequency(fls) = 0;
    else
        frequency(fls) = 1e03*length(pks)/DTA;
    end
    
    % sort data
    load stimNirenON;
    stim                            = stim(deadTime+1:end);
    BinInput                        = stim > 0;
    clear j stim S;
    i = windowB+1;
    while i < length(BinInput)-windowA
        if BinInput(i) == 1
            rng              = V(i-windowB:i+windowA-1);
            container(idx,:) = rng;
            [pks,locs]       = findpeaks(rng,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',3e-03/dt);
            buffer           = zeros(windowB+windowA,1);
            buffer(locs)     = 1;
            
            % count only the first AP before and after the input
            k                = windowB+1;
            while buffer(k) == 0 && k < length(buffer)
                k = k+1;
            end
            buffer(k+1:end)  = 0;
            k                = windowB;
            while buffer(k) == 0 && k > 1
                k = k-1;
            end
            buffer(1:k-2)    = 0;
            
            psth             = psth+buffer;
            idx              = idx+1;
            i                = i+10; % to account for the fact that spikes are square pulses in the input spike train
        end
        i = i+1;
    end
    V_avg = mean(container);
    
    % output to the terminal
    fprintf(' Done\n');
    
    % clear variables at the end of the loop
    clear time S gAMPA gNMDA gTot pks locs m h n m_iA h_iA m_iT h_iT Ca...
        O P1 OL V iInj iNa ih;
    
end

% downsample the psth and the voltage
j = 1;
for i = 1:binarization/dt:length(psth)
    rng             = i:min((i+binarization/dt),length(psth));
    PSTH_bin(j)     = sum(psth(rng));
    V_bin(j)        = mean(V_avg(rng));
    j               = j+1;
end
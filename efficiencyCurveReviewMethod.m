
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
dt          = 0.05e-03;         % sec
binarization= 3e-03;            % sec
windowBase  = 5;                % number of binarized steps

% load conductances
load('ampa conductance');                               % nS
load('nmda conductance');                               % nS
ampa = ampa(deadTime+1:end);
nmda = nmda(deadTime+1:end);

% load files and start the processing
for fls = 1:1:nFiles
    
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
    idx = fls;
    
    % find APs
    [pks,locs]      = findpeaks(S(:,end),'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',3e-03/dt);
    if length(pks) == 0
        frequency(fls) = 0;
    else
        frequency(fls) = 1e03*length(pks)/DTA;
    end
    
    % adjust the size of the window
    window = min(windowBase,floor(1/(2*frequency(fls)*binarization)));
    
    % binarize and downsample the spike trains
    voltage                         = V;
    binarizedSpikeTrain             = zeros(length(voltage),1);
    binarizedSpikeTrain(locs)       = 1;
    load stimNirenON;
    stim                            = stim(deadTime+1:end);
    BinInput                        = zeros((length(stim) - mod(length(stim),binarization/dt))/(binarization/dt),1);
    BinOutput                       = BinInput;
    j                               = 1;
    for i = 1:binarization/dt:length(stim)
        rng             = i:min((i+binarization/dt),length(stim));
        BinInput(j)     = sum(stim(rng));
        BinOutput(j)    = sum(binarizedSpikeTrain(rng));
        j = j+1;
    end
    BinInput    = BinInput > 50;
    BinOutput   = BinOutput > 0;
    clear j stim binarizedSpikeTrain voltage S;
%     BinInput    = (rand(size(BinInput)) < sum(BinInput)/length(BinInput));
    
    % calculate probabilities
    p1 = 0; p2 = 0; p3 = 0; p4 = 0; p5 = 0;
    i = 1;
    while i < length(BinInput)
        if BinInput(i) == 1
            if sum(BinOutput(i:min(length(BinOutput),i+window))) > 0
                p1 = p1+1;
                i  = i+window;
            else
                p2 = p2+1;
            end
        else
            if sum(BinOutput(i:min(length(BinOutput),i+window))) > 0
                p4 = p4+1;
                i  = i+window;
            else
                p5 = p5+1;
            end
        end
        i = i+1;
    end
    pOutput = sum(BinOutput)/length(BinOutput);
    
    % calculate entropies
    p           = sum(BinInput)/length(BinInput);
    Htotal      = -p*log2(p)-(1-p)*log2(1-p);
    HT(idx)     = Htotal/binarization;
    terms       = [p1/(p1+p4) p4/(p1+p4) p2/(p2+p5) p5/(p2+p5)];
    STR(fls,:)  = [terms pOutput];
    terms       = terms.*log2(terms);
    terms(isnan(terms)) = 0;
    Hcond       = pOutput*sum(terms(1:2))+(1-pOutput)*sum(terms(3:4));
    HN(idx)     = Hcond/binarization;
    HH(idx)     = Htotal/binarization+Hcond/binarization;
    fprintf(' Done\n');

    % clear variables at the end of the loop
    clear time S gAMPA gNMDA gTot pks locs m h n m_iA h_iA m_iT h_iT Ca...
        O P1 OL V iInj iNa ih;
    
end
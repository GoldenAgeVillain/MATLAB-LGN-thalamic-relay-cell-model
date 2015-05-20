%
% meta script to run all the simulations for Julia's Neuron paper
%
% author:       RJ
% created:      2014.05.04
% last update:  2014.05.25
% update:       added interleaved loops to generate data in a more useful
%               way allowing faster assessment of the final results
%

% try-catch for gain = 1
gain = 1;
try
    [time,S] = BazhenovModel(gain,1,0,'CTRL',[]);
catch
end
try
    [time,S] = BazhenovModel(gain,1,0,'NMDA',[]);
catch
end

% loop for part #1
gain = [0.5 1.5 0.8 2 7];
for i = 1:1:4
    try
        [time,S] = BazhenovModel(gain(i),1,0,'CTRL',[]);
    catch
    end
    try
        [time,S] = BazhenovModel(gain(i),1,0,'NMDA',[]);
    catch
    end
end

% loop for part #2
gain = [0.7 0.3 1.4 0.9 9 1.1 0.6 12 0.1 0.2 0.4 5 1.2 1.3 3];
for i = 1:1:length(gain)
    try
        [time,S] = BazhenovModel(gain(i),1,0,'CTRL',[]);
    catch
    end
    try
        [time,S] = BazhenovModel(gain(i),1,0,'NMDA',[]);
    catch
    end
end
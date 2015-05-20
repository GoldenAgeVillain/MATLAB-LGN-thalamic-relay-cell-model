%  
% 	Simple test of Mickey London's implementation of the CTW compression
% 	algorithm. Compares the theoretical predictions based on both the input
% 	parameters and on the effectively measured probabilities and verifies
% 	for the effects of repetitions in the input sequence. Adds the
% 	possibility to test the effects of a shift between input and ouput
% 	sequences. In here, the lag is the same across the 5 repetitions.
% 
% 	Started:                    2015.02.10
%  	Written by: 				Renaud Jolivet
%  	Last modification:			2015.02.10
%  	

% Adjust path
tic, addpath('/Users/renaud/Documents/programming/MATLAB/MCTW/matlab/');

% Parameters
bins    = 3;        % ms
freq    = 20;       % Hz
p       = freq/1000*bins;
dur     = 50000;    % ms
depth   = 20;
shift   = 15;       % ms // see also lines 34-35 below

% Input train
X       = dur/5*sort(rand(floor(p*(dur/5)/bins),1));   
X       = [X; X+dur*1/5; X+dur*2/5; X+dur*3/5; X+dur*4/5];
idx     = 1;        % generic index
del     = [0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99];

% Output trains
for pdel = del
    deletions   = (rand(floor(p*(dur/5)/bins),1) < 1-pdel);
    actualshift = 2*randn(size(deletions))+shift;                                       % adds random shift to the  
    actualshift = [actualshift; actualshift; actualshift; actualshift; actualshift];    % output spike train 
    deletions   = [deletions; deletions; deletions; deletions; deletions];                                                       
    Y           = (X+actualshift).*deletions;
    Y(Y == 0)   = [];   
    [mi(idx),h(idx),ch(idx),chshuffle(idx)] = mutual_information_ctw(X,Y,0,dur,bins,depth,shift);
    py(idx)     = length(Y)/dur*bins;
    H(idx)      = (-py(idx)*log2(py(idx))-(1-py(idx))*log2(1-py(idx)))/(bins*1e-03);    
    p1m(idx)    = sum(deletions);
    p2m(idx)    = sum(1-deletions);
    p4m(idx)    = 0;
    p5m(idx)    = dur/bins-length(X);
    idx         = idx+1;
end

% Plot CTW results
plot(del,mi,'-r'); hold on;
plot(del,H,'-b');
xlabel('Fraction of deleted APs');
ylabel('Information (bits/sec)');

% Theory
p1      = p*(1-del);
p2      = p*del;
p4      = 0;
p5      = (1-p);
terms   = [p1./(p1+p4); p4./(p1+p4); p2./(p2+p5); p5./(p2+p5)];
terms   = terms.*log2(terms);
terms(isnan(terms)) = 0;
HInput  = -p*log2(p)-(1-p)*log2(1-p);
HNoise  = py.*sum(terms(1:2,:))+(1-py).*sum(terms(3:4,:));

% Plot theoretical predictions
plot(del, (HInput+HNoise)/(bins*1e-03), '-g');

% Theory (measured)
terms   = [p1m./(p1m+p4m); p4m./(p1m+p4m); p2m./(p2m+p5m); p5m./(p2m+p5m)];
terms   = terms.*log2(terms);
terms(isnan(terms)) = 0;
HNoise  = py.*sum(terms(1:2,:))+(1-py).*sum(terms(3:4,:));

% Plot theoretical predictions (measured)
plot(del, (HInput+HNoise)/(bins*1e-03), '--g');
legend({'MI (bits/sec)','output channel carrying capacity (bits/sec)','Theoretical prediction (bits/sec)'});
ylim([0 110]); toc
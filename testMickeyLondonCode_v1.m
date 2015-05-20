%  
% 	Simple test of Mickey London's implementation of the CTW compression
% 	algorithm.
% 
%
%
%
%
% 	Started:                    2015.02.03
%  	Written by: 				Renaud Jolivet
%  	Last modification:			2015.02.10
%  	

% Adjust path
addpath('/Users/renaud/Documents/programming/MATLAB/MCTW/matlab/');

% Parameters
bins    = 3;        % ms
freq    = 20;       % Hz
p       = freq/1000*bins;
dur     = 50000;    % ms
depth   = 20;

% Input train
X       = dur*sort(rand(floor(p*dur/bins),1));   
idx     = 1;        % generic index
del     = [0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99];

% Output trains
for pdel = del
    Y           = X.*(rand(floor(p*dur/bins),1) < 1-pdel);
    Y(Y == 0)   = [];
    [mi(idx),h(idx),ch(idx),chshuffle(idx)] = mutual_information_ctw(X,Y,0,dur,bins,depth,0);
    py(idx)     = length(Y)/dur*bins;
    H(idx)      = (-py(idx)*log2(py(idx))-(1-py(idx))*log2(1-py(idx)))/(bins*1e-03);    
    idx         = idx+1;
end

% Plot results
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
legend({'MI (bits/sec)','output channel carrying capacity (bits/sec)','Theoretical prediction (bits/sec)'});
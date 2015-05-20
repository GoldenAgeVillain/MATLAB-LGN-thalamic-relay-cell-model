%  
% 	Simple test of Jiantao Jiao's implementation of the CTW compression
% 	algorithm. Adds the possibility to randomly shift the output spike
% 	train
% 
% 	Started:                    2015.02.05
%  	Written by: 				Renaud Jolivet
%  	Last modification:			2015.02.10
%  	

% Adjust path
tic, addpath('/Users/renaud/Documents/programming/MATLAB/DI_code_release/');

% Parameters
bins    = 3;        % ms
freq    = 20;       % Hz
p       = freq/1000*bins;
dur     = 45000;    % ms
depth   = 6;

% Input train
X       = rand(floor(dur/bins),1) < p;   
idx     = 1;        % generic index
del     = [0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99];

% Output trains
for pdel = del
    
    % Generate output spike train
    Y           = X.*(rand(floor(dur/bins),1) < 1-pdel);  
    
    % Adds a random shift to the output spike train
    i = 1;
    while i < length(Y)
        if Y(i) == 1
            actualshift     = floor((2*randn(1)+15)/bins); 
            Y(i)            = 0;
            Y(i+actualshift)= 1;
            i               = i+actualshift;
        end
        i = i+1;
    end
    try
        Y(length(X)+1:end) = [];
    catch
    end
    
    % Returns to normal operations
    [MI,DI,rev_DI] = compute_DI_MI(X',Y',2,depth,'E1',0,[],0);
    mi(idx)     = MI(end)/length(MI);
    di(idx)     = DI(end)/length(DI);
    rev_di(idx) = rev_DI(end)/length(rev_DI);
    py(idx)     = sum(Y)/dur*bins;
    H(idx)      = (-py(idx)*log2(py(idx))-(1-py(idx))*log2(1-py(idx)))/(bins*1e-03);    
    idx         = idx+1;
end

% Plot results
plot(del,mi/(bins*1e-03),'-r'); hold on;
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
legend({'MI (bits/sec)','Output channel carrying capacity (bits/sec)','Theoretical prediction (bits/sec)'}); 
ylim([0 120]); toc
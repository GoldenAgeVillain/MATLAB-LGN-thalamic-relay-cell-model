% NMDA data (from an Origin OPJ file)
data = [...
    -74	0.195495495495496	0.200877692939899   1
    -64	0.447715792708333   0.46004188129276    1
    -54	0.973206594691358   1                   1 
    -55 1                   NaN                 1e09    % added fake point... 
    -44	1.45164666233766    1.49161202796621    1       % here as well
    -34	2.98683752745098    3.06906831883751    1
    -24	4.44821014722222    4.57067407011656    1
    -14	7.08340533095238    7.278418960158      1
    -4	8.26679272666667    8.49438626059495    1
    46	9.36489211956522    9.62271748943007    1];                                          
                                                % added weights to force...
                                                % the fit through [-55,1]

% Fit
coeffNames      = {'a','b','c'};
myfun   = fittype(...
    'a/(1+b*exp(-c*V))',...
    'independent','V',...
    'coefficients',coeffNames);
options = fitoptions(...
    'method','NonLinearLeastSquares',...
    'StartPoint',[9 0.2 0.01],...
    'MaxFunEvals',5000,...
    'TolFun',1e-07,...
    'TolX',1e-07,...
    'Lower',[-Inf -Inf],...
    'Upper',[+Inf +Inf],...
    'Weights',data(:,4));
[cfun,gof] = fit(data(:,1),data(:,2),myfun,options);
p.data  = plot(data(:,1),data(:,2),'or'); hold('all');
rng     = -100:1:100;
p.fit   = plot(rng,cfun(rng),'-k','LineWidth',1);
set(p.data,'MarkerFaceColor','r');
xlabel('voltage (mV)');
ylabel('scaled NMDA conductance (a.u.)');
fprintf(['\nscaling at -55 mV = ' num2str(cfun(-55),'%2.2f') '\n\n']);
fprintf([...
    'parameters are:\n'...
    '\ta = ' num2str(cfun.a,'%2.2f') '\n'...
    '\tb = ' num2str(cfun.b,'%2.3f') '\n'...
    '\tc = ' num2str(cfun.c,'%2.3f') '\n\n']);
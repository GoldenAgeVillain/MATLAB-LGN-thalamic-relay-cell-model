%
% Extract steady-state IVs from Julia's recordings
%

filelist = dir('IVs/*.fig');

X = [];
Y = [];
Xs= [];
Ys= [];

Ri= [99.6 117.0 54.8 561.7 165.4 103.2 152.0 83.2 183.8 105.0 50.2 108.8 48.6 145.4 122.4 108.0 296.7 176.9];
Cm= [176.1 192.8 190.0 115.6 203.2 129.9 121.0 143.1 143.1 182.0 170.2 103.3 190.2 166.2 84.5 138.0 152.6 137.4];
Cs= Cm*1e-06;

for fls = 1:1:length(filelist)
    
    uiopen(['IVs/' filelist(fls).name],1)
    
    figHandles  = get(0,'Children');
    axesObjs    = get(figHandles(2), 'Children');
    dataObjs    = get(axesObjs, 'Children');
    objTypes    = get(dataObjs, 'Type');
    xdata       = get(dataObjs, 'XData');
    ydata       = get(dataObjs, 'YData');
    
    X           = cat(1,X,xdata{end}');
    Y           = cat(1,Y,ydata{end}');
    Xs          = cat(1,Xs,xdata{end}');
    Ys          = cat(1,Ys,ydata{end}'*1e-06/Cs(fls));
    
    fclose all, close all;
    
end

plot(Xs,Ys,'or'); hold on;
plot([-120 120],[0 0],'-k');

output = [Xs,Ys];
save('averageIVs','output');

% parameters
dt          = 0.05e-03;             % sec
dataFld     = 'data AP5';
cells       = dir(dataFld);
cells(1:3)  = [];
Ncells      = length(cells);
holdpot     = {'30','55','75'};
cond        = {'control','AP5','wash'};
for i = 1:1:3       % holding potential
    for j = 1:1:3   % condition
        dataSorted{i,j} = [];
        dataRaw{i,j}    = [];
    end
end

% main loop
for fld = 1:1:Ncells
    fprintf(['Now processing folder ' num2str(fld) ' out of ' num2str(Ncells) '\n']);
    buffer  = dir([dataFld '/' cells(fld).name '/*detrended*.mat']);
    buffer  = length(buffer);    
    fl      = 1;
    for i = 1:1:3       % holding potential
        for j = 1:1:3   % condition
            files = dir([dataFld '/' cells(fld).name '/*detrended*' holdpot{i} '*' cond{j} '*.mat']);
            for k = 1:1:length(files)                

                % load data
                fprintf(['\t file ' num2str(fl) ' out of ' num2str(buffer) '\n']);                
                load([dataFld '/' cells(fld).name '/' files(k).name]);
                if exist('vClamp')
                    data = -vClamp;
                end
                DAT                     = data/str2num(holdpot{i});         % divide current by driving force
                dataRaw{i,j}            = cat(2,dataRaw{i,j},DAT);          % store raw data
                DAT(DAT < 0)            = 0;                                % remove negative values
                DAT(DAT > 7*std(DAT))   = 0;                                % remove outliers
                F                       = smooth(DAT,51,'sgolay',2);        % filter traces            
                dataSorted{i,j}         = cat(2,dataSorted{i,j},F);         % store filtered traces                               
                fl                      = fl+1;                             % update file index
            end
        end
    end    
end
clear Ncells cells dataFld files fld i j k vClamp data;

% considered traces
total_131119c2  = mean(dataRaw{2,1}(:,1),2)-0.458;
ampa_131119c2   = mean(dataRaw{2,2}(:,1:2),2)-0.435;
nmda_131119c2   = total_131119c2-ampa_131119c2;

total_131120c1  = mean(dataRaw{2,1}(:,2:3),2)-0.415;
ampa_131120c1   = mean(dataRaw{2,2}(:,3:4),2)-0.147;
nmda_131120c1   = total_131120c1-ampa_131120c1;

total_140123c1  = mean(dataRaw{2,1}(:,4:5),2)-0.0935;
ampa_140123c1   = mean(dataRaw{2,2}(:,5),2)-0.0984;
nmda_140123c1   = total_140123c1-ampa_140123c1;

total   = mean([total_131119c2 total_131120c1 total_140123c1],2);
ampa    = mean([ampa_131119c2 ampa_131120c1 ampa_140123c1],2);
nmda    = mean([nmda_131119c2 nmda_131120c1 nmda_140123c1],2);

total(total > 7*std(total)) = 0;
ampa(ampa > 7*std(ampa))    = 0;
nmda(nmda > 7*std(nmda))    = 0;

total               = smooth(total,101,'sgolay',2); 
ampa                = smooth(ampa,101,'sgolay',2); 
nmda                = smooth(nmda,101,'sgolay',2); 

total(total < 0)    = 0;
ampa(ampa < 0)      = 0;
nmda(nmda < 0)      = 0;

timeAxis = dt*(1:1:length(dataRaw{1,1}))';
plot(timeAxis,total,'-k'), hold on;
plot(timeAxis,ampa,'-r');
plot(timeAxis,nmda,'-b');

save('ampa conductance','ampa');
save('nmda conductance','nmda');


% timeAxis = dt*(1:1:length(dataSorted{1,1}))';
% for i = 1:1:3
%     for j = 1:1:3
%         subplot(3,3,(i-1)*3+j);
%         plot(timeAxis,mean(dataRaw{i,j},2)), hold on;
%         plot(timeAxis,mean(dataSorted{i,j},2),'-r');
%         if i == 1
%             title(cond{j});
%         end
%         if i == 3
%             xlabel('time [sec]');
%         end
%         if j == 1
%             ylabel(['held @ -' holdpot{i} ' [mV]']);
%         end
%     end
% end
% 
% total   = mean(dataSorted{2,1},2);
% ampa    = mean(dataSorted{2,2},2);
% nmda    = total-ampa;
% figure;
% plot(timeAxis,total,'-k'), hold on;
% plot(timeAxis,ampa,'-r');
% plot(timeAxis,nmda,'-b');

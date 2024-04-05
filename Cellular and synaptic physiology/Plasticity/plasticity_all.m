clear
close all

alldirs=uipickfiles;
cd(char(alldirs));

% get directory name and path separately
[pathToDirName, dirName] = fileparts(char(alldirs));

%Load xml file to get info about the linescan protocol
fileName = [strcat(dirName) '.xml']; %transform into structure

[s] = xml2struct(fileName);

%Get number of cycles
ncycle= size(s.PVScan.Sequence,2);

% load .csv data 
csvFiles = dir('**/*.csv').';

% pre-allocate cell arrays
csvNames = cell(0);
Vtimes = cell(0);
voltageTraces = cell(0);
currentTraces = cell(0);

%get all voltage traces
% iterate over files, collecting voltage and current data and csvNames
for csvFile = csvFiles
    try
        data = csvread(csvFile.name, 1, 0);
    catch
        %disp('Error reading .csv file. File may be empty. Skipped.');
        continue
    end

    if size(data, 2) == 3
        csvNames{end + 1} = csvFile.name(1:end-25);
        Vtimes{end + 1} = data(:, 1);
        voltageTraces{end + 1} = movmedian(data(:,2)*100,10); %in mV 
        currentTraces{end + 1} = data(:,3)*10E3; %in pA
    else
        %disp('.csv file not in correct format. Skipped.');
    end
    
    %clear('data');    
end

f7 = figure; hold on;
for i = 1:length(voltageTraces)    
    plot(Vtimes{1,i},voltageTraces{1,i})
    %pause, clf
end

f6 = figure; hold on;
for i = 1:length(currentTraces)
    plot(Vtimes{1,i},currentTraces{1,i})
    curr_inj(i) = max(currentTraces{i});
end
for i = 1:length(currentTraces)
    scatter(i, mean(currentTraces{1,i})), hold on
end

saveas(f7, 'alltraces.png')
saveas(f6, 'currinj.png')

deltaI = 0; %size of step, in nA
Vmin = 0;  %mV minimum, steady state, and baseline for hyperpol step
Vbase = 0; 
Vss = 0;
Ri = 0; %input resistance, in MO
Vpeak = [0;0]; %peak EPSP
%Vslope = 0; %first 2ms of EPSP slope (Feldman, '00)
Vm_avg = 0; %avg resting membrane potential before stim
%artifact_size = 0; %stim artifact size
time_stamp = {0}; %realtime
xx = {};
ppr = 0;

t = find(diff(currentTraces{1,1})> 25 | diff(currentTraces{1,1})< -25);%find current pulse
dvdt = diff(voltageTraces{1,1}(5:end-5)./Vtimes{1}(2)); 
ii = find(dvdt <= -4); 
stim = [ii(find(diff(ii)>1)) ii(find(diff(ii)>1)+1)];%find stim artifact
%stim = 50/Vtimes{1}(2); %%stim location is hardcoded at 50ms
%%
for i = 1:length(voltageTraces)
    
    %get input resistance
    deltaI(i) = currentTraces{1,i}(t(1)) - currentTraces{1,i}(t(length(t))); 
    Vbase(i) = mean(voltageTraces{1,i}(t(1)-500:t(1)));
    Vmin(i) = min(voltageTraces{1,i}(t(1):t(length(t))));
    Vss(i) = mean(voltageTraces{1,i}(t(length(t))-1000:t(length(t))));
    Ri(i) = ((Vbase(i)-Vss(i))/50)*1000;
    
    %EPSP stats
    %clip to just EPSP, with baseline and following stim (arbitrary)
    x = voltageTraces{1,i}(stim(1)-49/.05:stim(2)+100/.05); 
    xx{i} = x; %store trace 
    
    Vm_avg(i) = mean(voltageTraces{1,i}(stim(1)-49/0.05:stim(1))); %baseline before epsp
    peak1(i) =  max(voltageTraces{1,i}(stim(1):stim(1)+200))-Vbase(i); 
    peak2(i) =  max(voltageTraces{1,i}(stim(2):stim(2)+200))-Vbase(i); 
    %peak = mean([peak1(i),peak2(i)]); %take avg of both peaks?
    Vpeak(:,i) = [peak1(i) peak2(i)];
    
    ppr(i) = peak2(i)/peak1(i);
   % Vslope(i) = mean(dvdt(stim+3/.05:stim+3/.05+40)); %slope after stim, 2ms

   %artifact_size(i) = min(voltageTraces{1,i}(stim)); 
   
   d = s.PVScan.Sequence{1,i}.Attributes.time; 
   [y, m, d, h, mn, ss] = datevec(d,'HH:MM:SS.FFF');
   time_stamp{i} = [h mn ss];
end

%%
%get timestamps
time_stamp = {0}; %realtime
for i = 1:length(voltageTraces)
    d = s.PVScan.Sequence{1,i}.Attributes.time; 
   [y, m, d, h, mn, ss] = datevec(d,'HH:MM:SS.FFF');
   time_stamp{i} = [h mn ss];
end

%%
spikes = 0; %eliminate spiking trials
f8 = figure; hold on, 
for i = 1:length(xx)
    if max(xx{i})< 0
        plot(xx{i})
    else
        spikes = [spikes i];
    end
end

spikes = spikes(2:end);

all_traces = cell2mat(xx);

%spikes= [13];
all_traces(:,spikes) = [];
Vpeak(:,spikes) = [];
%Vslope(spikes) = [];
Ri(spikes) = [];
Vm_avg(spikes) = [];
time_stamp(spikes) = [];
mean_trace = mean(all_traces,2); %mean epsp trace

f4 = figure, subplot(5,1,1), plot(Vpeak(1,:)), title('EPSP peak')
hold on, plot(Vpeak(2,:))
%subplot(4,1,2), plot(Vslope),title('EPSP slope')
subplot(5,1,2), plot(Ri),title('Ri')
subplot(5,1,3),plot(Vm_avg),title('Vm resting')
subplot(5,1,4),plot(ppr),title('PPR')
subplot(5,1,5),plot(curr_inj), title('Current inj')
savefig(f4,'subplots')

%%
%close all
save('ts2')
%saveas(f8,'ts5_allepsps.png')
%%
clear
load 20_hz_drugs_pop.mat
%%
clear V_mat
%% 5reps group, 0ms group *07251 is too short, removed 9263 (weak/runup)
Vpeak = {Vpeak_07112 Vpeak_07251 Vpeak_07281 Vpeak_07287 Vpeak_08017 Vpeak_09131 Vpeak_09263 Vpeak_09292 Vpeak_01061};
Vm = {Vm_07112 Vm_07251 Vm_07281 Vm_07287 Vm_08017 Vm_09131 Vm_09263 Vm_09292 Vm_01061};
Ri = {Ri_07112 Ri_07251 Ri_07281 Ri_07287 Ri_08017 Ri_09131 Ri_09263 Ri_09292 Ri_01061};
time = {time_07112 time_07251 time_07281 time_07287 time_08017 time_09131 time_09263 time_09292 time_01061};
%% all -1250
Vpeak = {Vpeak_07121 Vpeak_08083 Vpeak_08183 Vpeak_08194 Vpeak_08235 Vpeak_09222 Vpeak_09296}; %Vpeak_09222?
Vm = {Vm_07121 Vm_08083 Vm_08183 Vm_08194 Vm_08235 Vm_09222 Vm_09296};
Ri = {Ri_07121 Ri_08083 Ri_08183 Ri_08194 Ri_08235 Ri_09222 Ri_09296};
time = {time_07121 time_08083 time_08183 time_08194 time_08235 time_09222 time_09296};

%% obliques, 0ms, 5 reps, gabazine

Vpeak = {Vpeak_11141, Vpeak_11143, Vpeak_11162, Vpeak_11163, Vpeak_11171, Vpeak_11172, Vpeak_11182, Vpeak_11184};
Vm = {Vm_11141, Vm_11143, Vm_11162, Vm_11163, Vm_11171, Vm_11172, Vm_11182, Vm_11184};
Ri = {Ri_11141, Ri_11143, Ri_11162, Ri_11163, Ri_11171, Ri_11172, Ri_11182, Ri_11184};
time = {time_11141, time_11143, time_11162, time_11163, time_11171, time_11172, time_11182, time_11184};
%% all -2250s 
Vpeak = {Vpeak_11142 Vpeak_12205 Vpeak_01068 Vpeak_01105 Vpeak_06301 Vpeak_07202};
Vm = {Vm_11142 Vm_12205 Vm_01068 Vm_01105 Vm_06301 Vm_07202};
Ri = {Ri_11142 Ri_12205 Ri_01068 Ri_01105 Ri_06301 Ri_07202};
time = {time_11142 time_12205 time_01068 time_01105 time_06301 time_07202};
%% all 0.55s
Vpeak = {Vpeak_01248 Vpeak_01278 Vpeak_02027 Vpeak_02065 Vpeak_02244 Vpeak_07244}; %Vpeak_02215
Vm = {Vm_01248 Vm_01278 Vm_02027 Vm_02065 Vm_02244 Vm_07244};
Ri = {Ri_01248 Ri_01278 Ri_02027 Ri_02065 Ri_02244 Ri_07244};
time = {time_01248; time_01278; time_02027; time_02065; time_02244; time_07244};
%% all -0.55s
Vpeak = {Vpeak_04065 Vpeak_05265 Vpeak_05313 Vpeak_06153 Vpeak_06296 Vpeak_07102 Vpeak_06261 Vpeak_07123 Vpeak_07133}; 
Vm = {Vm_04065 Vm_05265 Vm_05313 Vm_06153 Vm_06296 Vm_07102 Vm_06261 Vm_07123 Vm_07133}; 
Ri = {Ri_04065 Ri_05265 Ri_05313 Ri_06153 Ri_06296 Ri_07102 Ri_06261 Ri_07123 Ri_07133};
time = {time_04065 time_05265 time_05313 time_06153 time_06296 time_07102  time_06261 time_07123 time_07133};
%% post only
Vpeak = {Vpeak_02037 Vpeak_02064 Vpeak_03131 Vpeak_03133 Vpeak_03134 Vpeak_04111 Vpeak_04124 Vpeak_04174}; %Vpeak_02036 too much rundown
Vm = {Vm_02037 Vm_02064 Vm_03131 Vm_03133 Vm_03134 Vm_04111 Vm_04124 Vm_04174};
Ri = {Ri_02037 Ri_02064 Ri_03131 Ri_03133 Ri_03134 Ri_04111 Ri_04124 Ri_04174};
time = {time_02037 time_02064 time_03131 time_03133 time_03134 time_04111 time_04124 time_04174};

%% pre only
Vpeak = {Vpeak_11216, Vpeak_11222, Vpeak_11225, Vpeak_12015, Vpeak_01123 Vpeak_04187};
Vm = {Vm_11216, Vm_11222, Vm_11225, Vm_12015, Vm_01123 Vm_04187};
Ri = {Ri_11216, Ri_11222, Ri_11225, Ri_12015, Ri_01123 Ri_04187};
time = {time_11216, time_11222, time_11225, time_12015, time_01123 time_04187};

%% NMDA block
Vpeak = {Vpeak_12122, Vpeak_12124, Vpeak_12131, Vpeak_12133 Vpeak_07251, Vpeak_07252};
Vm = {Vm_12122, Vm_12124, Vm_12131, Vm_12133 Vm_07251, Vm_07252};
Ri = {Ri_12122, Ri_12124, Ri_12131, Ri_12133 Ri_07251 Ri_07252};
time = {time_12122, time_12124, time_12131, time_12133 time_ind_07251, time_ind_07252};
%% 60hz APs 
Vpeak = {Vpeak_12144, Vpeak_12151, Vpeak_12152, Vpeak_01133, Vpeak_01137 Vpeak_04194};
Vm = {Vm_12144, Vm_12151, Vm_12152, Vm_01133, Vm_01137 Vm_04194};
Ri = {Ri_12144, Ri_12151, Ri_12152, Ri_01133, Ri_01137 Ri_04194};
time = {time_12144, time_12151, time_12152, time_01133, time_01137 time_04194};
%% ppr
% ppr = {ppr_04142 ppr_05123 ppr_05174 ppr_05175 ppr_05185 ppr_06074 ppr_06293 ppr_07154 ppr_07197 ppr_07204 ppr_08012};
%%

avg_pair = cellfun(@nanmean, Vpeak, 'UniformOutput', false);

%%
figure, hold on
for i = 1:length(avg_pair)
    plot(time{i}(1:3:end)',avg_pair{i}(1:3:end),'-o')
end
% %%
% figure, hold on
% for i = 1:length(Vpeak)
%     plot(time{i}(1:3:end),Vpeak{i}(:,1:3:end),'-o')
% end
% 
% figure, hold on
% for i = 1:length(Vm)
%     plot(time{i}(1:3:end),Vm{i}(:,1:3:end),'-o')
% end
% 
% figure, hold on
% for i = 1:length(Ri)
%     plot(time{i}(1:3:end),Ri{i}(:,1:3:end),'-o')
% end
%% extend NaNs to see all the data in matrix for alignment
for i = 1:length(Vpeak) 
   Vpeak{i} =  [Vpeak{i} NaN(2,(199-length(Vpeak{i})))]; %was 199
   Ri{i} = [Ri{i} NaN(1,(199-length(Ri{i})))];
   Vm{i} = [Vm{i} NaN(1,(199-length(Vm{i})))];
   time{i} = [time{i}; NaN((199-length(time{i})),1)];
end
% %%
% for i = 1:length(ppr) 
%    ppr{i} =  [ppr{i} NaN(1,(199-length(ppr{i})))];
% end
%%
Vm = cell2mat(Vm);
Vm = reshape(Vm,199,i);

Ri = cell2mat(Ri);
Ri = reshape(Ri,199,i);

time = cell2mat(time);
time = reshape(time,199,i);
%%
clear V_mat
ii = 1:2:length(Vpeak)*2;
for i = 1:length(Vpeak)
    V_mat{ii(i)} = Vpeak{i}(1,:);
    V_mat{ii(i)+1} = Vpeak{i}(2,:);
end

V_mat = cell2mat(V_mat);
V_mat = reshape(V_mat,199,i*2);
% %%
% ppr = cell2mat(ppr);
% ppr = reshape(ppr,199,i);

%%
V_mat = V_mat(1:2:end,:);
Vm = Vm(1:2:end,:);
Ri = Ri(1:2:end,:);
time = time(1:2:end,:);

% %%
% ppr = ppr(1:3:end,:);

%% manually adjust... align to induction period%%
V_mat(V_mat==0) = NaN
%%
Vm_mean = nanmean(Vm,2);
Vm_sem = nanstd(Vm, [], 2)./ sqrt(size(Vm,2));   

Ri_mean = nanmean(Ri,2);
Ri_sem = nanstd(Ri, [], 2)./ sqrt(size(Ri,2));   

Vpeak_mean = nanmean(V_mat,2);
Vpeak_sem = nanstd(V_mat, [], 2)./ sqrt(size(V_mat,2));   

%avg_time = nanmean(time,2);
avg_time = [0:mode(mode(diff(time))):max(max(time))];
%avg_time = avg_time(:,1:length(Vpeak_mean));
% %%
% ppr_mean = nanmean(ppr,2);
% ppr_sem = nanstd(ppr, [], 2)./ sqrt(size(ppr,2));   
% %%
% ppr(ppr==0) = NaN;
%%
figure, errorbar(avg_time(1:40),Vpeak_mean(1:40),Vpeak_sem(1:40), '-o','MarkerSize',10, 'CapSize',10, 'MarkerFaceColor',[0.00,0.45,0.74],'LineWidth',1)
hold on, line([1 max(avg_time(1:40))], [1 1],'LineStyle','--','Color',[0.50,0.50,0.50])
set(gca,'TickDir','out');
set(gca,'box','off');
axis([0 45 0.6 2]);
%%
%% pre-post boxplots for Ri, Vm
Ri(Ri==0) = NaN;
Vm(Vm==0) = NaN;

%%
Ri_pre = nanmean(Ri(1:14,:));
Ri_post = nanmean(Ri(14:end,:));
Ri_pre_mean = mean(Ri_pre);
Ri_post_mean = nanmean(Ri_post);

figure, plot(1,Ri_pre,'o','Color','k', 'MarkerSize',10), hold on
plot(1.2,Ri_post,'o','Color','b','MarkerSize',10)
line([ones(length(Ri_pre),1),1.2*ones(length(Ri_pre),1)]', [Ri_pre' Ri_post']','Color',[0.5 0.5 0.5])
plot(1,Ri_pre_mean,'.','MarkerSize',50,'Color','k')
plot(1.2,Ri_post_mean,'.','MarkerSize',50,'Color','b')
set(gca,'TickDir','out');
set(gca,'box','off');
axis([0.75 1.5 0 150]);

Vm_pre = nanmean(Vm(1:14,:));
Vm_post = nanmean(Vm(14:end,:));
Vm_pre_mean = mean(Vm_pre);
Vm_post_mean = mean(Vm_post);

figure, plot(1,Vm_pre,'o','Color','k', 'MarkerSize',10), hold on
plot(1.2,Vm_post,'o','Color','b', 'MarkerSize',10)
line([ones(length(Vm_pre),1),1.2*ones(length(Vm_pre),1)]', [Vm_pre' Vm_post']','Color',[0.5 0.5 0.5])
plot(1,Vm_pre_mean,'.','MarkerSize',50,'Color','k')
plot(1.2,Vm_post_mean,'.','MarkerSize',50,'Color','b')
set(gca,'TickDir','out');
set(gca,'box','off');
axis([0.75 1.5 -70 -55]);

figure, plot(1,ppr_all(:,1),'o','Color','k','MarkerSize',10), hold on
plot(1.2,ppr_all(:,2),'o','Color','b','MarkerSize',10)
line([ones(length(ppr_all),1),1.2*ones(length(ppr_all),1)]', [ppr_all(:,1) ppr_all(:,2)]','Color',[0.5 0.5 0.5])
plot(1,mean(ppr_all(:,1)),'.','MarkerSize',50,'Color','k')
plot(1.2,mean(ppr_all(:,2)),'.','MarkerSize',50,'Color','b')
set(gca,'TickDir','out');
set(gca,'box','off');
axis([0.75 1.5 0.5 2.5]);

%%
sem_2250 = nanstd(Vpeak_post_2250ms, [], 2)./ sqrt(size(Vpeak_post_2250ms,2));   
sem_1250 = nanstd(Vpeak_post_1250ms, [], 2)./ sqrt(size(Vpeak_post_1250ms,2));   
sem_n550 = nanstd(Vpeak_post_n550ms, [], 2)./ sqrt(size(Vpeak_post_n550ms,2));   
sem_550 = nanstd(Vpeak_post_550ms, [], 2)./ sqrt(size(Vpeak_post_550ms,2));   
sem_0 = nanstd(Vpeak_post_0ms, [], 2)./ sqrt(size(Vpeak_post_0ms,2)); 
sem_pre = nanstd(Vpeak_post_preonly, [], 2)./ sqrt(size(Vpeak_post_preonly,2)); 
sem_post = nanstd(Vpeak_post_postonly, [], 2)./ sqrt(size(Vpeak_post_postonly,2)); 

sem_60hz = nanstd(Vpeak_post_60hz, [], 2)./ sqrt(size(Vpeak_post_60hz,2));
sem_AP5 = nanstd(Vpeak_post_AP5, [], 2)./ sqrt(size(Vpeak_post_AP5,2)); 
sem_L4 = nanstd(Vpeak_post_L4, [], 2)./ sqrt(size(Vpeak_post_L4,2)); 
%%
figure, hold on
plot(0, Vpeak_post_0ms,'o','Color','r','MarkerSize',10)
plot(0.2,mean(Vpeak_post_0ms),'o','Color','r','MarkerSize',10),
errorbar(0.2,mean(Vpeak_post_0ms),sem_0,'LineWidth',1,'Color','r')

plot(-0.5, Vpeak_post_n550ms,'o','Color','r','MarkerSize',10)
plot(-0.3,mean(Vpeak_post_n550ms),'o','Color','r','MarkerSize',10),
errorbar(-0.3,mean(Vpeak_post_n550ms),sem_n550,'LineWidth',1,'Color','r')

plot(-1, Vpeak_post_1250ms,'o','Color',[0.8500 0.3250 0.0980],'MarkerSize',10)
plot(-0.8,mean(Vpeak_post_1250ms),'o','Color',[0.8500 0.3250 0.0980],'MarkerSize',10),
errorbar(-0.8,mean(Vpeak_post_1250ms),sem_1250,'LineWidth',1,'Color',[0.8500 0.3250 0.0980])

plot(-2, Vpeak_post_2250ms,'o','Color',[0.9290 0.6940 0.1250],'MarkerSize',10)
plot(-1.8,mean(Vpeak_post_2250ms),'o','Color',[0.9290 0.6940 0.1250],'MarkerSize',10),
errorbar(-1.8,mean(Vpeak_post_2250ms),sem_2250,'LineWidth',1,'Color',[0.9290 0.6940 0.1250])

plot(0.5, Vpeak_post_550ms,'o','Color',[0.6350 0.0780 0.1840],'MarkerSize',10)
plot(0.7,mean(Vpeak_post_550ms),'o','Color',[0.6350 0.0780 0.1840],'MarkerSize',10),
errorbar(0.7,mean(Vpeak_post_550ms),sem_550,'LineWidth',1,'Color',[0.6350 0.0780 0.1840])

plot(-2.5, Vpeak_post_postonly,'o','Color',[0.5 0.5 0.5],'MarkerSize',10)
plot(-2.3,mean(Vpeak_post_postonly),'o','Color',[0.5 0.5 0.5],'MarkerSize',10),
errorbar(-2.3,mean(Vpeak_post_postonly),sem_post,'LineWidth',1,'Color',[0.5 0.5 0.5])

plot(-3, Vpeak_post_preonly,'o','Color','k','MarkerSize',10)
plot(-2.8,mean(Vpeak_post_preonly),'o','Color','k','MarkerSize',10),
errorbar(-2.8,mean(Vpeak_post_preonly),sem_pre,'LineWidth',1,'Color','k')

hold on, line([-3.5 1], [1 1],'LineStyle','--','Color',[0.50,0.50,0.50])
set(gca,'box','off');
set(gca,'TickDir','out');
axis([-3.5 1 0.8 2.5 ]);
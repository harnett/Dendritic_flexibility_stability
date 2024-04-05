addpath 'D:\Courtney codes'
%% 
clear
%close all

alldirs=uipickfiles;
cd(char(alldirs));

% get directory name and path separately
[pathToDirName, dirName] = fileparts(char(alldirs));

%Load xml file to get info about the linescan protocol
fileName = [strcat(dirName) '.xml'];
%transform into structure
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
        voltageTraces{end + 1} = smooth(data(:,2)*100,10); %in mV
        currentTraces{end + 1} = data(:,3)*1E4; %in pA %just to check amt of current inj
    else
        %disp('.csv file not in correct format. Skipped.');
    end
    clear('data');    
end

figure, plot(voltageTraces{1})

% chop to individual sweeps
Vtimes = cell2mat(Vtimes);
voltageTraces = cell2mat(voltageTraces);
currentTraces = cell2mat(currentTraces);
%
clear sweeps
clear current_inj

start_idx = 1:70000:1000000; %hardcoded for 15 current injections, 3s interval, cut to 1.5
end_idx = 1.5E4:70000:1000000;

for i = 1:length(start_idx) 
    sweeps(:,i) = voltageTraces(start_idx(i):end_idx(i), 1); 
    current_inj(:,i) = currentTraces(start_idx(i):end_idx(i), 1); 
end

time = Vtimes(1:length(sweeps));
% figure, plot(time, current_inj)
% xlabel('ms')
% ylabel('pA')
figure, plot(time,sweeps)
xlabel('ms')
ylabel('mV')
%%
% Vm, Ri
%baseline_pA = round(mean(current_inj(1:1000,1))); 
%sweep_0pA = find(round(mean(current_inj(1000:2000,:))) == baseline_pA); %sweep w current injection = 0

%resting membrane potential
Vm = mean(sweeps(:,5)) 
%mean(sweeps(sweep_0pA,:)) 

%conversion factors
mV_conv = 1/1000; 
pA_conv = 1/1E6;
%input resistance
Ri = ((min(sweeps(5000:10000,5))-min(sweeps(5000:10000,1)))/(min(current_inj(5000:10000,5))-min(current_inj(5000:10000,1))))*mV_conv/pA_conv 

samplerate = Vtimes(2);

% AP properties
clear pks locs w p isi ap_widths

for i = 1:size(sweeps,2)
    [pks{i},locs{i},w{i},p{i}] = findpeaks(sweeps(:,i),'MinPeakHeight',0,'WidthReference','halfprom'); 
    try
        spike_rate(i) = 1000/(median(diff(locs{i}*samplerate))); %per step
    end
end

%hardcoded for same step, 15, @500 pA
ap_height = max(p{15}) %height of first AP
  
isis =  diff(locs{15}*samplerate); %isi, ms
min_isi = min(isis)
spike_adaptation_ratio = (1000/median(isis))/(1000/isis(1)) %reciprical for freq

half = (ap_height/2)+Vm;
hwhf = median(pulsewidth(sweeps(:,15),time,'StateLevels',[half-2 half+2]))

for i = 1:size(sweeps,2)
    %[dvdt(i),idx(i)] = max(diff(sweeps(:,i))/samplerate);
    [dv2dt2(i),idx2(i)] = max(diff(diff(sweeps(:,i)))/samplerate);
    ap_threshold(i) = sweeps(idx2(i),i);
end

%dvdt = dvdt(find(ap_widths_last));
ap_threshold = median(ap_threshold(9:15))

%sag @ -200pA
Vpeak = min(sweeps(:,1)) - Vm;
Vss = max(sweeps(9900:11000,1)) - Vm;
sag = (Vpeak-Vss)/Vpeak*100 %sag, percentage

%%
save(csvNames{1},"Vm","Ri","sag","spike_rate","spike_adaptation_ratio","ap_threshold","ap_height","hwhf","min_isi")
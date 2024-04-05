%%Glutamate uncaging analysis pipeline
%designed for experiments with grouping strategy for expected values
%assumes consecutive points (i.e. 1-5, 6-10, not necessarily groups of 5)
%assumes groups are in order (group 1 is first, group 2 next, etc).

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
% fluoNames = cell(0);
Vtimes = cell(0);
voltageTraces = cell(0);
currentTraces = cell(0);
fluoTimes = cell(0);
fluoTraces = cell(0);

%get all voltage and calcium traces
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
        %voltageTraces{end + 1} = data(:,2)*100; %in mV
        voltageTraces{end + 1} = smooth(data(:,2)*100,10); %in mV
        currentTraces{end + 1} = data(:,3)*1E4; %in pA %just to check amt of current inj

    elseif size(data,2) == 2
        %fluoNames{end + 1} = csvFile.name(1:end-25);
        fluoTimes{end + 1} = data(:,1);
        fluoTraces{end + 1} = data(:,2);
    else
        %disp('.csv file not in correct format. Skipped.');
    end
    
    clear('data');    
end

% normalize voltage Traces
for i = 1:ncycle
    base=mean(voltageTraces{i}([500:950])); %%%
    voltage_norm{i} = voltageTraces{i}-base;
end
        
%Get Ca data/parameters
%find maxima based on zscore (boundaries) for each cycle
%plot check
%average and store Ca signal

B = zeros(ncycle,2); %B holds boundaries for each cycle
% figure
P = {}; %linescan period for each cycle, in ms
F = {}; %average (mean) fluorescent trace, including artifact
Uncaging = {};

for i = 1:ncycle
    image = imread(strcat(csvNames{i}, '_Ch1_000001.ome.tif'));
    P{i}.scanLinePeriod = 1000*str2double(s.PVScan.Sequence{1,i}.Frame.PVStateShard.PVStateValue{1,4}.Attributes.value); %in ms
    Uncaging{i} = markpoints2([strcat(csvNames{i},'_MarkPoints') '.xml']);
    
    if size(image,1) > 10
       y1 = double(image);
       y = mean(zscore(double(image)'),2);
       A = find(round(y) == 1);
       if isempty(A)
           z = find(B(:,1)~=0);
           A = round(median(B(z,:)));
       end
       B(i,1) = min(A);
       B(i,2) = max(A);

%        clf
%        subplot(2,1,1)
%        imshow(image',[]);
%        hold on
%         line([0 size(image,1)],[min(A) min(A)],'Color','red');
%         line([0 size(image,1)],[max(A) max(A)],'Color','red');
%         title(sprintf('%d',i))
%         subplot(2,1,2)
%         plot(y1(:,min(A):max(A)))
%         hold on
%         plot(mean(y1(:,min(A):max(A)),2),'Color','r','LineWidth',1);
        
        F{i} = mean(y1(:,min(A):max(A)),2);

%        pause

    else
           F{i} = {};
           continue
    end
   
end

% close

% remove artifact and compute df/f
artifact = zeros(ncycle,2); %store artifact, not sure if useful
dF = {};

%for filtfilt
sampfreq=abs(1000/mean(fluoTraces{1}(3:end-1,1)-fluoTraces{1}(2:end-2,1)));
d = designfilt('lowpassfir', 'FilterOrder',4, ...
        'CutoffFrequency', 50, 'SampleRate', sampfreq);


for i = 1:ncycle
    if length(fluoTraces{i})>1
    limit1= floor((Uncaging{i}.InitialDelay)/(P{i}.scanLinePeriod));
    limit2= ceil(Uncaging{i}.InitialDelay/(P{i}.scanLinePeriod) + Uncaging{i}.TotalDuration/(P{i}.scanLinePeriod))+2; %%%add a few frames
        artifact(i,1) = limit1;
        artifact(i,2) = limit2;
    Fbaseline = mean(fluoTraces{i}(1:limit1));

    dF{i} = [fluoTraces{i}(1:limit1); fluoTraces{i}(limit2:length(fluoTraces{i}(:,1)))];
    dFtimes{i} = [fluoTimes{i}(1:limit1); fluoTimes{i}(limit2:length(fluoTimes{i}))];
    baseline = mean(dF{i}(15:limit1));
    dF{i} = filtfilt(d,(dF{i}-baseline)./baseline); 
    else 
        dF{i} = 0;
        continue
    end
end

% %%
%   %% plot ephys and dF/F with location of uncaging, sanity check
% figure
% for i=1:ncycle
%     clf
%     
%     if length(fluoTraces{i})>1
%     %Plot voltage
%     subplot(2,1,1)
%     %h=plot(Vtimes{i}, voltage_norm{i});
%     h=plot(Vtimes{i}, voltageTraces{i});
%     hold on;
%     h.LineWidth=1.5;
%     h.Color='b';
%     xlabel('ms')
%     ylabel('mV')
%     title('Voltage recording')
% 
%     %Add uncaging info on plot
%     dim = [0.65 0.6 0.3 0.3];
%     str = {strcat('Cycle = ',num2str(i)),strcat('Laser Power = ',num2str((Uncaging{i}.UncagingLaserPower))),strcat('Number of points = ',num2str((Uncaging{i}.nPoints))),strcat('Duration = ',num2str((Uncaging{i}.Duration))),strcat('InterPointDelay = ',num2str((Uncaging{i}.InterPointDelay)))};
%     a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
%     a.FontSize=8;
% 
%     %Add arrows indicating uncaging times
%     y=min(voltage_norm{i}); 
%     str1 = strcat('\uparrow ');%, num2str(j));
%     t=text([Uncaging{i}.InitialDelay,(Uncaging{i}.InitialDelay+Uncaging{i}.TotalDuration)],[y,y],str1,'FontSize',10);
% 
%     %Plot df/f
%     subplot(2,1,2)
%     h=plot(dFtimes{i},dF{i});
%     h.LineWidth=1.5;
%     h.Color='g';
%     xlabel('ms')
%     ylabel('df/f')   
%     title('Fluorescence')
% 
%     %Add arrows indicating uncaging times
%     y=min(dF{i});
%     t=text([Uncaging{i}.InitialDelay,(Uncaging{i}.InitialDelay+Uncaging{i}.TotalDuration)],[y,y],str1,'FontSize',10);
% 
%     pause
%     
%     else
%         continue
%     end 
% end 
% close
%             
%% Compute EPSP amplitude and area, Ca amplitude
%assumes constant sampling rate
%matrix, ucycles: (1) = nPoints, (2) = EPSP max, (3) = Ca
%max, (4) = EPSP area, (5) Ca area, (6), laser power

ucycles = zeros(6,ncycle);

 for i = 1:ncycle
    ucycles(1,i) = Uncaging{i}.nPoints;
    temp = mean(voltage_norm{i}(1:floor(Uncaging{2}.InitialDelay/(0.05)))); %baseline, prior to uncaging %may not need this with traces normalized
    %Get max & area of EPSP, after stimmulus is presented 
    timelimits = Uncaging{i}.InitialDelay:Uncaging{i}.InitialDelay+50; %%arbitrary window? Shorter linescans may be affected
    ucycles(2,i) = max(voltage_norm{i}(timelimits./.05))-temp;  
    ucycles(4,i)=trapz(Vtimes{i}(timelimits./.05),(voltage_norm{i}(timelimits./.05)-temp));

    %df/f max & area
    if length(dFtimes{i})>1
        ucycles(3,i)=max(dF{i}(ceil(timelimits(1:end-7)./P{i}.scanLinePeriod)));
        ucycles(5,i)= trapz(dFtimes{i}(ceil(timelimits(1:end-7)./P{i}.scanLinePeriod)),dF{i}(ceil(timelimits(1:end-7)./P{i}.scanLinePeriod)));
    else
       ucycles(3,i) = 0; 
       ucycles(5,i) = 0;
    end
    
    %laser powa
    ucycles(6,i)= Uncaging{i}.UncagingLaserPower;
 end

if mean(size(unique(ucycles(6,:))))~=1
    disp('Ya un laser powa change pendant la manip, fait gaffe a ce que tu fais bordel')
    ucycles(6,:)
else 
    disp('un laser powa ne change pas, dansez!')
end
% Sort and select groups; find cycles with single groups or multiple groups 
%lsequels: enter matrix of chosen cycles ex [2,4,5] or [2:4]
indx = cell(2,ncycle);
for i = 1:ncycle
    a = reshape(char(Uncaging{i}.GroupName),1,[]);
    indx{1,i} = str2double(regexp(a,'\d','Match'));  %group names; indexing this doesn't work for multiple groups because it counts each double...
    indx{2,i} = size(indx{1,i},2); %count how many groups 
    indx{3,i} = Uncaging{i}.nPoints; %how may points
end

A = [];
A(1,:) =  find([indx{2,:}] == 1 & [indx{3,:}]~=max(ucycles(1,:))); %single trial cycles, take out any with all points 
A(2,:) = [indx{1,A}]; %single trial groups
aa = unique(A(2,:)); %name of single groups
onegrup = {};
%group same single trials by index and group number/name
for i = 1:length(aa)
    b = find(A(2,:) == aa(i));
    onegrup{1,i} = A(1,b); %cycle index
    onegrup{2,i} = aa(i); %group name
    onegrup{3,i} = indx{3,onegrup{1,i}}; %no. of points
end

A = [];
B = {};
aa = [];
b = [];
A = find([indx{2,:}] > 1 | [indx{3,:}]== max(ucycles(1,:))); %find multiple groups or with all points
B = indx(1,A); %multiple group names
[a1,b,c] = unique(cellfun(@char,B,'un',0));
aa = B(b); %unique group names

multgrup = {};
%group same multiple trials by index and group number/name
for i = 1:length(b)
    d = find(c == i);
    multgrup{1,i} = A(d); %cycle index
    multgrup{2,i} = cell2mat(aa(i)); %groups run
    multgrup{3,i} = indx{3,multgrup{1,i}}; %no. of points
end
%% quick run

select_onegrup = onegrup;
select_multgrups = multgrup;
%%
%probably a better way to do this, but
xx = {};
if length(find(cell2mat(multgrup(3,:))== max(cell2mat(multgrup(3,:)))))>1 %check for groups with all points run as two different groups
    x = find(cell2mat(multgrup(3,:))== max(cell2mat(multgrup(3,:))));
    xx{1,1} = cell2mat(multgrup(1,x));
    [x1,y1] = max(cellfun('size',multgrup(2,:),2));
    xx{2,1} = cell2mat(multgrup(2,y1));
    xx{3,1} = max(cell2mat(multgrup(3,:)));
    multgrup(:,x) = [];
    multgrup = [multgrup xx];
end

select_onegrup = {};
j = .1:.1:3; %arbitrary units, just for text placememt, could be better but eh
for h = 1:size(onegrup,2)
    figure, hold on
    for ii = 1:length(onegrup{1,h})   
        f(ii) = plot(Vtimes{onegrup{1,h}(ii)}, voltage_norm{onegrup{1,h}(ii)});
        prut = get(f(ii), 'Color');
        idx1 = max(voltage_norm{onegrup{1,h}(ii)});
        text(100,idx1+j(ii),num2str(onegrup{1,h}(ii)),'Color',prut,'FontSize',20)
        text(150,idx1+j(ii),num2str(ucycles(6,onegrup{1,h}(ii))),'Color',prut,'FontSize',12)
    end   
    if length(onegrup{1,h})== 1
        select_onegrup{1,h} = onegrup{1,h}; 
        select_onegrup{2,h} = onegrup{2,h}; 
        select_onegrup{3,h} = onegrup{3,h}; 
    else
        select_onegrup{1,h} = input('lesquels')
        select_onegrup{2,h} = onegrup{2,h}; 
        select_onegrup{3,h} = onegrup{3,h}; 
    end
    close
end

select_multgrups = {};
for h = 1:size(multgrup,2)
    figure, hold on
    for ii = 1:length(multgrup{1,h})   
        f(ii) = plot(Vtimes{multgrup{1,h}(ii)}, voltage_norm{multgrup{1,h}(ii)});
        prut = get(f(ii), 'Color');
        idx1 = max(voltage_norm{multgrup{1,h}(ii)});
        text(100,idx1+j(ii),num2str(multgrup{1,h}(ii)),'Color',prut,'FontSize',20)
        text(150,idx1+j(ii),num2str(ucycles(6,multgrup{1,h}(ii))),'Color',prut,'FontSize',12)
    end    
    if length(multgrup{1,h})== 1
        select_multgrups{1,h} = multgrup{1,h}; 
        select_multgrups{2,h} = multgrup{2,h}; 
        select_multgrups{3,h} = multgrup{3,h}; 
    else
        select_multgrups{1,h} = input('lesquels')
        select_multgrups{2,h} = multgrup{2,h}; 
        select_multgrups{3,h} = multgrup{3,h}; 
    end
    close
end
%% Get measured and expected data
%measured, combined average of each trace

onegrup_results = zeros(6,size(select_onegrup,2));
for i = 1:size(select_onegrup,2)
    onegrup_results(1,i) = mean(ucycles(2,cell2mat(select_onegrup(1,i)))); %EPSP amp
    onegrup_results(2,i) = mean(ucycles(3,cell2mat(select_onegrup(1,i)))); %Ca amp
    onegrup_results(3,i) = mean(ucycles(4,cell2mat(select_onegrup(1,i)))); %EPSP area
    onegrup_results(4,i) = mean(ucycles(5,cell2mat(select_onegrup(1,i)))); %Ca area
    onegrup_results(5,i) = cell2mat(select_onegrup(2,i)); %grp name
    onegrup_results(6,i) = cell2mat(select_onegrup(3,i)); %no. of points
end

%expected, mV
N = cellfun(@numel, voltage_norm); %would be nice to only run this if the recording lengths are uneven
N = N(cell2mat(select_onegrup(1,:)));
len = min(N(find(N>1)));
singles = cell2mat(select_onegrup(1,:));
for i = 1:length(singles) %make recordings same length
    cut = voltage_norm{singles(i)}(1:len);
    cuttime = Vtimes{singles(i)}(1:len);
    voltage_norm{singles(i)} = cut;
    Vtimes{singles(i)} = cuttime;
end

expected_traces = [];
grup_avg = [];
duration = [];
for i = 1:size(select_onegrup,2)
        grup_avg{i} = mean([voltage_norm{select_onegrup{1,i}}],2); %make a trace from cycle avgs
        duration(i) = Uncaging{select_onegrup{1,i}(1,1)}.TotalDuration;
end

dFdelay = [];
for i = 1:size(select_onegrup,2)
    dFdelay(i)  = sum(duration(1:i));
end
delay = [0 dFdelay./.05];

for i = 1:size(grup_avg,2)
    %grup_avg{i} = delayseq(grup_avg{i},delay(i)); 
    grup_avgmax(i) = max(grup_avg{i}(timelimits./.05));
end

expected_traces = {};
expected = zeros(4,size(grup_avg,2));
expected_traces{1} = mean([grup_avg{:}],2); %average of all single groups, w/o respect to no. of points 
expected(1,1) = mean(grup_avgmax); %first point, average amp of compiled traces
expected(3,1) = trapz(expected_traces{1}(timelimits./.05));
for i = 2:size(grup_avg,2)
    expected_traces{i} = sum([grup_avg{1:i}],2);
    expected(1,i) = max(expected_traces{i}(timelimits./.05));
    expected(3,i) = trapz(expected_traces{i}(timelimits./.05));
end

% expected, dF %not sure if this is really necessary... 
N = cellfun(@numel, dF);
N = N(cell2mat(select_onegrup(1,:)));
len = min(N(find(N>1)));
singles = cell2mat(select_onegrup(1,:));

for i = 1:length(singles) 
    cut = dF{singles(i)}(1:len);
    cuttime = dFtimes{singles(i)}(1:len);
    dF{singles(i)} = cut;
    dFtimes{singles(i)} = cuttime;
end

dFexpected_traces = [];
dFgrup_avg = [];
for i = 1:size(select_onegrup,2)
    dFgrup_avg{i} = mean([dF{select_onegrup{1,i}}],2);
end

dFdelay = [0 dFdelay];

for i = 1:size(dFgrup_avg,2)
    %dFgrup_avg{i} = delayseq(dFgrup_avg{i},delay(i)); 
    dFgrup_avgmax(i) = max(dFgrup_avg{i}(timelimits));
end

dFexpected_traces = {};
dFexpected_traces{1} = mean([dFgrup_avg{:}],2); %average of all single groups, w/o respect to no. of points 
expected(2,1) = mean(dFgrup_avgmax); %first point
expected(4,1) = trapz(dFexpected_traces{1}(timelimits));
for i = 2:size(dFgrup_avg,2)
    dFexpected_traces{i} = sum([dFgrup_avg{1:i}],2);
    expected(2,i) = max(dFexpected_traces{i}(timelimits));
    expected(4,i) = trapz(dFexpected_traces{i}(timelimits));
end


%measured, both mV and dF
measured = [];
for i = 1:size(select_multgrups,2)
    measured(1,i) = mean(ucycles(2,cell2mat(select_multgrups(1,i)))); %EPSP amp
    measured(2,i) = mean(ucycles(3,cell2mat(select_multgrups(1,i)))); %Ca amp
    measured(3,i) = mean(ucycles(4,cell2mat(select_multgrups(1,i)))); %EPSP area
    measured(4,i) = mean(ucycles(5,cell2mat(select_multgrups(1,i)))); %Ca area
end

% Plot expected vs measured by no. of inputs

EPSPamp = [];
dFamp = [];
EPSParea = [];
dFarea = [];
gain = zeros(4, size(expected,2));
%EPSP amp
fp = mean(onegrup_results(1,1)); %get the expected value for one group 
no_inputs = [round(mean(onegrup_results(6,:))) cell2mat(select_multgrups(3,:))]; %avg no. of inputs per group
m = [fp measured(1,:)];
EPSPamp = [expected(1,:); m; no_inputs]; %expected, measured, and no of inputs    
gain(1,:) = EPSPamp(2,:)./EPSPamp(1,:)
mV1 = figure('Name', 'mV amp');
plot(EPSPamp(3,:),EPSPamp(2,:)) %measured
hold on, plot(EPSPamp(3,:),EPSPamp(1,:)) %expected
xlabel('Number of inputs')
ylabel('mV')
axis([min(no_inputs) max(no_inputs) 0 20])

%dF amp
fp = mean(onegrup_results(2,:));
m = [fp measured(2,:)];
dFamp = [expected(2,:); m; no_inputs];
gain(2,:) = dFamp(2,:)./dFamp(1,:);
dF1 = figure('Name', 'dF amp');
plot(dFamp(3,:),dFamp(2,:))
hold on, plot(dFamp(3,:),dFamp(1,:)) 
xlabel('Number of inputs')
ylabel('dF/F')

% %EPSP area
% fp = mean(onegrup_results(3,:)); 
% m = [fp measured(3,:)];
% EPSParea = [expected(3,:); m; no_inputs];
% gain(3,:) = EPSParea(2,:)./EPSParea(1,:);
% mV2 = figure('Name', 'mV area')
% plot(EPSParea(3,:),EPSParea(2,:)) 
% hold on, plot(EPSParea(3,:),EPSParea(1,:)) 
% xlabel('Number of inputs')
% ylabel('mV*s')
% 
% %dF area
% fp = mean(onegrup_results(4,:)); 
% m = [fp measured(4,:)];
% dFarea = [expected(4,:); m; no_inputs];
% gain(4,:) = dFarea(2,:)./dFarea(1,:);
% dF2 = figure('Name', 'dF area')
% plot(dFarea(3,:),dFarea(2,:)) 
% hold on, plot(dFarea(3,:),dFarea(1,:))
% xlabel('Number of inputs')
% ylabel('dF/F*s')


g1 = figure; hold on;
plot(EPSPamp(3,:), gain(1,:))
xlabel('Number of inputs')
ylabel ('gain')
axis([min(no_inputs) max(no_inputs) 0 2])
savefig(g1,'gain')
saveas(g1,'gain.png')
%% Pretty traces
%plot expected traces
r = zeros(size(expected_traces,2),3);
q = figure('Name','exp_mV'), hold on
title('Expected mV')
for h = 1:size(expected_traces,2)
    r(h,:) = 0 + 1*rand(1,3);
    plot(Vtimes{2}, expected_traces{h},'LineWidth',1,'Color',r(h,:))
end
xlabel('ms')
ylabel('mV')

%plot measured EPSP and dF traces, color by group
q2 = figure('Name','mV'); hold on
title('Measured mV')
plot(Vtimes{select_onegrup{1,2}(1,1)}, voltage_norm{select_onegrup{1,2}(1,1)},'LineWidth',1,'Color',r(1,:))
for h = 1:size(select_multgrups,2)
    for ii = 1:length(select_multgrups{1,h})   
        plot(Vtimes{select_multgrups{1,h}(ii)}, voltage_norm{select_multgrups{1,h}(ii)},'LineWidth',1,'Color',r(h+1,:))
    end
end
xlabel('ms')
ylabel('mV')

q3 = figure('Name','dF'); hold on
title('Measured dF')
plot(dFtimes{select_onegrup{1,1}(1,1)}, dF{select_onegrup{1,1}(1,1)},'LineWidth',1,'Color',r(1,:))
for h = 1:size(select_multgrups,2)
    for ii = 1:length(select_multgrups{1,h})   
        plot(dFtimes{select_multgrups{1,h}(ii)}, dF{select_multgrups{1,h}(ii)},'LineWidth',1,'Color',r(h+1,:))
    end
end
xlabel('ms')
ylabel('dF/F')

%%
save(dirName)%,'gain', 'EPSPamp', 'EPSParea', 'dFamp', 'dFarea') %others?
savefig(q, 'exp_mV')
savefig(q2, 'mV')
savefig(q3, 'dF')
savefig(mV1, 'mVvsInput')
%savefig(mV2, 'mVsvsInput')
savefig(dF1, 'dFvsInput')
%savefig(dF2, 'dFsvsInput')

saveas(q, 'exp_mV.png')
saveas(q2, 'mV.png')
saveas(q3, 'dF.png')
saveas(mV1, 'mVvsInput.png')
%%

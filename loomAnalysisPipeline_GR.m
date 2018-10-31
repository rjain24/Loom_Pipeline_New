%%
% Gokul Rajan, Institut Curie.

% August 2018: free-swim script adapted for analysing tap-induced escapes

% working implementations:
    % single fish analysis
    % find dish center and cut out edges - CHANGE HERE (currently dish center is
    % hard-coded @505/505) - NOT REQ IN TAP
    % replace nan frames
    % orientation correction -- verified on 17/05/18
    % bout-wise analysis
    % tap stim
    
    % #8 and 17 are important for tap associated Vm and Wm resp
    % peak translational velo and peak angular velocity (avg of 10 ms window around the peak velocity)

% Sept 2018- to be adapted for analysing loom-induced escapes

%% Figures on/off

set(0, 'DefaultFigureVisible', 'on');

%%

%% loading data, cleaning it and reshaping

% Loom file Data key
% 1 frame #
% 2 camera time in microsec (for data transfer purposes microseconds instead of nanoseconds used)
% 3 system time in ms
% 4 xpos, start at 0
% 5 ypos
% 6 angle in rad
% 7 stimulus (0=off, 1=on)
clearvars; clc;


ff=1;

disp('load your Loom file');
[FileName, PathName] = uigetfile('*.bin'); %read fish params (loom file)
%cd(PathName);

disp('load your StimFile');
[StimName, stimPathName] = uigetfile('*.bin'); %read stim params (stimFile)

% extract parameters from filename
tmp_str = strsplit(FileName, '_');

% save parameters in strings - FORMAT: Date_Time_ExperimentName_Animal#_Remark_Trial#_.bin
acquis_date = tmp_str{1, 1}; acquis_time = tmp_str{1, 2}; exp_type    = tmp_str{1, 3}; fish_num    = tmp_str{1, 4}; trial_num   = tmp_str{1, 6}; trial_comment=tmp_str{1, 5};

% acquisition parameters
window_size         = 0;
num_data_categories = 7+window_size^2; %6 for swim / 7 for loom-only / 7 for tap-only
camscale_px_per_mm  = 20.2; % px/mm changed on 180926 from 20.6
datarate_Hz         = 750;  % Hertz
NumVar = 7; %6 for swim / 7 for loom-only / 7 for tap-only

%%

% Read and reorganize the "Loom" .bin file
h        = fopen([PathName, FileName]);
tmp_data = fread(h, inf, 'float');
fclose(h);

%StimFile data key
% 1) frameCount
% 2) loom R
% 3) stim timer
% [6 adtl variables - these were removed 10/05/18 RJ]

% Read and reorganize the stim file
lam        = fopen([stimPathName, StimName]);
stm_data = fread(lam, inf, 'float');
fclose(lam);

tmp_data = tmp_data(1:(end-mod(size(tmp_data,1), num_data_categories)), 1); % cuts when parts of the data categories are missing at the end
tmp_data = reshape(tmp_data, [num_data_categories, size(tmp_data, 1)/num_data_categories])';

freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.datarate_Hz = datarate_Hz;

freeSwim.boutAnalysis.boutLength = 300; % considering the bout length to be 200ms / 750Hz acq
freeSwim.boutAnalysis.VthreshON = 1;
freeSwim.boutAnalysis.BthreshOFF = 2;

%%
% CHECK FOR TIMING PROBLEMS AND LOST FRAMES

% TIMER COUNTER: 
% time difference between frames in microseconds, based on the cameras 32bit time stamp counter (25Mhz, 40ns)

time_diff      = [0; diff(tmp_data(:, 2))];                      

% linearize the "saw-tooth function" for the timecounter
idx            = time_diff <= -2^32/1000 + 1.5*median(time_diff); % find the frames when 32bit timecounter was reset to zero
time_diff(idx) = median(time_diff);                               % reset the time difference between these frames to the median in microseconds

% camera frame length in microseconds calculated from the camera timecounter
frame_length_calc_ms = median(time_diff)/1000;

% aquisition datarate in Hertz calculated from the camera timecounter
datarate_calc_Hz = (1/(median(time_diff)).*10^6);  

% CHECK for timing problems (e.g. frames that took unusually long or are
% unusually shorter than what is expected from the used datarate)

idx_time       = abs(time_diff)/1000 >= 1.01*frame_length_calc_ms;  % searches for time differences between frames that are +-1% of the expected frame duration


% FRAME COUNTER:
% index difference between frames, based on the cameras 24bit frame counter

frame_diff = [0; diff(tmp_data(:, 1))]; 

% linearize the "saw-tooth function" for the frame counter (should not
% happen at low datarates) 
idx             = frame_diff <= -2^24 + 1.5*median(frame_diff); % find the frames when 24bit framecounter was reset to zero
frame_diff(idx) = 1;                                            % reset the frame difference between these frames to 1

% CHECK for missing frames
idx_frame  = frame_diff > 1;                         % index of missing frames
idx_lost   = find(idx_frame == 1);                   % first frame in the block of missed frames

% checks if missed timestamps coincide with missed frames, is 0 if inconsistent timestamps outside of missed frames 
isTime = isequal(idx_time, idx_frame);

% calculate total duration of video

duration = tmp_data(end,1)/datarate_Hz/60;

% implememnt this for loom
% idx_tap = tmp_data(:,7)==1;
% idx_tap = find(idx_tap==1);

loomEscape(ff).datarate_Hz = datarate_Hz;

loomEscape(ff).boutAnalysis.boutLength = 300; % considering the bout length to be 200ms / 750Hz acq
loomEscape(ff).boutAnalysis.VthreshON = 1;
loomEscape(ff).boutAnalysis.BthreshOFF = 2;

% prints the above calculated values

fprintf('\nacquisition time: %s %s', acquis_date, acquis_time); 
fprintf('\nexperiment type: %s', exp_type); 
fprintf('\nfish number: %s', fish_num); 
fprintf('\ntrial number: %s', trial_num);
fprintf('\nvideo duration: %2.2f min', duration);
fprintf('\n\nfirst frame in the block of missed frames : number of frames lost\n');
fprintf('\n %d: %d',  [idx_lost, frame_diff(idx_lost)].');
fprintf('\n\ntiming flawed (outside of lost frames):  %d  \n', ~isTime  );

% %% dish center
% 
% xpos = tmp_data(:, 4);
% ypos = tmp_data(:, 5);
% 
% snap = imread('petriplate.jpg');
% dish_center = determine_dish_centre(snap,430);
% tmp_radialloc = sqrt((xpos - dish_center(1) ).^2 + (ypos - dish_center(2)).^2);
% tmp_inmiddle  = tmp_radialloc < 430; %depends on pixel resolution/ plate diameter
% idx_edge = find(tmp_inmiddle==0);
% tmp_data(idx_edge,2:end) = nan;



%%
% % INSERT nans for lost frames... 

% define anonymous function that inserts (nxm) blocks into (oxm) matrices
insert_blocks = @(insert_block, matrix, n) cat(1,  matrix(1:n-1,:), insert_block, matrix(n:end,:) );

data_raw = tmp_data;

for ii = nnz(idx_frame):-1:1 % starts from the last row in the matrix to keep the indizes for block-insertion
 
    nan_block       = nan(frame_diff(idx_lost(ii)) - 1, num_data_categories);
    nan_block(:, 1) = tmp_data(idx_lost(ii)-1, 1)+1: tmp_data(idx_lost(ii)-1, 1) + frame_diff(idx_lost(ii))-1; % fill the first column of the Nan blocks with frame numbers that were missing
    
    tmp_data        = insert_blocks(nan_block, tmp_data, idx_lost(ii));
    
end

tmp_data(:,1) = tmp_data(:,1) - tmp_data(1,1) + 1; % framecounter starts at 1


% % read the fish images
% 
% % Read the 3600px associated with each frame and store it row-wise
% tmp_CXP = NaN(size(tmp_data,1),size(tmp_data,2)-NumVar);
% for gg = 1:size(tmp_data,1)
%     pp = 1;
%     for rr = 7:((size(tmp_data,2))-NumVar)
%         tmp_CXP(gg,pp) = tmp_data(gg,rr);
%         pp=pp+1;
%     end
% end
% 
% tmp_CXP = tmp_CXP'; %Transpose to column-wise
% 
% % reshape each 3600 column into a 60x60 px matrix and join each new consecutive
% % frame onto it
% CXPimg = zeros(100,100,(size(tmp_CXP,2)));
% for tt = 1:size(tmp_CXP,2)
%     res = reshape(tmp_CXP(:,tt),100,[])';
%     CXPimg(:,:,tt) = res;
%     f(tt)=im2frame(uint8(CXPimg(:,:,tt)),gray(256)); % create frame for video writing
% end
% figure(3);imshow3D(CXPimg);
% 
% % video conversion
% 
% vidObj = VideoWriter('awesomeDanionella.avi');
% vidObj.FrameRate=20;
% vidObj.Quality=100;
% open(vidObj)
% writeVideo(vidObj,f);
% close(vidObj);


%% to do: get the timestamps of stim onset (and also grab the corresp R)
% Identify onset times of looming stimuli [RJ 10-18-18]

tmp_data_nans = find(isnan(tmp_data(:,7)));
tmp_data(tmp_data_nans,:)=0;
raw_loom=tmp_data(:,7);

%extracts out the booleans for stimuli
loomwin = movsum(raw_loom, [0 19]); %calc sliding window of 20 frames that "real" stimuli must be on for
loomon = find(raw_loom); %locate all frames where stimulus was on
loomstart=[];
for cloom = 1:length(loomon)
    if cloom==1
    if loomon(cloom)==1 && loomwin(loomon(cloom))==20 %if loom starts at 1 and is on for at least 20 frames
        loomstart = loomon(cloom); %store loomon(cloom) as the first value in loomstart
    end
    elseif cloom>1
    if (loomon(cloom)-loomon(cloom-1))>1 && loomwin(loomon(cloom))==20 %if there's a 2sec stretch with no stim before the stim started
        loomstart = [loomstart; loomon(cloom)]; %store loomon(cloom+1) as the next value in loomstart
    end
    end
end

tmp_data(tmp_data_nans,:)=NaN;

%%
%-------------------------------------------------------------------------------------
%
% IDENTIFICATION OF SWIM BOUTS (bout speeds, IBIs etc)
%
%-------------------------------------------------------------------------------------

xpos = tmp_data(:, 4);
ypos = tmp_data(:, 5);
tmp_ori_deg = rad2deg(tmp_data(:,6));

% plot fish position 
xx = [xpos'; xpos']; tt = [ypos'; ypos'];
z  = 1:size(xpos', 2); zz = [z; z];       % create frame-vector

figure
hs = surf(xx, tt, zz, zz, 'EdgeColor', 'interp', 'LineWidth', 2);
colormap('parula');
view([0 90 0]); 
xlabel('x-position');
ylabel('y-position');
zlabel('frames');

plot(xpos,ypos);

dxV   = [0; diff(xpos)]; % distance between two consecutive x-coordinates
dy   = [0; diff(ypos)]; % distance between two consecutive y-coordinates

tmp_dist_unfilt           = sqrt(dxV.^2 + dy.^2);
tmp_dist_unfilt           = tmp_dist_unfilt./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate_Hz;         % convert to velocity in mm/s

%% 
idx_nan     = isnan(dxV);
dxV(idx_nan) = 0; % for filtering nan values need to be removed 
dy(idx_nan) = 0;

% filters used in the analysis
filterB = ones(1,100)./100; %for event detection 100 ms
filterV = ones(1,30)./30; %for more precise onset detection 30 ms
filterF = ones(1,4)./4; %for escapes 4.5 ms

freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.filters.B = filterB;
freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.filters.V = filterV;
freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.filters.F = filterF;

%% orientation

tmp_delta_ori = [nan; diff(tmp_ori_deg)];

% correction of discontinuties when fish turns from 0 to 2*pi or vice versa

for kk = 1: length(tmp_delta_ori)
    
    if tmp_delta_ori(kk) > 180 % right turn
        tmp_delta_ori(kk) =  tmp_delta_ori(kk) - 2*180;
        
    elseif tmp_delta_ori(kk) < -180 % left turn
        tmp_delta_ori(kk) =  2*180 + tmp_delta_ori(kk);
    end
    
end

%% Event detection
dxB        = filtfilt(filterB, 1, dxV);
dyB        = filtfilt(filterB, 1, dy);

tmp_dist_fB = sqrt(dxB.^2 + dyB.^2);     % distance moved between iterations, in pixels
tmp_dist_fB = tmp_dist_fB./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fB  = tmp_dist_fB.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fB(idx_nan) = nan; % re-insert the nan values

loomEscape(ff).boutAnalysis.acquis_date.exp_type.fish_num.filters.B = filterB;
loomEscape(ff).boutAnalysis.acquis_date.exp_type.fish_num.filters.V = filterV;
loomEscape(ff).boutAnalysis.acquis_date.exp_type.fish_num.filters.F = filterF;

% plot(tmp_vel_fB);hold on;vline(loomstart);hold off;
% gcf;

%% Onset detection
dxV        = filtfilt(filterV, 1, dxV);
dyV        = filtfilt(filterV, 1, dy);

tmp_dist_fV = sqrt(dxV.^2 + dyV.^2);     % distance moved between iterations, in pixels
tmp_dist_fV = tmp_dist_fV./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fV  = tmp_dist_fV.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fV(idx_nan) = nan; % re-insert the nan values

%% Peak detection
dxF        = filtfilt(filterF, 1, dxV);
dyF        = filtfilt(filterF, 1, dy);

tmp_dist_fF = sqrt(dxF.^2 + dyF.^2);     % distance moved between iterations, in pixels
tmp_dist_fF = tmp_dist_fF./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fF  = tmp_dist_fF.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fF(idx_nan) = nan; % re-insert the nan values

%% orientation filter

% removing nans for filtering orientation
tmp_delta_ori(idx_nan) = 0;

tmp_delta_ori_filtered = tmp_delta_ori;
tmp_delta_ori_filtered(isnan(tmp_delta_ori_filtered))=0;
windowWidth = 25; %larger window leads to more averaging --> depends on your signal
polynomialOrder = 3; %larger order --> less smotthing
tmp_delta_ori_filtered = sgolayfilt(tmp_delta_ori_filtered, polynomialOrder, windowWidth);
plot(tmp_delta_ori_filtered);
hold on;
% another small window filtering to further smoothen the signal, useful for
% finding local maximas and minimas.
windowWidth2 = 10;
polynomialOrder2 = 3;
tmp_delta_ori_filtered2=sgolayfilt(tmp_delta_ori_filtered, polynomialOrder, windowWidth);
plot(tmp_delta_ori_filtered2);title('check filtered output before proceeding');
hold off;

%% new - 12 june
tmp_ang_vel = (tmp_delta_ori_filtered2*datarate_Hz); % per sec
tmp_ang_vel2 = sgolayfilt(tmp_ang_vel, polynomialOrder, windowWidth);

%%
%re-inserting nans
tmp_delta_ori(idx_nan) = NaN;
tmp_delta_ori_filtered(idx_nan) = NaN;
tmp_delta_ori_filtered2(idx_nan) = NaN;
tmp_ang_vel(idx_nan) = NaN;

plot(tmp_ang_vel-10); title('ang velo compare');
hold on;
plot(tmp_vel_fB);
hold off;

%% new - 20180829 - to ID tap idx and look for peaks in those regions.
% this was not used in the end to compute stim-related events
% id the tap signals
idx_tap=loomstart;

% create ROI for taps - not in use
for ullu=1:length(idx_tap)
tapRegion(:,ullu)=idx_tap(ullu):(idx_tap(ullu)+2000);
end

% tapRegion_Master =[];
% 
% for gkl=1:size(tapRegion,2)
% tapRegion_Master=[tapRegion_Master; tapRegion(:,gkl)];
% end

plot(tmp_vel_fF);hold on;vline(loomstart);hold off;
gcf;

%% find peaks in the velocity
[pks,locsz] = findpeaks(tmp_vel_fB,'MinPeakProminence',6,'MinPeakDistance',100); %minPeakProminence & minPeakDist as bout interval and min Velo resp

escapeLocs=locsz(pks>30);
escapePks=pks(pks>30);

%% new - 20180829 - verify if locs are associated with taps and operate acc

% create a draft RT matrix - (all taps - all escapes in tmp_difftap)
for kkhh=1:length(idx_tap)
tmp_diffTap(:,kkhh)=idx_tap(kkhh)-locsz; 
end

% now filter out the matrix to retain only the best (as described below) escape 
% velocity associated with each tap in idx_realEscape (Pn)

klm=1;
for sss=1:size(tmp_diffTap,2)
    for aaa=1:size(tmp_diffTap,1)
        if (tmp_diffTap(aaa,sss)<=0 && abs(tmp_diffTap(aaa,sss))==min(abs(tmp_diffTap(:,sss))) && abs(tmp_diffTap(aaa,sss))<=3000)
            idx_realEscape(klm)=aaa;
            klm=klm+1;
        end
    end
    idx_realEscape(klm)=inf; %inf after every column and use this info to find non-resp fish
    klm=klm+1;
end

idx_realEscape(1:1:end)

% remove the infs
countTap=1;
for iiff=1:length(idx_realEscape)
   if idx_realEscape(iiff)==inf
       if(iiff>1)
       idx_corrTap(countTap)=idx_realEscape(iiff-1);
       countTap=countTap+1;   
       elseif(iiff==1)
       idx_corrTap(countTap)=inf;
       countTap=countTap+1;
       end
   end
end


% insert NaN for missed taps (no rxn from animal)
idx_corrTap(find(isinf(idx_corrTap)))=NaN;
idx_corrTap2=idx_corrTap;
idx_corrTap2(find(isnan(idx_corrTap)))=[];

idx_escape=locsz(idx_corrTap2);
fast_locs=locsz(idx_corrTap2);

%tmp_diffTap=sort(tmp_diffTap,'descend');

% hojaye=1;
% locs=[];
% 
% for jugaad=1:length(locsz)
%     if any(tapRegion_Master)==locs(jugaad)
%         locs(hojaye)=locsz(jugaad);
%         hojaye=hojaye+1;
%     end
% end

%% to-do: find all the bouts and then use the above escape bout info to extract our escape-related kinematics
% done on 5th oct 2018
%%

if locsz>=1
    
for qq = 1:size(locsz)
    
    if (locsz(qq) - 3*loomEscape(ff).boutAnalysis.boutLength <=0)... % get rid of any half bout in the begining
                ||(locsz(qq) + 3*loomEscape(ff).boutAnalysis.boutLength >=size(tmp_vel_fB,1))... % half bout in the end
                ||(any(isnan(locsz(qq):locsz(qq)+3*loomEscape(ff).boutAnalysis.boutLength)))... % remove nan-area
                ||(any(isnan(locsz(qq):locsz(qq)-3*loomEscape(ff).boutAnalysis.boutLength)))... % nan-area again
                ||tmp_vel_fB(locsz(qq))<= loomEscape(ff).boutAnalysis.VthreshON... %check w/ Vthresh
                ||tmp_vel_fB(locsz(qq))<= loomEscape(ff).boutAnalysis.BthreshOFF %check w/ Bthresh
            
            locsz(qq)= NaN;
    end
end


%% ID swim bouts and operate on them

%bout key
%1 start frame
%2 end frame
%3 bout duration in ms
%4 orientation before bout in deg
%5 orientation after bout in deg
%6 turn in deg

%7 mean velo
%8 peak velo
%9 total distance
%10 total L yaw in deg
%11 total R yaw in deg
%12 angular velocity in deg/sec

%13 inter-bout-interval in ms

%14 local bend minima? [RJ]
%15 local bend maxima? [RJ]
%16 total half bends? [RJ]
%17 peak angular velocity? Is this the peak or the average angular velocity? [RJ]
%18 ??? Tail Beat Frequency associated velocity? [RJ]
%19 ??? Reaction time? [RJ]

locsz(isnan(locsz)) = []; %imp - esle, error w/ for loop - sunbscript indices must be real integers or logicals

tmp_swim_bouts = zeros(size(locsz,1),19); %Changed 15 to 19 [10-11-18 RJ]

for mm = 1:size(locsz,1)
   
    %find start frame
    sss = 0;
    while locsz(mm)-sss>=2*loomEscape(ff).boutAnalysis.boutLength...
            && tmp_vel_fV(locsz(mm)-sss) >= loomEscape(ff).boutAnalysis.VthreshON % using narrow filter for onset
        sss = sss + 1;
    end   
    
    tmp_swim_bouts(mm,1) = locsz(mm)-sss+1; %bout start frame
    
    
    %find end frame
    eee = 0;
    while locsz(mm)<=size(tmp_vel_fB,1)...
            && tmp_vel_fB(locsz(mm)+eee) >= loomEscape(ff).boutAnalysis.BthreshOFF % using broad filter for offset
        eee = eee + 1;
    end
    
    tmp_swim_bouts(mm,2) = locsz(mm)+eee-1; %bout end frame
    
    %bout duration in ms
    tmp_swim_bouts(mm,3) = (tmp_swim_bouts(mm,2) - tmp_swim_bouts(mm,1))*frame_length_calc_ms;
    
    %orientation before bout /avg of 10 frames or 13.3 ms
    tmp_swim_bouts(mm,4) = nanmean(tmp_ori_deg((tmp_swim_bouts(mm,1)-11):(tmp_swim_bouts(mm,1)-1)));
    
    %orientation after bout /avg of 10 frames or 13.3 ms
    tmp_swim_bouts(mm,5) = nanmean(tmp_ori_deg((tmp_swim_bouts(mm,2)+1):(tmp_swim_bouts(mm,2)+11)));
        
    %TURN (ori after - ori before)
    tmp_swim_bouts(mm,6) = tmp_swim_bouts(mm,5)-tmp_swim_bouts(mm,4);
    
    %correct for change in orientation
    if tmp_swim_bouts(mm,6) > 180 % right turn 
        tmp_swim_bouts(mm,6) =  tmp_swim_bouts(mm,6) - 2*180;
        
    elseif tmp_swim_bouts(mm,6) < -180 % left turn
        tmp_swim_bouts(mm,6) =  2*180 + tmp_swim_bouts(mm,6);
    end
    
    % mean velocity during bout
    tmp_swim_bouts(mm,7) = nanmean(tmp_vel_fV(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2)));
    
    % peak velocity during bout
    tmp_swim_bouts(mm,8) = nanmax(tmp_vel_fV(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2)));
    
    % total distance covered during the bout
    tmp_swim_bouts(mm,9) = nansum(tmp_dist_fV(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2)));
    
    %head yaw during bouts
    yaws = tmp_delta_ori_filtered2(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2));
    tmp_swim_bouts(mm,10)= nansum(yaws(yaws>0)); %summation of total left yaws
    tmp_swim_bouts(mm,11)= nansum(yaws(yaws<0)); %summation of total right yaws
    
    %avg angular velocity
    tmp_swim_bouts(mm,12) = nansum(abs(yaws))/tmp_swim_bouts(mm,3); %deg/sec
    
    %tail half bends
    %tmp_swim_bouts(mm,14)= nansum(islocalmin(tmp_delta_ori_filtered2((tmp_swim_bouts(mm,1)):(tmp_swim_bouts(mm,2))),'MinSeparation',2)); %changed from 20
    %tmp_swim_bouts(mm,15)= nansum(islocalmax(tmp_delta_ori_filtered2((tmp_swim_bouts(mm,1)):(tmp_swim_bouts(mm,2))),'MinSeparation',2)); %changed from 20
    %tmp_swim_bouts(mm,16)= tmp_swim_bouts(mm,14)+tmp_swim_bouts(mm,15);
        
    %minmax(mm).min = tmp_swim_bouts(mm,1)+find(islocalmin(tmp_delta_ori_filtered2((tmp_swim_bouts(mm,1)):(tmp_swim_bouts(mm,2))))==1);
    %minmax(mm).max = tmp_swim_bouts(mm,1)+find(islocalmax(tmp_delta_ori_filtered2((tmp_swim_bouts(mm,1)):(tmp_swim_bouts(mm,2))))==1);
    
    %tail half bends - Calc'd by Roshan without using islocalmin/max
    locmin= findpeaks(-tmp_delta_ori_filtered2((tmp_swim_bouts(mm,1)):(tmp_swim_bouts(mm,2))),'MinPeakDistance',2); %changed from 20 to 2
    locmax= findpeaks(tmp_delta_ori_filtered2((tmp_swim_bouts(mm,1)):(tmp_swim_bouts(mm,2))),'MinPeakDistance',2); %changed from 20
    tmp_swim_bouts(mm,14)= nansum(locmin); 
    tmp_swim_bouts(mm,15)= nansum(locmax); 
    tmp_swim_bouts(mm,16)= tmp_swim_bouts(mm,14)+tmp_swim_bouts(mm,15);
        
    minmax(mm).min = tmp_swim_bouts(mm,1)+find(locmin==1);
    minmax(mm).max = tmp_swim_bouts(mm,1)+find(locmax==1);
    
    
    %delay - corrected 19th May 2018 - take Vp as the point to calculate
    %delay not bout start; to take care of cases where the animal was
    %already engaged in a bout when the tap was triggered
    tmp_swim_bouts(mm,13) = (locsz(mm)-idx_tap(mm))*frame_length_calc_ms; %if there is an error --> error in tap/escape ID

     %peak angular velocity
     tmp_swim_bouts(mm,17) = nansum(abs(tmp_ang_vel(locsz(mm,1):locsz(mm,2))))/((locsz(mm,2)-locsz(mm,1))*1.33)/1000; %in deg/sec already
     
     %tmp_swim_bouts(mm,18) = (nansum(islocalmax(tmp_delta_ori_filtered2((locsz(mm,2):locsz(mm,1))),'MinSeparation',2)))/((locsz(mm,2)-locsz(mm,1))*1.33)/1000;
     %TBF associated
     
    %rxn time - to-do after stim onset is extracted
    
    %tmp_swim_bouts(mm,19)=min(abs(idx_escape(mm)-idx_tap));
end

%% grab the escape bouts here
for cFast=1:length(escapeLocs)
    for  cBouts=1:size(tmp_swim_bouts,1)
    escapeTest=sum((tmp_swim_bouts(cBouts,1):tmp_swim_bouts(cBouts,2))==escapeLocs(cFast));
    if escapeTest==1
        tmp_escape_ID(cFast)=cBouts; % store the bout which corresponds to an escape

    end    
    end
end

for zzs=1:length(tmp_escape_ID)
    if tmp_escape_ID(zzs)>0 % 10-17-18 there is a problem where 0's are being placed into the beginning
                            % of tmp_escape_ID, which is coming from cBouts
    tmp_escape_bouts(zzs)=tmp_swim_bouts(tmp_escape_ID(zzs));
    end
end

%% - commented out since IBI is not imp for stim-induced escape analysis

% %% IBI section
% tmp_vel_IBI_fV = tmp_vel_fV;
% 
% % ID overlapping bouts and NaN them
% % very poorly written; had an indexing error and then just worked
% % around making it poorer and poorer, but it works!
% % changed to more efficient one on 22/01/18
% 
% for oo = 1:size(tmp_swim_bouts,1)
%       if (oo>=2 && tmp_swim_bouts(oo,1)<=tmp_swim_bouts(oo-1,2))...
%            || (oo<=size(tmp_swim_bouts,1)-1 && tmp_swim_bouts(oo,2)>=tmp_swim_bouts(oo+1,1))
%        
%        tmp_swim_bouts(oo, 3:end) = NaN; % convert all bout parameters corresp to this bout to NaN;
%        locs(oo) = NaN;
%        tmp_vel_IBI_fV(tmp_swim_bouts(oo,1):tmp_swim_bouts(oo,2))= NaN; % the velo vector to NaN as well - imp to ID NaN values in IBI cals  
%        
%       end
% end
% 
%             tmp_swim_bouts(isnan(tmp_swim_bouts(:,3)), : ) = [];
%             locs(isnan(locs))                              = [];
%                        
% % ID consecutive NaNs in velocity vector - imp for IBI identification
% tmp_vel_CC = bwconncomp(isnan(tmp_vel_IBI_fV));
% pixel_Idx = cellfun('prodofsize',tmp_vel_CC.PixelIdxList);
% tmp_vel_NaNs = zeros(size(tmp_vel_IBI_fV));
% 
% for zzz = 1:tmp_vel_CC.NumObjects
%     tmp_vel_NaNs(tmp_vel_CC.PixelIdxList{zzz})=pixel_Idx(zzz);
% end 
% 
% % Calculate the IBI
% tmp_swim_bouts(:,13) = NaN;
% 
% for ibi=2:size(tmp_swim_bouts,1)
%    tmp_swim_bouts(ibi,13)= (tmp_swim_bouts(ibi,1)-tmp_swim_bouts(ibi-1,2))*frame_length_calc_ms;
%    
%    if tmp_vel_NaNs(ibi)>=1.5*freeSwim.boutAnalysis.boutLength % >300 ms
%        tmp_swim_bouts(ibi,13)=NaN;
%    end
%     
% end

freeSwim.tmp_swim_bouts=tmp_swim_bouts;
freeSwim.tmp_escape_bouts=tmp_escape_bouts;

end

%% responsiveness index per fish

% We want to calculate: 
% A) the escape latency from stimulus onset to escape onset
% B) r x L (?) of stimulus at the initiation of escape for each bout
% C) Max velocity for each escape
% D) Max Angle of initial turn of each escape (head angle change)
% E) Direction of initial escape relative to stimulus
% F) Responsiveness of each fish (#escapes/#stim) - this might need to get
% adjusted in cases of multiple trials for the same fish

% Escape Latency: Find Escape initiation, calc EscInit - Loomstart

% Max Angle: For Escape Bout, Find first dOrient from 0 to when it changes
% sign, and get sum of these instantaneous dOrient's

% Direction: if Max angle is +: Right, if Max Angle is -: Left (?)

%% some plots
%fig1 = figure;
%hold on; 
%plot(tmp_vel_unfilt); 
%plot(tmp_vel_fV,'LineWidth', 3);
%plot(tmp_vel_fB,'LineWidth', 3);
%plot(tmp_vel_fF,'LineWidth', 3);
%plot(tmp_vel_fF,'LineWidth', 3);                

%plot(20*pks);
%plot(20*dxB-10);
%plot(dxV-10);
%plot(20*dyB-20);
%plot(dyV-20);
    % vline(idx_tap);
    % plot(10*tmp_delta_ori-1);
    % plot(10*tmp_delta_ori_filtered2-1.5, 'b-', 'LineWidth', 2);
%vline(tmp_swim_bouts(:,1),'m');
%vline(tmp_swim_bouts(:,2),'m');
%vline(idx_tap+1880,'k');
%vline(locsz,'c');
%hold off;

% fig2 = figure; histfit(tmp_swim_bouts(:,3),50);
% fig3 = figure; histogram(tmp_swim_bouts(:,7),50);
% fig4 = figure; histogram(tmp_swim_bouts(:,8),50);
% fig5 = figure; histogram(tmp_swim_bouts(:,12),50);
% fig6 = figure; histogram(tmp_swim_bouts(:,13),50);
% fig7 = figure; histfit(tmp_swim_bouts(:,6),50);
% fig8 = figure; histfit(tmp_swim_bouts(:,12),50);

% dcm1 = datacursormode(fig1);
% set(dcm1, 'UpdateFcn', @Data_Cursor_precision, 'Enable', 'on');


%savePath = strcat(PathName,'Danio_TAP_16June2018.mat');
%save(savePath,'freeSwim');

%% test on 12th June - for ang velo check

% fig2 = figure;
% yyaxis left
% plot(tmp_vel_fB);
% yyaxis right
% plot(tmp_delta_ori_filtered2);
% 
% hold on;
% 
% vline(tmp_swim_bouts(:,1));
% vline(tmp_swim_bouts(:,2));
% hold off;
% 
% dcm1 = datacursormode(fig2);
% set(dcm1, 'UpdateFcn', @Data_Cursor_precision, 'Enable', 'on');

%% sanity check on 09/09/2018
% plot(tmp_vel_fB);
% hold on;
% vline(idx_tap);
% vline(idx_escape);
% hold off;
% save(freeSwim,'test.mat');

%% added for loom
%function h = circle(x,y,r)
%hold on
%th = 0:pi/50:2*pi;
%xunit = r * cos(th) + x;
%yunit = r * sin(th) + y;
%h = plot(xunit, yunit);
%hold off
%end

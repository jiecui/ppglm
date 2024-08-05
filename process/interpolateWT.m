% Authors: Maxime Baud, Timothée Proix
% License: GPL-3.0-only

function [invWT] = interpolateWT(ts, mu)
% interpolateWT linearly interpolates the gap in the time series 
% and obtain the inverse wavelet transform with gaps filled when short enough
% compare to the period of the wavelet

% INPUTS:
% - ts: timeseries sampled hourly (fs = 24h per day), with NaNs replacing gaps. 
% - per: periods to evaluate
% - mu: bands of interest

% OUTPUTS
% - invWT: filled inverse wavelet transform 



%% Parameters
ts=ts(:)';
t =length(ts);
fs = 24;                                                                % 24h/day, so period units will therefore be days
per = [.1:.05:1.3 , 1.4:.1:2, 2.2:.2:4 , 4.5:.5:10 , 11:1:45];
fourier_factor = 1.0330;
scales = per/fourier_factor;                                            % period = fourier_factor*scale;
dt = 1/fs;                                                              % fs = 24 h/d Hz -> sampling period dt = 1/24 day;

%% FIND GAPS 
ON = find(diff(isnan(ts))==1)+1;                                        % Gap ON
OFF = find(diff(isnan(ts))==-1);                                        % Gap OFF
if isnan(ts(1)), OFF(1) = []; end                                      % No need to take care of first gap                            
if isnan(ts(end)), ON(end) = []; end                                   % No need to take care of final gap
gapsize=(OFF-ON) + 1;      % -1 changed to +1 here                                              % Length of gaps
maxgap = floor(max(per)*24*.2);                                        % 20% of longest period is OK to cross
gaps2fill = find(gapsize<maxgap);
g_ON = ON(gaps2fill);
g_OFF = OFF(gaps2fill);

%% FIND BLOCKS of DATA
gaps = find(gapsize >= maxgap);
ts_ON = find(~isnan(ts),1,'first');
ts_OFF = find(~isnan(ts),1,'last');
ts_ON = unique([ts_ON ; OFF(gaps)'+1]);                                  % Gap ON is ts OFF
ts_OFF = unique([ts_OFF ; ON(gaps)'-1]);                                 % Unique also makes sure to sort
blocksize = ts_OFF-ts_ON;                                               % Block duration
minblock = floor(max(per)*24);                                          % Take blocks of minimum 3 continuous months
ts_ON = ts_ON(blocksize > minblock);
ts_OFF = ts_OFF(blocksize > minblock);

%% LINEAR INTERPOLATION                                                  % First path as to be able to calculate WT
nfill = numel(g_ON);
if ~isempty(gaps2fill)
    for g=1:nfill
        fillsize = (g_OFF(g)-g_ON(g))+1;
        bef=ts(max(1,g_ON(g)-ceil(gapsize(g))):g_ON(g)-1);            % BEFORE GAP
        aft=ts(g_OFF(g)+1:min(t,g_OFF(g)+ceil(gapsize(g))));                 % AFTER GAP
        rgap=randn(1,fillsize)/3;                                       % FILL IN random data for gap
        v=2*std([bef aft]);                                             % VARIANCE base on pre- and post-gap
        trnd=linspace(mean(bef),mean(aft),fillsize);                    % TREND in data from BEFORE to AFTER
        gfill=rgap*v+trnd;                                       % FILL material
        ts(g_ON(g):g_OFF(g))=gfill;
    end
end

%% WAVELET TRANSFORM
invWT = nan(numel(mu),t);
nblocks = numel(ts_ON);
for g=1:nblocks
    for imu=1:length(mu)
        block = ts(ts_ON(g):ts_OFF(g));
        WT = cwtft({block,dt},'scales',scales,'wavelet','morl');            % Morlet wavelet
        bdw = 1/3*mu(imu);
        bd = per > mu(imu) - bdw & per < mu(imu) + bdw;
        WT.cfs(~bd,:)=0;
        invWT(imu, ts_ON(g):ts_OFF(g)) = icwtft(WT)-mean(icwtft(WT));         % inverse Morlet wavelet
    end
end

%% REMOVE CONE OF INFLUENCE
for imu = 1:numel(mu)
    fullmu = mu(imu)*24;                                                  % full period
    partmu = fullmu*0.2;                                                  % fraction (20%) of each period that can be interpolated
    coi=ceil(fullmu*fourier_factor);
    
    for g=1:nblocks
        invWT(imu,ts_ON(g):ts_ON(g)+coi) = NaN;                                % beginning
        invWT(imu,ts_OFF(g)-coi:ts_OFF(g)) = NaN;                             % end
    end
    
    for g=1:numel(gaps2fill)
        if gapsize(g) > partmu                                           % if not acceptable to interpolate
            coi_idx = max(1,g_ON(g)-coi):min(t,g_OFF(g)+coi);           % bound
            invWT(imu,coi_idx) = NaN;                                         % REMOVE CONE OF INFLUENCE from WT
        end
    end
end
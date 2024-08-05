% Author: Timothée Proix
% License: GPL-3.0-only

clear all; close all;
format LONGG

showPlot = 0;
saveData = 1;
preProcessed = 1; % set to 1 to use the data provided on Zenodo
filterType = 'acausal4'; 
maxMultiRhythmPeriod = 50; % in days. Maximal period of the multidien rhythms considered
gapInterpolationLimit = 10*24; % in hours. Bigger gaps are not interpolated
timeDayStart = 12;
timeDayEnd = 11;
lengthTrainParam = [480*24 0.60]; % [hours, percentage] minimal duration of train data, either in absolute value or percentage

% Set your own path roots
rootData = 'YourPath'
rootProcess = 'YourPath'
datasetc = 'UCSF_patients'
pathData = [rootData '/forecasting_over_days/' datasetc '/patients/'];
pathProcess = [rootProcess '/forecasting_over_days/' datasetc '/'];

%% extract patients data from D structure
% Full list
patientList = readmatrix([pathData 'patients.csv'], 'OutputType', 'char');
stats = cell2table(cell(length(patientList), 8), 'VariableNames', {'name', 'peakMulti', 'gapLengths', ...
                   'lengthtDisc', 'lengthBlocks', 'nbSzBlocks', 'plvs', 'sz_days'});

for idxPatient = 1:length(patientList)

    patient = patientList{idxPatient};
    disp(patient);

    if ~preProcessed
        stats{idxPatient, 'name'} = {patient};

        load([pathData 'D_' patient '.mat'])
        % to deal with different names in UCSF dataset for validated and not validated seizures
        if ~isfield(D.Selection.Hourly, 'C_sz')
            D.Selection.Hourly.C_sz = D.Selection.Hourly.E_sz;
        end
        if isfield(D.Selection, 'block2')
            disp('Two blocks for this patient');
        end
        if isfield(D.Selection, 'block3')
            error('the code does not deal with more than two blocks in the data');
        end

        %% Get times series
        if size(strfind(patient, '_block2'), 1)==0
            [~, idxBlockHourlyt1] = min(abs(D.Selection.Hourly.Time-D.Selection.block1(1)));
            [~, idxBlockHourlyt2] = min(abs(D.Selection.Hourly.Time-D.Selection.block1(2)));
        else
            [~, idxBlockHourlyt1] = min(abs(D.Selection.Hourly.Time-D.Selection.block2(1)));
            [~, idxBlockHourlyt2] = min(abs(D.Selection.Hourly.Time-D.Selection.block2(2)));
        end

        % IEA time series
        IEAAll = transpose(D.Selection.Hourly.IEA1((idxBlockHourlyt1:idxBlockHourlyt2)));
        % if there is a second IEA time series
        boolTwoIEAs = isfield(D.Selection.Hourly, 'IEA2') && ...
                      sum(~isnan(D.Selection.Hourly.IEA2(idxBlockHourlyt1:idxBlockHourlyt2)))>0 && ...
                      sum(~isnan(D.Selection.Hourly.IEA1(idxBlockHourlyt1:idxBlockHourlyt2)))>0; 
        if boolTwoIEAs
            IEAAll = [IEAAll transpose(D.Selection.Hourly.IEA2(idxBlockHourlyt1:idxBlockHourlyt2))];
            assert(all(isnan(D.Selection.Hourly.IEA1(idxBlockHourlyt1:idxBlockHourlyt2))==isnan(D.Selection.Hourly.IEA2(idxBlockHourlyt1:idxBlockHourlyt2))), ...
                   'NaNs are not the same in the two IEA time series');
        end

        % Seizure times series
        SzAll = transpose(D.Selection.Hourly.C_sz(idxBlockHourlyt1:idxBlockHourlyt2));

        % Time times series
        TimeAll = transpose(D.Selection.Hourly.Time(idxBlockHourlyt1:idxBlockHourlyt2));

        % select multidien peaks between 2 and 50 days
        Pks = D.Selection.pks_IEA1;
        peaksMultiAll = {Pks(Pks>=2 & Pks<maxMultiRhythmPeriod)};
        if boolTwoIEAs
            Pks = D.Selection.pks_IEA2;
            peaksMultiAll(2) = {Pks(Pks>=2 & Pks<maxMultiRhythmPeriod)};
        end
        nbPeaks = sum(cellfun(@(x) size(x, 2), peaksMultiAll));
        if strcmp(filterType, 'lumped')
            nbPeaks = 1;
        end
        stats{idxPatient, 'peakMulti'} = {cell2mat(peaksMultiAll)};
          
        %% Divide data in blocks around large gaps
        % cut initial or trailing NaNs
        idxNotNan = find(~(isnan(IEAAll(:, 1)) | isnan(SzAll)));
        IEAAll =  IEAAll(idxNotNan(1):idxNotNan(end), :);  
        TimeAll = TimeAll(idxNotNan(1):idxNotNan(end));
        SzAll = SzAll(idxNotNan(1):idxNotNan(end));
        assert(~isnan(IEAAll(1, 1)), 'initial NaNs');
        assert(~isnan(IEAAll(end, 1)), 'trailing NaNs');

        % find large gaps
        idxGap = transpose(find(isnan(IEAAll(:, 1)) | isnan(SzAll))); %NaNs are the same in both IEA times series
        if length(idxGap)==0
            startBlock = [1];
            endBlock = [size(IEAAll, 1)];
            stats{idxPatient, 'gapLengths'} = {[]};
        else
            endPosGap = [idxGap(diff(idxGap)>1) idxGap(end)];
            startPosGap = [idxGap(1) idxGap(find(diff(idxGap)>1)+1)];
            gapLength = endPosGap-startPosGap;
            stats{idxPatient, 'gapLengths'} = {gapLength};
            endBlock = [startPosGap(find(gapLength>gapInterpolationLimit))-1 size(IEAAll, 1)];
            startBlock = [1 endPosGap(find(gapLength>gapInterpolationLimit))+1];
        end

        tDiscAll = [];
        SzDiscFilledAll = [];
        IEADiscFilledAll = [];
        for iBlock = 1:length(startBlock)
            disp([datetime(TimeAll(startBlock(iBlock)), 'ConvertFrom', 'datenum') datetime(TimeAll(endBlock(iBlock)), 'ConvertFrom', 'datenum')]);

            iMarker = 0;
            for iIEA = 1:size(IEAAll, 2)
                IEA = IEAAll(startBlock(iBlock):endBlock(iBlock), iIEA);
                Sz = SzAll(startBlock(iBlock):endBlock(iBlock));
                Time = TimeAll(startBlock(iBlock):endBlock(iBlock));
                peaksMulti = peaksMultiAll{iIEA};

                assert(length(IEA)>50*24, 'block too small - check your data')

                if showPlot==1
                    figure();
                    sgtitle([patient ' raw data'])
                    subplot(211)
                    plot(24*(Time-Time(1)), IEA);
                    xlim([0 length(Time)]);
                    subplot(212);
                    plot(24*(Time-Time(1)), Sz);
                    xlim([0 length(Time)]);
                end
                
                %% Z-scoring the IEA data: already done in the data
                % tvisit = D.Settings.Visit;
                % indVisits = find(ismember(Time, tvisit));
                % tLast = 1;
                % zIEA = [];
                % for it=1:size(indVisits, 2)
                %     currDat = IEA(1, tLast:indVisits(it));
                %     zIEA = [zIEA (currDat - nanmean(currDat))/nanstd(currDat)];
                %     tLast = indVisits(it)+1;
                % end
                % currDat = IEA(1, tLast:end);
                % zIEA = [zIEA (currDat - nanmean(currDat))/nanstd(currDat)];
                zIEA = IEA;

                if showPlot==2
                    figure();
                    sgtitle([patient ' normalized data']) %2018b
                    subplot(311)
                    plot(24*(Time-Time(1)), IEA);
                    xlim([0 length(Time)]);
                    subplot(312)
                    plot(24*(Time-Time(1)), zIEA);
                    xlim([0 length(Time)]);        
                    subplot(313);
                    plot(24*(Time-Time(1)), Sz);
                    xlim([0 length(Time)]);
                end
                
                %% We start and end all blocks at the time timeDayStart and timeDayEnd
                hStart = hour(dateshift(datetime(Time(1),'ConvertFrom','datenum'), 'start', 'hour', 'nearest'));
                if hStart<=timeDayStart
                    tStart = timeDayStart - hStart + 1;
                elseif hStart>timeDayStart
                    tStart = 24 + timeDayStart - hStart + 1;
                end
                hEnd = hour(dateshift(datetime(Time(end),'ConvertFrom','datenum'), 'start', 'hour', 'nearest'));
                if hEnd<timeDayEnd
                    tEnd = length(Time) - (timeDayEnd + hEnd) - 2;
                elseif hEnd>=timeDayEnd
                    tEnd = length(Time) - (hEnd - timeDayEnd);
                end
                % check
                assert(hour(dateshift(datetime(Time(tStart),'ConvertFrom','datenum'), 'start', 'hour', 'nearest'))==timeDayStart);
                assert(hour(dateshift(datetime(Time(tEnd),'ConvertFrom','datenum'), 'start', 'hour', 'nearest'))==timeDayEnd);

                if mod(tEnd-tStart + 1, 24)~=0
                    disp(['Warning, some time points are missing in the times series, ' ...
                         'cutting so that the times series is a multiple of 24, '...
                         'the times series will not finish at the time asked for']);
                    tEnd = tEnd - mod(tEnd-tStart + 1, 24);
                end

                IEADisc = zIEA(tStart:tEnd);
                SzDisc = Sz(tStart:tEnd);           
                tDisc = Time(tStart:tEnd);
                if iBlock==1
                    startAll = tDisc(1); % initial start time
                end
                tvisit = D.Settings.Visit;
                indVisitsDisc = find(ismember(tDisc, tvisit))';
                tVisitUsedDisc = tDisc(ismember(tDisc, tvisit))';
                stats{idxPatient, 'lengthtDisc'} = {length(tDisc)};
                
                if showPlot==3
                    figure();
                    sgtitle([patient ' IEA discarded']) %2018b
                    subplot(211)
                    plot(24*(tDisc-tDisc(1)), IEADisc);
                    xlim([0 length(tDisc)]);
                    subplot(212);
                    plot(24*(tDisc-tDisc(1)), SzDisc);
                    xlim([0 length(tDisc)]);
                end

                %% Filling small IEA gaps
                [IEADiscFilled, SzDiscFilled] = fillGaps(tDisc, IEADisc, SzDisc, peaksMulti, indVisitsDisc, tVisitUsedDisc, patient, showPlot);

                % stats     
                if iIEA==1
                    if iBlock==1
                        stats{idxPatient, 'lengthBlocks'} = {length(IEADiscFilled)};
                        stats{idxPatient, 'nbSzBlocks'} = {sum(SzDiscFilled>0)};
                    else
                        stats{idxPatient, 'lengthBlocks'} = {[stats{idxPatient, 'lengthBlocks'} length(IEADiscFilled)]};
                        stats{idxPatient, 'nbSzBlocks'} = {[stats{idxPatient, 'nbSzBlocks'} sum(SzDiscFilled>0)]};
                    end
                end

                if showPlot==6
                    figure();
                    sgtitle([patient ' filled seizures'])
                    subplot(211)
                    plot(24*(tDisc-tDisc(1)), IEADiscFilled); hold on;
                    plot(24*(tDisc-tDisc(1)), IEADisc);
                    subplot(212);
                    plot(24*(tDisc-tDisc(1)), SzDiscFilled); hold on;
                    plot(24*(tDisc-tDisc(1)), SzDisc);
                    
                end
                
                if iIEA==1
                    IEADiscFilledBoth = IEADiscFilled;
                else
                    IEADiscFilledBoth = [IEADiscFilledBoth IEADiscFilled];
                end
            end
            tDiscAll = [tDiscAll; tDisc];
            SzDiscFilledAll = [SzDiscFilledAll; SzDiscFilled];
            IEADiscFilledAll = [IEADiscFilledAll; IEADiscFilledBoth];
        end
    else
        datPreprocessed = readtable([pathData 'preprocessed_' patient '.csv']);
        tDiscAll = datPreprocessed.Time;
        SzDiscFilledAll = datPreprocessed.Seizures;
        if sum(strcmp('IEA_2', datPreprocessed.Properties.VariableNames))
            IEADiscFilledAll = [datPreprocessed.IEA_1, datPreprocessed.IEA_2];
        else
            IEADiscFilledAll = datPreprocessed.IEA;
        end
        load([pathData '/peaksMultiAll_' patient])
    end

    %% Now reestimate the rhythms with the filled data for the covariates:
    lengthTrainData = min([lengthTrainParam(1), floor(lengthTrainParam(2) * length(tDiscAll))]);

    if strcmp(filterType, 'acausal4')   

        [ampC, phaseC, filtSigC] = acausal4(IEADiscFilledAll, repmat({[1]}, 1, size(IEADiscFilledAll, 2)));
        [ampM, phaseM, filtSigM] = acausal4(IEADiscFilledAll, peaksMultiAll);
    
    else

        error('asked filter does not exist');
    end


    if showPlot==7
        figure()
        plot(IEADiscFilledAll); hold on;
        plot(filtSigC); 
        yyaxis right
        sgtitle([patient ' circadian filtered IEA'])
        figure()
        plot(IEADiscFilledAll); hold on;
        plot(filtSigM); 
        yyaxis right
        sgtitle([patient ' multidien filtered IEA'])
    end

    if showPlot==8
        figure()
        plot(IEADiscFilledAll); hold on;
        plot(cos(phaseC));
        plot(ampC); 
        yyaxis right
        sgtitle([patient ' circadian abs and cos IEA'])
        figure()
        plot(IEADiscFilledAll); hold on;
        plot(cos(phaseM));
        plot(ampM);
        yyaxis right
        sgtitle([patient ' multidien abs and cos IEA'])
    end

    procDatHour = table();
    procDatDay = table();
    procDatHour.Time = tDiscAll;
    procDatDay.Time = mean24(tDiscAll);
    procDatHour.Seizures = SzDiscFilledAll;
    procDatDay.Seizures = mean24(SzDiscFilledAll)*24;
    procDatHour.IEA = IEADiscFilledAll;
    procDatDay.IEA = median24(IEADiscFilledAll);
    for iIEA = 1:size(filtSigC, 2)
        procDatHour.(['CircadianAbs' num2str(iIEA)]) = ampC(:, iIEA);
        procDatDay.(['CircadianAbs' num2str(iIEA)]) = median24(ampC(:, iIEA));
        procDatHour.(['CircadianAngle' num2str(iIEA)]) = phaseC(:, iIEA);
        procDatDay.(['CircadianAngle' num2str(iIEA)]) = median24(phaseC(:, iIEA));
        procDatHour.(['CircadianCos' num2str(iIEA)]) = cos(phaseC(:, iIEA));
        procDatDay.(['CircadianCos' num2str(iIEA)]) = cos(median24(phaseC(:, iIEA)));
        procDatHour.(['CircadianSin' num2str(iIEA)]) = sin(phaseC(:, iIEA));
        procDatDay.(['CircadianSin' num2str(iIEA)]) = sin(median24(phaseC(:, iIEA)));
        procDatHour.(['CircadianFiltSig' num2str(iIEA)]) = filtSigC(:, iIEA);
        procDatDay.(['CircadianFiltSig' num2str(iIEA)]) = median24(filtSigC(:, iIEA));
    end
    for iRhythm = 1:size(filtSigM, 2)
        procDatHour.(['MultidienAbs' num2str(iRhythm)]) = ampM(:, iRhythm);
        procDatDay.(['MultidienAbs' num2str(iRhythm)]) = median24(ampM(:, iRhythm));
        procDatHour.(['MultidienAngle' num2str(iRhythm)]) = phaseM(:, iRhythm);
        procDatDay.(['MultidienAngle' num2str(iRhythm)]) = median24(phaseM(:, iRhythm));
        procDatHour.(['MultidienCos' num2str(iRhythm)]) = cos(phaseM(:, iRhythm));
        procDatDay.(['MultidienCos' num2str(iRhythm)]) = cos(median24(phaseM(:, iRhythm)));
        procDatHour.(['MultidienSin' num2str(iRhythm)]) = sin(phaseM(:, iRhythm));
        procDatDay.(['MultidienSin' num2str(iRhythm)]) = sin(median24(phaseM(:, iRhythm)));
        procDatHour.(['MultidienFiltSig' num2str(iRhythm)]) = filtSigM(:, iRhythm);
        procDatDay.(['MultidienFiltSig' num2str(iRhythm)]) = median24(filtSigM(:, iRhythm));
    end

        
    %% Time series for seizure circadian distribution
    % only use the train data to estimate the distribution
    meanSzHourOfDay = zeros(24, 1);
    stdSzHourOfDay = zeros(24, 1);
    for hourOfDay = 0:23
        idxHour = find(round(mod(tDiscAll(1:lengthTrainData)*24, 24))==hourOfDay);
        meanSzHourOfDay(hourOfDay+1) = nanmean(SzDiscFilledAll(idxHour));
        stdSzHourOfDay(hourOfDay+1) = nanstd(SzDiscFilledAll(idxHour));
    end
    
    % Build the Circadian time series Sz
    tsMeanSzCirc = nan(length(SzDiscFilledAll), 1);
    tsStdSzCirc = nan(length(SzDiscFilledAll), 1);
    for iSzDisc = 1:length(SzDiscFilledAll)
        tsMeanSzCirc(iSzDisc) = meanSzHourOfDay(round(mod(tDiscAll(iSzDisc)*24, 24))+1);
        tsStdSzCirc(iSzDisc) = stdSzHourOfDay(round(mod(tDiscAll(iSzDisc)*24, 24))+1);
    end

    if find(isnan(tsMeanSzCirc))
        disp("*******************NaN values left in circ Sz!!! ****************")
    end     

    if showPlot==9
        h = figure();
        sgtitle([patient ' circadian seizures'])
        subplot(411)
        plot(meanSzHourOfDay);
        subplot(412)
        plot(stdSzHourOfDay);
        subplot(413)
        plot(tsMeanSzCirc); hold on;
        plot(SzDiscFilledAll/max(SzDiscFilledAll(:)));
        subplot(414)
        plot(tsStdSzCirc); hold on;
        plot(SzDiscFilledAll/max(SzDiscFilledAll(:)));
    end

    % add circadian rhythms, no circadian rhythms for daily
    procDatHour.meanDailySeizures = tsMeanSzCirc;
    procDatHour.stdDailySeizures = tsStdSzCirc;

    %% Weekly distribution for Sz
    % only use the train data to estimate the distribution
    meanSzDayOfWeek = zeros(7, 1);
    stdSzDayOfWeek = zeros(7, 1);
    for dayOfWeek = 0:6
        idxDay = find(floor(mod(tDiscAll(1:lengthTrainData), 7))==dayOfWeek);
        meanSzDayOfWeek(dayOfWeek+1, 1) = nanmean(SzDiscFilledAll(idxDay));
        stdSzDayOfWeek(dayOfWeek+1, 1) = nanstd(SzDiscFilledAll(idxDay));
    end

    % Build the weekly time series Sz
    tsMeanSzWeek = nan(length(SzDiscFilledAll), 1);
    tsStdSzWeek = nan(length(SzDiscFilledAll), 1);
    for iSzDisc = 1:length(SzDiscFilledAll)
        tsMeanSzWeek(iSzDisc) = meanSzDayOfWeek(floor(mod(tDiscAll(iSzDisc), 7))+1);
        tsStdSzWeek(iSzDisc) = stdSzDayOfWeek(floor(mod(tDiscAll(iSzDisc), 7))+1);
    end


    if find(isnan(tsMeanSzWeek))
        disp("*******************NaN values left in week Sz!!! ****************")
    end     

    if showPlot==10
        h = figure();
        sgtitle([patient ' Weekly seizures'])
        subplot(411)
        plot(meanSzDayOfWeek);
        subplot(412)
        plot(stdSzDayOfWeek);
        subplot(413)
        plot(tsMeanSzWeek); hold on;
        plot(SzDiscFilledAll/max(SzDiscFilledAll(:)));
        subplot(414)
        plot(tsStdSzWeek); hold on;
        plot(SzDiscFilledAll/max(SzDiscFilledAll(:)));
    end

    % add weekly rhythms, only for daily
    procDatDay.meanWeeklySeizures = mean24(tsMeanSzWeek);
    procDatDay.stdWeeklySeizures = mean24(tsStdSzWeek);

    % plot concatenated data
    if showPlot==11
        figure();
        subplot(211);
        plot(procDatHour.Time, procDatHour.MultidienCos1); hold on;
        for isz = find(procDatHour.Seizures)
            scatter(procDatHour.Time(isz), procDatHour.MultidienCos1(isz), 'r', 'filled');
        end
        xlim([min(procDatHour.Time) max(procDatHour.Time)])
        subplot(212);
        plot(procDatHour.Time, procDatHour.IEA); hold on;
        scatter(procDatHour.Time(isz), procDatHour.IEA(isz), 'r', 'filled');
        sgtitle([patient ' concatenated'])
        xlim([min(procDatHour.Time) max(procDatHour.Time)])
    end

    if ~preProcessed
        % Finally collect some stats
        if boolTwoIEAs
            stats{idxPatient, 'plvs'} = {[max(D.Selection.PLV_IEA1), max(D.Selection.PLV_IEA2)]};
        else
            stats{idxPatient, 'plvs'} = {max(D.Selection.PLV_IEA1)};
        end
        stats{idxPatient, 'sz_days'} = {D.Selection.Daily.sz_days};
    end

    if saveData
        if ~exist([pathProcess 'features'], 'dir')
            mkdir([pathProcess 'features']);
        end
        if ~exist([pathProcess 'stats'], 'dir')
            mkdir([pathProcess 'stats']);
        end
        writetable(procDatHour, [pathProcess 'features/features_hour_' patient '_' filterType '.csv']);
        writetable(procDatDay, [pathProcess 'features/features_day_' patient '_' filterType '.csv']);
        if ~preProcessed
            writetable(stats, [pathProcess 'stats/stats_processing.csv']);
            save([pathProcess 'stats/stats_processing.mat'], 'stats');
        end
    end 
end
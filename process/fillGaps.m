% Author: Timothée Proix
% License: GPL-3.0-only

function [IEADiscFilled, SzDiscFilled] = fillGaps(tDisc, IEADisc, SzDisc, peaksMulti, indVisitsDisc, tVisitUsedDisc, patient, showPlot)

    % First interpolation (linear) to obtain the long multidian rhythms
    invWT = interpolateWT(IEADisc, peaksMulti);
    multiDien = nansum(invWT, 1);

    if showPlot==4
        figure();
        sgtitle([patient ' interpolated inverse wavelet transform'])
        subplot(211)
        plot(24*(tDisc-tDisc(1)), IEADisc); hold on;
        plot(24*(tDisc-tDisc(1)), invWT);
        xlim([0 length(tDisc)]);
        subplot(212);
        plot(24*(tDisc-tDisc(1)), SzDisc);
        xlim([0 length(tDisc)]);
    end

    % Circadian distribution for IEA
    indVisitsAll = [1 indVisitsDisc length(IEADisc)];
    meanIEAHourOfDay = zeros(24, length(indVisitsAll)-1);
    stdIEAHourOfDay = zeros(24, length(indVisitsAll)-1);
    % init
    initHourOfDay = mod(tDisc(1)*24,24);
    meanIEAHourOfDay(initHourOfDay+1, 1) = IEADisc(1);
    stdIEAHourOfDay(initHourOfDay+1, 1) = 0;
    for iVisit = 2:length(indVisitsAll)
        for hourOfDay = 0:23
            idxHour = find(mod(tDisc*24, 24)==hourOfDay);
            idxHourVisit = idxHour(idxHour>indVisitsAll(iVisit-1) & idxHour<=indVisitsAll(iVisit));
            if ~isempty(idxHourVisit)
                if ~isnan(nanmean(IEADisc(idxHourVisit)))
                    meanIEAHourOfDay(hourOfDay+1, iVisit-1) = nanmean(IEADisc(idxHourVisit));
                    stdIEAHourOfDay(hourOfDay+1, iVisit-1) = nanstd(IEADisc(idxHourVisit));
                end
            end
        end
    end

    % take care of first point if no other points in that visit range
    initHourOfDay = mod(tDisc(1)*24,24);
    if isnan(meanIEAHourOfDay(initHourOfDay+1, 1))
        meanIEAHourOfDay(initHourOfDay+1, 1) = IEADisc(1);
        stdIEAHourOfDay(initHourOfDay+1, 1) = 0;
    end
    
    % Build the Circadian time series IEA
    tsMeanIEACirc = nan(length(IEADisc), 1);
    tsStdIEACirc = nan(length(IEADisc), 1);
    for iIEADisc = 1:length(IEADisc)-1
        tsMeanIEACirc(iIEADisc) = meanIEAHourOfDay(mod(tDisc(iIEADisc)*24, 24)+1, find(iIEADisc<=indVisitsAll(2:end), 1, 'first'));
        tsStdIEACirc(iIEADisc) = stdIEAHourOfDay(mod(tDisc(iIEADisc)*24, 24)+1, find(iIEADisc<=indVisitsAll(2:end), 1, 'first'));
    end
    nanmean(IEADisc(idxHour(idxHour>indVisitsAll(iVisit-1) & idxHour<=indVisitsAll(iVisit))));
    tsMeanIEACirc(end) =  meanIEAHourOfDay(mod(tDisc(end)*24, 24)+1, end);
    tsStdIEACirc(end) =  stdIEAHourOfDay(mod(tDisc(end)*24, 24)+1, end);

    if find(isnan(tsMeanIEACirc))
        disp("*******************NaN values left in circ IEA!!! ****************")
        find(isnan(IEADiscFilled))
    end

    if showPlot==5
        figure();
        sgtitle([patient ' circadian distribution IEA'])
        subplot(411)
        plot(meanIEAHourOfDay);
        subplot(412)
        plot(stdIEAHourOfDay);
        subplot(413)
        plot(tsMeanIEACirc); hold on;
        plot(IEADisc/max(IEADisc(:)));
        subplot(414)
        plot(tsStdIEACirc); hold on;
        plot(IEADisc/max(IEADisc(:)));    
    end

    % find and fill small IEA gaps
    IEADiscFilled = IEADisc;
    idxGapCurrIEA = transpose(find(isnan(IEADisc)));
    for iGap = idxGapCurrIEA
        IEADiscFilled(iGap) = multiDien(iGap) + tsMeanIEACirc(iGap) + tsStdIEACirc(iGap) .* randn(1, 1);

        firstIndexMean = indVisitsDisc(find(tVisitUsedDisc<tDisc(iGap), 1, 'last'));
        if isempty(firstIndexMean)
            firstIndexMean = 1;
        end
        lastIndexMean = indVisitsDisc(find(tVisitUsedDisc>tDisc(iGap), 1, 'first'));
        if isempty(lastIndexMean)
            lastIndexMean = size(tDisc, 2);
        end
        iMin = nanmin(IEADisc(firstIndexMean:lastIndexMean));
        if IEADiscFilled(iGap)<iMin
            IEADiscFilled(iGap) = iMin;
        end
        iMax = nanmax(IEADisc(firstIndexMean:lastIndexMean));
        if IEADiscFilled(iGap)>iMax
            IEADiscFilled(iGap) = iMax;
        end
    end
    
    % Filling small Sz gaps
    idxGapCurrSz = find(isnan(SzDisc));
    SzDiscFilled = SzDisc;
    SzDiscFilled(idxGapCurrSz) = 0;

    % check
    if find(isnan(IEADiscFilled))
        disp("*******************NaN values left in IEA!!! ****************")
        find(isnan(IEADiscFilled))
    end
    if find(isnan(SzDiscFilled))
        disp("*******************NaN values left in Sz!!! ****************")
    end   
end
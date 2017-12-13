clear all

%Subject number in first column, block number in 2nd column, go RT in 3rd
%column, SSD in 4th column, sig-resp RT in 5th column. Replace non-numbers in columns 3, 4, and 5 with a number
%below the SSD min (I usually use -500). This kludge seems to work
[SubjectSeq Block GoRTSeq SSDSeq SigRespRT] = textread('B&L20112up1downSigResp&PrecedingRT.txt', '%f%f%f%f%f');
%Subject Numbers in first column, mean go RT excluding omissions in second
%column
[SubjectNum MeanRT] = textread('B&L20112up1downMeanRT.txt', '%f%f');

SigRespCountCutoff = 5; %Threshold for the number of signal-respond trials at a specific SSD for that subject to be computed
MinimumSubjectsForAverage = 5; %Threshold for the number of subjects that pass the SigRespCountCutoff at that SSD for that SSD to be included in the group average

%Could hardcode SSDMin or SSDMax if you only want to evaluate a subset of
%the SSD distribution
SSDMin = min(SSDSeq(SSDSeq > -500));
SSDMax = max(SSDSeq);
NumberOfSSDs = size(unique(SSDSeq), 1)-1; %last -1 is to account for the -500 kludge mentioned above 

for a=1:(size(SubjectNum, 1))
    SubjectNumber = SubjectNum(a);
    d = 1; 
    for c=SSDMin:50:SSDMax
        e = 1; 
        for b=1:size(SubjectSeq)-1
            if(GoRTSeq(b) > 0 && SigRespRT(b+1) > 0 && SSDSeq(b+1) == c && SubjectSeq(b) == SubjectNumber)
                NoStopOutput(e, d, a) = GoRTSeq(b);
                SigRespOutput(e, d, a) = SigRespRT(b+1);
                e = e + 1;
            end
        end
        d = d + 1; 
    end
end

NoStopOutput(NoStopOutput==0) = NaN;
SigRespOutput(SigRespOutput==0) = NaN;

%Finds the mean sig-resp RT and mean preceding no-stop RT for each subject at each SSD in which they have >= SigRespCountCutoff number of signal-respond trials
for f=1:(size(SubjectNum, 1))
    SubjectNumber2 = SubjectNum(f);
    h = 1; 
    for g=SSDMin:50:SSDMax
        if(isnan(NoStopOutput(SigRespCountCutoff, h, f))) == 0;
            SigRespScatterplotY(h, f) = (nanmean(SigRespOutput(:, h, f)));
            SigRespScatterplotX(h, f) = g;
            NoStopScatterplotY(h, f) = (nanmean(NoStopOutput(:, h, f)));
            NoStopScatterplotX(h, f) = g;
        else
            SigRespScatterplotY(h, f) = NaN;
            SigRespScatterplotX(h, f) = NaN;
            NoStopScatterplotY(h, f) = NaN;
            NoStopScatterplotX(h, f) = NaN;            
        end
    h = h +1; 
    end
end

%Excludes SSDs from the grand mean averages for which there are fewer than
%MinimumSubjectsForAverage numbers of subjects at that SSD that had >=
%SigrespCountCutoff number of signal-respond trials. 
for l=1:NumberOfSSDs
    FullSSDList(l, 1) = SSDMin+(l-1)*50;
    if sum(isnan(NoStopScatterplotY(l, :))) <= size(SubjectNum, 1) - MinimumSubjectsForAverage
        OnlySSDsUsedForAverage(l, 1) = SSDMin+(l-1)*50; 
        meanNoStopScatterplotY(l, 1) = nanmean(NoStopScatterplotY(l, :));
        meanSigRespScatterplotY(l, 1) = nanmean(SigRespScatterplotY(l, :));
        NoSigMinusSigResp = meanNoStopScatterplotY - meanSigRespScatterplotY;
    end
end

meanNoStopScatterplotY(meanNoStopScatterplotY==0) = NaN;
meanSigRespScatterplotY(meanSigRespScatterplotY==0) = NaN;
NoSigMinusSigResp(NoSigMinusSigResp==0) = NaN;

figure;
for f=1:(size(SubjectNum, 1))
    plot(NoStopScatterplotX(:, f), (NoStopScatterplotY(:, f) - (SigRespScatterplotY(:, f))), 'b')
    hold on; 
end

plot(OnlySSDsUsedForAverage(:, 1), NoSigMinusSigResp, 'r', 'LineWidth', 4)
xlabel('SSD')
ylabel('PrecedingNoStopRT-StopFailRT (red=average)')

figure; 
for f=1:(size(SubjectNum, 1))
    plot(SigRespScatterplotX(:, f), SigRespScatterplotY(:, f), 'b')
    hold on; 
end

plot(OnlySSDsUsedForAverage(:, 1), meanSigRespScatterplotY, 'r', 'LineWidth', 4)
hold on; 
plot(OnlySSDsUsedForAverage(:, 1), meanNoStopScatterplotY, 'g', 'LineWidth', 4)
xlabel('SSD')
ylabel('RT (blue=individualSF, red=meanSF, green=precedingNS)')

%Create a cumulative frequency distribution for each subject at each SSD
%that includes sig-resp RT and no-stop-signal RT that precedes sig-resp. 

SortedSigRespOutput = sort(SigRespOutput, 1);
SortedNoStopOutput = sort(NoStopOutput, 1); 

for i=1:NumberOfSSDs
    for f=1:(size(SubjectNum))
        k = 0; 
        for j=1:size(SigRespOutput, 1)
            if(isnan(SigRespOutput(j, i, f))) == 0;
                k = k + 1;
            else
                break; 
            end  
        end
        figure
        plot(SortedSigRespOutput(1:k, i, f), 1:k, 'r')  
        hold on;
        plot(SortedNoStopOutput(1:k, i, f), 1:k, 'g')  
    end
end
                
    
%Produce a scatterplot of GoRT and SSD over time for each subject, only
%using the 3 most frequent SSDs and the goRTs that immediately precede them
for ff=1:size(SubjectNum)
    SubjectNumber2 = SubjectNum(ff);
    for gg=1:size(SubjectSeq)
        if (mod(gg, TrialCount) == 1) && SubjectSeq(gg) == SubjectNumber2
            for hh=1:TrialCount
                if SSDSeq(gg+hh-1) ==  (SSDMin + (SortedSSDIndices(1, ff)*50-50)) || SSDSeq(gg+hh-1) ==  (SSDMin + (SortedSSDIndices(2, ff)*50-50)) || SSDSeq(gg+hh-1) ==  (SSDMin + (SortedSSDIndices(3, ff)*50-50))
                    SSDScatterplotY(hh, ff) = SSDSeq(gg+hh-1);
                    SSDScatterplotX(hh, ff) = hh;
                else
                    SSDScatterplotY(hh, ff) = NaN;
                    SSDScatterplotX(hh, ff) = NaN;
                end
                if hh > 1    
                    if GoRTSeq(gg+hh-2) > -500 && (SSDSeq(gg+hh-1) ==  (SSDMin + (SortedSSDIndices(1, ff)*50-50)) || SSDSeq(gg+hh-1) ==  (SSDMin + (SortedSSDIndices(2, ff)*50-50)) || SSDSeq(gg+hh-1) ==  (SSDMin + (SortedSSDIndices(3, ff)*50-50)))
                        GoRTScatterplotY(hh, ff) = GoRTSeq(gg+hh-2);
                        GoRTScatterplotX(hh, ff) = hh;
                    else
                        GoRTScatterplotY(hh, ff) = NaN;
                        GoRTScatterplotX(hh, ff) = NaN;
                    end
                else
                    GoRTScatterplotY(hh, ff) = NaN;
                    GoRTScatterplotX(hh, ff) = NaN;
                end
            end
        end
    end
    figure
    scatter(SSDScatterplotX(:, ff), SSDScatterplotY(:, ff), 'b')
    hold on;
    scatter(GoRTScatterplotX(:, ff), GoRTScatterplotY(:, ff), 'g')
end

% Create a SSDCount x SSD x subj matrix of go RTs that immediately precede
% stop trials and a sig-resp count x SSD x subj matrix of signal-respond
% trials
 for t=1:(size(SubjectNum))
    w = SSDMin;
    for u=1:(size(SSDCount,1))
        x=1;
        cc=1;
        for v=2:(size(SubjectSeq))
            if(w == SSDSeq(v)) && GoRTSeq(v-1) > 0 && Block(v) == Block(v-1) && SubjectSeq(v) == SubjectNum(t)
                GoRTMinus1(x, u, t) = GoRTSeq(v-1);
                x = x + 1;
            end
            if(w == SSDSeq(v)) && SigRespRT(v) > 0 && SubjectSeq(v) == SubjectNum(t)
                SigRespRTList(cc, u, t) = SigRespRT(v);
                cc = cc + 1;
            end
        end
    w = w + 50;
    end
 end

% Create stop trial number x N Matrix
for b=1:(size(SubjectNum))
    SubjectNumber = SubjectNum(b);
    for c = 1:(size(SSD))
        if (SubjectNumber == SubjectSSD(c))
            SSDOut(a, b) = SSD(c);
            a = a + 1;
        end
    end
    a = 1;
end

%Sort stop trial number x N Matrix of SSD by descending SSDs
SSDFinal = sort(SSDOut, 1);

% Creates a SSDRange x N Matrix that counts the number of trials at each
% SSD for each subjects. Named SSDCount
for e=1:(size(SubjectNum))
    for g=SSDMin:50:SSDMax
        for f=1:(size(SSDFinal, 1))
            if(g == SSDFinal(f, e));
                h = h + 1;
            end
        end
        SSDCount(i, e) = h;
        i = i + 1;
        h = 0; 
    end
    i=1;
end

%Sorts SSDCount Matrix from most to least frequent SSDs for each subject (Sorted SSDs).
%Index for the corresponding SSD in SortedSSDIndeces. Then finds the
%observed mean sig resp, observed median signal resp, and observed
%p(resp|signal) at the NSSD most frequent SSDs
for j=1:size(SubjectNum)
    [SortedSSD, SortedSSDIndex] = sort(SSDCount(:, j), 'descend');
    SortedSSDs(:, j) = SortedSSD;
    SortedSSDIndices(:, j) = SortedSSDIndex;
    for r=1:NSSD
        ObservedMeanMax(j, r) = ObservedMeanSigResp(j, SortedSSDIndices(r, j));
        ObservedMedianMax(j, r) = ObservedMedianSigResp(j, SortedSSDIndices(r, j));
        ObservedPRespMax(j, r) = 1 - ObservedStopAcc(j, SortedSSDIndices(r, j));
        ObservedSigRespCountMax(j, r) = ObservedSigRespCount(j, SortedSSDIndices(r, j));
    end
end

% Create histograms of SSD distributions for each subject
for bb=1:(size(SubjectNum))
    figure
    SSDBins = SSDMin:50:SSDMax;
    hist(SSDFinal(:, bb), SSDBins);   
    xlim([plottedSSDMin, plottedSSDMax])
    yLimit = max(max(SortedSSDs));
    ylim([0, yLimit])
    pbaspect([2 1 1])
    [nelements,centers] = hist(SSDFinal(:, bb), SSDBins)
end


skew = skewness(SSDFinal);
SD = std(SSDFinal);
InterQuartileRange = iqr(SSDFinal);  

% Create a SSDCount x SSD x subj matrix of go RTs that immediately precede
% stop trials and a sig-resp count x SSD x subj matrix of signal-respond
% trials
 for t=1:(size(SubjectNum))
    w = SSDMin;
    for u=1:(size(SSDCount,1))
        x=1;
        cc=1;
        for v=2:(size(SubjectSeq))
            if(w == SSDSeq(v)) && GoRTSeq(v-1) > 0 && Block(v) == Block(v-1) && SubjectSeq(v) == SubjectNum(t)
                GoRTMinus1(x, u, t) = GoRTSeq(v-1);
                x = x + 1;
            end
            if(w == SSDSeq(v)) && SigRespRT(v) > 0 && SubjectSeq(v) == SubjectNum(t)
                SigRespRTList(cc, u, t) = SigRespRT(v);
                cc = cc + 1;
            end
        end
    w = w + 50;
    end
 end

 
GoRTMinus1(GoRTMinus1==0) = NaN;
SigRespRTList(SigRespRTList==0) = NaN;
SortedGoRTMinus1 = sort(GoRTMinus1);
SortedSigRespRTList = sort(SigRespRTList);

% Get some measures of the RTs that immediately precede a stop trial
for y=1:size(SubjectNum)
    for z=1:NSSD
        ObservedMeanMaxRTMinus1(y, z) = nanmean(SortedGoRTMinus1(:, SortedSSDIndices(z, y), y));
        ObservedMedianMaxRTMinus1(y, z) = nanmedian(SortedGoRTMinus1(:, SortedSSDIndices(z, y), y));
        ObservedRTMinus1Count(y, z) = numel(SortedGoRTMinus1(:, SortedSSDIndices(z, y), y)) - sum(isnan(SortedGoRTMinus1(:, SortedSSDIndices(z, y), y)));
    end
end  
 
%Get the Indices for after sorting the NSSD Max values for each subject by
%their stop-signal delay
[ae, af] = sort(SortedSSDIndices(1:NSSD, :));
RankOrderedSSDs = transpose(af);

bootstrapIterations = 5000;
for dd=1:size(SubjectNum)
    for ee=1:NSSD
        if SigRespRTList(5, SortedSSDIndices(ee, dd), dd) > 0
            [bootstat, bootsam] = bootstrp(bootstrapIterations, @nanmedian, SigRespRTList(:, SortedSSDIndices(ee, dd), dd));
            bootstatOutput(dd, ee, :) = bootstat;
            ConfidenceInterval = bootci(bootstrapIterations, {@nanmedian, SigRespRTList(:, SortedSSDIndices(ee, dd), dd)},'alpha', 0.05, 'type', 'per');
            CIOutput(dd, ee, :) = ConfidenceInterval;
        end
    end
end  

CIOutputLower(:, :) = CIOutput(:, :, 1);
CIOutputHigher(:, :) = CIOutput(:, :, 2);

%Produce a scatterplot of GoRT and SSD over time for each subject. 
for ff=1:size(SubjectNum)
    SubjectNumber2 = SubjectNum(ff);
    for gg=1:size(SubjectSeq)
        if (mod(gg, TrialCount) == 1) && SubjectSeq(gg) == SubjectNumber2
            for hh=1:TrialCount
                if SSDSeq(gg+hh-1) > -500 
                    SSDScatterplotY(hh, ff) = SSDSeq(gg+hh-1);
                    SSDScatterplotX(hh, ff) = hh;
                else
                    SSDScatterplotY(hh, ff) = NaN;
                    SSDScatterplotX(hh, ff) = NaN;
                end
                if GoRTSeq(gg+hh-1) > -500 
                    GoRTScatterplotY(hh, ff) = GoRTSeq(gg+hh-1);
                    GoRTScatterplotX(hh, ff) = hh;
                else
                    GoRTScatterplotY(hh, ff) = NaN;
                    GoRTScatterplotX(hh, ff) = NaN;
                end
            end
        end
    end
    figure
    scatter(SSDScatterplotX(:, ff), SSDScatterplotY(:, ff), 'b')
    hold on;
    scatter(GoRTScatterplotX(:, ff), GoRTScatterplotY(:, ff), 'g')
end

%Generate matrix of SSD changes up or down for each subject
for ii = 1:size(SubjectNum)
    for jj = 2:length(SSDOut)
        if SSDOut(jj, ii) > SSDOut(jj-1, ii);
            UpOrDown(jj-1, ii) = 1;
        else
            UpOrDown(jj-1, ii) = 0;
        end
    end
end
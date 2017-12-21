clear all

%Subject number in first column, block number in 2nd column, go RT in 3rd
%column, SSD in 4th column, sig-resp RT in 5th column. Replace non-numbers in columns 3, 4, and 5 with a number
%below the SSD min (I usually use -500). This kludge seems to work
[SubjectSeq Block GoRTSeq SSDSeq SigRespRT] = textread('SequentialRTsSSDsTurkSimpleHigh.txt', '%f%f%f%f%f');

SigRespCountCutoff = 3; %Threshold for the number of signal-respond trials at a specific SSD for that subject to be computed
MinimumSubjectsForAverage = 1; %Threshold for the number of subjects that pass the SigRespCountCutoff at that SSD for that SSD to be included in the group average

%Could hardcode SSDMin or SSDMax if you only want to evaluate a subset of
%the SSD distribution
SSDMin = min(SSDSeq(SSDSeq > -500)); %change to a number below the minimum SSD
SSDMax = max(SSDSeq);
NumberOfSSDs = size(unique(SSDSeq), 1)-1; %last -1 is to account for the -500 kludge mentioned above 
%SubjectNum = unique(SubjectSeq); %create a list of unique subject numbers
[SubjectNum] = textread('TurkN339.txt', '%f'); % can load in a subset of subject numbers here, as is necessary for our
% Turk data
SSDRequired = [150 200 250 300]; % use this to look only at subjects who pass the SigRespCountCutoff on these SSDs

NoStopOutput = zeros(1, NumberOfSSDs, size(SubjectNum, 1));
SigRespOutput = zeros(1, NumberOfSSDs, size(SubjectNum, 1));

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
    h = h + 1; 
    end
end

SigRespScatterplotYRequired = zeros(NumberOfSSDs, size(SubjectNum, 1));
SigRespScatterplotXRequired = zeros(NumberOfSSDs, size(SubjectNum, 1));
NoStopScatterplotYRequired = zeros(NumberOfSSDs, size(SubjectNum, 1));
NoStopScatterplotXRequired = zeros(NumberOfSSDs, size(SubjectNum, 1));

for j=1:(size(SubjectNum, 1))
    if(mean(ismember(SSDRequired(:), NoStopScatterplotX(:, j))) == 1);
        SigRespScatterplotYRequired(:, j) = SigRespScatterplotY(:, j);
        SigRespScatterplotXRequired(:, j) = SigRespScatterplotX(:, j);
        NoStopScatterplotYRequired(:, j) = NoStopScatterplotY(:, j);
        NoStopScatterplotXRequired(:, j) = NoStopScatterplotX(:, j);
    else
        SigRespScatterplotYRequired(:, j) = NaN;
        SigRespScatterplotXRequired(:, j) = NaN;
        NoStopScatterplotYRequired(:, j) = NaN;
        NoStopScatterplotXRequired(:, j) = NaN;           
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
    scatter(NoStopScatterplotX(:, f), (NoStopScatterplotY(:, f) - (SigRespScatterplotY(:, f))), 'b')
    hold on; 
end

plot(OnlySSDsUsedForAverage(:, 1), NoSigMinusSigResp, 'r', 'LineWidth', 4)
xlabel('SSD')
ylabel('PrecedingNoStopRT-StopFailRT (red=average)')

figure; 
for f=1:(size(SubjectNum, 1))
    scatter(SigRespScatterplotX(:, f), SigRespScatterplotY(:, f), 'b')
    hold on; 
end

plot(OnlySSDsUsedForAverage(:, 1), meanSigRespScatterplotY, 'r', 'LineWidth', 4)
hold on; 
plot(OnlySSDsUsedForAverage(:, 1), meanNoStopScatterplotY, 'g', 'LineWidth', 4)
xlabel('SSD')
ylabel('RT (blue=individualSF, red=meanSF, green=precedingNS)')

figure;
for m=1:(size(SubjectNum, 1))
    if(ismember(SSDRequired(1), NoStopScatterplotXRequired(:, m)));
        scatter(NoStopScatterplotXRequired(:, m), (NoStopScatterplotYRequired(:, m) - (SigRespScatterplotYRequired(:, m))), 'b')
    hold on; 
    end
end

n = 1; 
for o=min(SSDRequired):50:max(SSDRequired)
    meanNoStopScatterplotYRequired(o) = nanmean(NoStopScatterplotYRequired(n, :));
    meanNoStopScatterplotXRequired(o) = nanmean(NoStopScatterplotXRequired(n, :));
    meanSigRespScatterplotYRequired(o) = nanmean(SigRespScatterplotYRequired(n, :));
    meanSigRespScatterplotXRequired(o) = nanmean(SigRespScatterplotXRequired(n, :));
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
n = n + 1; 
end


plot(SSDRequired, mean(NoStopScatterploYRequired(NoSigMinusSigResp, 'r', 'LineWidth', 4)
xlabel('SSD')
ylabel('PrecedingNoStopRT-StopFailRT (red=average)')

SortedSigRespOutput = sort(SigRespOutput, 1);
SortedNoStopOutput = sort(NoStopOutput, 1); 
clear all

%Subject number in first column, block number in 2nd column, go RT in 3rd
%column, SSD in 4th column, sig-resp RT in 5th column. Replace non-numbers in columns 3, 4, and 5 with a number
%below the SSD min (I usually use -500). This kludge seems to work
[SubjectSeq Block GoRTSeq SSDSeq SigRespRT] = textread('SequentialRTsSSDsTurkMotorSelec.txt', '%f%f%f%f%f');

SigRespCountCutoff = 2; %Threshold for the number of signal-respond trials at a specific SSD for that subject to be computed
MinimumSubjectsForAverage = 5; %Threshold for the number of subjects that pass the SigRespCountCutoff at that SSD for that SSD to be included in the group average

%Could hardcode SSDMin or SSDMax if you only want to evaluate a subset of
%the SSD distribution
SSDMin = min(SSDSeq(SSDSeq > -500)); %change -500 to a number below the minimum SSD
SSDMax = max(SSDSeq);
SSDIncrement = 50; %incremental difference between each SSD

NumberOfSSDs = size(SSDMin:SSDIncrement:SSDMax, 2);
%SubjectNum = unique(SubjectSeq); %create a list of unique subject numbers
[SubjectNum] = textread('TurkN339.txt', '%f'); % can load in a subset of subject numbers here, as is necessary for our
% Turk data
SSDRequired = [100 150 200]; % use this to look only at subjects who pass the SigRespCountCutoff on these SSDs


v = 1; 
for t=SSDMin:SSDIncrement:SSDMax
    w = 1; 
    for u=min(SSDRequired):SSDIncrement:max(SSDRequired)
        if(u==t)
            SSDRequiredIndices(w) = v; 
        end
        w = w + 1; 
    end
    v = v + 1;
end
 
NoStopOutput = zeros(1, NumberOfSSDs, size(SubjectNum, 1));
SigRespOutput = zeros(1, NumberOfSSDs, size(SubjectNum, 1));
SubjectLevelViolation = zeros(1, NumberOfSSDs, size(SubjectNum, 1));

for a=1:(size(SubjectNum, 1))
    SubjectNumber = SubjectNum(a);
    d = 1;
    for c=SSDMin:SSDIncrement:SSDMax 
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
    for g=SSDMin:SSDIncrement:SSDMax
        if(isnan(NoStopOutput(SigRespCountCutoff, h, f))) == 0;
            SigRespScatterplotY(h, f) = (nanmean(SigRespOutput(:, h, f)));
            SigRespScatterplotX(h, f) = g;
            NoStopScatterplotY(h, f) = (nanmean(NoStopOutput(:, h, f)));
            NoStopScatterplotX(h, f) = g;
            if(NoStopScatterplotY(h, f) < SigRespScatterplotY(h, f));
                Violation(h, f) = 1;
            else
                Violation(h, f) = 0; 
            end
        else
            SigRespScatterplotY(h, f) = NaN;
            SigRespScatterplotX(h, f) = NaN;
            NoStopScatterplotY(h, f) = NaN;
            NoStopScatterplotX(h, f) = NaN;
            Violation(h, f) = NaN;
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
    FullSSDList(l, 1) = SSDMin+(l-1)*SSDIncrement;
    if sum(isnan(NoStopScatterplotY(l, :))) <= size(SubjectNum, 1) - MinimumSubjectsForAverage
        OnlySSDsUsedForAverage(l, 1) = SSDMin+(l-1)*SSDIncrement; 
        meanNoStopScatterplotY(l, 1) = nanmean(NoStopScatterplotY(l, :));
        meanSigRespScatterplotY(l, 1) = nanmean(SigRespScatterplotY(l, :));
        NoSigMinusSigResp = meanNoStopScatterplotY - meanSigRespScatterplotY;
        ViolationRateAverage(1, l) = nanmean(Violation(l, :)); 
    end
end

meanNoStopScatterplotY(meanNoStopScatterplotY==0) = NaN;
meanSigRespScatterplotY(meanSigRespScatterplotY==0) = NaN;
NoSigMinusSigResp(NoSigMinusSigResp==0) = NaN;
%ViolationRateAverage(ViolationRateAverage==0) = NaN;

DiffIncompleteCasesX = transpose(OnlySSDsUsedForAverage);
DiffIncompleteCasesY = transpose(NoSigMinusSigResp);

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
for f=1:(size(SubjectNum, 1))
    if(ismember(SSDRequired(1), NoStopScatterplotXRequired(:, f)));
        scatter(NoStopScatterplotXRequired(:, f), (NoStopScatterplotYRequired(:, f) - (SigRespScatterplotYRequired(:, f))), 'b')
    hold on; 
    end
end

plot((nanmean(SigRespScatterplotXRequired(SSDRequiredIndices, :), 2)), ((nanmean(NoStopScatterplotYRequired(SSDRequiredIndices, :), 2)) - (nanmean(SigRespScatterplotYRequired(SSDRequiredIndices, :), 2))), 'r', 'LineWidth', 4);
DiffXRequired = transpose((nanmean(SigRespScatterplotXRequired(SSDRequiredIndices, :), 2)));
DiffYRequired = transpose((nanmean(NoStopScatterplotYRequired(SSDRequiredIndices, :), 2)) - (nanmean(SigRespScatterplotYRequired(SSDRequiredIndices, :), 2)));

SortedSigRespOutput = sort(SigRespOutput, 1);
SortedNoStopOutput = sort(NoStopOutput, 1); 

DifferenceFullMatrix = NoStopScatterplotY - SigRespScatterplotY;

figure;
plot(OnlySSDsUsedForAverage(:, 1), ViolationRateAverage, 'r', 'LineWidth', 4)
xlabel('SSD')
ylabel('Violation Rate')
axis([0 inf 0 1])
%{

Displays results from several experiments combined together.
version 180XXX v0.1

S. Herbert sherbert@pasteur.fr

%}

function displayOverall()
% Display overall results based on combined experiments

clear

% PARAMS
PARAMS.scriptVersion = 'v0.1';

% Load individual experimental files
filePathAna = uipickfiles('Prompt', 'Please select the path to the individual analysis files');

% Initialize global analysis output files
globalInterMinTimes = []; % list of time between minima

% Initialize output files
localMinsVarOutList = {'time' 'prominence' 'speed' 'delayToNextMin' 'frameInterval' 'experimentName'}; % all variables of interest
allLocalMinsList = cell2table(cell(0,length(localMinsVarOutList)), 'VariableNames', localMinsVarOutList);

% do the same for the interMin times

% For each file
PARAMS.md.timeUnits = '';
expName = {};
for dataFile = 1:length(filePathAna)
    % Parse and concatenate results of interest
    %     [~,fileName,~] = fileparts(filePathAna{dataFile});
    indivResults = load(filePathAna{dataFile});
    
    % Check that time units is always the same
    if dataFile == 1
        PARAMS.md.timeUnits = indivResults.outputAnalysis.PARAMS.md.timeUnits;
    elseif ~strcmp(PARAMS.md.timeUnits, indivResults.outputAnalysis.PARAMS.md.timeUnits)
        error('Time units are not always the same! Stopping the analysis!')
    end
    
    % Create structure to save individual results (could use table)
    localMins{dataFile}.expName = indivResults.outputAnalysis.PARAMS.fileName;
    localMins{dataFile}.localMinsP = indivResults.outputAnalysis.dataTable.localMinsP;
    localMins{dataFile}.timeCourse = indivResults.outputAnalysis.dataTable.timeCourse;
    localMins{dataFile}.speed = indivResults.outputAnalysis.dataTable.speedSmoothed;
    
    % Calculation of interminima durations
    localMins{dataFile}.interMinTimes = ...
        diff(localMins{dataFile}.timeCourse(localMins{dataFile}.localMinsP~=0));
    
    % Merge into a single array the individual analyses
    globalInterMinTimes = [globalInterMinTimes; localMins{dataFile}.interMinTimes];
    
    % Merge into a single output table
    temp = array2table([localMins{dataFile}.timeCourse(localMins{dataFile}.localMinsP~=0)  ...
        localMins{dataFile}.localMinsP(localMins{dataFile}.localMinsP~=0)  ...
        localMins{dataFile}.speed(localMins{dataFile}.localMinsP~=0) ...
        [localMins{dataFile}.interMinTimes;nan]]); % Create a table with the experimental results
    [expNameTable{1:height(temp)}] = deal(localMins{dataFile}.expName);  % Create a maching table with the experiment name
    [expFrameIntervalTable{1:height(temp)}] = deal(indivResults.outputAnalysis.PARAMS.md.frameInterval); % % idem for frame interval
    temp = [temp expFrameIntervalTable' expNameTable']; % Combine the 3 tables    
    temp.Properties.VariableNames = localMinsVarOutList; % Adapt the variable names
    
    allLocalMinsList = [allLocalMinsList ; temp]; % Merge with previous runs
    
    clear temp expNameTable expFrameIntervalTable
end

% Additionnal calculations


% Ask user where to put the overall analysis
PARAMS.filePathOutput = uigetdir([],'Please select the path to save the analysis output');

% Display figures
displayInterMins(globalInterMinTimes, localMins, PARAMS);

% Save all data
writetable(allLocalMinsList, [PARAMS.filePathOutput filesep 'allLocalMinsList'])

end


function displayInterMins(globalInterMinTimes, localMins, PARAMS)
% Display boxplot of the global and individual results for the inter minimum times
figure
hold on

colorExp = lines(length(localMins));

% Display individual measures
for exp = 1:length(localMins)
    h{exp} = scatter(ones(1,length(localMins{exp}.interMinTimes))+...
        randn(1,length(localMins{exp}.interMinTimes))/20,... % adds jiggling
        localMins{exp}.interMinTimes,...
        [], colorExp(exp,:));
    legs{exp} = localMins{exp}.expName;
end

% Display overall boxplot
boxplot(globalInterMinTimes);

title('Duration distribution between successive minima');
ylabel( sprintf('Duration (%s)', PARAMS.md.timeUnits) );
legend(legs,'Location','southoutside');
ylim([0 max(globalInterMinTimes)+min(globalInterMinTimes)]);

saveas(gcf,sprintf('%scycleDurationDistri',...
    [PARAMS.filePathOutput filesep]));
saveas(gcf,sprintf('%scycleDurationDistri.png',...
    [PARAMS.filePathOutput filesep]));
clear h

% Display histogram
figure
h{1} = histogram(globalInterMinTimes,20);

medianExp = median(globalInterMinTimes);
perExp75 = prctile(globalInterMinTimes,75);
perExp25 = prctile(globalInterMinTimes,25);
x = [perExp25, perExp75, perExp75, perExp25];
y = [0 0 max(h{1}.Values) max(h{1}.Values)];
h{2} = line([medianExp, medianExp], [0 max(h{1}.Values)],  ylim, 'Color', 'r','LineStyle','-');
h{3} = line([perExp25, perExp75 ; perExp25, perExp75], [0 0 ;max(h{1}.Values) max(h{1}.Values)],...
    ylim, 'Color', 'r','LineStyle','--');
patch(x, y, 'r', 'EdgeColor','None','FaceAlpha',0.1);
uistack(h{1},'top');

title('Duration distribution between successive minima');
xlabel(sprintf('Duration (%s)', PARAMS.md.timeUnits)); ylabel('N');

legend([h{1} h{2} h{3}(1,:)],...
    {'Distribution' sprintf('median = %0.1f µm',medianExp)...
    sprintf('exp iqr [%0.1f-%0.1f] µm', perExp25, perExp75)});

saveas(gcf,sprintf('%shisto_cycleDurationDistri',...
    [PARAMS.filePathOutput filesep]));
saveas(gcf,sprintf('%shisto_cycleDurationDistri.png',...
    [PARAMS.filePathOutput filesep]));

end





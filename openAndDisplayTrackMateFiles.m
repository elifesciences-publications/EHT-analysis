%{ 

Tracking of A-P luminal distance 
Input is the xml output of the trackMate plugin in Fiji. It expects only 2
tracks of the 2 extremities.

version 180605 v1.3
- saves the output of the trajectories into a .mat file

version 180612 v1.3.1
S. Herbert sherbert@pasteur.fr 

%}



function openAndDisplayTrackMateFiles()

clear
% close all

% PARAMS
PARAMS.scriptVersion = 'v1.3';

% Display and save the analysis
PARAMS.doDisplays = 0; % Display the figures
PARAMS.doSave = 1; % Save the analysis file


% Use a smoothing factor in the display => Always keep 1 as the first value
% to keep the raw data
PARAMS.smoothFact = [1 5 11]; % to smooth the openings display and analyses
PARAMS.smoothFactAdv = 11; % to apply and advanced smoothing => Change to the number of steps to smooth onto

% Minimum prominence value before local minimum rejection 
PARAMS.promThresh = 0.01; % Minimum prominence value before local minimum rejection 


%% Double check some parameters of the acquisition
prompt = {'Enter frame interval:', 'Enter time unit:', 'Enter space unit:'};
dlg_title = 'Check inputs';
num_lines = 1;
defaultans = {'2','min','micron'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

PARAMS.md.frameInterval = str2double(answer{1});
PARAMS.md.timeUnits = answer{2};
PARAMS.md.spaceUnits = answer{3};

% Select xml files of interest
filePathSpotFeat = uipickfiles('Prompt', 'Please select the path to the .xml track files (keep same imaging conditions!)');
% filePathSpotFeat = ...
%     {'/media/sherbert/Data/Projects/OG_projects/Project4_ML/movies/160328_projected/280316_extremities1and2.xml'};

localMins = {};
for dataFile = 1:length(filePathSpotFeat)
    
    % Do the analysis for each datafile
    [path,PARAMS.fileName,~] = fileparts(filePathSpotFeat{dataFile});
    PARAMS.fullPath = [path filesep PARAMS.fileName];
    mkdir(PARAMS.fullPath);
    cd(PARAMS.fullPath);
    fprintf('Processing xml file: %s\n', PARAMS.fileName);
    
    % Import the data table associated
    [ spot_table, spot_ID_map ] = trackmateSpots( filePathSpotFeat{dataFile} );
    edge_map = trackmateEdges( filePathSpotFeat{dataFile} );
    track_names = edge_map.keys;
       
    % Do the analysis of this experiment
    analysisTracks = analyseTracks(PARAMS, spot_table, spot_ID_map, track_names, edge_map);
   
    % Display the analysis of this experiment
    if PARAMS.doDisplays
        displayAnalysis(PARAMS, analysisTracks, track_names);
        close all
    end
    
    % Create and save output structure
    if PARAMS.doSave % Save the analysis file
        
        outputAnalysis = {};
        outputAnalysis.dataTable = table(analysisTracks.timeCourse,...
            analysisTracks.openingRTadv.distance(:,3), ...
            analysisTracks.openingRTadv.speed(:,3),...
            analysisTracks.localMinIdx);

        outputAnalysis.dataTable.Properties.VariableNames = {'timeCourse' 'distanceSmoothed'...
            'speedSmoothed' 'localMinsP'};
        outputAnalysis.PARAMS = PARAMS;
        
        save('outputAnalysis','outputAnalysis');
    end 
    
end



end




function analysisTracks = analyseTracks(PARAMS, spot_table, spot_ID_map, track_names, edge_map)

% Check that there are only 2 tracks
n_tracks = numel(track_names);
if n_tracks>2
    error('WARNING: %d trajectories detected. Stopping the analysis.\n',...
        n_tracks);
end

% Recreate the tracks by spot ID
track_spot_IDs = recreate_IDs(n_tracks, track_names, edge_map);

% Reshape tracks into simpleTracks for simple handling
clipZ = true;
simpleTracks = reshapeTracks(n_tracks, track_spot_IDs,...
    spot_ID_map, spot_table, clipZ, PARAMS);

% Reshape simpleTracks to incorporate the empty positions => and make it
% Real Time ? 
simpleTracksRT = rtTracks(n_tracks, simpleTracks, PARAMS);

% Apply smoothing?

%% %%%%%%%%%%%%%%%% 2D only %%%%%%%%%%%%%%%%

% Calculate distance between the 2 extremities
openingRT = InterExtremDist(simpleTracksRT,PARAMS.smoothFact);

% Calculate the opening change
openingRTspeed = [-diff(openingRT) ; NaN(1,numel(PARAMS.smoothFact))]/PARAMS.md.frameInterval;

% Calculate advanced filtering
openingRTadv = advancedFilter(simpleTracksRT, PARAMS.smoothFactAdv, PARAMS);
% => first field = distance ; second field = speed;

% Set timecourse
timeCourse = simpleTracksRT{1}.time;

% Automated detection of local minima
[localMinIdx, localMinIdxRej] = findLocalMin(openingRTadv.speed(:,3), PARAMS);


% Packing analysis into a single structure
analysisTracks = {};
analysisTracks.simpleTracksRT = simpleTracksRT;
analysisTracks.openingRT = openingRT;
analysisTracks.openingRTspeed = openingRTspeed;
analysisTracks.openingRTadv = openingRTadv;
analysisTracks.localMinIdx = localMinIdx;
analysisTracks.localMinIdxRej = localMinIdxRej;
analysisTracks.timeCourse = timeCourse;

end


function displayAnalysis(PARAMS, analysisTracks, track_names)

% unpacking analysis into a separate structures
simpleTracksRT = analysisTracks.simpleTracksRT;
openingRT = analysisTracks.openingRT;
openingRTspeed = analysisTracks.openingRTspeed;
openingRTadv = analysisTracks.openingRTadv;
localMinIdx = analysisTracks.localMinIdx;
localMinIdxRej = analysisTracks.localMinIdxRej;
timeCourse = analysisTracks.timeCourse;

n_tracks = numel(track_names);

% Preping the legends
for smoothing = 1: numel(PARAMS.smoothFact)
    if PARAMS.smoothFact(smoothing)==1
        legs{smoothing} = 'No smooth';
    else
        legs{smoothing} = sprintf('Smooth over %d tp',PARAMS.smoothFact(smoothing));
    end
end

%% tracks and surface size
% Display the distance between the 2 extremities and associate tracks
figure;

subplot(2,2,1); % Display the opening size
dispOpeningDia(openingRT, timeCourse, PARAMS, legs); 

subplot(2,2,2); % Display the tracks
dispTracks(n_tracks, simpleTracksRT, track_names, PARAMS);

subplot(2,2,3);
dispOpeningSpeed(openingRTspeed, timeCourse, PARAMS, legs);

subplot(2,2,4);
nbins = 20;
dispOpeningHisto(openingRTspeed, PARAMS, legs, nbins);

saveas(gcf,sprintf('%s_morphoAnalysis',...
    [PARAMS.fullPath, filesep, PARAMS.fileName]));

%% Display the tracks intensities along time
figure;

subplot(2,2,1); % Maximum intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'maxInt',...
    'Max intensity', track_names, PARAMS);

subplot(2,2,2); % Total intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'totInt',...
    'Total intensity', track_names, PARAMS);

subplot(2,2,3); % Mean intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'meanInt',...
    'Mean intensity', track_names, PARAMS);

subplot(2,2,4); % Median intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'medianInt',...
    'Median intensity', track_names, PARAMS);

saveas(gcf,sprintf('%s_fluoAnalysis',...
    [PARAMS.fullPath, filesep, PARAMS.fileName]));

%% Display overlayed closing speed and closing distance
% figure;
%
% dispOverlayDistVsSpeed(timeCourse, openingRT, openingRTspeed);
%
% saveas(gcf,sprintf('%s_originalOverlay', [path, filesep, PARAMS.fileName], PARAMS.smoothFactAdv));
    
%% Display overlayed closing speed and closing distance for advanced filtering

figure;

dispOverlayAdvanced(timeCourse, openingRTadv, PARAMS.smoothFactAdv, PARAMS);

tempFig = gcf;
saveas(tempFig,sprintf('%s_overlayAdvFiltering_%ddt',...
    [PARAMS.fullPath, filesep, PARAMS.fileName], PARAMS.smoothFactAdv));
% saveas(gcf,sprintf('%s_overlayAdvFiltering_%ddt.png',[path, filesep, PARAMS.fileName], PARAMS.smoothFactAdv));
% set(tempFig,'PaperOrientation','landscape');
% print(tempFig, '-fillpage', '-dpdf','tempName')%, sprintf('%s_overlayAdvFiltering_%ddt.pdf',...
% %     [path, filesep, PARAMS.fileName], PARAMS.smoothFactAdv));

%% Display local minimums
dispLocalMin(timeCourse, openingRTadv.speed(:,3), localMinIdx, localMinIdxRej, PARAMS);
tempFig = gcf;
saveas(tempFig,sprintf('%s_localMinsP',...
    [PARAMS.fullPath, filesep, PARAMS.fileName]));
saveas(tempFig,sprintf('%s_localMinsP.png',...
    [PARAMS.fullPath, filesep, PARAMS.fileName]));


end


function [localMinIdx, localMinIdxRej] = findLocalMin(speed, PARAMS)
% Estimate automatically where the minima are in the speed curve

[~, localMinIdx] = islocalmin(speed);
% [~, localMinIdx] = islocalmin(speed,'MinProminence',PARAMS.promThresh); % use a
% treshold on minima prominence
localMinIdxRej = localMinIdx;
localMinIdxRej(localMinIdxRej>PARAMS.promThresh) = 0;
localMinIdx = localMinIdx - localMinIdxRej;

end

function dispLocalMin(timeCourse, speed, localMinIdx, localMinIdxRej, PARAMS)
% Display local minima in curve

localColors = lines(2);

% display minima
figure
hold on
h{1} = plot(timeCourse, speed, '.-', 'LineWidth',1 , 'Color', localColors(1,:)); % speed curve
h{2} = plot(timeCourse(localMinIdx~=0),speed(localMinIdx~=0)-0.01,...
    'v', 'MarkerFaceColor', [0,0.5,0], 'MarkerEdgeColor', [0,0.5,0]); % minima on the speed curve
h{3} = plot(timeCourse(localMinIdx~=0),localMinIdx(localMinIdx~=0),...
    'v', 'MarkerFaceColor', localColors(2,:), 'MarkerEdgeColor', localColors(2,:)); % prominence of the minima
h{4} = plot(timeCourse(localMinIdxRej~=0),speed(localMinIdxRej~=0)-0.01,...
    'v', 'MarkerEdgeColor', [0,0.5,0]); % minima on the speed curve
h{5} = plot(timeCourse(localMinIdxRej~=0),localMinIdx(localMinIdxRej~=0),...
    'v', 'MarkerEdgeColor', localColors(2,:)); % prominence of the minima
grid on
title(sprintf('Closing speed (%s)', PARAMS.fileName));
ylabel( sprintf('Closing speed (%s/%s)', PARAMS.md.spaceUnits, PARAMS.md.timeUnits) );
xlabel( [ 'Time (' PARAMS.md.timeUnits ')' ] );
legend({'Speed', 'Local minimum', 'Prominence', 'Rejected local min', 'Rejected Prominence'});

end



function dispIntensities(n_tracks, timeCourse, simpleTracksRT, fieldToPlot,...
    varName, track_names, PARAMS)

hold on;
for s = 1 : n_tracks
    %     track_name = track_names{s};
    plot(timeCourse, simpleTracksRT{s}.(fieldToPlot), '.-');    
end
title(sprintf('%s against time', varName));
xlabel( ['Time (' PARAMS.md.timeUnits ')'] )
ylabel(varName);
legend(track_names, 'Location', 'eastoutside');
end

function [simpleTracksRT, maxTime] = rtTracks(n_tracks, simpleTracks, PARAMS)
% Create and pad with nan the empty positions

simpleTracksRT = simpleTracks;
maxTimeTrack = zeros(n_tracks,1);
for tk = 1 : n_tracks
    maxTimeTrack(tk) = max(simpleTracks{tk}.frame);
end
maxTime = max(maxTimeTrack);

% Prepare the empty line to add at the tracks
TableFields = simpleTracks{1}.Properties.VariableNames;
emptyLine = num2cell(NaN(1,length(TableFields)));

for tk = 1 : n_tracks % for each track
    emptyFrameTable = table;
    frameExist = ismember(0:maxTime,simpleTracks{tk}.frame);
    for f = 1:length(frameExist) % for each frame
        if ~frameExist(f) % if the frame doesn't exist in the original track
            % adapt time and frame
            emptyLine{1} = f-1;
            emptyLine{2} = emptyLine{1}*PARAMS.md.frameInterval;
            emptyFrameTable = [emptyFrameTable ; emptyLine];
        end
    end
    if ~isempty(emptyFrameTable)    
        % Rename the fields of the emptyFrameTable
        emptyFrameTable.Properties.VariableNames = TableFields;
        % Merge tables of empty and filled frames
        simpleTracksRT{tk} = [simpleTracks{tk};emptyFrameTable];
    end
    % Resort the table based on frames
    simpleTracksRT{tk} = sortrows(simpleTracksRT{tk} ,{'frame'},{'ascend'});
end

end

function simpleTracks = reshapeTracks(n_tracks, track_spot_IDs,...
    spot_ID_map, spot_table, clipZ, PARAMS)
% Reshape the tracks into a table to handle them more easily

simpleTracks = {};
for s = 1 : n_tracks
    track_spot_ID = track_spot_IDs{ s };
    rows = cell2mat(spot_ID_map.values(num2cell(track_spot_ID)));
    frame = spot_table.FRAME(rows);
    time = spot_table.FRAME(rows) .* PARAMS.md.frameInterval;
    xPos = spot_table.POSITION_X(rows);
    yPos = spot_table.POSITION_Y(rows);
    if ~clipZ
        zPos = spot_table.POSITION_Z(rows);
    end
    maxInt = spot_table.MAX_INTENSITY(rows);
    totInt = spot_table.TOTAL_INTENSITY(rows);
    meanInt = spot_table.MEAN_INTENSITY(rows);
    medianInt = spot_table.MEDIAN_INTENSITY(rows);
    % Could also import the edge displacement but the 'row' order is not very
    % clear to me...
    
    % merge data into a table and sort the tracks based on the .frame info
    if clipZ
        tempTable = table(frame,time,xPos,yPos,maxInt,totInt,meanInt,medianInt);
    else
        tempTable = table(frame,time,xPos,yPos,xPos,maxInt,totInt,meanInt,medianInt);
    end
    tempTable = sortrows(tempTable,{'frame'},{'ascend'});
    simpleTracks{s} = tempTable;  
end

end

function dispTracks(n_tracks, simpleTracksRT, track_names, PARAMS)
% Display the tracks of the 2 extremities (2D only for the moment)
hold on
for s = 1 : n_tracks
    track_name = track_names{ s};
    x = simpleTracksRT{s}.xPos;
    y = simpleTracksRT{s}.yPos;
    % Plot the tracks by coloring spots.
    plot( x, y, '.-' , 'DisplayName', track_name)
end

title('Trajectories of the extremities');
ylabel( [ 'Y (' PARAMS.md.spaceUnits ')' ] )
xlabel( [ 'X (' PARAMS.md.spaceUnits ')' ] )
axis equal
legend toggle
end

function track_spot_IDs = recreate_IDs(n_tracks, track_names, edge_map)
track_spot_IDs = cell( n_tracks, 1 );
% Recreate the tracks by spot ID 

for s = 1 : n_tracks
    track_name = track_names{s};
    edge_table = edge_map( track_name );
    track_spot_IDs{ s } = unique( [...
        edge_table.SPOT_SOURCE_ID...
        edge_table.SPOT_TARGET_ID] );
end
end

function openingRTadv = advancedFilter(simpleTracksRT, smoothFactAdv, PARAMS)

% Calculate distance with every filters
diam = InterExtremDist(simpleTracksRT,[1 smoothFactAdv]); % => no smooth and single smooth
diam(:,3) = smooth(diam(:,2),smoothFactAdv); % => 2 smoothing steps

diam = diam .* ~isnan(diam(:,1));
diam(diam==0) = nan;

% Calculate speed with every filters
speed = ([-diff(diam) ; NaN(1,size(diam,2))]/PARAMS.md.frameInterval);

% merge into main structure
openingRTadv.distance = diam;
openingRTadv.speed = speed;

% Prepape legends
openingRTadv.dLegend = {'no smooth diameter',...
    sprintf('diameter smoothed once (%ddt)',smoothFactAdv-1),...
    sprintf('diameter smoothed twice (%ddt)',smoothFactAdv-1)};
openingRTadv.sLegend = {'no smooth speed',...
    sprintf('speed smoothed once (%ddt)',smoothFactAdv-1),...
    sprintf('speed smoothed twice (%ddt)',smoothFactAdv-1)};

end

function opening = InterExtremDist(simpleTracksRT,smoothFact)
% show the distance between the 2 extremeties of the cell surface
% 2D only

opening = zeros(length(simpleTracksRT{1}.frame),length(smoothFact));

track1 = [simpleTracksRT{1}.xPos simpleTracksRT{1}.yPos];
track2 = [simpleTracksRT{2}.xPos simpleTracksRT{2}.yPos];

nanTrack = ~isnan(simpleTracksRT{1}.xPos .* simpleTracksRT{2}.xPos);

for smoothing = 1:numel(smoothFact)
    opening(:,smoothing) = smooth(sqrt(sum((track2-track1).^2,2)),smoothFact(smoothing));
end

opening = opening.*nanTrack;
opening(opening==0) = nan;

end

function dispOpeningDia(openingRT, timeCourse, PARAMS, legs)
% Display the size of the surface diameter (between the extremities)

plot(timeCourse, openingRT, '.-');
title('Distance between extremities');
xlabel( ['Time (' PARAMS.md.timeUnits ')'] );
ylabel( ['Distance (' PARAMS.md.spaceUnits ')' ]);
legend(legs)

end

function dispOpeningSpeed(openingRTspeed, timeCourse, PARAMS, legs)
% Display the speed at which the surface is closing

hold on
plot(timeCourse, openingRTspeed, '.-');
title('Closing speed');
xlabel( sprintf('Time (%s)', PARAMS.md.timeUnits) ); 
ylabel( sprintf('Closing speed (%s/%s)', PARAMS.md.spaceUnits, PARAMS.md.timeUnits) );
legend(legs)

end

function dispOpeningHisto(openingRTspeed, PARAMS, legs, nbins)
% Display of the opening speeds as an histogram

lineColors = lines(numel(legs));

minmax = min(min(openingRTspeed)) : ...
    abs(min(min(openingRTspeed))/max(max(openingRTspeed)))/nbins : ...
    max(max(openingRTspeed));

h = zeros(numel(minmax)-1,numel(legs));

for smoothing = 1:numel(legs)
    h(:,smoothing) = histcounts(openingRTspeed(:,smoothing),minmax);
end
histoData = bar(minmax(1:end-1),h);

for pop = 1:numel(legs)
   histoData(pop).FaceColor = lineColors(pop,:);
   histoData(pop).EdgeColor = 'None';
end

title('Distribution of the closing speeds');
xlabel( sprintf('Closing speed (%s/%s)', PARAMS.md.spaceUnits, PARAMS.md.timeUnits) );
ylabel('N'); 
legend(legs);
legend boxoff;

end

function dispOverlayDistVsSpeed(timeCourse, openingRT, openingRTspeed)
% Display overlayed closing speed and closing distance

yyaxis left
plot(timeCourse,openingRT);

yyaxis right
plot(timeCourse,openingRTspeed);

end

function dispOverlayAdvanced(timeCourse, openingRTadv, smoothFactAdv, PARAMS)
% Display overlayed closing speed and closing distance for advanced filtering

hold on

yyaxis left
plot(timeCourse,openingRTadv.distance);
ylabel( sprintf('Diameter (%s)', PARAMS.md.spaceUnits) );

yyaxis right
plot(timeCourse,openingRTadv.speed);
ylabel( sprintf('Closing speed (%s/%s)', PARAMS.md.spaceUnits, PARAMS.md.timeUnits) );

title(sprintf('Advanced filtering (twice over %ddt <=> %0.1fmin)',...
    smoothFactAdv,smoothFactAdv*PARAMS.md.frameInterval));
xlabel( sprintf('Time (%s)', PARAMS.md.timeUnits) ); 

legend([openingRTadv.dLegend, openingRTadv.sLegend],'Location','EastOutside')

end







function openAndDisplayTrackMateFiles()

clear
close all

filePathTracks = uipickfiles('num',1 ,'Prompt', 'Please select the path to the tracks');
% filePathTracks = ...
%     {'/media/sherbert/Data/Projects/OG_projects/Project4_ML/movies/160328_projected/280316_extremities1and2_trackExport.xml'};

filePathSpotFeat = uipickfiles('num',1 ,'Prompt', 'Please select the path to the whole dataset');
% filePathSpotFeat = ...
%     {'/media/sherbert/Data/Projects/OG_projects/Project4_ML/movies/160328_projected/280316_extremities1and2.xml'};

% Use a smoothing factor in the display
smoothFact = [1 5 11]; % to smooth the openings display and analyses


% Import tracks
clipZ = true;
scaleT = false;
[~, md] = importTrackMateTracks(filePathTracks{1}, clipZ, scaleT);

% Until we figure out why there is an issue with the default time interval
temp = inputdlg('define temporal step','define temporal step (in sec)',1,{num2str(md.frameInterval)});
md.frameInterval = str2num(temp{1});

% Import the data table associated
[ spot_table, spot_ID_map ] = trackmateSpots( filePathSpotFeat{1} );
edge_map = trackmateEdges( filePathSpotFeat{1} );
track_names = edge_map.keys;

n_tracks = length(edge_map);
if n_tracks>2
    fprintf('WARNING: %d trajectories detected. Stopping the analysis.\n',...
        n_tracks);
    return
end


% Recreate the tracks by spot ID 
track_spot_IDs = recreate_IDs(n_tracks, track_names, edge_map);

% Reshape tracks into simpleTracks for simple handling
simpleTracks = reshapeTracks(n_tracks, track_spot_IDs,...
    spot_ID_map, spot_table, clipZ, md);

% Reshape simpleTracks to incorporate the empty positions => and make it
% Real Time ? 
simpleTracksRT = rtTracks(n_tracks, simpleTracks, md);

% Apply smoothing?

%% %%%%%%%%%%%%%%%% 2D only %%%%%%%%%%%%%%%%

% Calculate distance between the 2 extremities
openingRT = InterExtremDist(simpleTracksRT,smoothFact);

% Calculate the opening change
openingRTspeed = [-diff(openingRT) ; NaN(1,numel(smoothFact))]/md.frameInterval;

% Set timecourse
timeCourse = simpleTracksRT{1}.time;

% Preping the legends
for smoothing = 1: numel(smoothFact)
    if smoothFact(smoothing)==1
        legs{smoothing} = 'No smooth';
    else
        legs{smoothing} = sprintf('Smooth over %d tp',smoothFact(smoothing));
    end
end

%% tracks and surface size
% Display the distance between the 2 extremities and associate tracks
figure;

subplot(2,2,1); % Display the opening size
dispOpeningDia(openingRT, timeCourse, md, legs); 

subplot(2,2,2); % Display the tracks
dispTracks(n_tracks, simpleTracksRT, track_names, md);

subplot(2,2,3);
dispOpeningSpeed(openingRTspeed, timeCourse, md, legs);

subplot(2,2,4);
dispOpeningHisto(openingRTspeed, md, legs);

%% Display the tracks intensities along time
figure;

subplot(2,2,1); % Maximum intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'maxInt',...
    'Max intensity', track_names, md);

subplot(2,2,2); % Total intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'totInt',...
    'Total intensity', track_names, md);

subplot(2,2,3); % Mean intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'meanInt',...
    'Mean intensity', track_names, md);

subplot(2,2,4); % Median intensity 
dispIntensities(n_tracks, timeCourse, simpleTracksRT, 'medianInt',...
    'Median intensity', track_names, md);


    
end

function dispIntensities(n_tracks, timeCourse, simpleTracksRT, fieldToPlot,...
    varName, track_names, md)

hold on;
for s = 1 : n_tracks
    %     track_name = track_names{s};
    plot(timeCourse, simpleTracksRT{s}.(fieldToPlot), '.-');    
end
title(sprintf('%s against time', varName));
xlabel( ['Time (' md.timeUnits ')'] )
ylabel(varName);
legend(track_names, 'Location', 'eastoutside');
end




function [simpleTracksRT, maxTime] = rtTracks(n_tracks, simpleTracks, md)
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
            emptyLine{2} = emptyLine{1}*md.frameInterval;
            emptyFrameTable = [emptyFrameTable ; emptyLine];
        end
    end
    % Rename the fields of the emptyFrameTable
    emptyFrameTable.Properties.VariableNames = TableFields;
    % Merge tables of empty and filled frames
    simpleTracksRT{tk} = [simpleTracks{tk};emptyFrameTable];
    % Resort the table based on frames
    simpleTracksRT{tk} = sortrows(simpleTracksRT{tk} ,{'frame'},{'ascend'});
end

end

function simpleTracks = reshapeTracks(n_tracks, track_spot_IDs,...
    spot_ID_map, spot_table, clipZ, md)
% Reshape the tracks into a table to handle them more easily

simpleTracks = {};
for s = 1 : n_tracks
    track_spot_ID = track_spot_IDs{ s };
    rows = cell2mat(spot_ID_map.values(num2cell(track_spot_ID)));
    frame = spot_table.FRAME(rows);
    time = spot_table.FRAME(rows) .* md.frameInterval;
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


function dispTracks(n_tracks, simpleTracksRT, track_names, md)
% Display the tracks of the 2 extremities (2D only for the moment)
hold on
for s = 1 : n_tracks
    track_name = track_names{ s};
    x = simpleTracksRT{s}.xPos;
    y = simpleTracksRT{s}.yPos;
    % Plot the tracks by coloring spots.
    plot( x, y, '.-' , 'DisplayName', track_name)
end

units = md.spaceUnits;
title('Trajectories of the extremities');
ylabel( [ 'Y (' units ')' ] )
xlabel( [ 'X (' units ')' ] )
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


function dispOpeningDia(openingRT, timeCourse, md, legs)
% Display the size of the surface diameter (between the extremities)

plot(timeCourse, openingRT, '.-');
title('Distance between extremities');
xlabel( ['Time (' md.timeUnits ')'] );
ylabel( ['Distance (' md.spaceUnits ')' ]);
legend(legs)

end

function dispOpeningSpeed(openingRTspeed, timeCourse, md, legs)
% Display the speed at which the surface is closing

hold on
plot(timeCourse, openingRTspeed, '.-');
title('Closing speed');
xlabel( sprintf('Time (%s)', md.timeUnits) ); 
ylabel( sprintf('Closing speed (%s/%s)', md.spaceUnits, md.timeUnits) );
legend(legs)

end

function dispOpeningHisto(openingRTspeed, md, legs)
% Display of the opening speeds as an histogram

nbins = 20;

minmax = min(min(openingRTspeed)) : ...
    abs(min(min(openingRTspeed))/max(max(openingRTspeed)))/nbins : ...
    max(max(openingRTspeed));

h = zeros(numel(minmax)-1,numel(legs));

for smoothing = 1:numel(legs)
    h(:,smoothing) = histcounts(openingRTspeed(:,smoothing),minmax);
end
bar(minmax(1:end-1),h);

title('Distribution of the closing speeds');
xlabel( sprintf('Closing speed (%s/%s)', md.spaceUnits, md.timeUnits) );
ylabel('N'); 
legend(legs);
legend boxoff;

end




% Thesis: RSSI-based localization - Simulation
% Generate robot path based on inputed points
% Lai Quoc Bao - 01/04/2021


velocity = 1; % m/s              
sampTime = 0.01; % s
turnRadi = 2/pi; % m, turning circle radious

turnPoint      = [ 1 1; 9 1; 9 9; 1 9; 1 1 ]; % turning points/boundary coordinates of the AGV
% turnPoint      = [ 1 1; 8 1; 8 5; 13 5; 13 9; 6 9; 6 5; 1 5; 1 1 ]; % zigzag line
initPose       = [turnPoint(1,:), 0];
euDis          = zeros(1,length(turnPoint)-1);

%% Create straight line to location: Turnpoint-turnRadi
% calculate euclid distance from point to point
for ii = 1: length(turnPoint)-1
    euDis(ii) = sqrt(sum( (turnPoint(ii+1,:)-turnPoint(ii,:)).^2 ));
end

straightNum = round((euDis - 2*turnRadi)/(velocity*sampTime)); % initialize the number of straight points

% initialize zeros cell array of straight  points's coordination
straightCoord = cell(4,1);
for ii = 1:length(straightCoord)
    straightCoord{ii} = nan(straightNum(ii),2);
end


for ii = 1:length(straightNum)
    headAngle = atan2d(turnPoint(ii+1,2)-turnPoint(ii,2),turnPoint(ii+1,1)-turnPoint(ii,1));
    
    for jj = 1:straightNum(ii)
        straightCoord{ii}(jj,1) = turnPoint(ii,1)+ cosd(headAngle)*turnRadi + cosd(headAngle)*velocity*sampTime*jj; % update X
        straightCoord{ii}(jj,2) = turnPoint(ii,2)+ sind(headAngle)*turnRadi + sind(headAngle)*velocity*sampTime*jj; % update Y
    end
%     plot(straightCoord{ii}(:,1),straightCoord{ii}(:,2))
end

%% Create turning radious for each turnPoint
% calculate the turning angle
turnAngle = nan(1,length(turnPoint)-1);

for ii = 1 : length(turnPoint)-1
    if ii == length(turnPoint)-1
        turnAngle(ii) = atan2d(turnPoint(2,2)-turnPoint(1,2),turnPoint(2,1)-turnPoint(1,1))...
                    -atan2d(turnPoint(ii+1,2)-turnPoint(ii,2),turnPoint(ii+1,1)-turnPoint(ii,1));
        headAngle     = atan2d(turnPoint(ii+1,2)-turnPoint(ii,2),turnPoint(ii+1,1)-turnPoint(ii,1));
        % Make all the angle positive
        if headAngle < 0
            headAngle = 360 + headAngle;
        end
    else
        % head angle on current segment of path
        headAngle     = atan2d(turnPoint(ii+1,2)-turnPoint(ii,2),turnPoint(ii+1,1)-turnPoint(ii,1));
        % Make all the angle positive
        if headAngle < 0
            headAngle = 360 + headAngle;
        end
        % head angle on the next segment of path
        headAnglenext = atan2d(turnPoint(ii+2,2)-turnPoint(ii+1,2),turnPoint(ii+2,1)-turnPoint(ii+1,1));
        if headAnglenext < 0
            headAnglenext = 360 + headAnglenext;
        end
        turnAngle(ii) = headAnglenext - headAngle;
    end
end

turnArray = cell(1,length(turnAngle));
turnCoord = cell(1,length(turnAngle));
for ii = 1:length(turnArray)
    
    turnRes  = (abs(turnAngle(ii))/180*pi)/(velocity/turnRadi*sampTime); % turning Resolution
    
    for jj = 1:turnRes
        turnArray{ii}(jj)   = turnAngle(ii)/turnRes;
        if jj == 1
            turnCoord{ii}(jj,1) = straightCoord{ii}(end,1) + cosd(90+headAngle+turnArray{ii}(jj))*velocity*sampTime;
            turnCoord{ii}(jj,2) = straightCoord{ii}(end,2) + sind(90+headAngle+turnArray{ii}(jj))*velocity*sampTime;
        else
        turnCoord{ii}(jj,1) = turnCoord{ii}(jj-1,1) + cosd(90+headAngle+turnArray{ii}(jj))*velocity*sampTime;
        turnCoord{ii}(jj,2) = turnCoord{ii}(jj-1,2) + sind(90+headAngle+turnArray{ii}(jj))*velocity*sampTime;
        end
        headAngle = headAngle+turnArray{ii}(jj);
    end
%     plot(turnCoord{ii}(:,1),turnCoord{ii}(:,2));
end

%% Join all segments
% calculate total number of points
turnLength     = zeros(1,length(turnCoord));
straightLength = turnLength;
for ii = 1:length(turnCoord)
    turnLength(ii)     = length(turnCoord{ii});
    straightLength(ii) = length(straightCoord{ii});
end
pathLength = sum(turnLength) + sum(straightLength);
% Start joining the arrays into one
pathCoord = cell(1,4);
for ii = 1:length(turnCoord)
    pathCoord{ii} = cat(1,straightCoord{ii},turnCoord{ii});
end
pathCoord = cat(1,pathCoord{:});
pathCoord(end+1,:)=pathCoord(1,:);

%% Plot junk
figure(1);
clf
hold on;
grid on;
plot(pathCoord(:,1),pathCoord(:,2))
scatter(turnPoint(:,1),turnPoint(:,2),'r','x')
axis square

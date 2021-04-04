%% Thesis: RSSI-based localization - Simulation
%  Path generating function. Input velocity and sampling time
%  Output path points coordinates' and AGV's head angle
%  Lai Quoc Bao - 25/03/2021

%% Initialize all variables
%==========================================================================
velocity = 1; % m/s              
sampTime = 0.05; % s

turnPoint      = [ 1 1; 9 1; 9 9; 1 9; 1 1 ]; % turning points/boundary coordinates of the AGV
initPose      = [turnPoint(1,:), 0];
segmentLength  = zeros(length(turnPoint)-1,1); % initialize the segment's length vector
pointsNum      = zeros(length(segmentLength),1); % initialize the number of points vector
% End Initialize all variables
%--------------------------------------------------------------------------

%% Start the Path generation
%==========================================================================
for ii = 1: length(turnPoint)-1
    
    segmentLength(ii) = sqrt( sum( (turnPoint(ii+1,:) - turnPoint(ii,:)).^2));
    % calculate the number of points along each segment
    pointsNum(ii)      = segmentLength(ii)/(velocity*sampTime);
    
    if ii == 1
        startPose = initPose; % start pose = initial pose only in 1st iteration
    else
        startPose = lastPose;
    end
    
    % start the points' coordinate calculation
    for jj = 1: pointsNum(ii)
        
        if jj == 1
            pointX  = turnPoint(ii,1);
            pointY  = turnPoint(ii,2);
            oldPose = startPose;
            pointTheta1 = oldPose(3);
            pointTheta2 = 0;
        elseif jj == 160
            pointX     = turnPoint(ii+1,1);
            pointY     = turnPoint(ii+1,2);
            pointTheta2 = atan2d(pointY - oldPose(2), pointX - oldPose(1) );
            lastPose   = [pointX pointY pointTheta];
            
        else
            pointTheta1 = 0;
            pointX = oldPose(1) + cosd(pointTheta1)*velocity*sampTime;
            pointY = oldPose(2) + sind(pointTheta1)*velocity*sampTime;
            pointTheta2 = oldPose(3) + atan2d( pointY - oldPose(2), pointX - oldPose(1) );
        end
        newPose = [pointX pointY pointTheta];
        oldPose = newPose;
        robotPose(jj,:) = newPose;
    end
       
end






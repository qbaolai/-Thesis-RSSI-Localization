% Thesis: RSSI-based localization - Simulation
% Generate the WiFi node and Robot in simulation environment
% data has been filtered through Particle filter with error calculated
% with reference to the predetermined route, then to EKF localization algorithm
% Lai Quoc Bao - 09/04/2021

% Start simulation
%==========================================================================

clc;
clear;
close all; 

% Create Environment
%==========================================================================
% Environment Variables
A = 38;
n = 2;
len = 15;   % environment length (m)
wid = 15;   % environment width (m) 
% End Import Environment
%--------------------------------------------------------------------------

% Generate Nodes
%==========================================================================
% Node Params
numOfNodes = 5; % 3 for testing trialteration
nodeRes = 8;

% Pre select the node position
prenode(1).pos = [3 3];
prenode(2).pos = [7 3];
prenode(3).pos = [7 7]; 
prenode(4).pos = [3 7];
prenode(5).pos = [5 5];
for ii = 1:numOfNodes
    % Types can be Outdoors, Hallway
    [node(ii),A(ii),n(ii)] = generateNode(prenode(ii).pos,nodeRes, 3, 'Outdoors');
end

% End Generate Nodes
%--------------------------------------------------------------------------

% Generate Robot
%==========================================================================
% Robot Params
velocity = 0.5; % m/s              
sampTime = 0.5; % s
% robotPath: array of velocity and delta rotation
% robotPose: the X,Y and head angle of the robot, at ii = 1, it's initPose
% pathCoord: Coordination of calculated path points
[robotPath, robotPose,pathCoord] = generatePath(velocity,sampTime);

robotPath = [robotPath;robotPath;robotPath]; % repeat n times
pathCoord = [pathCoord;pathCoord;pathCoord];

errPer = 0.1;  % odometry error percentage
% End Generate Robot
%--------------------------------------------------------------------------

o2nhist = zeros(1,numOfNodes);  % prepare odometry history
% Simulation
%--------------------------------------------------------------------------
for ii = 1:length(robotPath)
    % Follow specified path
    %======================================================================
    % Rotate some amount
    oldRobotPose = robotPose;
    newRobotPose(3) = robotPose(3) + robotPath(ii,2);
    % Move forward X meters
    newRobotPose(1) = robotPose(1) + cosd(newRobotPose(3))*robotPath(ii,1);
    newRobotPose(2) = robotPose(2) + sind(newRobotPose(3))*robotPath(ii,1);
    robotPose = newRobotPose;
    locationHistory(:,ii) = robotPose;
    
    % End Follow Specified Path
    %----------------------------------------------------------------------
    
    % Get odometry data
    %======================================================================
    % Calculate odometry to introduce real errors that may be encountered
    % during actual implementation
    disp = (robotPose(1:2) - oldRobotPose(1:2));
    disp = sqrt(disp(1)^2 + disp(2)^2) + normrnd(0,robotPath(ii,1)*errPer);
    u = [disp, (robotPose(3) - oldRobotPose(3)) + normrnd(0,1)];
    
    % Get dead reckoning position
    if(ii == 1)
        dr_pose = oldRobotPose + [u(1)*cosd(oldRobotPose(3)+u(2)) u(1)*sind(oldRobotPose(3)+u(2)) u(2)];
        ohist = u(1);   % odometry history;
    else
        dr_pose = [dr_pose; (dr_pose(ii-1,1) + u(1)*cosd(dr_pose(ii-1,3)+u(2))), ...
            (dr_pose(ii-1,2) + u(1)*sind(dr_pose(ii-1,3)+u(2))), ...
            (dr_pose(ii-1,3) + u(2))];
        ohist = [ohist;u(1)];
        % Get dead reckoned distance to each node
        for jj = 1:numOfNodes
            temp(jj) = euclidDist(dr_pose(ii,1:2),node(jj).pos);
        end
        o2nhist = [o2nhist;temp];
    end
    %     u = [robotPath(ii,1),robotPath(ii,2)];
    % End odometry data
    %----------------------------------------------------------------------
    
    % Get RSSI Value and distance estimate to each node
    %======================================================================
    for jj = 1:numOfNodes
        rssi(jj) = getRSSI(node(jj),robotPose);
        if(rssi(jj)>70)
            rssi(jj) = oldrssi(jj);
        end
        oldrssi(jj) = rssi(jj);
        dist2node(jj) = euclidDist(node(jj).pos,robotPose);
        pathDist2node(jj) = euclidDist(node(jj).pos,pathCoord(ii,:));
        d(jj) =  10^((rssi(jj)-A(jj))/(10*n(jj))); % Log-Distance Model
%         if(dist2node(jj) > max(node(jj).dist))
%             warning('Distance to node exceeds maximum value.');
%         end
    end
    
    % Record history of rssi and distance values
    if(ii == 1)
        dhist = d;  % estimated distance history
        rhist = rssi;   % measured rssi history
        d2nhist = dist2node;    % distance to node history
    else
        dhist = [dhist;d];
        rhist = [rhist;rssi];
        d2nhist = [d2nhist;dist2node];
    end
    % End Get RSSI Value and distance estimate to each node
    %----------------------------------------------------------------------
    % Get change in distance relative to each node
    %======================================================================
    if(ii == 1)
        x = [oldRobotPose];
        for jj = 1:numOfNodes
            x((jj-1)*2 + 4) = node(jj).pos(1);
            x((jj-1)*2 + 5) = node(jj).pos(2);
            initNodePos(jj,:) = [x((jj-1)*2 + 4),x((jj-1)*2 + 5)];
        end
    end
    if(ii == 2)  % Will run filter after 2 measurements are made to integrate
        % Also to match vector lenght of rhist after velocity differentiation
        for jj = 1:numOfNodes
            d1 = euclidDist([x((jj-1)*2 + 4),x((jj-1)*2 + 5)],x(1:2));
            % get updated state variables from measured odometry
            odom_x = [x(1) + u(1)*cosd(x(3)+u(2)),x(2) + u(1)*sind(x(3)+u(2))];
            d2 = euclidDist([x((jj-1)*2 + 4),x((jj-1)*2 + 5)],odom_x(1:2));
            v_d2n(jj) = d2 - d1;
        end
        vhist = v_d2n;
    elseif(ii > 2)
        for jj = 1:numOfNodes
            d1 = euclidDist([x((jj-1)*2 + 4),x((jj-1)*2 + 5)],x(1:2));
            % get updated state variables from measured odometry
            odom_x = [x(1) + u(1)*cosd(x(3)+u(2)),x(2) + u(1)*sind(x(3)+u(2))];
            d2 = euclidDist([x((jj-1)*2 + 4),x((jj-1)*2 + 5)],odom_x(1:2));
            v_d2n(jj) = d2 - d1;
        end
        vhist = [vhist;v_d2n];
    end
    % End get change in distance relative to each node
    %----------------------------------------------------------------------
    
    % RSSI Particle Filter
    %======================================================================
    if(ii == 1)
        partDist = dist2node;
    else
        if(~exist('numParts'))
            numParts = 1000;
            for jj = 1:numParts
                % Actual random particle values
                % particles(jj,:) = rand([1 numOfNodes])*10 + normrnd(0,3,[1,numOfNodes]);
                % Initialize particles with estimated distances
%                 particles(jj,:) = d + normrnd(0,1,[1,numOfNodes]);
                % Cheese it with perfect initialization
                particles(jj,:) = dist2node + normrnd(0,0.1,[1,numOfNodes]);
            end
            newParticles = zeros(size(particles));
        end
        % Use odometry and estimated states to update particles
        for jj = 1:numOfNodes
            % get updated state variables from measured odometry
            odom_x = [x(1) + u(1)*cosd(x(3)+u(2)),x(2) + u(1)*sind(x(3)+u(2))];
            % Extract landmark positions from state vector
            lmx = x((jj-1)*2 + 4);
            lmy = x((jj-1)*2 + 5);
            % Prediction vector contains changes to each node
            prediction(jj) = euclidDist(odom_x,[lmx,lmy]) - ...
                euclidDist(x(1:2),[lmx,lmy]);
        end
        if(ii <= 2)  
            p_hist = prediction;
        else
            p_hist = [p_hist;prediction];
        end
        % Apply changes to particles
        for jj = 1:numParts
            particles(jj,:) = particles(jj,:) + prediction;
        end
        % Get fitness value for each particle
        err = zeros(size(particles));
        for jj = 1:numParts
            err(jj,:) = abs(pathDist2node-particles(jj,:));
        end
        for jj = 1:numParts
            err(jj,:) = abs(err(jj,:) - max(err));
        end
        for jj = 1:numParts
            err_dist(jj,:) = err(jj,:)./sum(err);
        end
        % Resample from distribution
        cdf = zeros(size(particles));
        for i=1:numParts
            cdf(i,:)=sum(err_dist(1:i,:));
        end
        
        for jj = 1:numParts
            for kk = 1:numOfNodes
                randval = min(find(cdf(:,kk) > rand(1)));
                newParticles(jj,kk) = particles(randval,kk);
            end
        end
        particles = newParticles;
        % Get most probable position
        partDist = sum(particles.*err_dist)./sum(err_dist);
        % Perturb particles
        particles = particles + normrnd(0,0.5,[numParts,numOfNodes]);
    end
    % Calculate R
    % R is calculated based on the difference between the chance in
    % estimated RSSI distance and the odometry distance.
    for jj = 1:numOfNodes
        if(ii > 1)
            % Calculate R using prediction vector from particle filter
%             R(jj) = abs(partDist(jj) - dhist(ii,jj))^2;
            R(jj) = ((partDist(jj)*abs(dhist(ii,jj) - dhist(ii-1,jj))));
            %                  R(jj) = 100;
            % Calculate R using linear odom value
            %                 lmx = x((jj-1)*2 + 4);
            %                 lmy = x((jj-1)*2 + 5);
            %                 % idx is the index of the observed node
            %                 e_d2n = sqrt((x(1) - lmx)^2 + (x(2) - lmy)^2);
            %                 R(jj) = abs(abs(d(jj) - e_d2n) - robotPath(ii,1))^2;
        else
            R(jj) = 100;
        end
    end
    % End RSSI Particle Filter
    %----------------------------------------------------------------------
    
    % RSSI EKF RO-SLAM
    %======================================================================
    % Initialize EKF Variables
    if(ii == 1)
        x = oldRobotPose;
        for jj = 1:numOfNodes
            x((jj-1)*2 + 4) = node(jj).pos(1);
            x((jj-1)*2 + 5) = node(jj).pos(2);
            initNodePos(jj,:) = [x((jj-1)*2 + 4),x((jj-1)*2 + 5)];
        end
        % Covariance Matrix
        P = eye(length(x)).*0.1; % Initializing with near-exact values
        P(1,1) = 0.1; P(2,2) = 0.1; P(3,3) = 0.1;
        % Measurement Noise
        R = ones(numOfNodes,1)*100; % Set initial uncertainty (noise) high
        % Process Noise
        C = 0.0;
        W = [u(1)*cosd(x(3)) u(1)*sind(x(3)) u(2)]';
        Q = zeros(size(P));
        Q(1:3,1:3) = W*C*W';
    end
    
    % Apply EKF for each landmark/node observation
    for jj = 1:numOfNodes
        if(jj > 1)
            u = [0, 0];
        end
        % Apply Range-Only EKF SLAM
        % Calculate Process Noise
        C = 0.1;
        W = [u(1)*cosd(x(3)) u(1)*sind(x(3)) u(2)]';
        Q = zeros(size(P));
        Q(1:3,1:3) = W*C*W';
        % Apply filter
        [x,P] = EKF_RO_SLAM(x,P,partDist(jj),u,jj,R(jj),Q);
%         [x,P] = EKF_RO_SLAM(x,P,dist2node(jj),u,jj,100,Q); % Test EKF
    end
    % record estimated position history
    if(~exist('xhist'))
        xhist = x(1:2);
    else
        xhist = [xhist;x(1:2)];
    end
    % End EKF SLAM
    %----------------------------------------------------------------------
    
    % Plot Movement
    %======================================================================
%     figure(1);
    clf;
    axis([-2 len+1 -2 wid+1]);
    hold on;
    grid on;
    scatter(robotPose(1),robotPose(2),'k','o');
    for jj = 1:numOfNodes
        scatter(node(jj).pos(1),node(jj).pos(2),'green','o');
    end
    % Plot dead reckoned position
    scatter(dr_pose(ii,1),dr_pose(ii,2),'blue','x');
    % Plot filtered position
    scatter(xhist(ii,1),xhist(ii,2),'red','*');
    pause(0.02); % to animate the robot
end
    % End Plot Movement
    %----------------------------------------------------------------------
    
% Calculate + Plot Localization Errors
    %==========================================================================
    for ii = 1:length(robotPath)
        % dead-reckoning error
        drError(ii)     = (euclidDist(dr_pose(ii,1:2),locationHistory(1:2,ii)'));
        % Filter error
        filterError(ii) = (euclidDist(xhist(ii,:),locationHistory(1:2,ii)'));
    end
    % Plot errors
    figure(2);
    hold on;
    plot(drError,'blue');
    plot(filterError,'red');
    grid on;
    legend(strcat('Dead-reackoning error - Mean: ',num2str(mean(drError))),...
           strcat('Filtered error - Mean: ',num2str(mean(filterError))));
    % End Calculate + Plot Localization Errors
    %----------------------------------------------------------------------
    
% Plot Results
    %==========================================================================
    figure(3);
    hold on;
    % true position path
    plot(locationHistory(1,:),locationHistory(2,:),'black');
    % dead reckoning path
    plot(dr_pose(:,1),dr_pose(:,2),'blue');
    title('Simulated Trajectory');
    % particle filtered path
    plot(xhist(:,1),xhist(:,2),'red');
    
    for jj = 1:numOfNodes
        scatter(node(jj).pos(1),node(jj).pos(2),'green','o'); 
    end
    scatter(x(1),x(2),'black','x');
    legend('Ground truth','dead-reckoned','Particle + EKF');
    grid on;
    
%     for jj = 1:numOfNodes
%         figure();
%         subplot(4,1,1:2); hold on;
%         plot(dhist(:,jj),'red');
%         plot(d2nhist(:,jj),'k','linew',1.5);
%         plot(o2nhist(:,jj),'blue');
%         legend('unfiltered','actual','dead-reckoned');
%         title(strcat('Distance to Node over Time: Node ',num2str(jj)));
%         v_rssi = vectorDerivative(dhist(:,jj));
%         subplot(4,1,3:4); hold on;
%         plot(rhist(:,jj),'red');
%         title(strcat('RSSI overtime: Node ',num2str(jj)));
%     end
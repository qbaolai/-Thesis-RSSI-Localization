% Thesis: RSSI-based localization - Simulation
% Generate the WiFi node and visualize the 3D radiation pattern
% Lai Quoc Bao - 17/03/2021

% Start simulation
%==========================================================================

clc;
clear;
close all;

% Import Environment
%==========================================================================
% Environment Variables
len = 5;   % environment length (m)
wid = 5;   % environment width (m)
% End Import Environment
%--------------------------------------------------------------------------

% Generate Nodes
%==========================================================================
% Node Params
nodeRes = 8;
node = generateNode([len/2,wid/2],nodeRes, 0, 'Hallway');
% End Generate Nodes
%--------------------------------------------------------------------------

% Simulation
%--------------------------------------------------------------------------
kk = 1;
for ii = 1:0.1:wid
    hh = 1;
    for jj = 1:0.1:len
        rssi(kk,hh) = getRSSI(node,[ii,jj]);
        hh = hh + 1;
    end
    kk = kk + 1;
end
% End Simulation
%--------------------------------------------------------------------------
figure();
surf(rssi);
title('Simulated Node Radiaton Pattern');
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('RSSI (-dBm)');
set(gca,'XTickLabel',1:1:wid+1 );
set(gca,'YTickLabel',1:1:len+1 );
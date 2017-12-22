function [ intensity ] = splineInt_threeD( sphereCoordinatesModel, sphereElectrodePosition, inputVector )
% SPLINEINT_THREED 
% The function is used to calculate 3D spline interpolation in 3D space.
%
% INPUTS:
%   sphereCoordinatesModel - projection of bsphereCoordinatesModelin model points on the sphere, M x 3, where M is number of model points
%   sphereElectrodePosition - projection of electrode positions points on the sphere, N x 3, where N is number of electrodes
%   inputVector - values of signal under electrodes
%
% OUTPUTS:
%   intensity - intensity of signal in each point of the bsphereCoordinatesModelin model
%
% EXAMPLE
%
% BRIEF EXPLANATION:
%
% SEE ALSO:
%
% Author: Vaclava Piorecka (vaclava.piorecka@fbmi.cvut.cz, vaclava.piorecka@nudz.cz)
% Date:   
% 2017-10-01    creation of function


%% Initialization of vasphereElectrodePositionables
numOfElecs = size(sphereElectrodePosition,1);
numOfPoints = size(sphereCoordinatesModel,1);
m = 3;      % order of spline interpolation

%% Computation of matsphereElectrodePositionx H (this matsphereElectrodePositionx reflects distatnces between elctrodes in space, depeneds on the order of spline interpolation)
H = zeros(numOfElecs);
% w = 0;        % smoothing pasphereCoordinatesModelmeter, not used yet

for NOE1 = 1 : 1 : numOfElecs
    for NOE = 1 : 1 : numOfElecs
        H(NOE1,NOE) = sqrt(((sphereElectrodePosition(NOE1,1)-sphereElectrodePosition(NOE,1))^2 + (sphereElectrodePosition(NOE1,2)-sphereElectrodePosition(NOE,2))^2 + (sphereElectrodePosition(NOE1,3)-sphereElectrodePosition(NOE,3))^2)^(2*m-3));
%         H(NOE1,NOE) = ((sphereElectrodePosition(NOE1,1)-sphereElectrodePosition(NOE,1))^2 + (sphereElectrodePosition(NOE1,2)-sphereElectrodePosition(NOE,2))^2 + (sphereElectrodePosition(NOE1,3)-sphereElectrodePosition(NOE,3))^2 + w^2)^2*log((sphereElectrodePosition(NOE1,1)-sphereElectrodePosition(NOE,1))^2 + (sphereElectrodePosition(NOE1,2)-sphereElectrodePosition(NOE,2))^2 + (sphereElectrodePosition(NOE1,3)-sphereElectrodePosition(NOE,3))^2 +w^2);
    end
end

%% Computation of matsphereElectrodePositionx F (this matsphereElectrodePositionx descsphereElectrodePositionbes relationship between coordinates of each electrode, depeneds on the order of spline interpolation)
F = zeros(numOfElecs,10);
for NOE = 1 : 1 : numOfElecs
    
    F(NOE,:) = [1 sphereElectrodePosition(NOE,1) sphereElectrodePosition(NOE,2) sphereElectrodePosition(NOE,3) sphereElectrodePosition(NOE,1)^2  sphereElectrodePosition(NOE,1)*sphereElectrodePosition(NOE,2) sphereElectrodePosition(NOE,1)*sphereElectrodePosition(NOE,3) sphereElectrodePosition(NOE,2)^2 sphereElectrodePosition(NOE,2)*sphereElectrodePosition(NOE,3) sphereElectrodePosition(NOE,3)^2];
    
end

%% Computation of matsphereElectrodePositionx Ft (this matsphereElectrodePositionx is tsphereCoordinatesModelnsposition of matsphereElectrodePositionx F)
Ft = transpose(F);

%% Help vectors
O1 = zeros(10);
O2 = zeros(1,10);

%% Computation of coefficients P,Q
V = inputVector;

A = [H F; Ft O1];
C = [V O2];

B = linsolve(A,C');
% B = pinv(A)*C';

P = B(1:numOfElecs,1);
Q = B((numOfElecs+1):end);

% test = A*B - C';      % This part test the solution.

%% Computation of intensity in each point of bsphereCoordinatesModelin model
Hm = zeros(numOfPoints,numOfElecs);
U1 = zeros(numOfPoints,1);
for NOP = 1 : 1 : numOfPoints
    for NOE = 1 : 1 : numOfElecs
        Hm(NOP,NOE) = sqrt(((sphereCoordinatesModel(NOP,1)-sphereElectrodePosition(NOE,1))^2 + (sphereCoordinatesModel(NOP,2)-sphereElectrodePosition(NOE,2))^2 + (sphereCoordinatesModel(NOP,3)-sphereElectrodePosition(NOE,3))^2)^(2*m-3));
    end
    
    U1(NOP,1) = sum(P'.*Hm(NOP,:));
end


Qm = zeros(numOfPoints,10);
U2 = zeros(numOfPoints,1);
for NOP = 1 : 1 : numOfPoints   
    Qm(NOP,:) = [1 sphereCoordinatesModel(NOP,1) sphereCoordinatesModel(NOP,2) sphereCoordinatesModel(NOP,3) sphereCoordinatesModel(NOP,1)^2 sphereCoordinatesModel(NOP,1)*sphereCoordinatesModel(NOP,2) sphereCoordinatesModel(NOP,1)*sphereCoordinatesModel(NOP,3) sphereCoordinatesModel(NOP,2)^2 sphereCoordinatesModel(NOP,2)*sphereCoordinatesModel(NOP,3) sphereCoordinatesModel(NOP,3)^2];
    U2(NOP,1) = sum(Q'.*Qm(NOP,:));
end

intensity = U1 + U2;



end


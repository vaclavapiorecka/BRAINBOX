function [ sphericalSplineCoefficients ] = splineInt_spherical( sphereCoordinatesModel, sphereElectrodePosition, m)
% SPLINEINT_SPHERICAL 
% The function is used to calculate spherical spline interpolation in 3D space.
%
% INPUTS:
%   sphereCoordinatesModel - projection of bsphereCoordinatesModelin model points on the sphere, M x 3, where M is number of model points
%   sphereElectrodePosition - projection of electrode positions points on the sphere, N x 3, where N is number of electrodes
%   m - order of spline interpolation
%
% OUTPUTS:
%   sphericalSplineCoefficients - 
%
% EXAMPLE
%
% BsphereElectrodePositionEF EXPLANATION:
%
% SEE ALSO:
%
% Author: Vaclava Piorecka (vaclava.piorecka@fbmi.cvut.cz, vaclava.piorecka@nudz.cz)
% Date:   
% 2017-10-01    creation of function


%% Initialization of vasphereElectrodePositionables
numOfElecs = size(sphereElectrodePosition,1);
numOfPoints = size(sphereCoordinatesModel,1);


%% Calculate C matrix, coefficients of electrodes
Cmatrix = zeros(numOfElecs,numOfElecs);
h = waitbar(0,'Computing spline matrix, coefficients of electrodes');
steps = numOfElecs;

delta = pdist2(sphereElectrodePosition,sphereElectrodePosition);
ei = ones(numOfElecs,numOfElecs)-delta

for i = 1 : 1 : numOfElecs
    for j = 1:numOfElecs
        Cmatrix(i,j) = calculategm(ei(i,j),3,7);
    end
    waitbar(i / steps)
end
% C = inv(Cmatrix);
sphericalSplineCoefficients.C = Cmatrix\eye(size(Cmatrix));
close(h);

%% Calculate Gx matrix, coefficients of model
Gmatrix = zeros(numOfPoints,numOfElecs);
h = waitbar(0,'Computing spline matrix, coefficients of model');
numOfSteps = numOfPoints;

clear ei delta
delta = pdist2(sphereCoordinatesModel,sphereElectrodePosition);
ei = ones(numOfPoints,numOfElecs)-delta;

for i = 1 : 1 : numOfPoints
    gm = zeros(1,numOfElecs);
    for j = 1:numOfElecs
        gm(j) = calculategm(ei(i,j),3,7);
    end
    Gmatrix(i,:) = gm;
    waitbar(i / numOfSteps)
end
sphericalSplineCoefficients.Gx = Gmatrix;

end


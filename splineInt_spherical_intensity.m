function [ intensity ] = splineInt_spherical_intensity( sphericalSplineCoefficients, inputVector, method )
% SPLINEINT_SPHERICAL_INTENSITY 
% The function is used to calculate intenstity in 3D brain model using spherical spline coefficients.
%
% INPUTS:
%   sphericalSplineCoefficients - include matrix of coefficients of model points (Gx) and matrix of coefficients of electrodes
%   inputVector - values of signal under electrodes
%   method - method of computation
%       'meanCorrection' - This correction method consists in calculating the mean value of the signal, subtracting it from the input vector, mapping itself and then adding the mean value back to it.
%       'normal' - Calculate intensity without correction. 
%
% OUTPUTS:
%   intensity - intensity of signal in each point of the sphereCoordinatesModel
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

    [numOfPoints,numOfElecs] = sphericalSplineCoefficients.Gx;
    
    switch method
        case 'meanCorrection'
            meanval = mean(inputVector); 
            inputVector = inputVector - meanval; 

            C = pinv([sphericalSplineCoefficients.C;ones(1,numOfElecs)]) * [inputVector(:);0];

            intensity = zeros(numOfPoints,1);
            for i = 1 : 1 : numOfPoints
                intensity(i) = dot(sphericalSplineCoefficients.C,sphericalSplineCoefficients.Gx(i,:)); 
            end
            intensity = intensity + meanval;

        case 'normal'
            C = pinv([sphericalSplineCoefficients.C;ones(1,numOfElecs)]) * [inputVector(:);0];

            intensity = zeros(numOfPoints,1);
            for i = 1 : 1 : numOfPoints
                intensity(i) = dot(sphericalSplineCoefficients.C,sphericalSplineCoefficients.Gx(i,:)); 
            end
    end

end


function [u0, windDirections, windDistribution, Nwd, Nws] = setWindParameters()


    u0=[8 12 17]; %units: m/s, u0 must be a row vector
    Nws = numel(u0);
    windDirections = [0 10 20];
    Nwd = numel(windDirections);

    %windDistribution must be Nwd by Nws matrix
    windDistribution = [ones(2,1)*[0.0039    0.0085    0.0115];
                        0.0039    0.0108    0.0135];

% 
%     u0=[12]; %units: m/s, u0 must be a row vector
%     Nws = numel(u0);
%     windDirections = [0 90];
%     Nwd = numel(windDirections);
% 
%     %windDistribution must be Nwd by Nws matrix
%     windDistribution = [1;1];

end


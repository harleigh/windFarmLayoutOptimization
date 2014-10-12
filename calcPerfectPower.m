%{
    Perfect Power is the power output of the wind farm under no wake effect
    whatsoever; i.e.: the largest power output possible for the wind farm.

    u0 is a row vector: 1 by Nws
%}
function [ perfectPwr ] = calcPerfectPower( Nwt, Nwd, u0, windDistribution, ...
                                            powerCurveDomain, powerCurve )

    perfectPwr=0;
    U0 = repmat( u0, Nwt,1); %onset wind at each turbine is the free wind speed
    for w=1:Nwd
        turbinePower = interp1(powerCurveDomain, powerCurve, U0,'linear');
        perfectPwr = perfectPwr + sum(sum(bsxfun(@times, windDistribution(w,:),turbinePower)));
    end

end


function [alpha, Ct, gamma, Nwt, Trad, D, L, dMin, powerCurve, powerCurveDomain] = setTurbineParameters()

    h = 60;                %hub height
    Tdiam = 40;            %diameter of turbine blades
    Trad = Tdiam/2;        %radius of turbine blades
    z0 = 0.3;              %roughness constant
    alpha = 0.5/log(h/z0); %wake spreading constant
    Ct = 0.88;             %turbine thrust coefficient
    a=(1-sqrt(1-Ct))/2;    %for calculation of D on next line
    D = Tdiam*sqrt((1-a)/(1-2*a)); %downstream rotor diameter
    gamma = 1-sqrt(1-Ct); %this is for the calculation of tildeU
    powerCurveDomain = [0:30];
    powerCurve = makeMosettiPowerCurve(powerCurveDomain);
    Nwt = 10;
    L = 1800;         %Effective size of a side of the square Windfarm
    dMin = 4*Tdiam;   %minimum safe (relative) distance between turbines
end


function [pwrCurve] = makeMosettiPowerCurve(domU)
    totalPts = numel(domU);
    pwrCurve = zeros(1,totalPts);
    for i=1:totalPts
        if( domU(i) < 2 )
            pwrCurve(i)=0;
        elseif( 2 <= domU(i) && domU(i)<12.8 )
            pwrCurve(i)=0.3*domU(i)^3;
        elseif( 12.8 <= domU(i) && domU(i)< 18 )
            pwrCurve(i) = 629.1;
        else
            pwrCurve(i)=0;
        end %--conditional block
    end%--for loop
end%--function

function [pwrCurve] = waspPowerCurve()
    pwrCurve = [0,0,0,0,2,97,255,459,726,1004,1330,1627,1772,1797,1802 * ones(1,10),1800,1800,zeros(1,5)];
end
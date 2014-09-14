%{
   Given a list of turbine sites on a (square) windfarm of width
   'farmWidth' (see code) plot the turbine locations and display the index
   location of the turbine (index within the turbineSites vector)

   Update: This function now has the features
               Able to draw the wake to a turbine
               Draw the minimum safe distance

   Inputs:
     turbineSites: column vector of Nwt turbine sites such that
         * turbineSite(1:Nwt) are the x positions of the turbines
         * turbineSite(Nwt+1:2*Nwt) are the y positions of the turbines
     eg: Given that there are 30 turbines, the first turbine is located at
         the (x,y) coordinate at (turbineSite(1),turbineSite(31))
     farmWidth: width of the square wind farm
     figTitle: Optional parameter, if unspecified this function sets an
     empty figure title, otherwise displays the title on the figure
%}
function [] = plotTurbineScheme( turbineSites, farmWidth, showWake, ...
                                 U, includeTurbineRadius, figTitle)

    %caller of this function does not have to supply a title to the figure;
    %in the case of no caller-supplied title, make the figure title empty
    if nargin < 3
        showWake=0;
        figTitle='';
    end
    if nargin < 4 %first three prameters defined: draw all wakes with a fixed length 
        U=[];
        figTitle='';
    end
    if nargin < 5 %user as supplied U, so the wake lengths are going to be variable
        allWakesFixedLength=false;
        includeTurbineRadius=false;
        figTitle='';
    end
    if nargin < 6
        figTitle='';
    end
    
    Nwt = numel(turbineSites)/2;
    figure;
    
    if showWake==1;
        hold on;
    end
    
    %First we draw the farm and the turbine positions
    indexNumber = 1:Nwt;
    strValues = strtrim(cellstr(num2str(indexNumber(:),'(%d)')));
    scatter( turbineSites(1:Nwt), turbineSites(Nwt+1:2*Nwt) );
    text(turbineSites(1:Nwt),turbineSites(Nwt+1:2*Nwt),strValues,'VerticalAlignment','bottom');
    grid on;
    xlim([0,farmWidth]);ylim([0,farmWidth]);
    title(figTitle);

    %Next, if desired, we draw the wake for each turbine
    if showWake == 1        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for the Vestas V80-1800 wind turbine for GRADY CASE I
        h = 60;                %hub height
        Tdiam = 40;            %diameter of turbine blades
        Trad = Tdiam/2;        %radius of turbine blades
        z0 = 0.3;              %roughness constant
        alpha = 0.5/log(h/z0); %wake spreading constant
        Ct = 0.88;             %turbine thrust coefficient
        a=(1-sqrt(1-Ct))/2;    %for calculation of D on next line; a is 'axial force'
        D = Tdiam*sqrt((1-a)/(1-2*a)); %downstream rotor diameter
        gamma = (1-sqrt(1-Ct));  %for use in wakeStrength
        
        %Parameters for drawing the wake: density and length
        WakeStep = 50;  %determines how fine we draw the wake
        %MaxWakeLength = 60*Tdiam; %Dr. Wang's code had 4*cellWidth, where cellWidth is 5*Tdiam
        TOL = 10^-1;
        %in this case, the user wants to see the wake, but did not specify
        %the onset winds to the turbines.  Hence all wakes will have a
        %fixed length
        if(any(U)==0)
            allWakesFixedLength = true;
        else
            allWakesFixedLength = false;
        end
        MaxWakeLength = 60*Tdiam;
        
        for n=1:Nwt  %iterate through each turbine
            xPos = turbineSites(n);
            yPos = turbineSites(n+Nwt);

            %For this turbine, Draw the wake as a family of straight lines, and
            %color them based on wakeStrength: the lighter the color, the weaker
            %the wake
            for downWindDist = 0:WakeStep:MaxWakeLength
                if includeTurbineRadius
                    wakeRad = Trad + alpha*downWindDist; %this draws the turbine radius
                else
                    wakeRad = alpha*downWindDist; %draws wake like Dr. Wang; without turbine radius
                end
                
                %this is the wake strength; also known as the velocity deficit:
                %if you multiplied by u(i) you would have uTilde
                velocityDeficit = gamma*(D/(D+2*alpha*downWindDist))^2;
                wakeStrength = 1 - velocityDeficit;
                
                %if we are drawing variable length wakes, stop drawing the
                %wake if the speed of the wake is close to the onset wind
                %to that turbine i.e.:the wake has returned to the original
                %speed hitting that turbine (so stop drawing the wake; it's
                %confusing)
                if(~allWakesFixedLength && abs(U(n)-U(n)*wakeStrength)<TOL)
                %if(~allWakesFixedLength && wakeStrength)
                    disp(strcat(num2str(n),num2str(downWindDist)));
                    break;
                end
                
                plot( [xPos-wakeRad xPos+wakeRad], ...  %horizontal line from xPos-wakeRad to xPos+wakeRad
                      (yPos-downWindDist).*[1 1], ...   %location on y-axis to draw the horizonal line
                      'color',[0 1 1].*wakeStrength^3, ...  %the weaker the wake, the more baby-blue the line
                      'lineWidth',0.5)
            end
            %now draw the min safe distance about this turbine
            drawCircle(xPos,yPos, 2*Tdiam);
        end
        hold off;   %all done drawing the wake
    end
    axis square;
end

%used to draw the minimum safe distance about each turbine
function drawCircle(centX,centY,circRad)
    tickStep=0.1;
    thetaAngle=0:tickStep:2*pi; 
    xPt=circRad*cos(thetaAngle);
    yPt=circRad*sin(thetaAngle);
    plot(centX+xPt,centY+yPt, 'b');
end
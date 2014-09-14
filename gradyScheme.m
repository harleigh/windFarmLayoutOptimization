%{
   A nice function to return the turbine positions from various Grady
   Papers.

   Inputs:
     gradyCase: Integer; The grady case 1 2 or 3 from their papers.
     ResultFrom: 'GradyPaper' for results from Grady, otherwise the
                 schemes from Serrano Gonzalezs are returned

   Output:
     turbineSites: column vector of turbine sites from the given paper
         * turbineSite(1:Nwt) are the x positions of the turbines
         * turbineSite(Nwt+1:2*Nwt) are the y positions of the turbines
     where Nwt is the number of turbines in the specified Grady Scheme.

%}
function [ turbineSites ] = gradyScheme( gradyCase, ResultFrom )

    %%%%%%%%%-------parameters from the grady papers: 
    h = 60;
    D = 40;
    z0 = 0.3;
    C_T = 0.88;
    farmWidth = 1800;
    CellWidth = 5*D;
    %%%%%%%%%-------end parameters
    
    switch gradyCase
        case 1
            % The results from Gradys paper (Renewable Energy 2005)
            PositionMatrix = [ones(1,10); zeros(4,10); ones(1,10); ...
                zeros(3,10); ones(1,10)];
        case 2
            switch ResultFrom
                case  'GradyPaper'
                    % The results from Gradys paper (Renewable Energy 2005)
                    PositionMatrix = [1 0 1 1 1 1 0 1 0 1
                        1 0 0 0 0 0 0 1 0 0
                        1 0 0 1 0 0 1 0 1 0
                        0 0 1 0 0 1 0 0 0 1
                        1 0 0 1 0 0 1 1 0 0
                        0 1 0 0 1 0 0 0 1 0
                        0 0 1 0 0 1 0 1 0 1
                        1 0 0 1 0 0 1 0 0 0
                        0 0 0 0 1 0 0 1 0 1
                        1 1 0 1 0 1 0 1 0 1];
                otherwise
                    % The results from Serrano Gonzalezs paper (Renewable Energy 2010)
                    PositionMatrix = [1 0 1 1 1 0 1 0 1 1
                        1 0 0 0 0 0 1 0 0 0
                        0 0 0 1 0 0 0 0 0 1
                        1 1 0 0 0 1 0 1 0 1
                        0 0 1 0 1 0 0 0 0 1
                        1 0 0 0 0 0 1 0 1 0
                        0 1 0 1 0 1 0 0 0 1
                        0 1 0 1 0 0 0 1 0 1
                        1 0 0 0 1 0 0 0 0 1
                        1 0 1 0 1 0 1 1 0 1];
            end
        case 3
            switch ResultFrom
                case 'GradyPaper'
                    % The results from Gradys paper (Renewable Energy 2005)
                    PositionMatrix = [1 0 1 1 1 0 1 1 0 1
                        1 0 0 0 0 0 0 0 1 0
                        1 0 0 0 1 1 0 0 0 1
                        1 0 1 0 0 0 0 1 0 0
                        1 0 0 0 1 1 0 0 0 1
                        1 0 0 0 1 0 0 0 1 0
                        1 0 1 0 0 0 1 0 0 1
                        1 0 0 0 1 0 0 0 0 1
                        0 0 1 0 0 0 1 0 1 0
                        1 0 1 0 1 0 1 1 0 1];
                otherwise
                    % The results from Serrano Gonzalezs paper (Renewable Energy 2010)
                    PositionMatrix = [1 1 1 1 1 1 0 0 1 1
                        1 0 0 0 0 0 1 1 0 0
                        1 0 0 0 1 0 0 0 0 1
                        1 0 0 0 0 0 0 1 0 1
                        1 0 0 1 0 1 0 0 0 0
                        1 0 0 0 0 0 0 0 1 1
                        1 0 0 0 1 0 1 0 0 0
                        0 1 0 1 0 0 0 0 1 1
                        1 0 0 0 0 1 0 0 0 0
                        1 0 1 0 1 1 0 1 1 1];
            end
        otherwise
            disp('------Error: Wrong gradyCase parameter-----');
            turbineSites=[];
            return;
    end

    [row, col] = size(PositionMatrix);
    l = 0;
    PosX = []; PosY = [];
    for i = 1:row
        for j = 1:col
            if PositionMatrix(i,j) == 1
                PosX = [PosX (j-1)*CellWidth];
                PosY = [PosY farmWidth-(i-1)*CellWidth];
            end
        end
    end
    
    turbineSites = [PosX PosY]';

end


%{
   Assuming the wind farm is a square, given a width and number of desired
   turbines, return an empirical turbine layot (as a column vector; see
   output note below)

   Input:
     farmWidth: width of the square farm
     Nwt: number of turbines to be on the site
     
   Output:
     turbineSite: column vector of Nwt turbine sites such that
         * turbineSite(1:Nwt) are the x positions of the turbines
         * turbineSite(Nwt+1:2*Nwt) are the y positions of the turbines
  
   Example:
     EmpiricalScheme(1800,30,1)  returns an empirical layout of turbines on
     a square of length 1800 with 30 turbines on it.  Also, the turbine
     positions are plotted on a figure :)

   Notes:
     This function gives great initial guesses, and is a great test bed for
     testing partial shadowing
%}
function turbineSite = empiricalScheme(farmWidth,Nwt)
    switch Nwt
        case 1
            PosX = farmWidth/2; PosY = farmWidth/2;
        case 2
            PosX = [0 farmWidth]; PosY = [0 farmWidth];
        case 3
            PosX = [0 farmWidth/2 farmWidth]; PosY = [0 farmWidth 0];
        case 4
            PosX = [0 0 farmWidth farmWidth]; PosY = [0 farmWidth 0 farmWidth];
        case 5
            PosX = [0 0 farmWidth/2 farmWidth farmWidth]; PosY = [0 farmWidth farmWidth/2 0 farmWidth];
        otherwise
            ColumnNum =round(sqrt(Nwt));
            RowNum = round(Nwt/ColumnNum);
            FirstRow = round((Nwt - ColumnNum*(RowNum-2))/2);
            LastRow = Nwt - ColumnNum*(RowNum-2) - FirstRow;
            %fprintf('%d \t  %d \t  %d \t  %d \t  %d \n',n,RowNum,ColumnNum,FirstRow,LastRow);
            PosX = (0:1:FirstRow-1) * (farmWidth/(FirstRow-1));
            PosY = zeros(1,FirstRow);
            if RowNum >=3
                for i = 2: (RowNum-1)
                    for j = 1:ColumnNum
                        if mod(i,2) == 0
                            PosX = [PosX (2*j-1)*farmWidth/(2*ColumnNum-1)];
                        else
                            PosX = [PosX 2*(j-1)*farmWidth/(2*ColumnNum-1)];
                        end
                        PosY = [PosY (i-1)*farmWidth/(RowNum-1)];
                    end
                end
            end
            LastRowPosX = (0:1:LastRow-1) * (farmWidth/(LastRow-1));
            LastRowPosY = farmWidth*ones(1,LastRow);
            PosX = [PosX LastRowPosX];
            PosY = [PosY LastRowPosY];
    end
    
    if size(PosX,2)~= Nwt
        disp('Error')
        Nwt
    end
    
    turbineSite = [PosX PosY]';

end

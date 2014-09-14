%{
    Given the indexed turbine positions on the farm we calculate the
    down wind distance between each turbine.

        (i,j)th entry of dwD is:
             = 0         if i>=j  Since turbine i is downwind of turbine j
                                  because of the indexing of the turbines
                                  and orientation of the grid.
             = y(i)-y(j) if i<j   By the indexing of the turbines and the
                                  orientation of the grid, we have that
                                  y(i) > y(j)

    Inputs:
         yPos: Column Vector, Nwt by 1 of y-Positions of the indexed
               turbine locations

    Assumptions:
       The turbines must be indexed such that 
           yPos(1) >= yPos(2) >= ... >= yPos(Nwt)

   Output: 
       dwD: An upper triangular matrix with zeros along the diagonal,
            where the (i,j)th entry of dwD represents the downstream
            distance between turbine i and turbine j. Also, ever element of
            dwD is >=0
%}
function [ dwD ] = calcDownWindDist( yPos )

    global Nwt
    dwD = triu(repmat(yPos, 1, Nwt) - repmat(yPos', Nwt,1),1);
     
end
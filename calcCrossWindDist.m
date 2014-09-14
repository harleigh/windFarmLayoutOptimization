%{
    Given the indexed turbine positions on the farm we calculate the
    cross wind distance between each turbine.

        (i,j)th entry of cwD is:
             = 0            if i=j
             = |x(i)-x(j)|  otherwise

    By our implementation of the code cwD is calculated as   
 
        (i,j)th entry of cwD is:
             = 0            if i>=j these values are not used in the code
             = |x(i)-x(j)|  otherwise

    Inputs:
         xPos: Column Vector, Nwt by 1 of x-Positions of the indexed
               turbine locations

   Notes: cwD is chopped into an upper triangular matrix; we do not need
          the lower portion of cwD in any computations.

   Output: 
       cwD: An upper triangular matrix with zeros along the diagonal,
            where the (i,j)th entry of cwD represents the cross wind
            distance between turbine i and turbine j.

examples:
   repmat([1 2 3]',1, 3)     
       ans =
             1     1     1
             2     2     2
             3     3     3

    repmat([1 2 3],3, 1)
        ans =
             1     2     3
             1     2     3
             1     2     3
%}
function [ cwD ] = calcCrossWindDist( xPos )

    global Nwt
    cwD = triu( repmat(xPos,1,Nwt) - repmat(xPos',Nwt,1) ,1 );
    cwD = abs(cwD);
end
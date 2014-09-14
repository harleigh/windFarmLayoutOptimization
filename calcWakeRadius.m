%{
   Calculates the radius of the wake produced by turbine i (Ti) at the
   location of turbine j (Tj)

   wakeRadius(i,j) is (see note below also)
      = 0 if i>=j since, due to indexing the turbines by downwind distance,
                  and by the grid orientation, we know that if i>j then Ti
                  is downwind of Tj.  Hence the wake at Ti does not act on
                  Tj. And, for i=j, we just set the value to zero, since we
                  are looking at the same turbine.
      = turbineRadius + (wakeSpreadingConstant)*dwdDij if 1<=i<j<=Nwt

   Note: Consider the case where Ti and Tj lie exactly on the same line.
                 eg:
                        Ti        Tj
                         *        *
         Hence dwD(i,j)=0.  Thus, we have that
         wakeRadius(i,j)= turbineRadius which makes sense since Ti and Tj
         are at the same level.  

   Inputs:
      dwD: Matrix representing the down wind distance between
           each turbine
      turbRad: constant parameter; the radius of the turbine
      alpha:   constant parameter; wake spreading constant

   Outputs: 
     wakeRadius: where the (i,j)th entry of wakeRadius represents the
                 radius of the wake of Ti at the location of Tj.  This
                 matrix is upper triangular with zeros along the diagonal.

   Where is this used?
     This matrix allows us to build the boolean matrix inWake, and we need
     wakeRadius to calculate the velocity deficit matrix Q.

%}
function [ wakeRadius ] = calcWakeRadius(dwD, turbRad, alpha )

    wakeRadius = turbRad*triu(ones(size(dwD)),1) + alpha*dwD;
end
%{
  old code: Whenever dwDij=0, we forced wakeRadius_ij to be zero.

       wakeRadius = turbRad*(dwD>0) + alpha*dwD;
%}

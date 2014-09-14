%{
   This function calculates the onset wind speed at each turbine
   in the wind farm, taking into account the wake effect.

   Let U(i) be the onset wind speed of Ti.  Now, wue to the indexing of the
   wind turbines, and the choice of the grid for the wind farm:
                   (0,L)    (L,L)
                     *-------*
                     |       |
                     |       |
                     *-------*
                   (0,0)    (0,L)
   We have had symmetry all throughout the problem.  Indeed we also have 
   that the calculation of the oncet wind speed simplifies to
     U(1)=u0, and for k=2,...,Nwt we have
     U(k)= u0[1-sqrt(sum i=1 to k-1 of (Q(i,k)[1- tildeU(i,k)/u0])^2)]

   Where tidleU(i,k) is the wake of Ti acting on Tk, and is (broken into
   two lines for readability on an 80 character monitor):
      tidleU(i,k) = U(i)*(1-(1-sqrt(1-Ct))*X )
      X= turbineDiameter/(turbineDiameter+ 2*alpha*dwD(i,k))

   Special Notes:
     Consider Beta(i,k)= [1- tildeU(i,k)/u0], then we will see how our
     implementation below uses the two norm. Now, explicitly, for the first
     say, five turbines, their onset wind speeds are
        U(1)=u0
        U(2)=u0(1-sqrt([Q(1,2)*Beta(1,2)]^2))
        U(3)=u0(1-sqrt([Q(1,3)*Beta(1,3)]^2+[Q(2,3)*Beta(2,3)]^2))
        U(4)=u0(1-sqrt([Q(1,4)*Beta(1,4)]^2+[Q(2,4)*Beta(2,4)]^2 +
                [Q(3,4)*Beta(3,4)]^2)
        U(5)=u0(1-sqrt([Q(1,5)*Beta(1,5)]^2+[Q(2,5)*Beta(2,5)]^2 +
                [Q(3,5)*Beta(3,5)]^2 + [Q(4,5)*Beta(4,5)]^2)
     Now, let's consider Q[k] a column vector from Q as 
     [ Q(1,k) Q(2,k) ... Q(k-1,k) ]' and Beta[k] as a column vector whose
     entries are [ Beta(1,k) Beta(2,k) ... Beta(k-1,k) ]'. By doing so, and
     considering '.*' as component-wise multiplication two vectors (i.e.:
     Q[k].*Beta[k] = [ Q(1,k)Beta(1,k) ... Q(k-1,k)Beta(k-1,k) ]') we now
     have that
          U(k) = u0(1-TwoNorm(Q[k].*Beta[k]))
     
   Inputs:
     Q:  Matrix whose (i,j)th entrie represents the velocity deficit due
         to Ti shadowing Tj
     dwD: Matrix whose (i,j)th entrie represents the downwind distance from
          Ti to Tj
     u0: Freedom wind speed.  Currently a scalar, but will be a vector in
         future releases
     D: Parameter; downstream rotor diameter
     Alpha: Parameter; wake spreading constant in the linear wake model
     Ct: Parameter; turbine thrust coefficient

   Outputs:
     * vector U (size Nwt by 1) where U(i) is the wind speed acting on Ti
       due to the initial wind speed (freedom wind speed u0) as well as the
       (linear) wake effect.
%}
function [ U ] = calcOnsetWind( Q, dwD, u0, D, alpha, Ct )


    Nwt = size(Q,1);
    U = zeros(Nwt,1);
    gamma = (1-sqrt(1-Ct)); %this is for the calculation of tildeU

    %by the way we indexed the turbines, we know that the onset wind speed
    %for T1 is u0, the free wind speed.
    U(1)=u0;
    
    %Build U; iterate over columns of Q, taking entries upto but not 
    %including the main diagonal; we start at k=2 as k=1 was done
    %above (i.e.: U(1)=u0)
    for k=2:Nwt
                            %everything past '(1/u0)*' is tildeU(i,k)
        Beta=ones(k-1,1)-(1/u0)*...
            (U(1:k-1,1).*(ones(k-1,1)-gamma*((D./(D+(2*alpha)*dwD(1:k-1,k))).^2)));
        %U(k) = u0*(1-norm(Q(1:k-1,k).*B));% Q should not be squared
        U(k) = u0*(1-sqrt(sum(Q(1:k-1,k).*Beta.*Beta)));
    end
   
end

%{
The following is equivalent to the above; here we use the formula where we
are not using (1-uTilde/u0) but rather (u0-uTilde)

Also, .*Beta.*Beta is replaced by Beta.^2

 B=u0*ones(k-1,1)-(U(1:k-1,1).*(ones(k-1,1)-gamma*((D./(D+(2*alpha)*dwD(1:k-1,k))).^2)));
 U(k) = u0-sqrt(sum(Q(1:k-1,k).*(B.^2)));
%}

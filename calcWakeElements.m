%{
   tidleU(i,k) is the wake of Ti acting on Tk, and is (broken into
   two lines for readability on an 80 character monitor):
      tidleU(i,k) = U(i)*(1-(1-sqrt(1-Ct))*X )
      X= turbineDiameter/(turbineDiameter+ 2*alpha*dwD(i,k))

%}
function [ Weighted_U_tilde,  ...
           Weighted_U_tilde_squared, ...
           D_Utilde_D_U ] = calcWakeElements( U, Q, dwD )

    global Nwt alpha gamma D u0
    
    Weighted_U_tilde = zeros(Nwt,Nwt);
    Weighted_U_tilde_squared = zeros(Nwt,Nwt);
    D_Utilde_D_U = zeros(Nwt,Nwt);%store partial derivative of \tilde U_{i,k} w.r.t. U_i

    for k=1:Nwt-1
        %all upper triangular, build by row
        D_Utilde_D_U(k,k+1:end) =  1-gamma*( ( D./(D+(2*alpha)*dwD(k,k+1:end)) ).^2);
        U_tilde = U(k)*D_Utilde_D_U(k,k+1:end);
        Weighted_U_tilde(k,k+1:end) = (u0-U_tilde).*Q(k,k+1:end);
        Weighted_U_tilde_squared(k,k+1:end) = (u0-U_tilde).^2.*Q(k,k+1:end);
        %Weighted_U_tilde_square(k,k+1:end) = Weighted_U_tilde(k,k+1:end).*(u0-U_tilde);
        % less accurate using Weighted_U_tilde to compute
        % Weight_U_tilde_dquare
        %Weighted_U_tilde_square(k,k+1:end) = ( (u0-U_tilde).*Q(k,k+1:end) ).^2;
        % q should not be squared. Typo in Wang's paper?
    end

end


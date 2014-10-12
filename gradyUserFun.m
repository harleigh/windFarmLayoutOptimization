%{
   Input
    x:  x contains all of the decision variables.  Our decision variables
        are ordered
            x(1:Nwt)   xPositions of the turbines
            x(Nwt+1:2*Nwt)  yPositions of the turbines
        Now for each of the Nwd wind directions, we have one onset wind
        constraint U_w (and each U_w has Nwd elements).  Hence
            x(2*Nwt+1:3*Nwt) is the first wind diretion (eg:from 0 degrees)
            x(3*Nwt+1:4*Nwt) is the second wind diretion (eg:from 30 degrees)
            ...
            x((1+w)*Nwt+1):(2+w)*Nwt) is the w-th wind diretion where w=1,...,Nwd
        

turbineSites(I,:) == indexedTurbines
indexedTurbines(rI,:) == turbineSites
Notes:
    reverseindex can be used to map the sorted turbines back to the original
    index. It may be useful for the computation of the jacobian.
    indexedTurbines = turbineProfiles(index,:);
    turbineProfiles = indexedTurbines(reverseindex,:);
%}
function [ F, G ] = gradyUserFun( x ) %x is a column vector

    global index4relD index4G scaleRelD scaleConstraints4U CostModification;
    global alpha gamma Nwt Trad D powerCurve powerCurveDomain;
    global movements;
    global u0 windDirections Nwd Nws windDistribution;
    
    %column k is the onset wind speed for the turbines in wind direction k
    U = zeros(Nwt,Nwd*Nws);
    Constraints_4_U = zeros(Nwt, Nwd*Nws);
    numRelDistConstraints = numel(index4relD);
    lengthF = 1 + numRelDistConstraints + Nwd*Nwt*Nws;
    fullG = NaN(lengthF,(Nwd*Nws+2)*Nwt);
    expectedPower=0;
    
    
    gradRelDwrtX = zeros(numRelDistConstraints,Nwt);
    gradRelDwrtY = zeros(numRelDistConstraints,Nwt);
    
    XLOC = 1;  %enumeration: location of xPositions in turbineSites is col 1
    YLOC = 2;  %enumeration: location of yPositions in turbineSites is col 2

    %%%%%%%----begin cost
    for w=1:Nwd

        %rotate the turbine a-y positions into the wind direction w
        rotatedTurbPos = [x(1:Nwt), x(Nwt+1:2*Nwt)]*getCWRotationMatrix(windDirections(w));

        %To keep the ordering of U (with respect to relative indexing, we
        %removed the onset wind from the turbine profiles, for by sorting,
        %no longer U(1)==u0 for each wind direction
        turbineProfiles = [ rotatedTurbPos(:,XLOC), rotatedTurbPos(:,YLOC)];

        %Now that the turbines are rotated into the wind direction, we index
        %them according to down-wind and cross wind.  Due to the chosen grid
        %orientation (square: (0,L) (L,L) (0,0) (L,O)--upLeft, upRight, lowLeft,
        %lowRight respectivly) we sort the turbines in Decending order wrt yPos
        %and then acending order wrt xPosition
        [indexedTurbines,index] = sortrows(turbineProfiles,[-YLOC,XLOC]);
        [~,reverseIndex]=sortrows(index,1);


        %Critticaly important: dwD and cwD are ordered w.r.t. the sorted turbines
        dwD = calcDownWindDist(indexedTurbines(:,YLOC));
        cwD = calcCrossWindDist(indexedTurbines(:,XLOC));

        wakeRadius = calcWakeRadius(dwD, Trad, alpha );
        inWake = calcInWake( cwD, wakeRadius, Trad );
        [Q,~] = calcVelocityDeficit( cwD, dwD, inWake, wakeRadius, Trad, CostModification);

        %Critticaly important: U is ordered with respect to the indexed
        %turbines.
        %Here we grab the onset wind for the current wind direction w from
        %x and using logical indexing to store it into U
        U((w-1)*Nwt*Nws+1:w*Nwt*Nws) = x((2+(w-1)*Nws)*Nwt+1 : (2+w*Nws)*Nwt);
        
        %Weighted_U_tilde is upper triangular with diagonal entries be zeros
        %Note, calcWake needs the indexed version of onset wind
%         [ Weighted_U_tilde, ...
%           Weighted_U_tilde_squared, ...
%           D_Utilde_D_U ] = calcWakeElements( U(:,(w-1)*Nws+1:w*Nws), Q, dwD ); 
      
        %the square root of the sum of the matrix Weighted_U_tilde_squared
        %is the expression under the square root in the formula for Ui in
        %the research papers.  Note that we transpose it to just make a
        %column vector, for the calculation of the constraints for onset
        %wind on the next line.
%         velocityLost = sqrt( sum(Weighted_U_tilde_squared,1)');

        velocityLost = zeros(Nwt,Nws);
        for k=2:Nwt
            
%             tA = Q(1:k-1,k);
%             tB = U(1:k-1,(w-1)*Nws+1:w*Nws);
%             tC = ones(k-1,1)-gamma*((D./(D+(2*alpha)*dwD(1:k-1,k))).^2);
%             tD = bsxfun(@times, tB, tC);
%             tE = bsxfun(@minus,u0, tD);
%             tF = bsxfun(@times, tA, tE.^2); %tA.*tE.^2; 
%             tH = sqrt(sum(tF,1));
%             tI(k,1:Nws) = tH;
            
            velocityLost(k,1:Nws) = sqrt(sum(bsxfun(@times, Q(1:k-1,k), bsxfun(@minus,u0, bsxfun(@times, U(1:k-1,(w-1)*Nws+1:w*Nws), ones(k-1,1)-gamma*((D./(D+(2*alpha)*dwD(1:k-1,k))).^2))).^2),1));

        end
        
        
        %Store the constraints for U along the current wind direction (it's
        %of size Nwt by Nws) In the end, constraints for U is Nwt by
        %Nwd*Nws
        %Constraints_4_U(:,w) = U(:,(w-1)*Nws+1:w*Nws) - bsxfun(@minus,u0,velocityLost); 
        Constraints_4_U(:,(w-1)*Nws+1:w*Nws) = U(:,(w-1)*Nws+1:w*Nws) - bsxfun(@minus,u0,velocityLost); 
        
        
        
        
        %Calculate the expected power for all turbines along this wind
        %direction.
        turbinePower = interp1(powerCurveDomain, powerCurve, U(:,(w-1)*Nws+1:w*Nws),'linear');
        expectedPower = expectedPower + sum(sum(bsxfun(@times, windDistribution(w,:),turbinePower)));
%{
        %Jacobian of Constraints_4_U w.r.t. U along the current direction,
        %build from Weighted_U_tilde, D_Utilde_D_U and velocity_lost which
        %are updated in calcWakeElements
        Grad_Cons_U = eye(Nwt);%lower triangular
        for k = 2:Nwt;
            %turn the column into a row:                          column                     column
            %Dvide by zero happens! Matlab gives NaN; careful when porting;
            %maybe add a check (currently covered below)
           Grad_Cons_U(k,1:k-1) = -(1/velocityLost(k))*(Weighted_U_tilde(1:k-1,k).*D_Utilde_D_U(1:k-1,k))';
        end
        Grad_Cons_U = scaleConstraints4U*Grad_Cons_U;
        %zero out the areas where Q is zero or NaN
        Grad_Cons_U((Q+eye(Nwt))'== 0) = 0; %zero out the areas where Q is zero
        %store the gradient of the onset wind constraint into the Jacobian
        fullG( (end-((Nwd-(w-1))*Nwt) +1):( (end-((Nwd-(w-1))*Nwt))+Nwt), ...
               (end-((Nwd-(w-1))*Nwt) +1):( (end-((Nwd-(w-1))*Nwt))+Nwt)) ...
                  = Grad_Cons_U(:,:); 
        %}
    %%%%%%%-----end cost
    end

    %{
       Postcondition:
         -- U is Nwt by Nwd matrix; col i is the onset wind in the ith
            direction
         -- Constraints_4_U is Nwt by Nwd matrix; col i is the constraints
            of the onset wind in the ith direction
    %}

    %The cost function: since we only consider wind at 12 m/s the power curve
    %                   is 0.3*U^3
% % %     expectedPower = sum(U(:).^3)*0.3;
    
    %Critticaly important: relative distance should be w.r.t. the original
    %index for the calculation of the gradient (as the sparsity was set to
    %the original index of the turbines)
    relDsquared = calcSquaredRelativeDist(cwD, dwD);
    relDsquared = triu( relDsquared(reverseIndex,reverseIndex)+...
                        relDsquared(reverseIndex,reverseIndex)', 1);
    
    %pack the constraints to the problem
    F = [ expectedPower;
          scaleRelD*relDsquared(index4relD); %we flatten the matrix to a column vec
          scaleConstraints4U*Constraints_4_U(:)]; %flatten Constraints_4_U matrix to a column vec
    
      
    %Jacobian with respect to cost function: sum(U(:).^3)*0.3;
% % %     fullG(1,1:2*Nwt) = 0;   %zero for x and y
% % %     fullG(1,2*Nwt+1:end) = 0.9*U(:).^2; %onset wind
    
    %Jacobian: relative distance constraints (x(i)-x(j))^2-(y(i)-y(j))^2
    absPos = [x(1:Nwt), x(Nwt+1:2*Nwt)];
    [sparsityRow, sparsityCol] = ind2sub([Nwt,Nwt],index4relD);
    for i=1:numRelDistConstraints
        gradRelDwrtX(i,sparsityRow(i))=2*(absPos(sparsityRow(i),XLOC) - absPos(sparsityCol(i),XLOC));
        gradRelDwrtX(i,sparsityCol(i))=-2*(absPos(sparsityRow(i),XLOC) - absPos(sparsityCol(i),XLOC));
        
        gradRelDwrtY(i,sparsityRow(i))=2*(absPos(sparsityRow(i),YLOC) - absPos(sparsityCol(i),YLOC));
        gradRelDwrtY(i,sparsityCol(i))=-2*(absPos(sparsityRow(i),YLOC) - absPos(sparsityCol(i),YLOC));
    end
    %store the wrt xPosition portion and then the y portion and lastly
    %onset wind portion
    fullG(2:(numRelDistConstraints+1), 1:Nwt) = scaleRelD*gradRelDwrtX;
    fullG(2:(numRelDistConstraints+1), (Nwt+1):2*Nwt) = scaleRelD*gradRelDwrtY;
    fullG(2:(numRelDistConstraints+1), 2*Nwt+1:end)=0;  %zero for onset wind

    
    G = fullG(index4G);

    
    %movements(:,end+1) = indexedTurbines(1:2*Nwt);%turbineSites(1:2*Nwt);
    
end

%
% end gradyUserFun.m
%
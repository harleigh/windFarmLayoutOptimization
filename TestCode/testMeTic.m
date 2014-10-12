Nwt = 30;
xPos = rand(Nwt,1);
NumberCalls = 50000;    %to simulate the number of calls to the nonlinear solver
runTimeBsxFun = 0;
runTimeRepMat = 0;

for i = 1:NumberCalls
    tic
    cwD = abs(triu(bsxfun(@minus, xPos,xPos'),1));
    runTimeBsxFun = runTimeBsxFun + toc;
end


for i = 1:NumberCalls
    tic
    cwD = abs(triu( repmat(xPos,1,Nwt) - repmat(xPos',Nwt,1) ,1 ));
    runTimeRepMat = runTimeRepMat + toc;
end
disp( strcat('Number of calls:', num2str(NumberCalls)) );
disp( strcat('Run time for bsxfun:', num2str(runTimeBsxFun/NumberCalls)) );
disp( strcat('Run time for repmat:', num2str(runTimeRepMat/NumberCalls)) );
clear all;close all;


WindDirectionAll = [0:30:330]';
WindSpeedAll = [4:1:30]; % Cut-in speed is 4 m/s
DirFre = [0.020, 0.044, 0.056, 0.076, 0.061, 0.053, 0.078, 0.083, 0.123, 0.158,0.167, 0.079]';
Scale_v = [6.4, 7.4, 8.8, 7.0, 7.3, 7.5, 8.0, 8.9, 9.8, 9.9, 9.3, 7.8]';
Shape_k = [1.79, 1.81, 1.65, 1.96, 1.83, 1.81, 1.89, 1.85, 1.96, 1.92, 1.93, 1.75]';
SpeedDistributionTable = zeros(size(DirFre,1),size(WindSpeedAll,2));
for i =1:size(DirFre,1)
    v = Scale_v(i); k = Shape_k(i);
    for j=1:size(WindSpeedAll,2)
        SpeedDistributionTable(i,j) = k/v*(WindSpeedAll(j)/v)^(k-1)...
            *exp(-(WindSpeedAll(j)/v)^k);
    end
end

DistributionTable = diag(DirFre) * SpeedDistributionTable;
save WindSpeedDistribution_WASP DistributionTable
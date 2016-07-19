

figure(3); clf; hold on; box on;

rLigandArray =[1:1:10 15 20 25 30:10:50];

POccluded = zeros(numel(rLigandArray),6);

for irLigand=1:numel(rLigandArray)
    
    M = dlmread(['../runs/varyKinaseSize4/run' num2str(rLigandArray(irLigand)) '.parse100']);
    
    NTCHECK = 1000;
    
    for iiSite=1:6
        POccluded(irLigand,iiSite) = mean(M(NTCHECK:end,10+iiSite)>1);
    end


end % loop through N

%%


figure(3);

subplot(2,1,1);
plot(0.3*rLigandArray,POccluded(:,1:2));

legend('Y72', 'Y83', 'Y111', 'Y123', 'Y142', 'Y153');

subplot(2,1,2);
plot([0 0.3*rLigandArray], [1 1; 1./(1-POccluded(:,1:2))]);
%set(gca, 'yscale', 'log')
set(gca, 'xtick', [0 3.2 3.5 5 10 15])


M = dlmread('../src/testOcclusion2');

%%

NTCHECK = 1;
reeEquil = M(NTCHECK:10:end,2);
rMEquil = M(NTCHECK:10:end,3);
rHEquil = M(NTCHECK:10:end,4);


figure(1); clf; hold on; box on;
ecdf(rHEquil)
ecdf(rMEquil)
ecdf(reeEquil)

%%

figure(2); clf; hold on; box on;

subplot(3,1,1); hold on; box on;
plot(M(:,10));

subplot(3,1,2); hold on; box on;
plot(M(:,2));

subplot(3,1,3); hold on; box on;
plot(M(:,11));

%%

figure(3); clf; hold on; box on;

POccluded = mean(M(NTCHECK:end,11)>1);

display(POccluded);
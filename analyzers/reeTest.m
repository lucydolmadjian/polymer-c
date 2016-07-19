

figure(1); clf; hold on; box on;

figure(2); clf; hold on; box on;


for N=3:100
    
    M = dlmread(['../runs/reeTest/testRee' num2str(N)]);

    NTCHECK = 1000;
    reeEquil = M(NTCHECK:10:end,2);

    figure(1);
    ecdf(reeEquil)
    
    
    reeMean = mean(reeEquil);
    reeRMS = sqrt(mean(reeEquil.^2));
    figure(2);
    plot(N,reeRMS, '+b');
    
end % loop through N

%%

figure(2); 
plot(1:100, sqrt(1:100), '-r');
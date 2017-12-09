%% Look at polymer location distributions for Electrostatics

clear all;
close all;

%% Initialize parameters

ZETA = 1;

if(ZETA)
    combinations=64;
    N=113;
    iSiteTotal = 6;
else
    combinations = 4;
    N = 57;
    iSiteTotal = 2;
end

NTCHECK = 200000;
NBINS = 3000;
KBT = 4.14;
BinSize = 2*N/NBINS;
bSiteTotal = 0;
model = 1;
save = 0;
    
folder = '~/Documents/pub/lclemens/polymer-c_runs';
%folder = '~/Documents/polymer-c_runs';

saveFolder = '~/Documents/polymer-c_data_and_figures/CD3Zeta_Electrostatics/20171207';
    
switch (model)
    
    case 1
        subfolder = '20171121CD3ZetaElectrostaticsSweep1/Dephos';
        filenamePrefix = 'CD3ZetaElectrostatics';
        
        PD = 10.^(-2:1:2);  % parabola depth
        PW = 10.^(-2:1:2);  % parabola width
        WK = 10.^(-2:1:2);  % wall parabola k
        ER = 10.^(-2:1:2);  % repulsive energy multiplier
        ZR = 10.^(-2:1:2);  % repulsive energy exponent multiplier

end
    
%% Read Files

fileDoesNotExist=0;
distributionData = zeros(length(PD),length(PW),iSiteTotal+1,NBINS);

for k = [1 3 5]
    for w = [1 3 5]
        for d = [1 3 5]
            for e = 1:length(ER)
                for z = 1:length(ZR)
                    
                    %% Open file and retrieve distribution data
                    switch (model)
                        case 1
                            filename = strcat(filenamePrefix,'.','PD','.',num2str(d),'PW','.',num2str(w),'WK','.',num2str(k),'ER','.',num2str(e),'ZR','.',num2str(z));      
                    end

                    if(exist(fullfile(folder,subfolder,filename))~=0)

                        M = dlmread(fullfile(folder,subfolder,filename));
                        if (ZETA)
                            ntTotal = M(1,12);
                        else
                            ntTotal = M(1,8);
                        end

                        for i=1:1:iSiteTotal+1
                            if (ZETA)
                                start = 7+9+4*iSiteTotal+3+bSiteTotal+(i-1)*NBINS+1;  
                            else
                                start = 3+9+4*iSiteTotal+3+bSiteTotal+(i-1)*NBINS+1;    
                            end
                            distributionData(d,w,i,:) = M(1,start:start+NBINS-1)./(ntTotal-NTCHECK);
                        end
                    else
                        fileDoesNotExist=fileDoesNotExist+1;
                        disp(filename);
                    end

                end
            end
        end
    end

    %% Plot histograms

    bins = 0:1:NBINS-1;
    xaxis = -N+bins.*BinSize;

    for d=1:1:3
        for w=1:1:3
            for i=1:iSiteTotal+1
                plotData = reshape(distributionData(d,w,i,:),[NBINS,1]);
                figure((k-1)*(iSiteTotal+1)+i); hold on;
                %subplot(length(PW),length(PD),d+length(PD)*(w-1));
                subplot(3,3,d+3*(w-1));
                plot(xaxis,plotData,'-r','LineWidth',1.2);
                xlim([-5,5]);
                xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5]);
                title1 = strcat('D: ', num2str(PD(d)), 'W: ', num2str(PW(w)));
                set(gcf,'units','centimeters','position',[1,4,60,46]);
                disp(i);
                disp((k-1)*(iSiteTotal+1)+i);
                xlabel1 = 'z-coordinate (Kuhn lengths)';
                ylabel1 = 'Probability';
                %maxData(j)=max([maxData(j);plotData(:)]);
    %         
    %             if(i==iSiteTotal+1)
    %                 title1 = strcat('Distribution of Polymer Tip');
    %             else
    %                 title1 = strcat('Distribution of Tyrosine: ', num2str(i));
    %             end

                xlabel(xlabel1,'FontName','Arial','FontSize',14);
                ylabel(ylabel1,'FontName','Arial','FontSize',14);
                title(title1,'FontName','Arial','FontSize',14);
            end
        end
    end

    %% Plot energy diagrams
    if (0)
        for d=1:length(PD)
            for w=1:length(PW)
                for i=2

                    plotData = reshape(distributionData(d,w,i,:),[NBINS,1]);
                    figure(((k-1)*(iSiteTotal+1)+i)*100); hold on;
                    subplot(length(PW),length(PD),d+length(PD)*(w-1));

                    plot(xaxis,KBT*log(plotData),'-r');
                    xlim([-20,20]);
                    ylabel('pNnm');
                    title(strcat('D: ', num2str(PD(d)), 'W: ', num2str(PW(w))));
                    set(gcf,'units','centimeters','position',[1,4,60,46]);

                end
            end
        end
    end

    %% Save Figures

    if (save)


        switch (model)
            case 1
                saveSubFolder = 'HarwallPiecewiseBasicsY';
        end

        for i=15:21
            figure(i); 
            saveName = strcat('DistributioniSite',num2str(i-14));
            saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
            saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
        end
        if (0)
        for i=11:17
            figure(i); 
            saveName = strcat('EnergyBarrieriSite',num2str(i-10));
            saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'fig');
            saveas(gcf,fullfile(saveFolder,saveSubFolder,saveName),'epsc');
        end
        end

    end
end


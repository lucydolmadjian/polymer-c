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
model = 5;
save = 1;
    
folder = '~/Documents/pub/lclemens/polymer-c_runs';
%folder = '~/Documents/polymer-c_runs';

saveFolder = '~/Documents/polymer-c_data_and_figures/CD3Zeta_ElectrostaticsExtended/Mini';
    
switch (model)
    
    case 1
        subfolder = 'Mar232017ElectrostaticsHardwallPiecewiseBasicsY';
        filenamePrefix = 'CD3ZetaElectrostaticsHardwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        depth = 10.^(-3:0.5:3);
        width = 10.^(-3:0.5:3); %ParabolaCenter
        wallParabolaK = 1;
        
    case 2
        subfolder = 'Mar232017ElectrostaticsHardwallPiecewiseYOnly';
        filenamePrefix = 'CD3ZetaElectrostaticsHardwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        depth = 10.^(-3:0.5:3);
        width = 10.^(-3:0.5:3); %ParabolaCenter
        wallParabolaK = 1;
        
    case 3
        subfolder = 'Mar232017CD3EpsilonElectrostaticsHardwallPiecewiseBasicsY';
        filenamePrefix = 'CD3EpsilonElectrostaticsHardwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        depth = 10.^(-3:0.5:3);
        width = 10.^(-3:0.5:3); %ParabolaCenter
        wallParabolaK = 1;
        
    case 4
        subfolder = 'Mar232017CD3EpsilonElectrostaticsHardwallPiecewiseYOnly';
        filenamePrefix = 'CD3EpsilonElectrostaticsHardwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        
        depth = 10.^(-3:0.5:3);
        width = 10.^(-3:0.5:3); %ParabolaCenter
        wallParabolaK = 1;
        
    case 5
        subfolder = 'Mar232017ElectrostaticsSoftwallPiecewiseBasicsY';
        filenamePrefix = 'CD3ZetaElectrostaticsSoftwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 10.^(0:1:5); %ParabolaDepth
        width = 10.^(0:1:5); %ParabolaCenter
        wallParabolaK = 10.^(-3:0.5:2);
        
    case 6
        subfolder = 'Mar232017ElectrostaticsSoftwallPiecewiseYOnly';
        filenamePrefix = 'CD3ZetaElectrostaticsSoftwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 10.^(0:1:5); %ParabolaDepth
        width = 10.^(0:1:5); %ParabolaCenter
        wallParabolaK = 10.^(-3:0.5:2);

    case 7
        subfolder = 'Mar232017CD3EpsilonElectrostaticsSoftwallPiecewiseBasicsY';
        filenamePrefix = 'CD3EpsilonElectrostaticsSoftwallPiecewiseBasicsY';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 10.^(0:1:5); %ParabolaDepth
        width = 10.^(0:1:5); %ParabolaCenter
        wallParabolaK = 10.^(-3:0.5:2);
        
    case 8
        subfolder = 'Mar232017CD3EpsilonElectrostaticsSoftwallPiecewiseYOnly';
        filenamePrefix = 'CD3EpsilonElectrostaticsSoftwallPiecewiseYOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 10.^(0:1:5); %ParabolaDepth
        width = 10.^(0:1:5); %ParabolaCenter  
        wallParabolaK = 10.^(-3:0.5:2);
        
      case 9
        subfolder = 'May022017CD3ZetaElectrostaticsSoftwallBasicsOnly';
        filenamePrefix = 'CD3ZetaElectrostaticsSoftwallBasicsOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 10.^(0:0.5:1.5); %ParabolaDepth
        %width = 1:1:10; %ParabolaCenter 
        width = 10.^(0:1:3); %ParabolaCenter 
        wallParabolaK = 10.^(-5:1:0);
        
      case 10
        subfolder = 'May022017CD3EpsilonElectrostaticsSoftwallBasicsOnly';
        filenamePrefix = 'CD3EpsilonElectrostaticsSoftwallBasicsOnly';
        depthName = 'ParabolaDepth';
        widthName = 'ParabolaWidth';
        wallParabolaKName = 'WallParabolaK';
        
        depth = 10.^(0:0.5:2); %ParabolaDepth
        width = 1:1:10; %ParabolaCenter  
        wallParabolaK = 10.^(-7:1:-5);

end
    



%% Read Files

for k = 3

distributionData = zeros(length(depth),length(width),iSiteTotal+1,NBINS);

fileDoesNotExist=0;
for w = 1:length(width)
    for d = 1:length(depth)
    
        
        switch (model)
            case 1
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w));
            case 2
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w));
            case 3
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w));
            case 4
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w));
            case 5
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w), '.', wallParabolaKName,'.', num2str(k));
            case 6
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w), '.', wallParabolaKName,'.', num2str(k));
            case 7
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w), '.', wallParabolaKName,'.', num2str(k));
            case 8
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w), '.', wallParabolaKName,'.', num2str(k));
            case 9
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w), '.', wallParabolaKName,'.', num2str(k));
            case 10
                filename = strcat(filenamePrefix,depthName,'.',num2str(d),'.',widthName,'.',num2str(w), '.', wallParabolaKName,'.', num2str(k));        
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

%% Plot histograms

bins = 0:1:NBINS-1;
xaxis = -N+bins.*BinSize;

            disp(k);
for d=1:1:3
    for w=1:1:3
        for i=1:iSiteTotal+1
            plotData = reshape(distributionData(d,w,i,:),[NBINS,1]);
            figure((k-1)*(iSiteTotal+1)+i); hold on;
            %subplot(length(width),length(depth),d+length(depth)*(w-1));
            subplot(3,3,d+3*(w-1));
            plot(xaxis,plotData,'-r','LineWidth',1.2);
            xlim([-5,5]);
            xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5]);
            title1 = strcat('D: ', num2str(depth(d)), 'W: ', num2str(width(w)));
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
for d=1:length(depth)
    for w=1:length(width)
        for i=2
            
            plotData = reshape(distributionData(d,w,i,:),[NBINS,1]);
            figure(((k-1)*(iSiteTotal+1)+i)*100); hold on;
            subplot(length(width),length(depth),d+length(depth)*(w-1));
            
            plot(xaxis,KBT*log(plotData),'-r');
            xlim([-20,20]);
            ylabel('pNnm');
            title(strcat('D: ', num2str(depth(d)), 'W: ', num2str(width(w))));
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
        case 2
            saveSubFolder = 'HardwallPiecewiseYOnly';
        case 3
            saveSubFolder = 'HardwallPiecewiseBasicsY';
        case 4
            saveSubFolder = 'HardwallPiecewiseYOnly';
        case 5
            saveSubFolder = 'SoftwallPiecewiseBasicsY';
        case 6
            saveSubFolder = 'SoftwallPiecewiseYOnly';
        case 7
            saveSubFolder = 'SoftwallPiecewiseBasicsY';
        case 8
            saveSubFolder = 'SoftwallPiecewiseYOnly';
        case 9
            saveSubFolder = 'BasicsOnly';
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


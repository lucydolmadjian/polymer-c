%% Analysis_EffectiveConcKernel

% Aka Figure 3 for Surface Effects 

%% Sigma(r) computations - Edited to remove parameter dependence
% Cell array version

clear all;
close all;

%%
NRODS = [25 64];
delta = 1;
rho = 0:1:50;
irLigand = [0 2 5 10];
cutoff = 1;
%stDevDex = (5.6/0.3)/(sqrt(3)); % divide by 0.3 to put in kuhn lengths, divide by sqrt(3) to get each individual dimension stdev
%set to zero to remove influence of dextran
stDevDex=0;
dataThin = 0;
NTCHECK = 200000;


save=1;
savefolder='~/Google Drive/polymer-c/polymer-c_data_and_figures/SurfaceEffects/EffectiveConcentrationKernel/UpdatedJune072017';

%%

for Nrods = NRODS
    folder = '~/Documents/polymer-c_runs/';
    
    subfolderprefix = 'SPROcclusion';
    
    %% Read in data

    % read in data
    dataFiles = SphereOcclusionReadFilesSHP1(Nrods,irLigand,folder,subfolderprefix,dataThin);

    % if N==30, average data
    if (0)
        dataFiles2 = SphereOcclusionReadFilesSHP1(Nrods,irLigand,folder,subfolderprefix,dataThin);
    end


    ligandBoundrData = dataFiles{1};
    ligandOffrData = dataFiles{2};
    ligandBoundcenterData = dataFiles{3};
    ligandOffcenterData = dataFiles{4};
    occlusionData = dataFiles{5};

    ratio = zeros(length(irLigand),length(rho));

    %% Shift data for ligandCenter2 by rho
    for membrane=0:1:1
            clear nearEnough ableToBind sigmaError BernoulliError;

            colors = parula(5);
            
        for irL = 1:length(irLigand)
            clear ligandCenter1 ligandCenter2 occlusion sigma fraction;
            
            ligandCenter1 = ligandBoundcenterData{1,irL,membrane+1};
            occlusion = occlusionData{1,irL,membrane+1};
            
            availability = ~occlusion;

                for r = 1:size(rho,2)
                    clear ligandCenter2shifted distance distanceClose bindingAvailability;

                    % decide if want same dextran shift per rho for
                    % consistency, or different

                    ligandCenter2shifted = ligandOffcenterData{1,irL,membrane+1};
                    ligandOffrDatashifted = ligandOffrData{1,irL,membrane+1};
                    
                    
                    % shift data by rho(r)
                    ligandCenter2shifted(:,1) = ligandCenter2shifted(:,1)+rho(r);
                    ligandOffrDatashifted(:,1) = ligandCenterrDatashifted(:,1)+rho(r);

                    % create random matrix of Dextran shifts - only use if in
                    % free-space
                    if (~membrane)
                        dextranShift1 = normrnd(0,stDevDex,[size(ligandCenter1,1),size(ligandCenter1,2)]);
                        dextranShift2 = normrnd(0,stDevDex,[size(ligandCenter2shifted,1),size(ligandCenter2shifted,2)]);
                        ligandCenter1 = ligandCenter1+dextranShift1;
                        ligandCenter2shifted = ligandCenter2shifted+dextranShift2;
                    end

                    distance(:) = sqrt((ligandCenter1(:,1)-ligandCenter2shifted(:,1)).^2 + (ligandCenter1(:,2)-ligandCenter2shifted(:,2)).^2 + (ligandCenter1(:,3)-ligandCenter2shifted(:,3)).^2);

                    distanceClose = (distance<cutoff);

                    nearEnough(irL,r) = size(find(distanceClose==1),2);

                    %bindingAvailability = distanceClose.*occlusion;
                    bindingAvailability = distanceClose.*availability;
                    %bindingAvailability = distanceClose;
                    ableToBind(irL,r) = size(find(bindingAvailability),2);

                    fraction(r) = size(find(bindingAvailability),2)/size(bindingAvailability,2);

                    BernoulliError(irL,r) = sqrt(fraction(r)*(1-fraction(r))/size(bindingAvailability,2));

                    %fraction = size(find(distanceClose),2)/size(distanceClose,2);

                    sigma(r) = fraction(r)/(4/3*pi*(cutoff)^3);

                    sigmaError(irL,r) = BernoulliError(irL,r)/(4/3*pi*(cutoff)^3);
                end



            % calculate magnitude of concentration ehancement from sigma at r = 0
            % where does the 2*pi come from?  I expect it should be 4*pi*N?
            cN = (4*pi/3)^(3/2)*(1/sigma(1));

            cN2 = (4*pi*delta^2*Nrods/3)^(3/2)*sigma(1);

            rhoArray = linspace(0,max(rho),200);
            sigmaTheory = (3/(4*pi*delta^2*Nrods))^(3/2)*exp(-3*rhoArray.^2/(4*Nrods*delta^2));

            %% Plot of 'able to bind' probability (i.e. close enough and open)
            %figure(25+membrane); hold on; box on;
            figure(25); hold on; box on;
            subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
            h = plot(rho, fraction, '-','LineWidth',1.2,'Color',colors(irL,:));
            errColor = get(h,'Color');
            errLine = errorbar(rho,fraction,BernoulliError(irL,:),'Color',errColor);
            set(get(get(errLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            xlabel1 = 'Separation Distance, r';
            ylabel1 = 'Probability of Access';
            title1 = strcat('N = ',num2str(Nrods),'     Membrane = ', num2str(membrane));
            set(gca,'FontName','Arial','FontSize',18);
            xlabel(xlabel1,'FontName','Arial','FontSize',18);
            ylabel(ylabel1,'FontName','Arial','FontSize',18);
            title(title1,'FontName','Arial','FontSize',18);
            ylim([0 0.01]);

            %% Plot figure of Sigma(r) Nondimensionalized

            figure(1); hold on; box on;
            subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
            plot(rho/(delta*sqrt(Nrods)),sigma*delta^3*Nrods^(3/2),'-*','LineWidth',1.2,'Color',colors(irL,:));
            xlabel1 = 'Separation Distance, r/l';
            ylabel1 = 'Concentration, \sigma*l^3';
            title1 = strcat('N = ',num2str(Nrods),'     Membrane = ', num2str(membrane));
            title(title1,'FontName','Arial','FontSize',18);
            %pos = get(gcf, 'position');
            %set(gcf,'units','centimeters','position',[pos(1:2),8,5.5]);
            set(gca,'FontName','Arial','FontSize',18);
            xlabel(xlabel1,'FontName','Arial','FontSize',18);
            ylabel(ylabel1,'FontName','Arial','FontSize',18);
            Fig1Max(irL,find(NRODS==Nrods),membrane+1) = max(sigma*delta^3*Nrods^(3/2));

            %% Plot figure of Sigma(r) Nondimensionalized

            figure(3); hold on; box on;
            subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
            plot(rho,sigma,'-*','LineWidth',1.2,'Color',colors(irL,:));
            xlabel1 = 'Separation Distance, r';
            ylabel1 = 'Concentration, \sigma';
            title1 = strcat('N = ',num2str(Nrods),'     Membrane = ', num2str(membrane));
            title(title1,'FontName','Arial','FontSize',18);
            %pos = get(gcf, 'position');
            %set(gcf,'units','centimeters','position',[pos(1:2),8,5.5]);
            set(gca,'FontName','Arial','FontSize',18);
            xlabel(xlabel1,'FontName','Arial','FontSize',18);
            ylabel(ylabel1,'FontName','Arial','FontSize',18);
            Fig3Max(irL,find(NRODS==Nrods),membrane+1) = max(sigma);



            %% Plot, normalized by magnitude of concentration

            % figure(2); hold on; box on;
            % plot(rho/(delta*sqrt(Nrods)),sigma/cN*delta^3*Nrods^(3/2),'-*');
            % plot(rhoArray/(delta*sqrt(Nrods)), sigmaTheory*delta^3*Nrods^(3/2),'-');
            % xlabel1 = 'Separation Distance, r/l';
            % ylabel1 = 'Concentration Normalized by Enhancement,\sigma/c_N*l^3';
            % %pos = get(gcf, 'position');
            % %set(gcf,'units','centimeters','position',[pos(1:2),8,5.5]);
            % set(gca,'FontName','Arial','FontSize',16);
            % xlabel(xlabel1,'FontName','Arial','FontSize',16);
            % ylabel(ylabel1,'FontName','Arial','FontSize',16);

            figure(20); hold on; box on;
            subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
            plot(rho/(delta*sqrt(Nrods)),sigma/cN2*delta^3*Nrods^(3/2),'-*','Color',colors(irL,:));
            xlabel1 = 'Separation Distance, r/l';
            ylabel1 = 'Concentration Normalized by Enhancement, \sigma/c_N*l^3';
            title1 = strcat('N = ',num2str(Nrods),'     Membrane = ', num2str(membrane));
            title(title1,'FontName','Arial','FontSize',18);
            %pos = get(gcf, 'position');
            %set(gcf,'units','centimeters','position',[pos(1:2),8,5.5]);
            set(gca,'FontName','Arial','FontSize',18);
            xlabel(xlabel1,'FontName','Arial','FontSize',18);
            ylabel(ylabel1,'FontName','Arial','FontSize',18);
            Fig20Max(irL,find(NRODS==Nrods),membrane+1) = max(sigma/cN2*delta^3*Nrods^(3/2));
            

            %% Ratio of Sigma(r) nondimensionalized without membrane vs with membrane
                if (membrane)
                    ratioMem(irL,:) = sigma*delta^3*Nrods^(3/2);
                    sigmaErrorMem = sigmaError*delta^3*Nrods^(3/2);
                else 
                    ratioNoMem(irL,:) = sigma*delta^3*Nrods^(3/2);
                    sigmaErrorNoMem = sigmaError*delta^3*Nrods^(3/2);
                end


    %         disp(nearEnough);
    %         disp(ableToBind);
        end

        figure(1);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        plot(rhoArray/(delta*sqrt(Nrods)), sigmaTheory*delta^3*Nrods^(3/2),'k--','LineWidth',1.5);
        legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5','Ligand Radius = 10','WLC Theory');
        pos = get(gcf, 'position');
        set(gcf,'units','centimeters','position',[[1 1],40,30]);

        figure(3);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        plot(rhoArray, sigmaTheory,'k--','LineWidth',1.5);
        legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5','Ligand Radius = 10','WLC Theory');
        pos = get(gcf, 'position');
        set(gcf,'units','centimeters','position',[[1 1],40,30]);

        figure(25);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        %plot(rhoArray/(delta*sqrt(Nrods)), sigmaTheory*delta^3*Nrods^(3/2),'-','LineWidth',1.5);
        legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5','Ligand Radius = 10','WLC Theory');
        pos = get(gcf, 'position');
        set(gcf,'units','centimeters','position',[[1 1],40,30]);


        %figure(1);
        %legend('N=30','Theory 30','N=100','Theory 100','N=250','Theory 250');
        %title1 = 'Nondimensionalized Half-Space Data Against Free Space Theory Curve';
        %title(title1,'FontName','Arial','FontSize',12); 

        %figure(2);
        %legend('N=30','Theory 30','N=100','Theory 100','N=250','Theory 250');


        figure(20+membrane);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        %plot(rhoArray/(delta*sqrt(Nrods)), sigmaTheory*delta^3*Nrods^(3/2),'-','LineWidth',1.5);
        legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5','Ligand Radius = 10','WLC Theory');   
        pos = get(gcf, 'position');
        set(gcf,'units','centimeters','position',[[1 1],40,30]);
    end
    %% Plot, cN vs N

    for irL = 1:length(irLigand)
        ratioError(irL,:) = abs(ratioMem(irL,:)./ratioNoMem(irL,:)).*sqrt((sigmaErrorNoMem(irL,:)./ratioNoMem(irL,:)).^2+(sigmaErrorMem(irL,:)./ratioNoMem(irL,:)).^2);
    end

    %figure(3); clf;
    % 
    % figure(3); hold on; box on;
    % plot(Nrods,cN2,'-*');
    % xlabel1 = 'Number of Rods, N';
    % ylabel1 = 'Concentration Enhancement, c_N';
    % %pos = get(gcf, 'position');
    % %set(gcf,'units','centimeters','position',[pos(1:2),8,5.5]);
    % set(gca,'FontName','Arial','FontSize',14);
    % xlabel(xlabel1,'FontName','Arial','FontSize',14);
    % ylabel(ylabel1,'FontName','Arial','FontSize',14);
    % axis([0 250 0 2]);



    % 
    % 
    % %end
    % 
%     figure(40); clf; hold on; box on;
%     colors = parula(20);
%     for irL = 1:12
%         h = plot(rho/(delta*sqrt(Nrods)), ratioMem(irL,:)./ratioNoMem(irL,:),'-*','LineWidth',1.2,'Color',colors(irL,:));
%         errLine = errorbar(rho/(delta*sqrt(Nrods)),ratioMem(irL,:)./ratioNoMem(irL,:),ratioError(irL,:),'Color',colors(irL,:));
%         set(get(get(errLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     end
%     xlabel1 = 'Separation Distance, r/l';
%     ylabel1 = 'Fold Change in Concentration';
%     xlabel(xlabel1,'FontName','Arial','FontSize',18);
%     ylabel(ylabel1,'FontName','Arial','FontSize',18);
%     legend('Ligand Radius = 1','Ligand Radius = 2','Ligand Radius = 3','Ligand Radius = 4','Ligand Radius = 5','Ligand Radius = 6','Ligand Radius = 7','Ligand Radius = 8','Ligand Radius = 9','Ligand Radius = 10','Ligand Radius = 11','Ligand Radius = 12');
%     ylim([0 5]);



    % figure(31); hold on;
    % for irL = irLigand
    %     plot(rho/(delta*sqrt(Nrods)), ratioNoMem(irL,:),'-*');
    % end
    % xlabel1 = 'Separation Distance, r/l';
    % ylabel1 = 'Fold Change in Concentration';
    % xlabel(xlabel1,'FontName','Arial','FontSize',18);
    % ylabel(ylabel1,'FontName','Arial','FontSize',18);
    % legend('Ligand Radius = 9','Ligand Radius = 10','Ligand Radius = 11','Ligand Radius = 12','Ligand Radius = 13','Ligand Radius = 14','Ligand Radius = 15','Ligand Radius = 16','Ligand Radius = 17','Ligand Radius = 18','Ligand Radius = 19','Ligand Radius = 20','Ligand Radius = 21','Ligand Radius = 22','Ligand Radius = 23','Ligand Radius = 24','Ligand Radius = 25','Ligand Radius = 26','Ligand Radius = 27','Ligand Radius = 28','Ligand Radius = 29');



end

%% Set Maximums on Figures

for Nrods = NRODS
    for membrane=0:1:1
        figure(1);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        ylim([0 max(max(max(Fig1Max)))]);

        figure(3);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        ylim([0 max(max(max(Fig3Max)))]);

        figure(20+membrane);
        subplot(2,2,membrane+(2*(find(NRODS==Nrods)-1)+1)); hold on;
        ylim([0 max(max(max(Fig20Max)))]);
        
    end
end





%% Save Figures

if (save)
       figure(25);
            savefilename = 'ProbabilityCloseAndOpen';
            saveas(gcf,fullfile(savefolder,savefilename),'fig');
            saveas(gcf,fullfile(savefolder,savefilename),'epsc'); 

       
        figure(1);
            savefilename = 'ConcentrationNondimVSSeparation';
            saveas(gcf,fullfile(savefolder,savefilename),'fig');
            saveas(gcf,fullfile(savefolder,savefilename),'epsc');     

            
        figure(3);
            savefilename = 'ConcentrationVSSeparation';
            saveas(gcf,fullfile(savefolder,savefilename),'fig');
            saveas(gcf,fullfile(savefolder,savefilename),'epsc'); 


        figure(20);
            savefilename = 'NormalizedConcentrationVSSeparation';
            saveas(gcf,fullfile(savefolder,savefilename),'fig');
            saveas(gcf,fullfile(savefolder,savefilename),'epsc');     

end



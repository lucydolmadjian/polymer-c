function [plot1,plot2] = Funct_EffConcPlots(saveFilename,saveTF,membrane)
%function [] = Funct_EffConcPlots(saveFilename,saveTF)
    
    load(saveFilename);

    colors = parula(5);

    %% Plot figure of Sigma(r) Nondimensionalized
    figure(1); clf; hold on; box on;
    for n = 1:length(NRODS)
        %subplot(2,2,membrane+(2*(n-1)+1)); hold on;
        for irL = 1:length(irLigand)
            sigmaPlot = reshape(sigma(n,irL,membrane+1,:),length(rho),1);
            plot(rho/(delta*sqrt(NRODS(n))),sigmaPlot*delta^3*NRODS(n)^(3/2),'-*','LineWidth',1.2,'Color',colors(irL,:));
        end
    end
    fig = plot(rhoArray/(delta*sqrt(NRODS(n))), sigmaTheory*delta^3*NRODS(n)^(3/2),'k--','LineWidth',1.5);
        legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5','Ligand Radius = 10','WLC Theory');
        pos = get(gcf, 'position');
        set(gcf,'units','centimeters','position',[[1 1],40,30]);
        xlabel1 = 'Separation Distance, r/l';
        ylabel1 = 'Concentration, \sigma*l^3';
        title1 = strcat('N = ',num2str(NRODS(n)),'     Membrane = ', num2str(membrane));
        title(title1,'FontName','Arial','FontSize',18);
        set(gca,'FontName','Arial','FontSize',18);
        xlabel(xlabel1,'FontName','Arial','FontSize',18);
        ylabel(ylabel1,'FontName','Arial','FontSize',18);
        set(gcf,'toolbar','figure');
        %Fig1Max(irL,n,membrane+1) = max(sigma(n,irL,membrane+1,:)*delta^3*NRODS(n)^(3/2));
        %plot1 = ancestor(fig,'figure');
    %% Plot figure of Sigma(r)
    figure(2); clf; hold on; box on;
    for n = 1:length(NRODS) 
        %subplot(2,2,membrane+(2*(n-1)+1)); hold on;
        for irL = 1:length(irLigand)
            sigmaPlot = reshape(sigma(n,irL,membrane+1,:),length(rho),1);
            plot(rho,sigmaPlot,'-*','LineWidth',1.2,'Color',colors(irL,:));
        end
    end
    fig2 = plot(rhoArray, sigmaTheory,'k--','LineWidth',1.5);
        legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5',...
               'Ligand Radius = 10','WLC Theory');
        pos = get(gcf, 'position');
        set(gcf,'units','centimeters','position',[[1 1],40,30]);
        xlabel1 = 'Separation Distance, r';
        ylabel1 = 'Concentration, \sigma';
        title1 = strcat('N = ',num2str(NRODS(n)),'     Membrane = ', num2str(membrane));
        title(title1,'FontName','Arial','FontSize',18);
        set(gca,'FontName','Arial','FontSize',18);
        xlabel(xlabel1,'FontName','Arial','FontSize',18);
        ylabel(ylabel1,'FontName','Arial','FontSize',18);
        %Fig3Max(irL,n,membrane+1) = max(sigma(n,irL,membrane+1,:));
        %plot2 = ancestor(fig2,'figure');

    %% Plot, normalized by magnitude of concentration

    % figure(4); hold on; box on;
    % plot(rho/(delta*sqrt(NRODS(n))),sigma/cN*delta^3*NRODS(n)^(3/2),'-*');
    % plot(rhoArray/(delta*sqrt(NRODS(n))), sigmaTheory*delta^3*NRODS(n)^(3/2),'-');
    % xlabel1 = 'Separation Distance, r/l';
    % ylabel1 = 'Concentration Normalized by Enhancement,\sigma/c_N*l^3';
    % %pos = get(gcf, 'position');
    % %set(gcf,'units','centimeters','position',[pos(1:2),8,5.5]);
    % set(gca,'FontName','Arial','FontSize',16);
    % xlabel(xlabel1,'FontName','Arial','FontSize',16);
    % ylabel(ylabel1,'FontName','Arial','FontSize',16);

%     figure(3); hold on; box on;
%     for n = 1:length(NRODS)
%         subplot(2,2,membrane+(2*(n-1)+1)); hold on;
%         plot(rho/(delta*sqrt(NRODS(n))),sigma/cN*delta^3*NRODS(n)^(3/2),'-*','Color',colors(irL,:));
%     end
%         xlabel1 = 'Separation Distance, r/l';
%         ylabel1 = 'Concentration Normalized by Enhancement, \sigma/c_N*l^3';
%         title1 = strcat('N = ',num2str(NRODS(n)),'     Membrane = ', num2str(membrane));
%         title(title1,'FontName','Arial','FontSize',18);
%         set(gca,'FontName','Arial','FontSize',18);
%         xlabel(xlabel1,'FontName','Arial','FontSize',18);
%         ylabel(ylabel1,'FontName','Arial','FontSize',18);
%         Fig20Max(irL,n,membrane+1) = max(sigma/cN*delta^3*NRODS(n)^(3/2));
%         legend('Ligand Radius = 0','Ligand Radius = 1','Ligand Radius = 5',...
%                'Ligand Radius = 10','WLC Theory');   
%         pos = get(gcf, 'position');
%         set(gcf,'units','centimeters','position',[[1 1],40,30]);
        
%% Set Maximums on Figures

%for n = 1:length(NRODS)
    %for membrane=0:1:1
        %figure(1);
        %subplot(2,2,membrane+(2*(n-1)+1)); hold on;
        %ylim([0 max(max(max(Fig1Max)))]);

        figure(1);
        title('Edits');
        %subplot(2,2,membrane+(2*(n-1)+1)); hold on;
        %ylim([0 max(max(max(Fig3Max)))]);

%         figure(3);
%         subplot(2,2,membrane+(2*(n-1)+1)); hold on;
%         ylim([0 max(max(max(Fig20Max)))]);

    %end
%end

%% Save Figures

if (saveTF)

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

end
%% Function to read in files for SHP1

function output = SphereOcclusionReadFiles(Nrods,irLigand,folder,subfolderprefix,dataThin)



ligandBoundrData = cell(size(Nrods,2),size(irLigand,2),2);
ligandOffrData = cell(size(Nrods,2),size(irLigand,2),2);
ligandBoundcenterData = cell(size(Nrods,2),size(irLigand,2),2);
ligandOffcenterData = cell(size(Nrods,2),size(irLigand,2),2);
occlusionData = cell(size(Nrods,2),size(irLigand,2),2);

for membrane=0:1:1

    if (membrane)
        subfoldersuffix = 'HalfSpace';
    else
        subfoldersuffix = 'FreeSpace';
    end
    
    foldername1 = strcat(folder,'June072017SigmaRSweep',subfoldersuffix,'LigandOn');
    foldername2 = strcat(folder,'June072017SigmaRSweep',subfoldersuffix,'LigandOff');

    subfolder1 = strcat(subfolderprefix,subfoldersuffix, 'LigandOn');
    subfolder2 = strcat(subfolderprefix,subfoldersuffix, 'LigandOff');

    
    for n = 1:length(Nrods)
        for irL = 1:length(irLigand)

            clear M1 M2 r1 r2 ligandCenter1 ligandCenter2 occlusion;

            %clearvars -except Nrods delta rho irLigand cutoff folder subfolder1 subfolder2 membrane irL ratio ratioMem ratioNoMem nearEnough available ableToBind;


            % read in data
            filename1 = strcat(subfolder1,'.',num2str(Nrods(n)),'.',num2str(irLigand(irL)));
            filename2 = strcat(subfolder2,'.',num2str(Nrods(n)),'.',num2str(irLigand(irL)));

            M1 = dlmread(fullfile(foldername1,subfolder1,filename1));
            M2 = dlmread(fullfile(foldername2,subfolder2,filename2));

            % truncate larger one by removing initial transient 
            % (is this going to cause problems since kinase bound should be bigger file?  Better to force the
            %shorter ones to run longer?  
            if (size(M1,1) > size(M2,1))
                if (dataThin)
                    MScrap = M1;
                    cutRows = randi([1 size(M1,1)],size(M2,1),1);
                    M1 = MScrap(cutRows(:),:);

                    disp(size(M1));
                    disp(size(M2));
                else
                    MScrap = M2;
                    addRows = randi([1 size(M2,1)],abs(size(M1,1)-size(M2,1)),1);
                    M2 = [M2; MScrap(addRows(:),:);];
                    
                    % double the size of the data
                    addRows1 = randi([1 size(M1,1)],size(M1,1),1);
                    addRows2 = randi([1 size(M2,1)],size(M2,1),1);
                    M1 = [M1; M1(addRows1(:),:);];
                    M2 = [M2; M2(addRows2(:),:);];

                    disp(size(M1));
                    disp(size(M2));
                end
            elseif (size(M1,1) < size(M2,1))
                if (dataThin)
                    MScrap = M2;
                    cutRows = randi([1 size(M2,1)],size(M1,1),1);
                    M2 = MScrap(cutRows(:),:);
                    
                    disp(size(M1));
                    disp(size(M2));
                else
                    MScrap = M1;
                    addRows = randi([1 size(M1,1)],abs(size(M2,1)-size(M1,1)),1);
                    M1 = [M1; MScrap(addRows(:),:);];
                    
                    % double the size of the data
                    addRows1 = randi([1 size(M1,1)],size(M1,1),1);
                    addRows2 = randi([1 size(M2,1)],size(M2,1),1);
                    M1 = [M1; M1(addRows1(:),:);];
                    M2 = [M2; M2(addRows2(:),:);];
                    
                    disp(size(M1));
                    disp(size(M2));
                end

            end

            % collect data for end point location and ligand center location
            % probably don't need end point location
            for i=1:3
                r1(:,i) = M1(:, 4+(i-1));
                r2(:,i) = M2(:,4+(i-1));
                ligandCenter1(:,i) = M1(:,7+3+(i-1));
                ligandCenter2(:,i) = M2(:,7+(i-1));
            end

%               for i=1:3
%                   r1(:,i) = M1(:, 12+(Nrods(n)-1)*3+(i-1)); % take end point from bound ligand sim
%                   r2(:,i) = M2(:,12+(Nrods(n)-1)*3+(i-1)); % take end point from unbound ligand sim
%                   ligandCenter1(:,i) = M1(:,7+(Nrods(n))*3+3+(i-1)); % take bound ligand center
%                   ligandCenter2(:,i) = M2(:,7+(Nrods(n))*3+(i-1)); % take unbound ligand center
%               end

            % create cell array of data

            ligandBoundrData{n,irL,membrane+1} = r1;
            ligandOffrData{n,irL,membrane+1} = r2;
            ligandBoundcenterData{n,irL,membrane+1} = ligandCenter1;
            ligandOffcenterData{n,irL,membrane+1} = ligandCenter2;

            % collect occlusion data from second data set
            occlusion(:) = M2(:,10);

            occlusionData{n,irL,membrane+1} = occlusion;
            %available(irL) = size(find(occlusion==0),2)
        end
    end

end

output = {ligandBoundrData,ligandOffrData,ligandBoundcenterData,ligandOffcenterData,occlusionData};


end
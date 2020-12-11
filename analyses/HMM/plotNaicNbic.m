function plotNaicNbic( aic, bic, nAic, nBic, resultSave )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    colors = getNumStateColor();
    maxK = length(nAic);
    
    dirName = resultSave.figuresDirName;
    if (~isdir(dirName))
        mkdir(dirName);
    end;
    
    %==================================
    % BIC
    %==================================
    hf = figure();
    set(hf,'visible','off');

    for iK = 1:maxK
        h = bar(iK,bic(iK));
        hold on;
        set(h, 'FaceColor', colors(iK,:));
    end;
   
    xlabel('Number of states', 'fontsize', resultSave.labelFontSize);
    ylabel('BIC', 'fontsize', resultSave.labelFontSize);
    margin = 0.1*(max(bic) - min(bic));
    ylim([min(bic)-margin max(bic)+margin]);
    xlim([0 maxK+1]);
    set(gca,'XTick',1:maxK);
    set(gca,'fontsize',resultSave.axisFontSize);
    set(gca, 'Color', 'none');
    box off;

    currDir = cd(dirName);
    fileName = 'bic.pdf';
    export_fig(fileName, '-transparent');
    cd(currDir);
    close(hf);
    
    
    
    %==================================
    % AIC
    %==================================
    hf = figure();
    set(hf,'visible','off');

    for iK = 1:maxK
        h = bar(iK,aic(iK));
        hold on;
        set(h, 'FaceColor', colors(iK,:));
    end;
    
    xlabel('Number of states', 'fontsize', resultSave.labelFontSize);
    ylabel('AIC', 'fontsize', resultSave.labelFontSize);
    margin = 0.1*(max(aic) - min(aic));
    ylim([min(aic)-margin max(aic)+margin]);
    xlim([0 maxK+1]);
    set(gca,'XTick',1:maxK);
    set(gca,'fontsize',resultSave.axisFontSize);
    set(gca, 'Color', 'none');
    box off;

    currDir = cd(dirName);
    fileName = 'aic.pdf';
    export_fig(fileName, '-transparent');
    cd(currDir);
    close(hf);

    

    %==================================
    % Proportion model selected by AIC
    %==================================
    hf = figure();
    set(hf,'visible','off');

    for iK = 1:maxK
        h = bar(iK,nAic(iK));
        hold on;
        set(h, 'FaceColor', colors(iK,:));
    end;

    xlabel('Number of states', 'fontsize', resultSave.labelFontSize);
    ylabel('Fraction selected, AIC', 'fontsize', resultSave.labelFontSize);
    ylim([0 1]);
    xlim([0 maxK+1]);
    set(gca,'XTick',1:maxK);
    set(gca,'fontsize',resultSave.axisFontSize);
    set(gca, 'Color', 'none');
    box off;

    currDir = cd(dirName);
    fileName = 'naic.pdf';
    export_fig(fileName, '-transparent');
    cd(currDir);
    close(hf);
    
    %==================================
    % Proportion model selected by BIC
    %==================================
    hf = figure();
    set(hf,'visible','off');

    for iK = 1:maxK
        h = bar(iK,nBic(iK));
        hold on;
        set(h, 'FaceColor', colors(iK,:));
    end;

    xlabel('Number of states', 'fontsize', resultSave.labelFontSize);
    ylabel('Fraction selected, BIC', 'fontsize', resultSave.labelFontSize);
    ylim([0 1]);
    xlim([0 maxK+1]);
    set(gca,'XTick',1:maxK);
    set(gca,'fontsize',resultSave.axisFontSize);
    set(gca, 'Color', 'none');
    box off;

    currDir = cd(dirName);
    fileName = 'nbic.pdf';
    export_fig(fileName, '-transparent');
    cd(currDir);
    close(hf);


end


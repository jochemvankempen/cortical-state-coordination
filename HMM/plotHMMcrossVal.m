function plotHMMcrossVal( data, resultSave )
% plots results of cross-validation tests of the HMM fit
%
% arguments:
% ---------
%   data        -   structure that contains data to plot, fields are cvError,
%                   mrError, spError1, spError2
%   resultSave  -   data strcuture that contains field figuresDirName, where
%                   the figures will be saved, and figures a boolean
%                   variable whether to save figrues

    numState = data.numState;
    
    colors = [ 148  0    211; ... 1  dark violet
               34   139  34; ...  2  forest green
               0    0    205; ... 3  medium blue
               233  150  122; ... 4  dark salmon
               0    0    205; ... 5  medium blue
             ]/255;
    

    dirName = [ resultSave.figuresDirName 'numState_' num2str(numState) '/crossVal/'];  
    if (resultSave.figures) 
        if (~isdir(dirName))
            mkdir(dirName);
        end;
    else
        return
    end;
    % switch to the figures directory
    currDir = cd(dirName);
    
    if (~isempty(data.R2))
        % plot variance explained
        hf = figure;
        set(hf,'visible','off');
        temp = mean(data.R2,3);
        plot(data.timeWindow, temp','-o', 'LineWidth', 3, 'MarkerSize', 4); 
        xlim([0 max(data.timeWindow)+0.1]);
        box off;
        set(gca,'fontsize',resultSave.axisFontSize);
        set(gca, 'Color', 'none');
        xlabel('Integration window (s)','fontsize', resultSave.labelFontSize);
        ylabel('Variance explained','fontsize', resultSave.labelFontSize);
        fileName = 'cvError.pdf';
        export_fig(fileName, '-transparent');
        close(hf);
    end;
    
    if (~isempty(data.veR2))
        % plot variance explained
        hf = figure;
        set(hf,'visible','off');
        temp = mean(data.veR2,3);
        plot(data.timeWindow, temp','-o', 'LineWidth', 3, 'MarkerSize', 4); 
        xlim([0 max(data.timeWindow)+0.1]);
        box off;
        set(gca,'fontsize',resultSave.axisFontSize);
        set(gca, 'Color', 'none');
        xlabel('Integration window (s)','fontsize', resultSave.labelFontSize);
        ylabel('Variance explained','fontsize', resultSave.labelFontSize);
        fileName = 'veError.pdf';
        export_fig(fileName, '-transparent');
        close(hf);
    end;
    
 
    %{
    % plot cvError and mrError
    hf = figure;
    set(hf,'visible','off');
    plot(data.mrError,'-o', 'Color', colors(4,:), 'MarkerFaceColor', colors(4,:), 'MarkerEdgeColor', colors(4,:), 'LineWidth', 3); 
    hold on
    plot(data.cvError,'-o', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', colors(3,:), 'LineWidth', 3);
    xlim([0 size(data.cvError,1)+1]);
    box off;
    set(gca,'fontsize',resultSave.axisFontSize);
    set(gca, 'Color', 'none');
    xlabel('Channel number','fontsize', resultSave.labelFontSize);
    ylabel('Normalized CVE','fontsize', resultSave.labelFontSize);
    fileName = 'cvError.pdf';
    export_fig(fileName, '-transparent');
    close(hf);
    
    % plot difference between mrErrror and cvError
    hf = figure;
    set(hf,'visible','off');
    plot(data.mrError-data.cvError,'-o', 'Color', colors(5,:), 'MarkerFaceColor', colors(5,:), 'MarkerEdgeColor', colors(5,:), 'LineWidth', 3);
    hold on;
    xlim([0 size(data.cvError,1)+1]);
    plot(xlim, [0 0], '--k');
    box off;
    set(gca,'fontsize',resultSave.axisFontSize);
    set(gca, 'Color', 'none');
    xlabel('Channel number','fontsize', resultSave.labelFontSize);
    ylabel('Reduction in CVE','fontsize', resultSave.labelFontSize);
    fileName = 'cvErrorReduc.pdf';
    export_fig(fileName, '-transparent');
    close(hf);
    %}
    
    % plot state prediction errors
    if (isfield(data,'spError1'))
        hf = figure;
        set(hf,'visible','off');
        plot(data.spError1,'-o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:), 'LineWidth', 3);
        hold on
        plot(data.spError2,'-o', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:), 'LineWidth', 3);
        xlim([0 size(data.cvError,1)+1]);
        plot(xlim, 1-[1 1]/numState, '--k');
        box off;
        set(gca,'fontsize',resultSave.axisFontSize);
        set(gca, 'Color', 'none');
        xlabel('Channel number','fontsize', resultSave.labelFontSize);
        ylabel('Fraction mismatch in decoded state','fontsize', resultSave.labelFontSize);
        fileName = 'spError.pdf';
        export_fig(fileName, '-transparent');
        close(hf);
    end;
    
    % switch back to the working directory
    cd(currDir);

end


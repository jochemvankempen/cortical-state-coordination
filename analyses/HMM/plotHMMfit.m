function plotHMMfit( mu, data, resultSave )
    
    % get the color pallet for states
    colors = getStateColor( );
  
    dirName = [ resultSave.figuresDirName 'numState_' num2str(data.numState) '/'];   
    % make figures
    if (resultSave.figures) 
        if (~exist(dirName,'dir'))
            mkdir(dirName);
        end;

        % switch to the figures directory
        currDir = cd(dirName);


        numChannel = size(mu,1);
        numTrial = min(resultSave.numTrialToPlot,size(mu,2));

        % determine unique states
        maxLength = max( cellfun(@(x) length(unique(x)), data.states ) );
        temp = cell2mat( cellfun(@(x) [unique(x), NaN(1,maxLength-length(unique(x))) ], data.states,'UniformOutput',0 ) );
        uniqueState = unique( temp(~isnan(temp)) );
        clear temp;
        numState = data.numState;

        if (numState>1)
            % plot individual trial examples
            for iTrial = 1:numTrial

                hf = figure;
                set(hf,'visible','off');
                defP = get(hf,'Position');
                set(hf, 'Position', [0 0 defP(3) defP(4)*1.5]);
                % combined mu activity
                subplot('Position',[0.12 0.77 0.8 0.12]);
                allSpike = vertcat( mu{:,iTrial} );
                trace = zeros(size(data.timeBins{iTrial}));
                for iSpike = 1:length(allSpike)
                    trace = trace + normpdf(data.timeBins{iTrial},allSpike(iSpike),0.001);
                end;
                trace = trace/numChannel;
                plot(data.timeBins{iTrial},trace,'-k','LineWidth',5);
                xlim(data.timeEpoch(iTrial,:));
                ylim([0 max(trace)]);
                box off;
                %set(gca, 'visible', 'off') ;
                set(gca,'fontsize',resultSave.axisFontSize);
                set(gca, 'Color', 'none');
                set(gca,'XTick',[]);
                ylabel('FR, Hz','fontsize', resultSave.labelFontSize);
                %
                % spike raster
                subplot('Position',[0.12 0.59 0.8 0.15]);
                for iChann = 1:numChannel
                    a = mu{iChann,iTrial};
                    a = a( a>data.timeEpoch(iTrial,1) & a<=data.timeEpoch(iTrial,2) );
                    if (~isempty(a))
                        plot(a, numChannel-iChann+1,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2 );
                    end;
                    hold on;
                end;
                xlim(data.timeEpoch(iTrial,:));
                ylim([0 numChannel+1]);
                box off;
                set(gca, 'visible', 'off') ;
                set(gca,'fontsize',resultSave.axisFontSize);
                set(gca, 'Color', 'none');
                %
                % state trajectory
                subplot('Position',[0.12 0.36 0.8 0.2])
                stairs(data.timeBins{iTrial},data.states{iTrial},'Color',[30 144 255]/255,'LineWidth',5);
                hold on
                box off;
                xlim(data.timeEpoch(iTrial,:));
                ylim([1 1.02*numState]);
                set(gca,'fontsize',resultSave.axisFontSize);
                set(gca, 'Color', 'none');
                set(gca,'XTick',[],'YTick',1:numState);
                ylabel('State','fontsize', resultSave.labelFontSize);
                %
                % probabilities of states
                subplot('Position',[0.12 0.13 0.8 0.2])
                for iState = 1:numState
                    plot(data.timeBins{iTrial},data.Pstates{iTrial}(iState,:),'Color', colors(iState,:),'LineWidth',5);
                    hold on
                end;
                plot(xlim,[0.8 0.8],'--k');
                box off;
                xlim(data.timeEpoch(iTrial,:));
                ylim([0 1.02]);
                set(gca,'fontsize',resultSave.axisFontSize);
                set(gca, 'Color', 'none');
                set(gca,'YTick',[0 1]);
                xlabel('Time, s','fontsize', resultSave.labelFontSize);
                ylabel('P(state)','fontsize', resultSave.labelFontSize);
                fileName = ['states_trial' num2str(iTrial) '.pdf'];
                export_fig(fileName, '-transparent');
                close(hf);  
            end;
        end;

        %==================================
        % emission matrix
        %==================================
        hf = figure;
        set(hf,'visible','off');
        for iState = 1:numState
            errorbar(1:numChannel, data.estEmis(iState,:)/data.binSize, (data.estEmis(iState,:) - squeeze(data.estEmisBcb(1,iState,:))')/data.binSize, ...
                (squeeze(data.estEmisBcb(2,iState,:))'-data.estEmis(iState,:))/data.binSize, '-o','Color', colors(iState,:),'LineWidth',3,'MarkerFaceColor', ...
                colors(iState,:),'MarkerEdgeColor', colors(iState,:));
            hold on;
        end;
        box off;
        ylim([0 1.05*max(max(data.estEmisBcb(2,:,:)))/data.binSize]);
        xlim([0 numChannel+1])
        set(gca,'fontsize',resultSave.axisFontSize);
        set(gca, 'Color', 'none');
        xlabel('Channel number','fontsize', resultSave.labelFontSize);
        ylabel('Firing rate, Hz','fontsize', resultSave.labelFontSize);
        tickMarks = 1:floor(size(data.estEmis,2)/4):size(data.estEmis,2);
        set(gca,'XTick', tickMarks);
        fileName = 'emissionMatrix.pdf';
        export_fig(fileName, '-transparent');
        close(hf); 


        %==================================
        % transition matrix
        %==================================
        if (numState>1)
            hf = figure;
            set(hf,'visible','off');
            t = imresize(data.estTrans,50,'nearest');
            imagesc(0.5+(1:numState*50)/50, 0.5+(1:numState*50)/50, t );
            colormap(flipud(gray(256)));
            caxis([0 1]);
            colorbar;

            box off;
            set(gca,'fontsize',resultSave.axisFontSize);
            set(gca, 'Color', 'none');
            xlabel('State label','fontsize', resultSave.labelFontSize);
            ylabel('State label','fontsize', resultSave.labelFontSize);
            set(gca,'XTick', 1:numState);
            set(gca,'YTick', 1:numState);
            fileName = 'transition.pdf';
            export_fig(fileName, '-transparent');
            close(hf);
        end;

        %==================================
        % biograph representation of the transition matrix
        %==================================
        if (numState>1)
            t = data.estTrans;
            t(data.estTrans<data.binSize) = 0;
            t = t & 1-eye(numState);
            if ( any(any(t)) ) % plot the graph only of there are edges
                labels  = cell(numState,1);
                for iState = 1:numState
                    labels{iState} = num2str(iState);
                end;
                bg = biograph(t,labels);
                set(bg.nodes,'FontSize',resultSave.labelFontSize,'Shape','circle','Linecolor', [0 0 0]);
                for iState = 1:numState
                    set(bg.nodes(iState),'Color',colors(iState,:));
                end;
                set(bg,'LayoutType','equilibrium','EdgeType','curved');
                set(bg.edges,'LineColor',[0 0 0],'LineWidth',2);

                view(bg);
                % hack around to get the figure handle
                child_handles = allchild(0);
                names = get(child_handles,'Name');
                k = strncmp('Biograph Viewer', names, 15);
                set(gca, 'Color', 'none');
                fileName = 'biograph_transition.pdf';
                export_fig(child_handles(k),fileName, '-transparent');
                close all hidden;
            end;
        end;


        %==================================
        % residence time distribution
        %==================================
        if (numState>1)
            for iState = 1:numState
                if (~isempty(data.residTimeCDF{iState}))
                    % cdf
                    hf = figure;
                    set(hf,'visible','off');

                    stairs(data.residTimeCDF{iState}(:,1), data.residTimeCDF{iState}(:,2), 'Color',colors(iState,:), 'LineWidth',5);
                    hold on

                    tMax = ceil( max(data.residTimeCDF{iState}(:,1)/data.binSize) );
                    bins = 0:tMax;
                    y = geocdf(bins, 1-data.estTrans(iState,iState));
                    stairs(bins*data.binSize, y, 'k', 'LineWidth',2);
                    xlim([0 tMax*data.binSize]);
                    box off;
                    set(gca,'fontsize',resultSave.axisFontSize);
                    set(gca, 'Color', 'none');
                    xlabel('Residence time, s','fontsize', resultSave.labelFontSize);
                    ylabel('Cumulative probability','fontsize', resultSave.labelFontSize);
                    fileName = ['residTimeCDF_state' num2str(iState) '.pdf'];
                    export_fig(fileName, '-transparent');
                    close(hf);

                    scale = (data.residTimePDF{iState}(2,1)-data.residTimePDF{iState}(1,1))/data.binSize;
                    % pdf log
                    hf = figure;
                    set(hf,'visible','off');

                    semilogy(data.residTimePDF{iState}(:,1), data.residTimePDF{iState}(:,2),'o', 'MarkerFaceColor',colors(iState,:),...
                        'MarkerEdgeColor',colors(iState,:), 'MarkerSize', 12);
                    hold on;

                    tMax = ceil( max(data.residTimePDF{iState}(:,1)/data.binSize) );
                    bins = 0:tMax;
                    y = geopdf(bins, 1-data.estTrans(iState,iState));
                    semilogy(bins*data.binSize, y*scale, '--k', 'LineWidth',3);
                    xlim([0 tMax*data.binSize]);
                    box off;
                    set(gca,'fontsize',resultSave.axisFontSize);
                    set(gca, 'Color', 'none');
                    xlabel('Residence time, s','fontsize', resultSave.labelFontSize);
                    ylabel('Probability density','fontsize', resultSave.labelFontSize);
                    fileName = ['residTimePDFlog_state' num2str(iState) '.pdf'];
                    export_fig(fileName, '-transparent');
                    close(hf);

                    % pdf
                    hf = figure;
                    set(hf,'visible','off');

                    shift = 0.5*(data.residTimePDF{iState}(2,1)-data.residTimePDF{iState}(1,1));
                    stairs(data.residTimePDF{iState}(:,1)-shift, data.residTimePDF{iState}(:,2), 'Color',colors(iState,:), 'LineWidth',5);
                    hold on;

                    tMax = ceil( max(data.residTimePDF{iState}(:,1)/data.binSize) );
                    bins = 0:tMax;
                    y = geopdf(bins, 1-data.estTrans(iState,iState));
                    plot(bins*data.binSize, y*scale, '--k', 'LineWidth',3);
                    xlim([0 tMax*data.binSize]);
                    box off;
                    set(gca,'fontsize',resultSave.axisFontSize);
                    set(gca, 'Color', 'none');
                    xlabel('Residence time, s','fontsize', resultSave.labelFontSize);
                    ylabel('Probability density','fontsize', resultSave.labelFontSize);
                    fileName = ['residTimePDF_state' num2str(iState) '.pdf'];
                    export_fig(fileName, '-transparent');
                    close(hf);

                end;
            end;
        end;

        %==================================
        % fraction of time spent in each state
        %==================================
        if (numState>1)
            hf = figure;
            set(hf,'visible','off');
            for iState = 1:numState
                h = bar(iState,data.fracTimeState(iState));
                hold on;
                set(h, 'FaceColor', colors(iState,:));
            end;

            xlabel('State label', 'fontsize', resultSave.labelFontSize);
            ylabel('Fraction of time occupied', 'fontsize', resultSave.labelFontSize);
            ylim([0 1]);
            xlim([0 numState+1]);
            set(gca,'XTick',1:numState);
            set(gca,'fontsize',resultSave.axisFontSize);
            set(gca, 'Color', 'none');
            box off;
            fileName = 'fracTimeState.pdf';
            export_fig(fileName, '-transparent');
            close(hf);
        end;


        % switch back to the working directory
        cd(currDir);
    
    end;

end















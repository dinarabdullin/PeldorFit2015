function PeldorFitErrorPlot(varargin)
% This function generates RMSD plots for the PeldorFit fitting parameters
% varargin: 
% 1st parameter: a maximal RMSD value
% 2nd parameter: a minimal RMSD value

    % The names of the PeldorFit data files
    fileName_r     = 'error_plot_1_2.dat  ';
    fileName_xi    = 'error_plot_3_4.dat  ';
    fileName_phi   = 'error_plot_5_6.dat  ';
    fileName_alpha = 'error_plot_7_8.dat  ';
    fileName_beta  = 'error_plot_9_10.dat ';
    fileName_gamma = 'error_plot_11_12.dat';
    fileNames = [fileName_r; fileName_xi; fileName_phi; ...
                 fileName_alpha; fileName_beta; fileName_gamma]; 
    % Figure settings 
    nPlots = 0;
    plots = []; 
    alignement = [1, 1; 1, 2; 1, 3; 2, 2; 2, 3; 2, 3];
    % Optional parameter
    numvarargs = length(varargin);
    if (numvarargs > 2)
        error('myfuns:somefun2Alt:TooManyInputs', ...
              'Too many input parameters');
    end
    optargs = {0, 0};
    optargs(1:numvarargs) = varargin;
    [zmax_user, zmin_user] = optargs{:};
    
    [fileName, pathName] = uigetfile('*.DAT','File Selector');
    if ~isequal(fileName,0)
        % Check which error plots are availbale
        for fn = 1:6
            filePath = strcat(pathName,fileNames(fn,:));
            if exist(filePath, 'file')
                plots = [plots fn];              
            end
        end
        nPlots = size(plots, 2);
        % Determine the max and min values of RMSD for all plots
        count = 0;
        zmax = 0; 
        zmin = 0;
        zmax1 = 0;
        zmin1 = 0;
        for fn = plots
            filePath = strcat(pathName,fileNames(fn,:));
            P = dlmread(filePath, '', 1, 0);
            if (count == 0)
                zmax = max(P(:,3));
                zmin = min(P(:,3)); 
            else
                zmax1 = max(P(:,3));
                zmin1 = min(P(:,3));
                if (zmax1 > zmax)
                    zmax = zmax1;
                end
                if (zmin1 < zmin)
                    zmin = zmin1;
                end
            end
            count = count + 1;
        end
        if (str2double(zmax_user) > 0)
            zmax = str2double(zmax_user);
        end
        class(zmax_user)
        if (str2double(zmin_user) > 0)
            zmin = str2double(zmin_user);
        end
        class(zmin_user)
        
        % Create a figure
        hfig = figure(1);
        set(gcf,'color','w'); 
        count = 1;
        for fn = plots
            % Read a data file
            filePath = strcat(pathName,fileNames(fn,:));
            P = dlmread(filePath, '', 1, 0);
            x = P(:,1);
            y = P(:,2);
            z = P(:,3);
            % Interpolate data on regular grid
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            xsteps = xmin:0.01*(xmax-xmin):xmax;
            ysteps = ymin:0.01*(ymax-ymin):ymax;
            [XI, YI] = meshgrid(xsteps, ysteps);
            ZI = griddata(x,y,z,XI,YI);
            % Make a plot
            ax(count) = subplot(alignement(nPlots,1),alignement(nPlots,2)+1,count);
            if count == alignement(nPlots,2)
                count = count + 2;
            else
                count = count + 1;
            end
            pcolor(XI, YI, ZI);
            colormap(flipud(jet))
            shading interp
            axis square
            set(gca,'clim',[zmin, zmax])
            set(gca,'FontSize',28,'linewidth',1,'TickDir','out','Box','on')
            % Calculate tick values
            if (fn == 1)
                xlow = 0.5*round(xmin/0.5);
                xhigh = 0.5*round(xmax/0.5);
                dx = 0.5*(xhigh-xlow);
                ylow = 0.1*round(ymin/0.1);
                yhigh = 0.1*round(ymax/0.1);
                dy = 0.5*(yhigh-ylow);
                xlim([xlow,xhigh]);
                ylim([ylow,yhigh]);
                set(gca,'XTick',[xlow:dx:xhigh], 'YTick',[ylow:dy:yhigh]);
            else
                xlow = 10*round(xmin/10);
                xhigh = 10*round(xmax/10);
                dx = 0.5*(xhigh-xlow);
                ylow = 10*round(ymin/10);
                yhigh = 10*round(ymax/10);
                dy = 0.5*(yhigh-ylow);
                xlim([xlow,xhigh]);
                ylim([ylow,yhigh]);
                set(gca,'XTick',[xlow:dx:xhigh], 'YTick',[ylow:dy:yhigh]);
            end
            % Axes labels
            if (fn == 1)
                xlabel('\itr\rm (nm)');
                ylabel('\rm\Delta\itr\rm (nm)');
            end
            if (fn == 2)
                xlabel('\it\xi\rm °');
                ylabel('\rm\Delta\it\xi\rm °');
            end
            if (fn == 3)
                xlabel('\it\phi\rm °');
                ylabel('\rm\Delta\it\phi\rm °');
            end
            if (fn==4)
                xlabel('\it\alpha\rm °');
                ylabel('\rm\Delta\it\alpha\rm °');
            end
            if (fn==5)
                xlabel('\it\beta\rm °');
                ylabel('\rm\Delta\it\beta\rm °');
            end
            if (fn==6)
                xlabel('\it\gamma\rm °');
                ylabel('\rm\Delta\it\gamma\rm °');
            end
        end
        hcolorbar = colorbar('FontSize',28);
        set(hcolorbar,'Position',[0.80 0.25 0.025 0.5]);
        set(get(hcolorbar,'title'),'string','\itRMSD','FontSize',28);
    end
end

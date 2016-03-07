function PeldorFitErrorPlot(varargin)
% This function generates the error plots for all geometric parameters 
% used by the program PeldorFit 
% varargin: Here the user can specify the maximal value for colorbar

    % The names of the PeldorFit data files
    fileName_r     = 'error_plot_1_2.dat  ';
    fileName_xi    = 'error_plot_3_4.dat  ';
    fileName_phi   = 'error_plot_5_6.dat  ';
    fileName_alpha = 'error_plot_7_8.dat  ';
    fileName_beta  = 'error_plot_9_10.dat ';
    fileName_gamma = 'error_plot_11_12.dat';
    fileNames = [fileName_r; fileName_xi; fileName_phi; fileName_alpha; fileName_beta; fileName_gamma]; 
    % Figure settings 
    nPlots = 0;
    plots = []; 
    alignement = [1, 1; 1, 2; 1, 3; 2, 2; 2, 3; 2, 3];
    % Optional parameter
    numvarargs = length(varargin);
    if (numvarargs > 1)
        error('myfuns:somefun2Alt:TooManyInputs', ...
              'Too many input parameters');
    end
    optargs = {0};
    optargs(1:numvarargs) = varargin;
    [zmax_user] = optargs{:};
    
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
            ax(count) = subplot(alignement(nPlots,1),alignement(nPlots,2),count);
            count = count + 1;
            pcolor(XI, YI, ZI);
            colormap(flipud(jet))
            shading interp
            axis square
            set(gca,'clim',[zmin, zmax])
            set(gca,'FontSize',24,'linewidth',1,'TickDir','out','Box','on')
            if (fn == 1)
                xlabel('\it\mu\rm (nm)');
                ylabel('\it\sigma\rm (nm)');
                xlim([0.5*round(xmin/0.5),0.5*round(xmax/0.5)]);
                ylim([0.1*round(ymin/0.1),0.1*round(ymax/0.1)]);
                set(gca,'XTick',[0.5*round(xmin/0.5):0.5:0.5*round(xmax/0.5)], ...
                        'YTick',[0.1*round(ymin/0.1):0.2:0.1*round(ymax/0.1)]);
            end
            if (fn == 2)
                xlabel('\it\xi\rm °');
                ylabel('\rm\Delta\it\xi\rm °');
                xlim([0,90]);
                ylim([0,90]);
                set(gca,'XTick',[0:30:90],'YTick',[0:30:90]);
            end
            if (fn == 3)
                xlabel('\it\phi\rm °');
                ylabel('\rm\Delta\it\phi\rm °');
                xlim([0,180]);
                ylim([0,180]);
                set(gca,'XTick',[0:60:180],'YTick',[0:60:180]);
            end
            if (fn==4)
                xlabel('\it\alpha\rm °');
                ylabel('\rm\Delta\it\alpha\rm °');
                xlim([0,180]);
                ylim([0,180]);
                set(gca,'XTick',[0:60:180],'YTick',[0:60:180]);
            end
            if (fn==5)
                xlabel('\it\beta\rm °');
                ylabel('\rm\Delta\it\beta\rm °');
                xlim([0,90]);
       
                ylim([0,90]);
                set(gca,'XTick',[0:30:90],'YTick',[0:30:90]);
            end
            if (fn==6)
                xlabel('\it\gamma\rm °');
                ylabel('\rm\Delta\it\gamma\rm °');
                xlim([0,180]);
                ylim([0,180]);
                set(gca,'XTick',[0:60:180],'YTick',[0:60:180]);
            end
        end
        hcolorbar = colorbar('FontSize',24);
        set(hcolorbar,'Position',[0.93 0.25 0.025 0.5]);
        set(get(hcolorbar,'title'),'string','\itRMSD','FontSize',24);
    end
end

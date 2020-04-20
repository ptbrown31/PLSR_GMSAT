close all
clear all

% load coast
% lat_coast = lat;
% lon_coast = long;

addpath(genpath('/Users/patrickbrown/Documents/MATLAB/'))

%% load

load('/Users/patrickbrown/Documents/MATLAB/PLS_decadal_prediction/interannual_SAT_GMST/CMIP5/MATS/hindcast_forecast_figs_all_obs_preproc_MLR.mat',...
    'land_frac_obs',...
    'obs_ann_time',...
    'obs_lon',...
    'obs_lat',...
    'forced_GMST_obs',...
    'obs_GMST_ann',...
    'obs_GMST_ann_unforced',...
    'forced_3d_ann_obs',...
    'obs_3d_ann_anom',...
    'obs_3d_ann_unforced',...
    'forced_GMST_obs_plus_forecast',...
    'obs_ann_time_plus_forecast')

%NOAA needs to be flipped

obs_names = {'BEST',...
             'GISTEMP',...
             'NOAA',...
             'HadCRUT4'};
         
%% douple check GISTEMP GMST lines up

GISTEMP_their_GMSAT = xlsread('/Users/patrickbrown/Documents/MATLAB/PLS_decadal_prediction/interannual_SAT_GMST/raw_obs/GISTEMP_GMSAT_1880_2017.xlsx','graph','B3:B140');
GISTEMP_their_GMSAT_year = 1880:2017;

%make sure its anomed the same
GISTEMP_their_GMSAT = GISTEMP_their_GMSAT - mean(GISTEMP_their_GMSAT(GISTEMP_their_GMSAT_year >= 1951 & GISTEMP_their_GMSAT_year <= 1980));

%% test plot

figure
hold on
    plot(GISTEMP_their_GMSAT_year,GISTEMP_their_GMSAT,'-k','LineWidth',2)
    scatter(GISTEMP_their_GMSAT_year,GISTEMP_their_GMSAT,25,'k','filled')
    plot(obs_ann_time{2},obs_GMST_ann{2},'-b','LineWidth',2)
    scatter(obs_ann_time{2},obs_GMST_ann{2},25,'b','filled')

%% set options

% shift

    num_lag_predictor = 2;

% Preproc option

    norm_predictors_flag = 1; % 0 = no, 1=yes,

% starting time for PREDICTEMP

    start_year = 1900;
    
% define forecast/backcast delimination
    crawling_forecast_begin_year = 2000;

    lasso_alpha = 0.5;
    kval_number = 10;
    mc_trials = 1;
    num_pls_comps = 4;
    
%% loop through all obs

xval_and_crawling_forecast_predictions_all_obs = {};
X_val_and_crawling_climatology_all_obs = {};
prediction_RMSEs_all_obs = {};
prediction_Rs_all_obs = {};
prediction_RMSSSs_all_obs = {};

for forecast_period_i = 1:4

        forecast_period = forecast_period_i;

    for obsi = 1:length(obs_names)

        disp('calculating for')
        obs_names{obsi}
        disp('lag')
        forecast_period
        
            optimum_1st_year_i = find(obs_ann_time{obsi} == start_year);

            obs_ann_time{obsi} = obs_ann_time{obsi}(optimum_1st_year_i:end);

            obs_3d_ann_anom{obsi} = obs_3d_ann_anom{obsi}(:,:,optimum_1st_year_i:end);
            forced_3d_ann_obs{obsi} = forced_3d_ann_obs{obsi}(:,:,optimum_1st_year_i:end);
            obs_3d_ann_unforced{obsi} = obs_3d_ann_unforced{obsi}(:,:,optimum_1st_year_i:end);

            obs_GMST_ann{obsi} = obs_GMST_ann{obsi}(optimum_1st_year_i:end);    

            forced_GMST_obs{obsi} = forced_GMST_obs{obsi}(optimum_1st_year_i:end);
            obs_GMST_ann_unforced{obsi} = obs_GMST_ann_unforced{obsi}(optimum_1st_year_i:end);

            obs_ann_time_plus_forecast{obsi} = obs_ann_time_plus_forecast{obsi}(optimum_1st_year_i:end);    
            forced_GMST_obs_plus_forecast{obsi} = forced_GMST_obs_plus_forecast{obsi}(optimum_1st_year_i:end);


        %% optionally normalize predictor feilds

        if norm_predictors_flag == 1        
            for lati = 1:length(obs_lat{obsi})
                for loni = 1:length(obs_lon{obsi})

                    %obs_3d_ann_unforced_norm{obsi}(lati,loni,:) = obs_3d_ann_unforced{obsi}(lati,loni,:)./std(squeeze(obs_3d_ann_unforced{obsi}(lati,loni,:)));
                    obs_3d_ann_unforced_norm{obsi}(lati,loni,:) = obs_3d_ann_anom{obsi}(lati,loni,:)./std(squeeze(obs_3d_ann_anom{obsi}(lati,loni,:)));

                end
            end
        end

        %% pick an obs dataset

        %Y_predictand = obs_GMST_ann_unforced{obsi}';
        Y_predictand = obs_GMST_ann{obsi};

        X_predictor_field = obs_3d_ann_unforced_norm{obsi};
        ann_time = obs_ann_time{obsi};

        crawling_forecast_ind = find(ann_time == crawling_forecast_begin_year);

        %% assemble predictor and predictand

                    Y_predictand_shift = Y_predictand(forecast_period+num_lag_predictor:end);

                % predictor part
                    if num_lag_predictor == 1

                        predictor_2D_space_lag_1 = X_predictor_field(:,:,num_lag_predictor:end-num_lag_predictor-forecast_period+1);

                        X_predictor_field_shift = cat(4,predictor_2D_space_lag_1);

                        X_predictor_matrix_shift = obs_GMST_ann_unforced{obsi}(num_lag_predictor:end-num_lag_predictor-forecast_period+1)';


                    end
                    if num_lag_predictor == 2

                        predictor_2D_space_lag_1 = X_predictor_field(:,:,num_lag_predictor:end-num_lag_predictor-forecast_period+2);
                        predictor_2D_space_lag_2 = X_predictor_field(:,:,num_lag_predictor-1:end-num_lag_predictor-forecast_period+1);

                        X_predictor_field_shift = cat(4,predictor_2D_space_lag_1,predictor_2D_space_lag_2);

                        X_predictor_matrix_shift = obs_GMST_ann_unforced{obsi}(num_lag_predictor:end-num_lag_predictor-forecast_period+2)';


                    end
                    if num_lag_predictor == 3

                        predictor_2D_space_lag_1 = X_predictor_field(:,:,num_lag_predictor:end-num_lag_predictor-forecast_period+3);
                        predictor_2D_space_lag_2 = X_predictor_field(:,:,num_lag_predictor-1:end-num_lag_predictor-forecast_period+2);
                        predictor_2D_space_lag_3 = X_predictor_field(:,:,num_lag_predictor-2:end-num_lag_predictor-forecast_period+1);

                        X_predictor_field_shift = cat(4,predictor_2D_space_lag_1,predictor_2D_space_lag_2,predictor_2D_space_lag_3);

                        X_predictor_matrix_shift = obs_GMST_ann_unforced{obsi}(num_lag_predictor:end-num_lag_predictor-forecast_period+3)';

                    end
                    if num_lag_predictor == 4

                        predictor_2D_space_lag_1 = X_predictor_field(:,:,num_lag_predictor:end-num_lag_predictor-forecast_period+4);
                        predictor_2D_space_lag_2 = X_predictor_field(:,:,num_lag_predictor-1:end-num_lag_predictor-forecast_period+3);
                        predictor_2D_space_lag_3 = X_predictor_field(:,:,num_lag_predictor-2:end-num_lag_predictor-forecast_period+2);
                        predictor_2D_space_lag_4 = X_predictor_field(:,:,num_lag_predictor-3:end-num_lag_predictor-forecast_period+1);

                        X_predictor_field_shift = cat(4,predictor_2D_space_lag_1,predictor_2D_space_lag_2,predictor_2D_space_lag_3,predictor_2D_space_lag_4);

                        X_predictor_matrix_shift = obs_GMST_ann_unforced{obsi}(num_lag_predictor:end-num_lag_predictor-forecast_period+4)';

                    end

        %% put into PLSR/Lasso function

        Y_predictand_shift_2D = NaN(1,1,length(Y_predictand_shift));
        Y_predictand_shift_2D(1,1,:) = Y_predictand_shift;

            [xval_and_crawling_forecast_predictions,...
            X_val_and_crawling_climatology,...
            prediction_RMSEs,...
            prediction_Rs,...
            prediction_RMSSSs] = ...
            PLSR_field_with_local_lasso_crawling_forecast(...
            Y_predictand_shift_2D,...
            X_predictor_matrix_shift,...
            X_predictor_field_shift,...
            crawling_forecast_ind,...
            lasso_alpha,...
            kval_number,...
            mc_trials,...
            num_pls_comps);

        xval_and_crawling_forecast_predictions_all_obs{obsi,forecast_period_i} = xval_and_crawling_forecast_predictions;
        X_val_and_crawling_climatology_all_obs{obsi,forecast_period_i} = X_val_and_crawling_climatology;
        prediction_RMSEs_all_obs{obsi,forecast_period_i} = prediction_RMSEs;
        prediction_Rs_all_obs{obsi,forecast_period_i} = prediction_Rs;
        prediction_RMSSSs_all_obs{obsi,forecast_period_i} = prediction_RMSSSs;        
        

    end
end

%% PLS forecasting

forecasts_foreward = NaN(4,length(obs_names));
forecasts_foreward_RMSE = NaN(4,length(obs_names));
frac_var_explained = NaN(num_pls_comps,4,length(obs_names));

XLoadings_maps_all_obs = {};
XScores_t_series_all_obs = {};

for obsi = 1:length(obs_names)
    
    % assign variables
        if norm_predictors_flag == 1        
            for lati = 1:length(obs_lat{obsi})
                for loni = 1:length(obs_lon{obsi})

                    %obs_3d_ann_unforced_norm{obsi}(lati,loni,:) = obs_3d_ann_unforced{obsi}(lati,loni,:)./std(squeeze(obs_3d_ann_unforced{obsi}(lati,loni,:)));
                    obs_3d_ann_unforced_norm{obsi}(lati,loni,:) = obs_3d_ann_anom{obsi}(lati,loni,:)./std(squeeze(obs_3d_ann_anom{obsi}(lati,loni,:)));

                end
            end
        end

        %Y_predictand = obs_GMST_ann_unforced{obsi}';
        Y_predictand = obs_GMST_ann{obsi};
        X_predictor_field = obs_3d_ann_unforced_norm{obsi};

    for forecast_periodi = 1:4

        [XLoadings,...
        YLoadings,...
        XScores,...
        YScores,...
        BETA,...
        PCTVAR,...
        RMS_prediction_error,...
        spread_ratio,...
        stats,...
        Y_foreward_forecast,...
        XLoadings_maps] = plsregress_leadlag_make_forecast(X_predictor_field,...
            Y_predictand,...
            num_pls_comps,...
            forecast_periodi,...
            num_lag_predictor);

        forecasts_foreward(forecast_periodi,obsi) = Y_foreward_forecast;
        forecasts_foreward_RMSE(forecast_periodi,obsi) = RMS_prediction_error;

        frac_var_explained(:,forecast_periodi,obsi) = PCTVAR(2,:)';
        
        if forecast_periodi == 1; XLoadings_maps_all_obs{obsi} = XLoadings_maps; end
        if forecast_periodi == 1; XScores_t_series_all_obs{obsi} = XScores; end

    end
end

%% make a persistence baseline (the above program uses a climatology baseline that is too easy to beat

persistence_rmse_all_obs = {};
persistence_r_all_obs = {};
persistence_hindcast_all_obs = {};
persistence_hindcast_time_all_obs = {};
persistence_hindcast_errs_all_obs = {};

forecast_periods = 1:4;
persistence_average_length = 5;

for obsi = 1:length(obs_names)
    
     predictand_t_series = obs_GMST_ann{obsi};
     
     optimum_1st_year_i = find(obs_ann_time{obsi} == start_year);
     full_time = obs_ann_time{obsi}(optimum_1st_year_i:end);

     [persistence_rmse,...
      persistence_r,...
      persistence_hindcast,...
      persistence_hindcast_time,...
      persistence_hindcast_errs] = persistence_hindcast_errors(full_time,...
      predictand_t_series,...
      persistence_average_length,...
      forecast_periods);
  
persistence_rmse_all_obs{obsi} = persistence_rmse;
persistence_r_all_obs{obsi} = persistence_r;
persistence_hindcast_all_obs{obsi} = persistence_hindcast;
persistence_hindcast_time_all_obs{obsi} = persistence_hindcast_time;
persistence_hindcast_errs_all_obs{obsi} = persistence_hindcast_errs;

end

%% make an extrapolate linear trend baseline

% Extrapolate linear trend    
    extrap_trend_rmse_all_obs = {};
    extrap_trend_r_all_obs = {};
    extrap_trend_hindcast_all_obs = {};
    extrap_trend_hindcast_time_all_obs = {};
    extrap_trend_hindcast_errs_all_obs = {};
    
    foreward_years_trend = 4;
    
    forecast_periods = 1:foreward_years_trend;
    trend_fit_length = 20;

    for obsi = 1:length(obs_names)

         predictand_t_series = obs_GMST_ann{obsi};

         optimum_1st_year_i = find(obs_ann_time{obsi} == start_year);
         full_time = obs_ann_time{obsi}(optimum_1st_year_i:end);

          [extrap_lin_trend_rmse,...
          extrap_lin_trend_r,...
          extrap_lin_trend_hindcast,...
          extrap_lin_trend_hindcast_time,...
          extrap_lin_trend_hindcast_errs] = extrap_lin_trend_hindcast_errors(full_time,...
          predictand_t_series,...
          trend_fit_length,...
          forecast_periods);

         extrap_trend_rmse_all_obs{obsi} = extrap_lin_trend_rmse;
         extrap_trend_r_all_obs{obsi} = extrap_lin_trend_r;
         extrap_trend_hindcast_all_obs{obsi} = extrap_lin_trend_hindcast;
         extrap_trend_hindcast_time_all_obs{obsi} = extrap_lin_trend_hindcast_time;
         extrap_trend_hindcast_errs_all_obs{obsi} = extrap_lin_trend_hindcast_errs;

    end

%% plot fig 3

close all

forecast_years = [2020 2021 2022 2023];

display_hindcast_forecast_mode_1 = 1;
display_hindcast_forecast_mode_2 = 2;

lead_colors = {'b','r','g','m'};
lead_sizes = [50 50 50 50];

for obsi = 1:length(obs_names)
        
    FigHandle = figure('Position', [100, 100, 1400, 800]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',14);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    
    legend_handles_plsr_forecasts = {};
            
        for forecast_period_i = 1:4
            
            subplot(2,2,forecast_period_i)
            hold on
            
            %Y_predictand_shift = Y_predictand(forecast_period_i+num_lag_predictor:end);
            Y_predictand_shift = obs_GMST_ann{obsi}(forecast_period_i+num_lag_predictor:end);

            ann_time_shift = ann_time(forecast_period_i+num_lag_predictor:end)';

            %plot hindcast/forecast cutoff

            lin_y_max = max(Y_predictand_shift)+0.5*std(Y_predictand_shift);
            lin_y_min = min(Y_predictand_shift)-0.5*std(Y_predictand_shift);

            plot([crawling_forecast_begin_year, crawling_forecast_begin_year],[lin_y_min lin_y_max],'--r','LineWidth',1)
            
            %plot obs
                plot(ann_time_shift,Y_predictand_shift,'-k','LineWidth',2);
%                 scatter(ann_time_shift,Y_predictand_shift,30,...
%                      'MarkerFaceColor','k','MarkerEdgeColor','k',...
%                      'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
             
            % plot persistence and errors
                        
%                 boundedline(persistence_hindcast_time_all_obs{obsi}{forecast_period_i},...
%                 persistence_hindcast_all_obs{obsi}{forecast_period_i},...
%                 persistence_rmse_all_obs{obsi}(forecast_period_i),...
%                 '-k','alpha','transparency',0.3)

%                 scatter(persistence_hindcast_time_all_obs{obsi}{forecast_period_i},...
%                 persistence_hindcast_all_obs{obsi}{forecast_period_i},lead_sizes(forecast_period_i),...
%                 'MarkerFaceColor','k','MarkerEdgeColor','k',...
%                  'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

            % plot extrapolated trend and errors
                
                xvals = extrap_trend_hindcast_time_all_obs{obsi}(forecast_period_i,:);
                yvals = extrap_trend_hindcast_all_obs{obsi}(forecast_period_i,:);
                good_inds = find(isnan(xvals) == 0);
            
%                 scatter(xvals(good_inds),...
%                 yvals(good_inds),...
%                 'MarkerFaceColor','k','MarkerEdgeColor','k',...
%                  'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%  
                boundedline(xvals(good_inds),...
                yvals(good_inds),...
                extrap_trend_rmse_all_obs{obsi}(forecast_period_i),...
                '-k','alpha','transparency',0.3)
                          
             % plot PLSR hindcasts
                boundedline(ann_time_shift,...
                    xval_and_crawling_forecast_predictions_all_obs{obsi,forecast_period_i},...
                    prediction_RMSEs_all_obs{obsi,forecast_period_i}(display_hindcast_forecast_mode_1),...
                    strcat('-',lead_colors{forecast_period_i}),'alpha','transparency',0.3)
            
                 plot(ann_time_shift,xval_and_crawling_forecast_predictions_all_obs{obsi,forecast_period_i},...
                     strcat('-',lead_colors{forecast_period_i}),...
                     'LineWidth',2)
                 
%                  scatter(ann_time_shift,xval_and_crawling_forecast_predictions_all_obs{obsi,forecast_period_i},lead_sizes(forecast_period_i),...
%                  'MarkerFaceColor',lead_colors{forecast_period_i},'MarkerEdgeColor',lead_colors{forecast_period_i},...
%                  'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
             
             %plot forecasts
             
                %plot the error bars

                    plot([forecast_years(forecast_period_i) forecast_years(forecast_period_i)],...
                        [forecasts_foreward(forecast_period_i,obsi)-2*forecasts_foreward_RMSE(forecast_period_i,obsi) forecasts_foreward(forecast_period_i,obsi)+2*forecasts_foreward_RMSE(forecast_period_i,obsi)],...
                        strcat('-',lead_colors{forecast_period_i}),'LineWidth',0.5)

                    plot([forecast_years(forecast_period_i) forecast_years(forecast_period_i)],...
                        [forecasts_foreward(forecast_period_i,obsi)-1*forecasts_foreward_RMSE(forecast_period_i,obsi) forecasts_foreward(forecast_period_i,obsi)+1*forecasts_foreward_RMSE(forecast_period_i,obsi)],...
                        strcat('-',lead_colors{forecast_period_i}),'LineWidth',1)

             %plot the central values
                scatter(forecast_years(forecast_period_i),...
                    forecasts_foreward(forecast_period_i,obsi),...
                    lead_sizes(forecast_period_i),lead_colors{forecast_period_i},'filled')

             
             xlim([min(ann_time) max(ann_time)+4])
             ylim([-0.7 1.4])
        
             %labels and such
                title(strcat('lead time = ',num2str(forecast_period_i),' year'))
                
                if forecast_period_i == 3 || forecast_period_i == 4; xlabel('Year'); end
                if forecast_period_i == 1 || forecast_period_i == 3; ylabel('Temperature Anomaly (C)'); end
                
             
        end  
end

%% plot fig 4

close all

display_hindcast_forecast_mode_1 = 1;
display_hindcast_forecast_mode_2 = 2;

lead_sizes = [50 30 10 5];

for obsi = 1:length(obs_names)
        
    FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
    set(gcf,'color',[1 1 1]);
    set(0, 'DefaultAxesFontSize',17);
    set(0,'defaultAxesFontName', 'helvetica')
    hold on
    
    legend_handles_plsr_forecasts = {};
        
        for forecast_period_i = 1:4
            
            %Y_predictand_shift = Y_predictand(forecast_period_i+num_lag_predictor:end);
            Y_predictand_shift = obs_GMST_ann{obsi}(forecast_period_i+num_lag_predictor:end);

            ann_time_shift = ann_time(forecast_period_i+num_lag_predictor:end)';

            %plot hindcast/forecast cutoff

            lin_y_max = max(Y_predictand_shift)+0.5*std(Y_predictand_shift);
            lin_y_min = min(Y_predictand_shift)-0.5*std(Y_predictand_shift);

            plot([crawling_forecast_begin_year, crawling_forecast_begin_year],[lin_y_min lin_y_max],'--r','LineWidth',1)
            
            % plot persistence and errors
            
%             if forecast_period_i == 1
%             
%                 boundedline(persistence_hindcast_time_all_obs{obsi}{forecast_period_i},...
%                 persistence_hindcast_all_obs{obsi}{forecast_period_i},...
%                 persistence_rmse_all_obs{obsi}(forecast_period_i),...
%                 '-k','alpha','transparency',0.3)
% 
%                 scatter(persistence_hindcast_time_all_obs{obsi}{forecast_period_i},...
%                 persistence_hindcast_all_obs{obsi}{forecast_period_i},lead_sizes(forecast_period_i),...
%                 'MarkerFaceColor','k','MarkerEdgeColor','k',...
%                  'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%              
%             end
            % plot extrapolation of linear trend
            
            if forecast_period_i == 1
                
                scatter(extrap_trend_hindcast_time_all_obs{obsi}(forecast_period_i,:),...
                extrap_trend_hindcast_all_obs{obsi}(forecast_period_i,:),...
                'MarkerFaceColor','k','MarkerEdgeColor','k',...
                 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
            
                boundedline(extrap_trend_hindcast_time_all_obs{obsi}(forecast_period_i,:),...
                extrap_trend_hindcast_all_obs{obsi}(forecast_period_i,:),...
                extrap_trend_rmse_all_obs{obsi}(forecast_period_i),...
                '-k','alpha','transparency',0.3)
             
            end
            
            %plot obs
            
%                 plot(ann_time_shift,Y_predictand_shift,'-k','LineWidth',2);
%                 h1 = scatter(ann_time_shift,Y_predictand_shift,50,...
%                      'MarkerFaceColor','k','MarkerEdgeColor','k',...
%                      'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);

                plot(ann_time,obs_GMST_ann{obsi},'-k','LineWidth',2);
                h1 = scatter(ann_time,obs_GMST_ann{obsi},50,...
                     'MarkerFaceColor','k','MarkerEdgeColor','k',...
                     'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
             
             
             % plot PLSR hindcasts
                 legend_handles_plsr_forecasts{forecast_period_i} = scatter(ann_time_shift,xval_and_crawling_forecast_predictions_all_obs{obsi,forecast_period_i},lead_sizes(forecast_period_i),...
                 'MarkerFaceColor',lead_colors{forecast_period_i},'MarkerEdgeColor',lead_colors{forecast_period_i},...
                 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
             
             %plot forecasts
                          
                %plot the error bars

                    plot([forecast_years(forecast_period_i) forecast_years(forecast_period_i)],...
                        [forecasts_foreward(forecast_period_i,obsi)-2*forecasts_foreward_RMSE(forecast_period_i,obsi) forecasts_foreward(forecast_period_i,obsi)+2*forecasts_foreward_RMSE(forecast_period_i,obsi)],...
                        strcat('-',lead_colors{forecast_period_i}),'LineWidth',0.5)

                    plot([forecast_years(forecast_period_i) forecast_years(forecast_period_i)],...
                        [forecasts_foreward(forecast_period_i,obsi)-1*forecasts_foreward_RMSE(forecast_period_i,obsi) forecasts_foreward(forecast_period_i,obsi)+1*forecasts_foreward_RMSE(forecast_period_i,obsi)],...
                        strcat('-',lead_colors{forecast_period_i}),'LineWidth',1)

             %plot the central values
                scatter(forecast_years(forecast_period_i),...
                    forecasts_foreward(forecast_period_i,obsi),...
                    lead_sizes(forecast_period_i),lead_colors{forecast_period_i},'filled')
             
        end
             
            xlim([min(ann_time) max(ann_time)+4])
            ylim([-0.7 1.3])
            
             %labels and such
                title(strcat(obs_names{obsi}))
                xlabel('Year')
                ylabel('Temperature Anomaly (C)');
                 
             % legend
                 h2 = legend_handles_plsr_forecasts{1};
                 h3 = legend_handles_plsr_forecasts{2};             
                 h4 = legend_handles_plsr_forecasts{3};
                 h5 = legend_handles_plsr_forecasts{4};
             
%              legend([h1 h2 h3 h4 h5],{'Observations',...
%              strcat('1 year lead, r=',...
%              num2str(round(prediction_Rs_all_obs{obsi,1}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_Rs_all_obs{obsi,1}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_r_all_obs{obsi}(1),2)),')',...
%              ', RMSE=',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,1}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,1}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_rmse_all_obs{obsi}(1),2)),')'),...
%              strcat('2 year lead, r=',...
%              num2str(round(prediction_Rs_all_obs{obsi,2}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_Rs_all_obs{obsi,2}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_r_all_obs{obsi}(2),2)),')',...
%              ', RMSE=',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,2}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,2}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_rmse_all_obs{obsi}(2),2)),')'),...
%              strcat('3 year lead, r=',...
%              num2str(round(prediction_Rs_all_obs{obsi,3}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_Rs_all_obs{obsi,3}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_r_all_obs{obsi}(3),2)),')',...
%              ', RMSE=',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,3}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,3}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_rmse_all_obs{obsi}(3),2)),')'),...
%              strcat('4 year lead, r=',...
%              num2str(round(prediction_Rs_all_obs{obsi,4}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_Rs_all_obs{obsi,4}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_r_all_obs{obsi}(4),2)),')',...
%              ', RMSE=',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,4}(display_hindcast_forecast_mode_1),2)),',',...
%              num2str(round(prediction_RMSEs_all_obs{obsi,4}(display_hindcast_forecast_mode_2),2)),...
%              ' (',...
%              num2str(round(persistence_rmse_all_obs{obsi}(4),2)),')')},...
%              'location','NorthWest')

             legend([h1 h2 h3 h4 h5],{'Observations',...
             strcat('1 year lead, RMSE=',...
             num2str(round(prediction_RMSEs_all_obs{obsi,1}(display_hindcast_forecast_mode_1),2)),',',...
             num2str(round(prediction_RMSEs_all_obs{obsi,1}(display_hindcast_forecast_mode_2),2)),...
             ' (',...
             num2str(round(persistence_rmse_all_obs{obsi}(1),2)),')'),...
             strcat('2 year lead, RMSE=',...
             num2str(round(prediction_RMSEs_all_obs{obsi,2}(display_hindcast_forecast_mode_1),2)),',',...
             num2str(round(prediction_RMSEs_all_obs{obsi,2}(display_hindcast_forecast_mode_2),2)),...
             ' (',...
             num2str(round(persistence_rmse_all_obs{obsi}(2),2)),')'),...
             strcat('3 year lead, RMSE=',...
             num2str(round(prediction_RMSEs_all_obs{obsi,3}(display_hindcast_forecast_mode_1),2)),',',...
             num2str(round(prediction_RMSEs_all_obs{obsi,3}(display_hindcast_forecast_mode_2),2)),...
             ' (',...
             num2str(round(persistence_rmse_all_obs{obsi}(3),2)),')'),...
             strcat('4 year lead, RMSE=',...
             num2str(round(prediction_RMSEs_all_obs{obsi,4}(display_hindcast_forecast_mode_1),2)),',',...
             num2str(round(prediction_RMSEs_all_obs{obsi,4}(display_hindcast_forecast_mode_2),2)),...
             ' (',...
             num2str(round(persistence_rmse_all_obs{obsi}(4),2)),')')},...
             'location','NorthWest')
             
end
    
%% plot RMSE plot

display_hindcast_forecast_mode = 2;

% forced CMIP5
    RMSEs_forced_CMIP5 = NaN(4,length(obs_names));
    
    for obsi = 1:length(obs_names)
        
        RMSEs_forced_CMIP5(:,obsi) = std(obs_GMST_ann_unforced{obsi});

    end

    %define alternatives
    
    alt_backcast_newman = [0.08595041 0.10661157 0.11322314 0.11652893]';
    alt_backcast_suckling_1 = [0.097879767 0.092800955 0.091863657 0.08293869]';
    alt_backcast_suckling_2 = [0.096400477 0.089842375 0.091863657 0.094477152]';
    alt_backcast_suckling_3 = [0.126577992 0.115877878 0.116715728 0.114595992]';
    alt_backcast_suckling_4 = [0.128648998 0.117948884 0.12085774 0.12080901]';
    alt_backcast_krueger_5 = [0.09483395 0.11343173 0.11933579 0.120369]';
    alt_backcast_laepple_6 = [0.10916823 0.112269 0.11456851 0.11526929]';
    
    alt_backcast_Sevellec_Drijfhout_ann_year_1_5 = 0.104;    
    
%     alt_backcasts = cat(2,alt_backcast_newman,...
%                           alt_backcast_suckling_1,...
%                           alt_backcast_suckling_2,...
%                           alt_backcast_suckling_3,...
%                           alt_backcast_suckling_4,...
%                           alt_backcast_krueger_5,...
%                           alt_backcast_laepple_6);
                      
    alt_backcasts = cat(2,alt_backcast_newman,...
                          alt_backcast_suckling_1,...
                          alt_backcast_suckling_3,...
                          alt_backcast_suckling_4,...
                          alt_backcast_krueger_5,...
                          alt_backcast_laepple_6);
                      
    benchmarks = NaN(4,4,3);
        for obsi = 1:length(obs_names)
            benchmarks(:,obsi,1) = persistence_rmse_all_obs{obsi}(:);
            benchmarks(:,obsi,2) = extrap_trend_rmse_all_obs{obsi}(:);
        end
            benchmarks(:,:,3) = RMSEs_forced_CMIP5;
                      
    %load the CMIP5 initilized w.r.t. each obs
    
     load('/Users/patrickbrown/Documents/MATLAB/PLS_decadal_prediction/interannual_SAT_GMST/CMIP5/MATS/bias_drift_correct_CMIP5_GMST_decadal_keep_ens_X_val.mat',...
     'CMIP5_decadal_backcasts_RMSE',...
     'CMIP5_decadal_backcasts_RMSE_ens',...
     'CMIP5_decadal_backcasts_RMSE_ens_X_val',...
     'model_names')
 
row_num = 2;
col_num = 2;

lin_dot_size = 55/2;
mean_dot_size_1 = 85/2;
mean_dot_size_2 = 40/2;
mean_dot_size_3 = 15/2;

FigHandle = figure('Position', [100, 100, 1100, 700]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName', 'helvetica')

    for obsi = 1:length(obs_names)
        
        subplot(row_num,col_num,obsi)
        hold on
        
        %plot the decadal hindcast experiments
        boundedline(1:4,nanmean(squeeze(CMIP5_decadal_backcasts_RMSE_ens_X_val(:,:,obsi)),2),nanstd(squeeze(CMIP5_decadal_backcasts_RMSE_ens_X_val(:,:,obsi)),[],2),'-r','alpha','transparency',0.1)
        
        for modeli = 1:length(model_names)
            plot(1:4,squeeze(CMIP5_decadal_backcasts_RMSE_ens_X_val(:,modeli,obsi)),'--r','LineWidth',0.1)
        end
        
           plot(1:4,nanmean(squeeze(CMIP5_decadal_backcasts_RMSE_ens_X_val(:,:,obsi)),2),'-r','LineWidth',3)
           scatter(1:4,nanmean(squeeze(CMIP5_decadal_backcasts_RMSE_ens_X_val(:,:,obsi)),2),lin_dot_size,'r','filled')
            
        %plot the benchmarks

        %boundedline(1:4,nanmean(squeeze(benchmarks(:,obsi,1)),2),nanstd(squeeze(benchmarks(:,obsi,1)),[],2),'-k','alpha','transparency',0.1)
        
        for bmarki = 1:2
            plot(1:4,squeeze(benchmarks(:,obsi,bmarki)),'-k','LineWidth',3)
            scatter(1:4,squeeze(benchmarks(:,obsi,bmarki)),lin_dot_size,'k','filled')
            %plot(1:4,squeeze(benchmarks(:,obsi,bmarki)),bmark_line_styles{bmarki},'LineWidth',2)
        end
           
            
        % plot the alternatives
%         boundedline(1:4,nanmean(alt_backcasts,2),nanstd(alt_backcasts,[],2),'-g','alpha','transparency',0.1)
%         
%         for alti = 1:size(alt_backcasts,2)
%             
%             plot(1:4,alt_backcasts(:,alti),'--g','LineWidth',0.1)
%             
%         end
%             plot(1:4,nanmean(alt_backcasts,2),'-g','LineWidth',3)        
%             scatter(1:4,nanmean(alt_backcasts,2),lin_dot_size,'g','filled')
            
        %plot PREDICTemp
        
        % prediction_RMSEs_all_obs{obsi,1}(display_hindcast_forecast_mode_1)
        
            Predict_RMSEs_all_obs_forecast_mode = NaN(4,1);
            Predict_RMSEs_all_obs_Xval_mode = NaN(4,1);
        
                for forecast_period_i = 1:4
                    
                    Predict_RMSEs_all_obs_Xval_mode(forecast_period_i) = prediction_RMSEs_all_obs{obsi,forecast_period_i}(1);                   
                    Predict_RMSEs_all_obs_forecast_mode(forecast_period_i) = prediction_RMSEs_all_obs{obsi,forecast_period_i}(2);

                end
                
            plot(1:4,Predict_RMSEs_all_obs_Xval_mode,'--b','LineWidth',3)
            scatter(1:4,Predict_RMSEs_all_obs_Xval_mode,lin_dot_size,'b','filled')         
            plot(1:4,Predict_RMSEs_all_obs_forecast_mode,'-b','LineWidth',3)
            scatter(1:4,Predict_RMSEs_all_obs_forecast_mode,lin_dot_size,'b','filled')            

        % --- plot the means --- 
        
            end_lead_to_mean_over = 4;
        
            %models
            
                for modeli = 1:length(model_names)
                    scatter(4.6,mean(CMIP5_decadal_backcasts_RMSE_ens_X_val(1:end_lead_to_mean_over,modeli,obsi),1),mean_dot_size_3,'r','filled')
                end
                
                    scatter(4.6,mean(mean(CMIP5_decadal_backcasts_RMSE_ens_X_val(1:end_lead_to_mean_over,:,obsi),1),2),mean_dot_size_1,'r','s','filled')
                    
            %benchmarks

                for bmarki = 1:2
                    scatter(4.4,mean(benchmarks(1:end_lead_to_mean_over,obsi,bmarki),1),mean_dot_size_3,'k','o','filled')
                end
                for bmarki = 1:2
                    scatter(4.4,mean(mean(benchmarks(1:end_lead_to_mean_over,obsi,bmarki),1),3),mean_dot_size_1,'k','s','filled')
                end
                
            %alternatives
                        
%                 for alti = 1:size(alt_backcasts,2)
%                     scatter(4.8,mean(alt_backcasts(:,alti),1),mean_dot_size_3,'g','o','filled')
%                 end           
%                     scatter(4.8,mean(mean(alt_backcasts(:,alti),1),2),mean_dot_size_1,'g','s','filled')
%                     
%                     scatter(4.8,alt_backcast_Sevellec_Drijfhout_ann_year_1_5,20,'m','s','filled')                    
            
            %PREDICTemp
            
                scatter(4.2,mean(Predict_RMSEs_all_obs_forecast_mode),mean_dot_size_1,'b','s','filled')
                scatter(4.2,mean(Predict_RMSEs_all_obs_Xval_mode),mean_dot_size_1,'b','s','filled')

        xlim([1 5.1])
        ylim([0.05 0.3])

%         if obsi == 3 || obsi == 4; xlabel('Lead time (years)'); end
%         if obsi == 1 || obsi == 3; ylabel('Root Mean Square Error (C)'); end
                                     xlabel('Lead time (years)');
                                     ylabel('Root Mean Square Error (C)');

        title(obs_names{obsi})
        
    end
    
%% make the map plots

lags_to_plot = 2;
PLS_comps_to_plot = 3;
row_num = PLS_comps_to_plot+1;
col_num = lags_to_plot;

% 1 2
% 3 4
% 5 6
% 7 8

panel_order = [2 1 4 3 6 5 8 7];

%output_file_letters = {'a','b','c','d','e','f','g','h','i'};
output_file_letters = {'a','b','c','d','e','f','g','h'};

    coast = load('coast');
    
%for obsi = 1:length(obs_names)
obsi = 2;
    
        obs_lat{obsi} = double(obs_lat{obsi});
        obs_lon{obsi} = double(obs_lon{obsi});

        FigHandle = figure('Position', [100, 100, 500, 900]); %[left bottom width height]
        %figure
        set(gcf,'color',[1 1 1]);
        set(0, 'DefaultAxesFontSize',12);
        set(0,'defaultAxesFontName', 'helvetica')
        colormap(redblue)

        lin_panel_count = 1;

        for pls_compsi = 1:PLS_comps_to_plot
            for lagi = 1:lags_to_plot

            subplot(row_num,col_num,panel_order(lin_panel_count))
            hold on

                axesm('robinson',...
                'Frame', 'on',...
                'Grid', 'off',...
                'maplatlim',[min(obs_lat{obsi}) max(obs_lat{obsi})],...
                'maplonlim',[min(obs_lon{obsi}) max(obs_lon{obsi})])
                 tightmap

                plot_map_vector = reshape(squeeze(XLoadings_maps_all_obs{obsi}(:,:,lagi,pls_compsi)),[length(obs_lat{obsi})*length(obs_lon{obsi}),1]);
                max_color_val = prctile(abs(plot_map_vector),95);

%                 min_color_range = -1.2*max_color_val;
%                 max_color_range = 1.2*max_color_val;
                 min_color_range = -7;
                 max_color_range = 7;

                num_contours = 20;

                contourfm(obs_lat{obsi},obs_lon{obsi},squeeze(XLoadings_maps_all_obs{obsi}(:,:,lagi,pls_compsi)),linspace(min_color_range,max_color_range,num_contours), 'LineStyle','none');
                caxis([min_color_range max_color_range])

                geoshow(coast.lat, coast.long, 'Color', 'black')

%                 t=contourcbar;
%                 set(get(t,'ylabel'),'String', 'loadings');
                
                if lagi > 1; title(strcat(num2str(lagi),' Years prior, PLSR # =',num2str(pls_compsi))); end
                if lagi == 1; title(strcat(num2str(lagi),' Year prior, PLSR # = ',num2str(pls_compsi),', VarEx=',num2str(round(100*frac_var_explained(pls_compsi,forecast_periodi,obsi),0)),'%')); end
                
                %write numbers to file
                
%                 csvwrite(strcat('/home/pbrown/Documents/MATLAB/PLS_decadal_prediction/interannual_SAT_GMST/fig_data/BC2018_Fig5_raw',...
%                     output_file_letters{lin_panel_count},'.dat'),squeeze(XLoadings_maps_all_obs{obsi}(:,:,lagi,pls_compsi)))
                
                
                lin_panel_count = lin_panel_count + 1;
                
            end
        end
        
        %plot the obs predictor forecast row
        
        forecast_lead_years = [2019 2018];
        
        for lagi = 1:lags_to_plot
            
            subplot(row_num,col_num,panel_order(lin_panel_count))
            hold on
            
                axesm('robinson',...
                'Frame', 'on',...
                'Grid', 'off',...
                'maplatlim',[min(obs_lat{obsi}) max(obs_lat{obsi})],...
                'maplonlim',[min(obs_lon{obsi}) max(obs_lon{obsi})])
                 tightmap
                 
                plot_map_vector = reshape(squeeze(obs_3d_ann_anom{obsi}(:,:,end-lagi+1)),[length(obs_lat{obsi})*length(obs_lon{obsi}),1]);
                max_color_val = prctile(abs(plot_map_vector),95);

%                 min_color_range = -1.2*max_color_val;
%                 max_color_range = 1.2*max_color_val;
                 min_color_range = -2.5;
                 max_color_range = 2.5;
                num_contours = 20;

                contourfm(obs_lat{obsi},obs_lon{obsi},squeeze(obs_3d_ann_anom{obsi}(:,:,end-lagi+1)),linspace(min_color_range,max_color_range,num_contours), 'LineStyle','none');
                caxis([min_color_range max_color_range])

                geoshow(coast.lat, coast.long, 'Color', 'black')

%                 t=contourcbar;
%                 set(get(t,'ylabel'),'String', 'Unforced SAT Anom (K)');

                title(strcat(num2str(forecast_lead_years(lagi)),' mean=',...
                    num2str(round(threed2oned_lat_lon_2(squeeze(obs_3d_ann_unforced{obsi}(:,:,end-lagi+1)),obs_lat{obsi},obs_lon{obsi}),2))));
                
                %write numbers to file
                
%                 csvwrite(strcat('/home/pbrown/Documents/MATLAB/PLS_decadal_prediction/interannual_SAT_GMST/fig_data/BC2018_Fig5_raw',...
%                     output_file_letters{lin_panel_count},'.dat'),squeeze(obs_3d_ann_unforced{obsi}(:,:,end-lagi+1)))               

                lin_panel_count = lin_panel_count + 1;

        end
    
%% plot the scores (like the principle component time series)

PLS_comp_colors = {'k','k','k'};

row_num = 3;
col_num = 1;

FigHandle = figure('Position', [100, 100, 600, 800]); %[left bottom width height]
set(gcf,'color',[1 1 1]);
set(0, 'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName', 'helvetica')

obsi = 2;
        
        %plot(obs_ann_time{obsi},obs_GMST_ann{obsi},'-k','LineWidth',2)
        
        for pls_compsi = 1:PLS_comps_to_plot
            
            subplot(row_num,col_num,pls_compsi)
            hold on
            
            plot(obs_ann_time{obsi}(num_lag_predictor+1:end),XScores_t_series_all_obs{obsi}(:,pls_compsi),strcat('-',PLS_comp_colors{pls_compsi}),'LineWidth',3)


            if pls_compsi == 3; xlabel('Year'); end
            ylabel('Scores (unitless)');
            
            if pls_compsi == 1; title(obs_names{obsi}); end

            xlim([start_year, 2024])
            %ylim([-0.35, 1.3])
            
        end
            
            %legend('PLSR Component 1','PLSR Component 2','PLSR Component 3','Location','NorthWest')
        
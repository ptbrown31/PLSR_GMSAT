function [persistence_rmse,...
          persistence_r,...
          persistence_hindcast,...
          persistence_hindcast_time,...
          persistence_hindcast_errs] = persistence_hindcast_errors(full_time,...
          predictand_t_series,...
          persistence_average_length,...
          forecast_periods)

%% load practice data

%         close all
%         clear all
%         
%         load('/Users/patrickbrown/Documents/MATLAB/PLS_decadal_prediction/interannual_SAT_GMST/CMIP5/MATS/hindcast_forecast_figs_all_obs_preproc_MLR.mat',...
%             'land_frac_obs',...
%             'obs_ann_time',...
%             'obs_lon',...
%             'obs_lat',...
%             'forced_GMST_obs',...
%             'obs_GMST_ann',...
%             'obs_GMST_ann_unforced',...
%             'forced_3d_ann_obs',...
%             'obs_3d_ann_anom',...
%             'obs_3d_ann_unforced',...
%             'forced_GMST_obs_plus_forecast',...
%             'obs_ann_time_plus_forecast')
%         
%         %NOAA needs to be flipped
%         
%         obs_names = {'BEST',...
%                      'GISTEMP',...
%                      'NOAA',...
%                      'HadCRUT4'};
%         
%         full_time = obs_ann_time{1};
%         predictand_t_series = obs_GMST_ann{1};
%         
%         persistence_average_length = 5;
%         forecast_periods = 1:4;
%         forecast_block_size = 1;

%% 1st just smooth series to create persistence forecast

    whole_predictand_smoothed = smooth(predictand_t_series,persistence_average_length);

    persistence_rmse = NaN(length(forecast_periods),1);
    persistence_r = NaN(length(forecast_periods),1);

    persistence_hindcast = {};
    persistence_hindcast_errs = {};
    persistence_hindcast_time = {};

    for forecast_period_i = 1:length(forecast_periods)

        forecast_period = forecast_periods(forecast_period_i);

        % Select forecast and validation indicies
            forecast_validation_inds = persistence_average_length + forecast_period_i : length(predictand_t_series);
            persistence_forecast_inds = 1:length(forecast_validation_inds);

        % Pull series values

            forecast_validation = predictand_t_series(forecast_validation_inds);

            persistence_hindcast{forecast_period_i} = whole_predictand_smoothed(persistence_forecast_inds);
            persistence_hindcast_errs{forecast_period_i} = whole_predictand_smoothed(persistence_forecast_inds) - forecast_validation;

            persistence_hindcast_time{forecast_period_i} = full_time(forecast_validation_inds);

        %calculate error metrics

            mean_square_errors = mean(persistence_hindcast_errs{forecast_period_i}.^2);

            persistence_rmse(forecast_period_i) = sqrt(mean_square_errors);

            % correlation

            corr_mat = corrcoef(whole_predictand_smoothed(persistence_forecast_inds),...
                                forecast_validation,...
                                'rows',...
                                'pairwise');

           persistence_r(forecast_period_i) = corr_mat(2,1);

           %test plot

%            figure
%            hold on
%     
%             plot(full_time,predictand_t_series,'-k','LineWidth',3)
%             plot(persistence_hindcast_time{forecast_period_i},forecast_validation,'-r','LineWidth',3)
%     
%             plot(persistence_hindcast_time{forecast_period_i},persistence_hindcast{forecast_period_i},'-b','LineWidth',3)
%             
%             title(strcat('forecast period = ',num2str(forecast_period),...
%                          'r = ',num2str(round(persistence_r(forecast_period_i),2)),...
%                          'rmse = ',num2str(round(persistence_rmse(forecast_period_i),2))))
%     
     end
end
 


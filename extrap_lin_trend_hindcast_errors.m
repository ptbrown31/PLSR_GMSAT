function [extrap_lin_trend_rmse,...
          extrap_lin_trend_r,...
          extrap_lin_trend_hindcast,...
          extrap_lin_trend_hindcast_time,...
          extrap_lin_trend_hindcast_errs] = extrap_lin_trend_hindcast_errors(full_time,...
          predictand_t_series,...
          trend_fit_length,...
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
%         trend_fit_length = 10;
%         forecast_periods = 1:4;
%         forecast_block_size = 1;

%% 1st just smooth series to create persistence forecast

    foreward_years_trend = max(forecast_periods);

    extrap_lin_trend_hindcast = [];
    extrap_lin_trend_predictand = [];
    extrap_lin_trend_hindcast_errs = [];
    extrap_lin_trend_hindcast_time = [];
    
%     figure
%     hold on
%     plot(full_time,predictand_t_series,'-k')

    %for year_i = 1+trend_fit_length:1:length(predictand_t_series)-length(forecast_periods)
    for year_i = 1+trend_fit_length:1:length(predictand_t_series)
        
            %pull out years to fit
                fit_years_correct_time_inds = year_i-trend_fit_length:year_i;
                fit_years_predictand = predictand_t_series(fit_years_correct_time_inds);
                %fit lm
                    mdl = fitlm(full_time(fit_years_correct_time_inds),fit_years_predictand,'linear');
                    %y = aX+b;
                    %y = slope*X+y_int
                    y_int = mdl.Coefficients.Estimate(1);
                    slope = mdl.Coefficients.Estimate(2);
                    forecast_now = NaN(foreward_years_trend,1);
                    forecast_now_time = NaN(foreward_years_trend,1);
                    predictand_now = NaN(foreward_years_trend,1);
                    errors_now = NaN(foreward_years_trend,1);
                    
                    if year_i+length(forecast_periods) <= length(predictand_t_series) %if you are not running out of room
                        for forecast_yearsi = 1:foreward_years_trend
                            xval = full_time(fit_years_correct_time_inds(end))+forecast_yearsi;
                            forecast_now(forecast_yearsi,1) = slope*xval + y_int;
                            predictand_now(forecast_yearsi,1) = predictand_t_series(fit_years_correct_time_inds(end)+forecast_yearsi);
                            errors_now(forecast_yearsi,1) = forecast_now(forecast_yearsi,1) - predictand_t_series(fit_years_correct_time_inds(end)+forecast_yearsi);
                        end
                        forecast_years_correct_time_inds = fit_years_correct_time_inds(end)+1:fit_years_correct_time_inds(end)+foreward_years_trend;
                        forecast_now_time(:,1) = full_time(forecast_years_correct_time_inds);
                    end
                    
                    if year_i+length(forecast_periods) > length(predictand_t_series) %if you are are running out of room
                        diff = year_i+length(forecast_periods) - length(predictand_t_series);
                        for forecast_yearsi = 1:foreward_years_trend-diff
                            xval = full_time(fit_years_correct_time_inds(end))+forecast_yearsi;
                            forecast_now(forecast_yearsi,1) = slope*xval + y_int;
                            predictand_now(forecast_yearsi,1) = predictand_t_series(fit_years_correct_time_inds(end)+forecast_yearsi);
                            errors_now(forecast_yearsi,1) = forecast_now(forecast_yearsi,1) - predictand_t_series(fit_years_correct_time_inds(end)+forecast_yearsi);
                        end
                        forecast_years_correct_time_inds = fit_years_correct_time_inds(end)+1:fit_years_correct_time_inds(end)+foreward_years_trend-diff;
                        forecast_now_time(1:length(forecast_years_correct_time_inds),1) = full_time(forecast_years_correct_time_inds);
                    end                    
                    
%                     plot(full_time(forecast_years_correct_time_inds),forecast_now(~isnan(forecast_now)),'-r')
%                     scatter(full_time(forecast_years_correct_time_inds),forecast_now(~isnan(forecast_now)),10,'r','filled')
%                     drawnow
                    
                    extrap_lin_trend_predictand = cat(2,extrap_lin_trend_predictand,predictand_now);
                    
                    extrap_lin_trend_hindcast = cat(2,extrap_lin_trend_hindcast,forecast_now);
                    extrap_lin_trend_hindcast_errs = cat(2,extrap_lin_trend_hindcast_errs,errors_now);
                    
                    
                    extrap_lin_trend_hindcast_time = cat(2,extrap_lin_trend_hindcast_time,forecast_now_time);
                                        
                        
    end
    
    extrap_lin_trend_se = extrap_lin_trend_hindcast_errs.^2;
    extrap_lin_trend_mse = nanmean(extrap_lin_trend_se,2);
    extrap_lin_trend_rmse = sqrt(extrap_lin_trend_mse);

    extrap_lin_trend_r = NaN(length(forecast_periods),1);
    
    for forecast_period_i = 1:length(forecast_periods)
        
        predictions_this_lead = extrap_lin_trend_hindcast(forecast_period_i,:);
        actuals_this_lead = extrap_lin_trend_predictand(forecast_period_i,:);
        
        corr_mat = corrcoef(predictions_this_lead,...
                    actuals_this_lead,...
                    'rows',...
                    'pairwise');

        extrap_lin_trend_r(forecast_period_i) = corr_mat(2,1);

    end
end


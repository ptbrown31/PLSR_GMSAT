function [XLoadings,...
    YLoadings,...
    XScores,...
    YScores,...
    BETA,...
    PCTVAR,...
    stats,...
    held_out_based_Y_prediction_1D_space_norm,...
    XLoadings_maps,...
    YLoadings_maps,...
    Y_forecast_2D_space_unnorm,...
    held_out_based_prediction_2D_space_unnorm,...
    cross_val_RMSE_map,...
    cross_val_corr_map,...
    cross_val_RMSSS_map] = PLSR_return_cross_val_backcasts_no_shift_Y2D_alpha(X_predictor_field,...
    X_predictor_field_for_forecast,...
    Y_predictand_field,...
    num_pls_comps)

%the predictor depth could be the number of lags or the number of different
%predictor variables or a combination

predictor_depth = size(X_predictor_field,4);
predictand_depth = size(Y_predictand_field,4); 

% fields should be 2D in spcae when they come into this function

%% normalize X_field Y_field

    X_predictor_field_means = NaN(size(X_predictor_field,1),size(X_predictor_field,2),1,predictor_depth);
    X_predictor_field_STDs = NaN(size(X_predictor_field,1),size(X_predictor_field,2),1,predictor_depth);
    X_predictor_field_norm = zeros(size(X_predictor_field));
    X_predictor_field_for_forecast_norm = zeros(size(X_predictor_field_for_forecast));

    for depthi = 1:predictor_depth
        for lati = 1:size(X_predictor_field,1)
            for loni = 1:size(X_predictor_field,2)

                    X_predictor_field_means(lati,loni,1,depthi) = mean(squeeze(X_predictor_field(lati,loni,:,depthi)));
                    X_predictor_field_STDs(lati,loni,1,depthi) = std(squeeze(X_predictor_field(lati,loni,:,depthi)));
                    
                    if X_predictor_field_STDs(lati,loni,1,depthi) ~= 0 %dont divide by zero and create a NaN if there is no data there (ice has this problem)

                        X_predictor_field_norm(lati,loni,:,depthi) = (squeeze(X_predictor_field(lati,loni,:,depthi)) - squeeze(X_predictor_field_means(lati,loni,1,depthi)))./X_predictor_field_STDs(lati,loni,1,depthi);
                        X_predictor_field_for_forecast_norm(lati,loni,1,depthi) = (squeeze(X_predictor_field_for_forecast(lati,loni,1,depthi)) - squeeze(X_predictor_field_means(lati,loni,1,depthi)))./X_predictor_field_STDs(lati,loni,1,depthi);
                        
                    end
            end
        end
    end

% normalize Y_field

    Y_predictand_field_means = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),1,predictand_depth);
    Y_predictand_field_STDs = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),1,predictand_depth);
    Y_predictand_field_norm = zeros(size(Y_predictand_field));

    for depthi = 1:predictand_depth
        for lati = 1:size(Y_predictand_field,1)
            for loni = 1:size(Y_predictand_field,2)

                    Y_predictand_field_means(lati,loni,1,depthi) = mean(squeeze(Y_predictand_field(lati,loni,:,depthi)));
                    Y_predictand_field_STDs(lati,loni,1,depthi) = std(squeeze(Y_predictand_field(lati,loni,:,depthi)));
                    
                    Y_predictand_field_norm(lati,loni,:,depthi) = (squeeze(Y_predictand_field(lati,loni,:,depthi)) - squeeze(Y_predictand_field_means(lati,loni,1,depthi)))./Y_predictand_field_STDs(lati,loni,1,depthi);

            end
        end
    end
    

%% make the X and Y fields 1D in space

    X_predictor_1D_space = NaN(size(X_predictor_field,3),size(X_predictor_field,1)*size(X_predictor_field,2)*predictor_depth);
    X_predictor_1D_space_forecast = NaN(1,size(X_predictor_field,1)*size(X_predictor_field,2)*predictor_depth);

    lin_count = 1;
    
    for depthi = 1:predictor_depth
            for lati = 1:size(X_predictor_field,1)
                for loni = 1:size(X_predictor_field,2)

                    X_predictor_1D_space(:,lin_count) = X_predictor_field_norm(lati,loni,:,depthi);
                    X_predictor_1D_space_forecast(lin_count) = X_predictor_field_for_forecast_norm(lati,loni,1,depthi);

                    lin_count = lin_count + 1;
                    
                end
            end
    end
    

% make the Y fields 1D in space

    Y_predictand_1D_space = NaN(size(Y_predictand_field,3),size(Y_predictand_field,1)*size(Y_predictand_field,2)*predictand_depth);

    lin_count = 1;
    
    for depthi = 1:predictand_depth
            for lati = 1:size(Y_predictand_field,1)
                for loni = 1:size(Y_predictand_field,2)

                    Y_predictand_1D_space(:,lin_count) = Y_predictand_field_norm(lati,loni,:,depthi);

                    lin_count = lin_count + 1;
                    
                end
            end
    end
    
%% go through and get rid of columns without full data

    % for X

        X_predictor_1D_space_pruned = [];
        X_predictor_1D_space_pruned_forecast = [];

        X_predictor_nan_column_inds = [];

        for columni = 1:size(X_predictor_1D_space,2)
            
            column_data = X_predictor_1D_space(:,columni);
            
            if sum(isnan(column_data)) ~= 0
                
                X_predictor_nan_column_inds = horzcat(X_predictor_nan_column_inds,columni);
                
            end
            
            if sum(isnan(column_data)) == 0
                
                X_predictor_1D_space_pruned = horzcat(X_predictor_1D_space_pruned,column_data);
                X_predictor_1D_space_pruned_forecast = horzcat(X_predictor_1D_space_pruned_forecast,X_predictor_1D_space_forecast(columni));
                
            end
        end

    % for Y

        Y_predictand_1D_space_pruned = [];
        Y_predictand_nan_column_inds = [];

        for columni = 1:size(Y_predictand_1D_space,2)
            
            column_data = Y_predictand_1D_space(:,columni);
            
            if sum(isnan(column_data)) ~= 0
                
                Y_predictand_nan_column_inds = horzcat(Y_predictand_nan_column_inds,columni);
                
            end
            
            if sum(isnan(column_data)) == 0
                
                Y_predictand_1D_space_pruned = horzcat(Y_predictand_1D_space_pruned,column_data);
                
            end
        end

%% do the cross validation/forecast part manually 

%%define empty arrays that dont come directly from the plsregress function

held_out_based_Y_prediction_1D_space_norm = NaN(size(Y_predictand_1D_space_pruned));

% at this point you have a predictor and predictand that are lined up and are the same size
% (in the important dimension)

    for held_out_i = 1:size(Y_predictand_field,3) %looping through years to hold out
        
        held_out_i/size(Y_predictand_field,3)

        %Seperate out 1 thing from the predictand (to test the prediction) from the rest of the
        %(to inform the prediction)

            indicies_informing_prediction = find((1:size(Y_predictand_field,3) ~= held_out_i));

            X_informing_prediction = X_predictor_1D_space_pruned(indicies_informing_prediction,:);
            Y_informing_prediction = Y_predictand_1D_space_pruned(indicies_informing_prediction,:);

            X_testing_prediction = X_predictor_1D_space_pruned(held_out_i,:);
            
            %do PLS regression
            
            [XLoadings_hold_out,YLoadings_hold_out,XScores_hold_out,YScores_hold_out,BETA_hold_out,PCTVAR_hold_out,MSE_hold_out] = plsregress(X_informing_prediction,Y_informing_prediction,num_pls_comps);
                    
             % predict based on held-out predictor
                    
             held_out_based_Y_prediction_1D_space_norm(held_out_i,:) = [1 X_testing_prediction]*BETA_hold_out;

    end
    
%% do regular PLS part using all data to make a forecast

% K_fold_cross_val = size(Y_predictand_1D_space_pruned,1);
% mc_num = 1;

% K_fold_cross_val = 4;
% mc_num = 100;

%size(X_predictor_1D_space_pruned)

%[XLoadings,YLoadings,XScores,YScores,BETA,PCTVAR,MSE,stats] = plsregress(X_predictor_1D_space_pruned,Y_predictand_1D_space_pruned,num_pls_comps,'cv',K_fold_cross_val,'mcreps',mc_num);
[XLoadings,YLoadings,XScores,YScores,BETA,PCTVAR,MSE,stats] = plsregress(X_predictor_1D_space_pruned,Y_predictand_1D_space_pruned,num_pls_comps);

%make prediction for new obs data
Y_forecast_1D_space_norm = [1 X_predictor_1D_space_pruned_forecast]*BETA;

%% rewrap X and Y Loading maps

    XLoadings_maps = NaN(size(X_predictor_field,1),size(X_predictor_field,2),predictor_depth,num_pls_comps);

    for ncompsi = 1:num_pls_comps
        
        lin_count = 1;
        
        %certain lat-lon-lag combos dont have vectors (bc they had NaNs)
        %these NaNs need to be put back in 
        
        for depthi = 1:predictor_depth
            for lati = 1:size(X_predictor_field,1)
                for loni = 1:size(X_predictor_field,2)
                    
                    %see if this time series had an NaN
                    
                    t_series_for_this_loc_depth_combo = squeeze(X_predictor_field(lati,loni,:,depthi));
                    
                    if sum(isnan(t_series_for_this_loc_depth_combo)) == 0 %if there are NO NaNs move over the data and increment
                        
                        XLoadings_maps(lati,loni,depthi,ncompsi) = XLoadings(lin_count,ncompsi);
                        
                        lin_count = lin_count + 1;
                        
                    end
                    
                end
            end
        end
    end

% rewrap Y loading maps

    YLoadings_maps = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),predictand_depth,num_pls_comps);

    for ncompsi = 1:num_pls_comps
        
        lin_count = 1;
        
        %certain lat-lon-lag combos dont have vectors (bc they had NaNs)
        %these NaNs need to be put back in 
        
        for depthi = 1:predictand_depth
            for lati = 1:size(Y_predictand_field,1)
                for loni = 1:size(Y_predictand_field,2)
                    
                    %see if this time series had an NaN
                    
                    t_series_for_this_loc_depth_combo = squeeze(Y_predictand_field(lati,loni,:,depthi));
                    
                    if sum(isnan(t_series_for_this_loc_depth_combo)) == 0 && any(t_series_for_this_loc_depth_combo) == 1 %if there are NO NaNs move over the data and increment
                        
                        YLoadings_maps(lati,loni,depthi,ncompsi) = YLoadings(lin_count,ncompsi);
                        
                        lin_count = lin_count + 1;
                        
                    end
                end
            end
        end
    end
    
%% rewrap the X-val Y predictions

        held_out_based_prediction_2D_space_norm = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),size(Y_predictand_field,3),predictand_depth);
        Y_forecast_2D_space_norm = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),1,predictand_depth);
        
        
        lin_count = 1;
        
        %certain lat-lon-lag combos dont have vectors (bc they had NaNs)
        %these NaNs need to be put back in 
        
        for depthi = 1:predictand_depth
            for lati = 1:size(Y_predictand_field,1)
                for loni = 1:size(Y_predictand_field,2)
                    
                    %see if this time series had an NaN
                    
                    t_series_for_this_loc_depth_combo = squeeze(Y_predictand_field(lati,loni,:,depthi));
                    
                    if sum(isnan(t_series_for_this_loc_depth_combo)) == 0 && any(t_series_for_this_loc_depth_combo) == 1 %if there are NO NaNs move over the data and increment
                        
                        held_out_based_prediction_2D_space_norm(lati,loni,:,depthi) = held_out_based_Y_prediction_1D_space_norm(:,lin_count);
                        Y_forecast_2D_space_norm(lati,loni,1,depthi) = Y_forecast_1D_space_norm(lin_count);
                        
                        lin_count = lin_count + 1;
                        
                    end
                end
            end
        end

%% un-norm the Y hindcast errors in order to compute RMSE

held_out_based_prediction_2D_space_unnorm = NaN(size(held_out_based_prediction_2D_space_norm));
Y_forecast_2D_space_unnorm = NaN(size(Y_forecast_2D_space_norm));

    for depthi = 1:predictand_depth
        for lati = 1:size(Y_predictand_field,1)
            for loni = 1:size(Y_predictand_field,2)

                    % see if this time series had an NaN
                    
                    t_series_for_this_loc_depth_combo = squeeze(Y_predictand_field(lati,loni,:,depthi));
                    
                    if sum(isnan(t_series_for_this_loc_depth_combo)) == 0  && any(t_series_for_this_loc_depth_combo) == 1 %if there are NO NaNs move over the data and increment

                        held_out_based_prediction_2D_space_unnorm(lati,loni,:,depthi) = held_out_based_prediction_2D_space_norm(lati,loni,:,depthi).*Y_predictand_field_STDs(lati,loni,1,depthi) + Y_predictand_field_means(lati,loni,1,depthi);

                        Y_forecast_2D_space_unnorm(lati,loni,1,depthi) = Y_forecast_2D_space_norm(lati,loni,1,depthi).*Y_predictand_field_STDs(lati,loni,1,depthi) + Y_predictand_field_means(lati,loni,1,depthi);
                        
                    end
            end
        end
    end
       
%% calulcate the error metrics

cross_val_RMSE_map = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),size(Y_predictand_field,4));
cross_val_corr_map = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),size(Y_predictand_field,4));
cross_val_RMSSS_map = NaN(size(Y_predictand_field,1),size(Y_predictand_field,2),size(Y_predictand_field,4));

    for depthi = 1:predictand_depth
        for lati = 1:size(Y_predictand_field,1)
            for loni = 1:size(Y_predictand_field,2)

                    % see if this time series had an NaN
                    
                    t_series_for_this_loc_depth_combo = squeeze(Y_predictand_field(lati,loni,:,depthi));
                    
                    if sum(isnan(t_series_for_this_loc_depth_combo)) == 0 && any(t_series_for_this_loc_depth_combo) == 1 %if there are NO NaNs move over the data and increment

                        hindcast_t_series = squeeze(held_out_based_prediction_2D_space_unnorm(lati,loni,:,depthi));
                        obs_t_series = t_series_for_this_loc_depth_combo;
                        
                        %calculations
                        
                                clim_prediction_RMSE = sqrt(mean((obs_t_series-mean(obs_t_series)).^2));
                                
                                PLSR_prediction_RMSE = sqrt(mean((hindcast_t_series-obs_t_series).^2));

                                rmsss = 1 - PLSR_prediction_RMSE./clim_prediction_RMSE;

                                corr_mat = corrcoef(hindcast_t_series,...
                                                    obs_t_series,...
                                                    'rows',...
                                                    'pairwise');

                                correlations = corr_mat(2,1);
                                
                                cross_val_RMSE_map(lati,loni,depthi) = PLSR_prediction_RMSE;
                                cross_val_corr_map(lati,loni,depthi) = correlations;
                                cross_val_RMSSS_map(lati,loni,depthi) = rmsss;


                    end
            end
        end
    end


end

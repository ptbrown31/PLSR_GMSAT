function [xval_and_crawling_forecast_predictions,...
          X_val_and_crawling_climatology,...
          prediction_RMSEs,...
          prediction_Rs,...
          prediction_RMSSSs] = ...
          PLSR_field_with_local_lasso_crawling_forecast(...
          Y_predictand_2D,...
          X_predictor_matrix,...
          X_predictor_fields,...
          crawling_forecast_ind,...
          lasso_alpha,...
          kval_number,...
          mc_trials,...
          num_pls_comps)
      
      Y_predictand = squeeze(Y_predictand_2D(1,1,:));

    %This function takes in
    
        %1) a predictand
        %2) a 2d matrix of predictors (time,var)
        %3) a 4d feild of predictors (lat,lon,time,var)
        
        %4) an index which represents the point in time for which the
        %crawling out-of-sample forecasts are performed
        
        %5) the number of PCAs to retain for the 4d feilds
        
        %6) alpha value for lasso
        
    %this function returns
    
        %1) predictions that are cross validated prior to the specified
        %index and out-of-sample forecasts after the cross validated index
        

%% make empty arrays
        
    xval_and_crawling_forecast_predictions = NaN(length(Y_predictand),1);
    X_val_and_crawling_climatology = NaN(length(Y_predictand),1);
    
    prediction_RMSEs = NaN(2,1); %2 elements, 1 for X val and 1 for crawling predictions
    prediction_Rs = NaN(2,1); %2 elements, 1 for X val and 1 for crawling predictions
    prediction_RMSSSs = NaN(2,1); %2 elements, 1 for X val and 1 for crawling predictions 
    
%% first thing - split into an X-val portion and a Crawling out-of-sample forecast portion

%xval portion
    X_predictor_matrix_xval = X_predictor_matrix(1:crawling_forecast_ind-1,:);
    X_predictor_fields_xval = X_predictor_fields(:,:,1:crawling_forecast_ind-1,:);
    
    Y_predictand_xval = Y_predictand(1:crawling_forecast_ind-1);
    Y_predictand_2D_xval = Y_predictand_2D(1,1,1:crawling_forecast_ind-1);
    
%Crawling forecast portion
    X_predictor_matrix_crawling_forecast = X_predictor_matrix(crawling_forecast_ind:end,:);
    X_predictor_fields_crawling_forecast = X_predictor_fields(:,:,crawling_forecast_ind:end,:);
    
    Y_predictand_2D_crawling_forecast = Y_predictand_2D(1,1,crawling_forecast_ind:end);
    Y_predictand_crawling_forecast = Y_predictand(crawling_forecast_ind:end);
    
%% do PLSR on xval portion

   [PLSR_XLoadings_xval,...
    PLSR_YLoadings_xval,...
    PLSR_XScores_xval,...
    PLSR_YScores_xval,...
    PLSR_BETA_xval,...
    PLSR_PCTVAR_xval,...
    PLSR_stats_xval,...
    PLSR_held_out_based_Y_prediction_1D_space_norm_xval,...
    PLSR_XLoadings_maps_xval,...
    PLSR_YLoadings_maps_xval,...
    PLSR_Y_forecast_2D_space_unnorm_xval,...
    PLSR_held_out_based_prediction_2D_space_unnorm_xval,...
    PLSR_cross_val_RMSE_map_xval,...
    PLSR_cross_val_corr_map_xval,...
    PLSR_cross_val_RMSSS_map_xval] = PLSR_return_cross_val_backcasts_no_shift_Y2D_alpha(...
    X_predictor_fields_xval,...
    X_predictor_fields_xval(:,:,end,:),...
    Y_predictand_2D_xval,...
    num_pls_comps);
    
%% prepare the lasso X val portion
    
    %Assemble large Predictor Matrix
    
    combined_X_predictor_matrix_xval = cat(2,X_predictor_matrix_xval,squeeze(PLSR_held_out_based_prediction_2D_space_unnorm_xval(1:length(X_predictor_matrix_xval))));
    
    %Norm large Predictor Matrix
    
    combined_X_predictor_matrix_xval_norm = normalize(combined_X_predictor_matrix_xval,1,'zscore');

    %% do lasso predictor selection in X val
    
        %https://www.mathworks.com/help/stats/lasso.html#mw_2715ae13-55dd-409c-aed8-49394e9c1964

        % Find the coefficients of a regularized linear regression model using 10-fold 
        % cross-validation and the elastic net method with Alpha = 0.75. Use the largest 
        % Lambda value such that the mean squared error (MSE) is within one standard error 
        % of the minimum MSE.

        for held_out_i = 1:length(Y_predictand_xval) %looping through models to hold out

            %Seperate out 1 thing from the predictand (to test the prediction) from the rest of the
            %(to inform the prediction)

                indicies_informing_xval_prediction = find((1:length(Y_predictand_xval) ~= held_out_i));

                X_informing_xval_prediction = combined_X_predictor_matrix_xval_norm(indicies_informing_xval_prediction,:);
                Y_informing_xval_prediction = Y_predictand_xval(indicies_informing_xval_prediction);

                X_testing_xval_prediction = combined_X_predictor_matrix_xval_norm(held_out_i,:);

                %do LASSO prediction
                
                %lasso predictor selection
    
                %https://www.mathworks.com/help/stats/lasso.html#mw_2715ae13-55dd-409c-aed8-49394e9c1964

                % Find the coefficients of a regularized linear regression model using 10-fold 
                % cross-validation and the elastic net method with Alpha = 0.75. Use the largest 
                % Lambda value such that the mean squared error (MSE) is within one standard error 
                % of the minimum MSE.

                    [B,FitInfo] = lasso(X_informing_xval_prediction,...
                        Y_informing_xval_prediction,...
                        'Alpha',...
                        lasso_alpha,...
                        'CV',...
                        kval_number,...
                        'MCReps',...
                        mc_trials);

                    idxLambdaMinMSE = FitInfo.IndexMinMSE;   
                    coef = B(:,idxLambdaMinMSE);
                    coef0 = FitInfo.Intercept(idxLambdaMinMSE);
                    
%                     idxLambda1SE = FitInfo.Index1SE;
%                      coef = B(:,idxLambda1SE);
%                      coef0 = FitInfo.Intercept(idxLambda1SE);

                    prediction = X_testing_xval_prediction*coef + coef0;
                
                %option to do MLR using only predictors that have non-zero
                %coeffiecients in LASSO
                                                        
%                     Mdl = fitlm(X_informing_xval_prediction(:,B(:,idxLambdaMinMSE)~=0),Y_informing_xval_prediction);
% 
%                     mlr_intercept = Mdl.Coefficients.Estimate(1);
%                     mlr_slopes = Mdl.Coefficients.Estimate(2:end);
% 
%                     %make a prediction using the held out point
% 
%                     prediction = mlr_intercept + sum(mlr_slopes'.*X_testing_xval_prediction(:,B(:,idxLambdaMinMSE)~=0));

                %compare prediction to real value

                xval_and_crawling_forecast_predictions(held_out_i) = prediction;


        end

            %calculations
            
                Y_xval_predictions = xval_and_crawling_forecast_predictions(1:crawling_forecast_ind-1);

                clim_prediction_RMSE = sqrt(mean((Y_predictand_xval-mean(Y_predictand_xval)).^2));

                prediction_RMSEs(1) = sqrt(mean((Y_xval_predictions-Y_predictand_xval).^2));

                prediction_RMSSSs(1) = 1 - prediction_RMSEs(1,1)./clim_prediction_RMSE;

                corr_mat = corrcoef(Y_xval_predictions,...
                                    Y_predictand_xval,...
                                    'rows',...
                                    'pairwise');

                prediction_Rs(1) = corr_mat(2,1);
                
%% crawling forecast portion                
    
num_crawling_inds = length(Y_predictand_2D_crawling_forecast);

crawling_climatology_forecast = NaN(size(Y_predictand_2D_crawling_forecast));

PLSR_forecasts_crawling = squeeze(PLSR_held_out_based_prediction_2D_space_unnorm_xval);

for crawling_ind = 1:num_crawling_inds
    
    universal_ind = crawling_ind + crawling_forecast_ind - 1;
    
    %add new predictor data 
    
    %pull new predictor data from the 2nd set of split X predictors
    
        X_predictor_matrix_to_add = X_predictor_matrix_crawling_forecast(1:crawling_ind,:);
        X_predictor_fields_to_add = X_predictor_fields_crawling_forecast(:,:,1:crawling_ind,:);
        
    %assemble the predictor and predictand for this year
    
        X_predictor_matrix_crawling = cat(1,X_predictor_matrix_xval,X_predictor_matrix_to_add);
        X_predictor_fields_crawling = cat(3,X_predictor_fields_xval,X_predictor_fields_to_add);
        
        Y_predictand_crawling = cat(1,Y_predictand_xval,Y_predictand_crawling_forecast(1:crawling_ind));
        
        Y_predictand_2D_crawling = NaN(1,1,length(Y_predictand_crawling));
        Y_predictand_2D_crawling(1,1,:) = Y_predictand_crawling;
        
        crawling_climatology_forecast(crawling_ind) = mean(Y_predictand_crawling);
                
        %DO PLSR forecast going forward
        
           [PLSR_forecast_XLoadings,...
            PLSR_forecast_YLoadings,...
            PLSR_forecast_XScores,...
            PLSR_forecast_YScores,...
            PLSR_forecast_BETA,...
            PLSR_forecast_PCTVAR,...
            PLSR_forecast_stats,...
            PLSR_forecast_Y_forecast_2D_space_unnorm] = PLSR_forecasts_noxval_no_shift_Y2D_alpha(...
            X_predictor_fields_crawling(:,:,1:end-1,:),...
            X_predictor_fields_crawling(:,:,end,:),...
            Y_predictand_2D_crawling(:,:,1:end-1,:),...
            num_pls_comps);
        
            PLSR_forecasts_crawling = vertcat(PLSR_forecasts_crawling,PLSR_forecast_Y_forecast_2D_space_unnorm);

            %Assemble large Predictor Matrix

            combined_X_predictor_matrix_crawling = cat(2,X_predictor_matrix_crawling,...
                                                         PLSR_forecasts_crawling);

            %Norm large Predictor Matrix

            combined_X_predictor_matrix_crawling_norm = normalize(combined_X_predictor_matrix_crawling,1,'zscore');
        
            % Do lasso predictor selection (it cant use the last datapoint
            % to train)

                    [B,FitInfo] = lasso(combined_X_predictor_matrix_crawling_norm(1:end-1,:),...
                                        Y_predictand_crawling(1:end-1),...
                                        'Alpha',...
                                        lasso_alpha,...
                                        'CV',...
                                        kval_number,...
                                        'MCReps',...
                                        mc_trials);

                    idxLambdaMinMSE = FitInfo.IndexMinMSE;        
                    coef = B(:,idxLambdaMinMSE);
                    coef0 = FitInfo.Intercept(idxLambdaMinMSE);
                    
%                     idxLambda1SE = FitInfo.Index1SE;
%                      coef = B(:,idxLambda1SE);
%                      coef0 = FitInfo.Intercept(idxLambda1SE);

                    %make a foreward looking prediction

                    prediction = combined_X_predictor_matrix_crawling_norm(end,:)*coef + coef0;
                
             %option to use MLR with only the non-zero lasso predictors
%                 
%                     Mdl = fitlm(combined_X_predictor_matrix_crawling_norm(1:end-1,B(:,idxLambdaMinMSE)~=0),Y_predictand_crawling(1:end-1));
% 
%                     mlr_intercept = Mdl.Coefficients.Estimate(1);
%                     mlr_slopes = Mdl.Coefficients.Estimate(2:end);
% 
%                     %make a prediction using the held out point
% 
%                     prediction = mlr_intercept + sum(mlr_slopes'.*combined_X_predictor_matrix_crawling_norm(end,B(:,idxLambdaMinMSE)~=0));
                
                %put into big array

                xval_and_crawling_forecast_predictions(universal_ind) = prediction;
        
end

            %calculations after all the crawling
            
                Y_crawling_predictions = xval_and_crawling_forecast_predictions(crawling_forecast_ind:end);

                clim_prediction_RMSE = sqrt(mean((squeeze(Y_predictand_2D_crawling_forecast)-squeeze(crawling_climatology_forecast)).^2));

                prediction_RMSEs(2) = sqrt(mean((Y_crawling_predictions-squeeze(Y_predictand_2D_crawling_forecast)).^2));

                prediction_RMSSSs(2) = 1 - prediction_RMSEs(1,1)./clim_prediction_RMSE;

                corr_mat = corrcoef(Y_crawling_predictions,...
                                    squeeze(Y_predictand_2D_crawling_forecast),...
                                    'rows',...
                                    'pairwise');

                prediction_Rs(2) = corr_mat(2,1);


%assemble the climatology series
    X_val_and_crawling_climatology(1:crawling_forecast_ind-1) = mean(Y_predictand_xval);
    X_val_and_crawling_climatology(crawling_forecast_ind:end) = crawling_climatology_forecast;
        
    
end


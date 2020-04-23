clear

% those are control parameters
dataset_size = 400;
% order = 2;

% generate random input with values +-1
u = 2*(rand(dataset_size,1)>0.5) - 1;

% preallocate and initialize all zeros for outputs
y = zeros(dataset_size,1);
y_2 = zeros(dataset_size, 1);
y_3 = zeros(dataset_size, 1);

% generate error sets with 4 and 9 variance
e_2 = 2*randn(dataset_size,1);
e_3 = 3*randn(dataset_size,1);

% calculate the uncorrupted and both corrupted outputs
% we start at 3, because first 2 values are 0(initial conditions)
for k = 3:dataset_size
    y(k) = 1.5*y(k-1) - 0.7*y(k-2) + u(k-1) + 0.5*u(k-2);
    y_2(k) = 1.5*y_2(k-1) - 0.7*y_2(k-2) + u(k-1) + 0.5*u(k-2)+ e_2(k);
    y_3(k) = 1.5*y_3(k-1) - 0.7*y_3(k-2) + u(k-1) + 0.5*u(k-2) + e_3(k);
end
%%
order = 1;
[current_variances, current_err_vect, current_thetas] = least_sqr_estimate(order, dataset_size, u, y_2, y_3);

%% 
order = 2;
[current_variances_2, current_err_vect_2, current_thetas_2] = least_sqr_estimate(order, dataset_size, u, y_2, y_3);

%%
order = 3;
[current_variances_3, current_err_vect_3, current_thetas_3] = least_sqr_estimate(order, dataset_size, u, y_2, y_3);
%%

% be careful with this part of the code
% it will generate as many plots, as many iterations you do
result = [];
for order = 1:50
    [current_variances_i, current_err_vect_i, current_thetas_i] = least_sqr_estimate(order, dataset_size, u, y_2, y_3);
    result = [result, current_err_vect_i(1)];
end
t = 1:50;
figure;
plot(t, result);

%%
% this function takes in the input, output vectors, the order of the system
% and the current index, and creates the current regression vector
function result = create_regression_vector(data_vector, input_vector, order, index)
    result = [];
    
    for i = 1:order
        result = [result; -data_vector(index - i)];
    end
    for i = 1:order
        result = [result; input_vector(index - i)];
    end
end

% this function takes in the order, the dataset_size, input and 2 outputs,
% and estimates parameters for all 2 cases, plots the prediction vs actual
% output and returns error values, variances of estimates and estimated
% parameters for all 3 cases
function [variances,  err_vect, thetas] = least_sqr_estimate(order, dataset_size, u, y_2, y_3)
    % allocating space for variances of the estimates for 2 cases
    variances = zeros(2*order, 2);

    % preallocate space for regression matrix for 2 cases
    PHI_2 = zeros(2*order,dataset_size);
    PHI_3 = zeros(2*order,dataset_size);

    % fill in the regression matrix for all 2 cases
    for k = order+1:dataset_size
        PHI_2(:,k) = create_regression_vector(y_2, u, order, k);
        PHI_3(:,k) = create_regression_vector(y_3, u, order, k);
    end

    % calculate covariance matrices for all 2 cases
    cov_2 = inv(PHI_2*transpose(PHI_2));
    cov_3 = inv(PHI_3*transpose(PHI_3));
    
    variances(:,1) = 4*diag(cov_2);
    variances(:,2) = 9*diag(cov_3);
    
    % use the formula derived in the class to estimate the parameters
    theta_2 = cov_2 * PHI_2 * y_2;
    theta_3 = cov_3 * PHI_3 * y_3;
    thetas = [theta_2, theta_3];

    % use estimated parameters and regression matrix to make a prediction
    pred_2 = transpose(PHI_2)*theta_2;
    pred_3 = transpose(PHI_3)*theta_3;

    % calculate the value of the cost function for all 2 cases
    cost_2 = sumsqr(pred_2 - y_2)/dataset_size;
    cost_3 = sumsqr(pred_3 - y_3)/dataset_size;
    err_2 = pred_2 - y_2;
    err_3 = pred_3 - y_3;
    err_vect = [cost_2; cost_3];

    % plot the actual and predicted value of the output for all 3 cases
    t = 1:dataset_size;
    figure;
    subplot(2, 2, 1);
    plot(t, [y_2,pred_2, err_2]);
    legend('Actual data,noise variance 4', 'Predicted value', "Prediction error, noise variance 4");
    title("Model order:" + order);

    subplot(2, 2, 2);
    plot(t, [y_3,pred_3, err_3]);
    legend('Actual data with error variance 9', 'Predicted value',  "Prediction error, noise variance 9");
    
    subplot(2,2,3);
    plot(t, err_2);
    legend("Prediction error, noise variance 4");
    
    subplot(2,2,4);
    plot(t, err_3);
    legend("Prediction error, noise variance 9");
    
    
end

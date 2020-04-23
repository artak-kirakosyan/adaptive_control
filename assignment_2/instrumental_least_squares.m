clear
 
% those are control parameters
dataset_size = 400;
order = 2;
 
% generate random input with values +-1
u = 2*(rand(dataset_size,1)>0.5) - 1;
% preallocate and initialize all zeros for outputs
y = zeros(dataset_size,1);
y_2 = zeros(dataset_size, 1);
y_3 = zeros(dataset_size, 1);
% generate error sets with 4 and 9 variance
e_2 = 2*randn(dataset_size,1);
e_3 = 3*randn(dataset_size,1);
% calculate the uncorrupted output
% we start at 3, because first 2 values are 0(initial conditions)
for k = 3:dataset_size
    y(k) = 1.5*y(k-1) - 0.7*y(k-2) + u(k-1);
end
 
% add the error after calculating uncorrupted output
y_2 = y + e_2;
y_3 = y + e_3;
%%
% use the function to estimate the thetas, their variances, get the
% prediction and cost function value
number_of_iterations = 3;
[theta_2, costs_2, variances_2, preds_2] = least_sqr_instrument(order, dataset_size, u, y_2, number_of_iterations);
[theta_3, costs_3, variances_3, preds_3] = least_sqr_instrument(order, dataset_size, u, y_3, number_of_iterations);
% plot the results and save to a file
figure_1 = plot_the_data(y, y_2, preds_2, costs_2, number_of_iterations, dataset_size, 4);
figure_2 = plot_the_data(y, y_3, preds_3, costs_3, number_of_iterations, dataset_size, 9);
% print(figure_1,'error_2.png','-dpng','-r1200');
% print(figure_2,'error_3.png','-dpng','-r1200');
 
%%
% this part does 10 iterations using dataset with error of variance of 4 to
% see the dynamics of the cost function
[theta, costs, variances, preds] = least_sqr_instrument(order, dataset_size, u, y_2, 10);
figure_3 = plot_the_data(y, y_2, preds, costs, 10, dataset_size, 4);
% print(figure_3,'error_2_10_iter.png','-dpng','-r1200');
 
%%
% this function takes in the order, the number of data points, input, 
% output and number of iterations and using instrumental method estimates
% the parameters and returns some characteristics
function [theta, costs, variances, predictions] = least_sqr_instrument(order, dataset_size, u, y, num_of_iterations)
    % preallocate and initialize the array for variances
    variances = zeros(2*order);
    
    % preallocate space to keep all predictions
    predictions = zeros(dataset_size, num_of_iterations+1);
    % preallocate space for regression matrix
    PHI = zeros(2*order,dataset_size);
    
    % fill in the regression matrix
    for k = order+1:dataset_size
        PHI(:,k) = create_regression_vector(y, u, order, k);
    end
 
    % calculate covariance matrix
    cov_1 = inv(PHI*transpose(PHI));
    %keep diagonal elements of the covariance matrix in the preallocated array 
    variances(:, 1) = diag(cov_1);
 
    % use the formula derived in the class to estimate the parameters
    theta = cov_1 * PHI * y;
    
    % make the initial prediction;
    pred = transpose(PHI) * theta;
    predictions(:,1) = pred;
    
    % calculate the cost function value for first estimation
    cost = sumsqr(pred - y)/dataset_size;
    
    % preallocate space for the cost values
    costs = zeros(num_of_iterations+1, 1);
    costs(1) = cost;
    
    for iteration_index = 1:num_of_iterations
        % preallocate space for new instrument and the new instrumental
        % matrix
        instrument = zeros(dataset_size, 1);
        eta = zeros(2*order, dataset_size);
       
        % fill in the instrumental matrix and create the instrument
        for k = order+1:dataset_size
            current_regression_vector = create_regression_vector(instrument, u, order, k);
            eta(:,k) = current_regression_vector;
            instrument(k) = transpose(current_regression_vector) * theta;
        end
        
        % calculate the covariance matrix  for current instrument set and
        % keep the values of parameter variances
        curr_cov = inv(eta * transpose(PHI));
        variances(:, iteration_index+1) = diag(curr_cov);
        %update the value of the parameter estimates
        theta = curr_cov * eta * y;
        
        % make a new prediction using updated parameter estimates
        pred = transpose(eta) * theta;
        predictions(:, iteration_index+1) = pred;
        % calculate updated cost function value and keep it
        cost = sumsqr(pred - y)/dataset_size;
        costs(iteration_index+1) = cost;
    end
end
 
% this function is used for plotting
function [fig] = plot_the_data(y, y_corrupted, predictions, costs, num_of_iterations, dataset_size, error_var)
    % create the legend for the plot
    leg = "Least squares";
    for k = 1:num_of_iterations
        leg = [leg, "Iteration: " + k];
    end
    % plot all predictions in one plot with actual measured data
    fig = figure;
    subplot(2, 1, 1);
    plot(1:dataset_size, [y, y_corrupted, predictions]);
    title("Predictions of the dataset with variance of " + error_var);
    legend(["Actual Data", "Measured data", leg]);
    
    subplot(2,1,2);
    plot(1:num_of_iterations+1, costs);
    title("Cost value dynamics over iterations");
end
 
% this function takes in the input, output vectors, the order of the system
% and the current index, and creates the current regression vector
function result = create_regression_vector(data_vector, input_vector, order, index)
    result = [];
    
    for i = 1:order
        result = [result; data_vector(index - i)];
    end
    for i = 1:order
        result = [result; input_vector(index - i)];
    end
end


clear
 
% those are control parameters
dataset_size = 400;
order = 2;
 
% generate random input with values +-1
u = 2*(rand(dataset_size,1)>0.5) - 1;
u_b = u;
u_b(100:300) = zeros(201,1);
u_c = u_b;
 
% preallocate and initialize all zeros for outputs
y = zeros(dataset_size,1);
y_b = zeros(dataset_size, 1);
y_c = zeros(dataset_size, 1);
 
% generate error sets with variance of 4
e = 2*randn(dataset_size,1);
 
% prepare data for all parts
for k = 3:dataset_size
    y(k) = 1.5*y(k-1) - 0.7*y(k-2) + u(k-1) + e(k);
    y_b(k) = 1.5*y_b(k-1) - 0.7*y_b(k-2) + u_b(k-1) + e(k);
    if k <= 250
        y_c(k) = 1.5*y_c(k-1) - 0.7*y_c(k-2) + u_c(k-1) + e(k);
    else
        y_c(k) = 1.7*y_c(k-1) - 0.9*y_c(k-2) + 5*u_c(k-1) + e(k);
    end
end
 
%%
% use the function to estimate the thetas, the trace of the covariance
% matrix and the final prediction
 
% Part A
forgetting_factor = 1;
[thetas, traces] = recursive_least_squares(order, dataset_size, u, y, forgetting_factor);
fig1 = plot_the_data(thetas, traces, y, "Lambda=1");
forgetting_factor = 0.95;
[thetas, traces] = recursive_least_squares(order, dataset_size, u, y, forgetting_factor);
fig2 = plot_the_data(thetas, traces, y, "Lambda=0.95");
 
% Part B
[thetas_b, traces_b] = recursive_least_squares(order, dataset_size, u_b, y_b, forgetting_factor);
fig3 = plot_the_data(thetas_b, traces_b, y_b, "Input to 0 from 100~300");
% Part C
[thetas_c, traces_c] = recursive_least_squares(order, dataset_size, u_c, y_c, forgetting_factor);
fig4 = plot_the_data(thetas_c, traces_c, y_c, "Change parameters from 250~");
 
 % % save figures to a file
% print(fig1,'part_a_lambda_1.png','-dpng','-r1200');
% print(fig2,'part_a_lambda_0_95.png','-dpng','-r1200');
% print(fig3,'part_b.png','-dpng','-r1200');
% print(fig4,'part_c.png','-dpng','-r1200');

% this function takes in the order, the number of data points, input, 
% output and number of iterations and using instrumental method estimates
% the parameters, calculates their variances, makes a prediction and
% calculates the cost function value
function [thetas, traces, cov_matrices] = recursive_least_squares(order, dataset_size, u, y, forgetting_factor)
 
    traces = zeros(dataset_size, 1);
    % initialize everything as zeros except diagonal elements of the first
    % covariance matrix
    thetas = zeros(2*order, dataset_size);
    cov_matrices = zeros(2*order, 2*order, dataset_size);
    % first covariance matrix is a diagonal matrix with huge positive
    % values on diagonal
    cov_matrices(:,:,1) = diag(repelem(100, 2*order));
 
    % initialize PHI as zeros
    PHI = zeros(2*order,dataset_size);
    
    % fill in the regression matrix
    for k = order+1:dataset_size
        PHI(:,k) = create_regression_vector(y, u, order, k);
    end
    
    for k = 2:dataset_size
        % use formulas derived in class to update covariance matrices and
        % new estimates
        traces(k) = trace(cov_matrices(:,:,k-1));
        numerator = cov_matrices(:,:,k-1)*PHI(:,k)*transpose(PHI(:,k)) * cov_matrices(:,:,k-1);
        denominator = forgetting_factor + transpose(PHI(:,k)) * cov_matrices(:,:,k-1) * PHI(:,k);
        
        cov_matrices(:,:,k) = (cov_matrices(:,:,k-1) - numerator/denominator)/forgetting_factor;
        thetas(:,k) = thetas(:,k-1) + cov_matrices(:,:,k) * PHI(:,k) * (y(k) - transpose(PHI(:,k)) * thetas(:,k-1));
    end
end
 
 
function [fig] = plot_the_data(thetas, traces, y, plt_title)
    dataset_size = size(y);
    dataset_size = dataset_size(1);
    disp(thetas(:,dataset_size));
    fig = figure;
    subplot(2,1,1);
    plot(1:dataset_size, thetas)
    title(plt_title + ", first plot: estimated parameters over time");
    legend(["y(k-1)", "y(k-2)","u(k-1)", "u(k-2)"], 'Location','southwest');
    legend("boxoff");
    subplot(2,1,2);
    plot(1:dataset_size, traces);
    title("Plot of the trace of covariance matrix over time");
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


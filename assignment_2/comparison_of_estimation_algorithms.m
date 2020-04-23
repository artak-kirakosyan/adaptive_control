clear
 
% those are control parameters
dataset_size = 400;
order = 1;
% generate random input with values +-1
x = randn(dataset_size, 1);
d = zeros(dataset_size, 1);
e_2 = 0*randn(dataset_size,1);
% calculate the uncorrupted output
for k = 2:dataset_size
    d(k) = 0.8*d(k-1) + 0.2*x(k)+ e_2(k);
end
 
%LMS estimation
alfa_lms = 0.4;
[thetas_lms, cost, pred] = least_mean_sqr(d, x, order, alfa_lms);
fig1 = plot_the_data_lms(thetas_lms, d, pred, order);
 
% RLS estimation
forgetting_factor = 1;
[thetas_rls, traces] = recursive_least_squares(order, dataset_size, x, d, forgetting_factor);
fig2 = plot_the_data_rls(thetas_rls, traces, d, "Recursive with alfa=1");
 
% STEEPEST DESCENT
%number of iterations for steepest descent
iter_num = 400;
% calculation of the R and p mmatrices is given in the report
R = [1/9, 0; 0, 1];
p = [4/45; 0.2];
% this is being used for cost calculation
sigma_sqr = R(1,1);
 
% calculate parameters using steep descent algorithm for 2 different cases
alfa_steep = 0.05;
initials = [0,0];
[thetas_steep] = steepest_descent(R, p, alfa_steep, iter_num, order, initials);
 
fig3 = figure;
plot(1:iter_num, thetas_steep);
legend("a", "b");
title("Parameter evolution with alfa 0.05 and initial values [0,0]");
 
% create first plot of contour lines in a small range close to true values
a = linspace(0, 1, 1001);
b = linspace(0, 1, 1001);
fig4 = plot_contours(a, b, sigma_sqr, p, R);
 
% create another plot for a wide range of values
a = linspace(-10, 10, 1001);
b = linspace(-10, 10, 1001);
fig5 = plot_contours(a, b, sigma_sqr, p, R);
%% 
% save figures to files
% print(fig1,'lms_plot.png','-dpng','-r1200');
% print(fig2,'rls_plot.png','-dpng','-r1200');
% print(fig3,'parameter_evolution.png','-dpng','-r1200');
% print(fig4,'contours_small.png','-dpng','-r1200');
% print(fig5,'contours_big.png','-dpng','-r1200');
%%
 
function [fig] = plot_contours(a, b, sigma_sqr, p, R)
    % create a meshgrid using a and b
    [A,B] = meshgrid(a,b);
 
    % preallocate space for cost values
    costs_steep = zeros(size(B));
 
    % calculate the cost value for all pairs of a and b
    for i = 1:length(a)
        for j = 1:length(b)
            costs_steep(i,j) = sigma_sqr - 2*[a(i),b(j)]*p + [a(i),b(j)]*R*[a(i);b(j)];
        end
    end
 
    fig = figure;
    contour(A, B, costs_steep, 10, 'ShowText','on');
end
 
 
% this function uses steepest descent algorithm to calculate the parameters
function [thetas] = steepest_descent(R, p, alfa, iter_num, order, initials)
 
    thetas = zeros(2*order, iter_num);
    thetas(:,1) = initials;
    for k = 2:iter_num
        thetas(:,k) = (eye(2*order) - alfa*R)*thetas(:,k-1) + alfa*p;
    end
end
 
% this function uses least mean square algorithm to estimate the systems
% parameters.
function [thetas, cost, pred] = least_mean_sqr(d, x, order, alfa)
    % preallocate space for all necessary variables
    dataset_size = size(x);
    dataset_size = dataset_size(1);
    PHI = zeros(2*order+1, dataset_size);
    thetas = zeros(2*order+1,dataset_size);
    
    % prepare the regression vector for later use
    for k = order+1:dataset_size
        PHI(:,k) = create_regression_vector(d, x, order, k);
    end
 
    for k = order+1:dataset_size
        % dirty estimate of information matrix and p vector
        p = d(k-1) * PHI(:,k-1);
        R = PHI(:,k-1)*transpose(PHI(:,k-1));
        % update the estimate
        thetas(:,k) = thetas(:,k-1) + alfa * p - alfa * R * thetas(:,k-1);
    end
    % make a prediction using latest predicitons
    pred = transpose(PHI) * thetas(:,dataset_size);
    cost = sumsqr(pred - d)/dataset_size;
end
 
% this function uses recursive least squares algorithm to estimate systems
% parameters
function [thetas, traces, cov_matrices] = recursive_least_squares(order, dataset_size, u, y, forgetting_factor)
 
    traces = zeros(dataset_size, 1);
    % initialize everything as zeros except diagonal elements of the first
    % covariance matrix
    thetas = zeros(2*order+1, dataset_size);
    cov_matrices = zeros(2*order+1, 2*order+1, dataset_size);
    % first covariance matrix is a diagonal matrix with huge positive
    % values on diagonal
    cov_matrices(:,:,1) = diag(repelem(100, 2*order+1));
    % initialize PHI as zeros
    PHI = zeros(2*order+1,dataset_size);
    
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
 
function fig = plot_the_data_lms(thetas, d, pred, order)
    leg = [];
    for k = 1:2*order+1
        leg = [leg, "Coefficient: " + k];
    end
    
    dataset_size = size(d);
    dataset_size = dataset_size(1);
    fig = figure;
    subplot(2,1,1);
    plot(1:dataset_size, thetas);
    title("Plot of estimates over time")
    legend(leg);
    subplot(2,1,2);
    plot(1:dataset_size, [d, pred]);
    title("Plot of the actual and predicted values over time");
    legend("Actual", "Predicted");
end
 
function fig = plot_the_data_rls(thetas, traces, y, plt_title)
    dataset_size = size(y);
    dataset_size = dataset_size(1);
    disp(thetas(:,dataset_size));
    fig = figure;
    subplot(2,1,1);
    plot(1:dataset_size, thetas)
    title(plt_title + ", first plot: estimated parameters over time");
    legend(["y(k-1)", "u(k)","u(k-1)"], 'Location','southwest');
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
    for i = 1:order+1
        result = [result; input_vector(index - i + 1)];
    end
end


clear;

% control parameters
number_of_states = 9;
number_of_actions = 4;
number_of_iterations = 2000;
gamma = 0.99;
theta = 0.0001;

V = zeros(number_of_states, 1);
% keep all results to plot later
Vs = zeros(number_of_states, number_of_iterations);
% iterations
for k = 1:number_of_iterations
    
    % save previous values
    V_prev = V;
    % update value function for each state
    for s = 1:number_of_states
        V(s) = 0;
        % for each state, iterate of all actions
        for a = 1:number_of_actions
            temp = 0;
            % for each action, iterate of next state and sum up all
            % immediate and discounted next reward
            for s_pr = 1:number_of_states
                temp = temp + P(s, s_pr, a) * (R(s, s_pr, a) + gamma*V_prev(s_pr)); 
            end
            % sum up for each aciton
            V(s) = V(s) + pi(s, a) * temp;
        end  
    end
    % save the current value function
    Vs(:, k) = V;
    % if there are not big changes for any of the states, stop the
    % iteration
    if abs( V_prev - V) < theta
        break;
    end
end

%plot the results
disp(reshape(V, [3,3]));
plot(1:k, Vs(:, 1:k));
legend(string(1:9), "Location", "Best");
title("Value evolution over iterations");

% this function returns the expected reward while taking the action a at
% state s and transfering to stat s_pr
function r = R(s, s_pr, a)
    % make a three digit number out of input arguments for easy
    % implementation.
    index = s*100 + s_pr * 10 + a;
    switch index
        case {191 192 193 194}
            r = 10;
        case {221 331 334 443 664 772 773 882 992 994}
            r = -1;
        otherwise
            r = 0;
    end
end

% return the proability of transfering to s_pr from s if action a is taken
function p = P(s, s_pr, a)
    % make a three digit number out of input arguments for easy
    % implementation.
    index = s*100 + s_pr * 10 + a;
    switch index
        case {191 192 193 194 
              221 252 213 234
              331 362 323 334
              411 472 443 454
              521 582 543 564
              631 692 653 664
              741 772 773 784
              851 882 873 894
              961 992 983 994}
            p = 1;
        otherwise
            p = 0;
    end
            
end

% this function immitates the policy. it return the probability of taking
% action a while being at state s
function p = pi(s, a)
    % make a three digit number out of input arguments for easy
    % implementation.
    index = 10 * s + a;
    switch index
        otherwise p = 0.25;
    end
end
clear

% control parameters
number_of_states = 9;
number_of_actions = 4;
number_of_iterations = 1000;
gamma = 0.9;
% generate initial policy: SxA matrix, where the row shows the state, and
% column shows the action
pi = zeros(number_of_states, number_of_actions) + (1/number_of_actions);
% find value function for initial policy
V = iterate_value(pi, gamma);
previous_pi = pi;

for k = 1:number_of_iterations
    % do for each state
    for s = 1:number_of_states
        % act_values keeps the state actions values for current
        % state:Q(s,a)
        act_values = zeros(number_of_actions, 1);
        % calculate the Q(s, a) for each action
        for a = 1:number_of_actions
            % iterate over the next states and find the Q(s, a)
            temp = 0;
            % for each action, iterate of next state and sum up all
            % immediate and discounted next reward
            for s_pr = 1:number_of_states
                temp = temp + P(s, s_pr, a) * (R(s, s_pr, a) + gamma*V(s_pr)); 
            end
            act_values(a) = temp;
        end
        % find all indexes for which the value is equal to the maximum
        % among them
        act_new = find(act_values == max(act_values));
        % create a list of zeros and change the ones who are maximum
        % distribute the probability among all maximum valued acitons
        temp = zeros(1, number_of_actions);
        temp(act_new) = 1/length(act_new);
        pi(s, :) = temp;
        % if the policy changed, find the value function of the current
        % policy
        if previous_pi(s, :) ~= pi(s, :)
            V = iterate_value(pi, gamma);
        end
    end
    % if the policy didnt change after iteration over all states, than we
    % have the optimal policy
    if previous_pi == pi
        break
    end
    % save the previous policy
    previous_pi = pi;
end
disp(reshape(V, [3, 3]));
disp(pi);

% this function takes the pi and the current gamma, and does the value
% iteration. Its almost the same code as in part 3
function V = iterate_value(pi, gamma)
    sz = size(pi);
    number_of_states = sz(1);
    number_of_actions = sz(2);
    V = zeros(number_of_states, 1);
    number_of_iterations = 2000;
    theta = 0.001;
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
        % if there are not big changes for any of the states, stop the
        % iteration
        if abs( V_prev - V) < theta
            break;
        end
    end

end

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
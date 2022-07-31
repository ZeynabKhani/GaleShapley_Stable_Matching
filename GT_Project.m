clc; clear; close all;

surgery_types = struct('A', 1, 'B', 2, 'C', 3); hospitals = struct('a', 1, 'b', 2, 'c', 3); dismissed = 4;
mus = [2 4 6; 15 5 10; 20 30 10]; surgery_costs = [2 7 5; 5 2 7; 7 5 2]; hospitals_preference_relations = [3 1 2; 2 3 1; 1 2 3];
slot = 8;

mean_surgery_times_k = []; mean_request_completion_times_k = []; mean_costs_k = []; max_rates_k = []; 
% simulating for different k's (maximum number of slots an entity waits in
% the queue, hoping to be matched to its best choice of hospitals)               
for k = 0:1:5
    clock = 0; rate = 1; surgery_index = 1; max_rate = 1;
    mean_surgery_times = []; mean_request_completion_times = []; mean_costs = []; max_rates = [];
    simulation_times = [];
    hospital_a_busy = false; hospital_b_busy = false; hospital_c_busy = false;
    queue = []; arrival_times = [];  queue_k = [];
    served_surgeries = []; server_hospitals = []; serving_times = []; serving_clock = [];
    serving_surgeries_indexes_in_queue = []; hospitals_busy_times = zeros(1, 3);

    for simulation_time = 0:100:2000
        simulation_times = [simulation_times simulation_time];
        max_rate = 1;
        while clock <= simulation_time
            surgeries_preference_relations = []; surgery_durations = [];
            if mod(clock, slot) == 0 && clock ~= 0
                % adding entities to the queue
                [arrival_times, queue, queue_k] = create_queue(arrival_times, queue, rate, clock - slot, queue_k, k);

                % calculating the surgeries preference relation of free hospitals
                % and surgery durations
                free_hospitals = [];
                if ~hospital_a_busy
                    free_hospitals = [hospitals.a];
                end
                if ~hospital_b_busy
                    free_hospitals = [free_hospitals hospitals.b];
                end
                if ~hospital_c_busy
                    free_hospitals = [free_hospitals hospitals.c];
                end
                [surgeries_preference_relations, surgery_durations, best_hospitals] = find_preference(mus, free_hospitals);

                % choosing the first 3 unserved entities of the queue (if any)
                if ~isempty(queue) && ~isempty(surgeries_preference_relations)
                    % considering k
                    [queue, queue_k] = wait_in_queue(surgeries_preference_relations, queue_k, queue, surgery_index, best_hospitals, free_hospitals);
                end

                requests = [];
                requests_indexes = [];
                requests_arrival_times = [];
                requests_k = [];
                unserved_requests = [];
                l = length(queue);
                if surgery_index == l
                    requests = queue(1,surgery_index);
                    requests_indexes = surgery_index;
                    requests_arrival_times = arrival_times(1, surgery_index);
                    requests_k = queue_k(1, surgery_index);
                elseif surgery_index+1 == l
                    requests = queue(1, surgery_index:surgery_index+1);
                    requests_indexes = [surgery_index, surgery_index + 1];
                    requests_arrival_times = arrival_times(1, surgery_index:surgery_index+1);
                    requests_k = queue_k(1, surgery_index:surgery_index+1);
                elseif surgery_index+2 <= l
                    requests = queue(1, surgery_index:surgery_index+2);
                    requests_indexes = [surgery_index, surgery_index+1, surgery_index+2];
                    requests_arrival_times = arrival_times(1, surgery_index:surgery_index+2);
                    requests_k = queue_k(1, surgery_index:surgery_index+2);
                end

                % Gale-Shaplely algorithm for free hospitals and surgeries in the queue (FIFO)
                if ~isempty(requests) && ~isempty(free_hospitals)
                    [served_surgeries, server_hospitals, unserved_requests, hospitals_busy_times, serving_times, serving_surgeries_indexes_in_queue, serving_clock] = Gale_Shapley_Alg(requests, requests_indexes, serving_surgeries_indexes_in_queue, free_hospitals, surgery_durations, surgeries_preference_relations, hospitals_preference_relations, served_surgeries, server_hospitals, hospitals_busy_times, serving_times, serving_clock, clock);
                    % taking unserved requests back to the queue
                    if any(unserved_requests)
                        [queue, queue_k, arrival_times, surgery_index] = queue_rearrangement(queue, queue_k, arrival_times, requests, requests_arrival_times, requests_k, unserved_requests, surgery_index);
                    end
                end

                % updating the hospitals' status
                for i =1:1:3
                    if hospitals_busy_times(i) <= slot
                        hospitals_busy_times(i) = 0;
                    else
                        hospitals_busy_times(i) = hospitals_busy_times(i) - slot;
                    end
                end

                if hospitals_busy_times(1) ~= 0
                    hospital_a_busy = true;
                else
                    hospital_a_busy = false;
                end
                if hospitals_busy_times(2) ~= 0
                    hospital_b_busy = true;
                else
                    hospital_b_busy = false;
                end
                if hospitals_busy_times(3) ~= 0
                    hospital_c_busy = true;
                else
                    hospital_c_busy = false;
                end

                % changing the requests rate to make a stable queue
                if length(queue) - length(served_surgeries) < length(free_hospitals)
                    rate = rate*1.1;
                    if rate > max_rate
                        max_rate = rate;
                    end
                elseif length(queue) - length(served_surgeries) == length(free_hospitals)
                    rate = 1;
                elseif length(queue) - length(served_surgeries) > length(free_hospitals)
                    rate = rate*0.9;
                end
            end
            clock = clock + 1;
        end

        % mean surgeries
        a = []; b = []; c = [];
        for i = 1:1:length(served_surgeries)
            if served_surgeries(i) == surgery_types.A
                a = [a serving_times(i)];
            elseif served_surgeries(i) == surgery_types.B
                b = [b serving_times(i)];
            else
                c = [c serving_times(i)];
            end
        end
        mean_surgery_times = [mean_surgery_times [mean(a); mean(b); mean(c)]];

        % mean request completion time
        a = []; b = []; c = [];
        for i = 1:1:length(served_surgeries)
            if served_surgeries(i) == surgery_types.A
                a = [a serving_times(i)+serving_clock(i)-arrival_times(serving_surgeries_indexes_in_queue(i))];
            elseif served_surgeries(i) == surgery_types.B
                b = [b serving_times(i)+serving_clock(i)-arrival_times(serving_surgeries_indexes_in_queue(i))];
            else
                c = [c serving_times(i)+serving_clock(i)-arrival_times(serving_surgeries_indexes_in_queue(i))];
            end
        end
        mean_request_completion_times = [mean_request_completion_times [mean(a); mean(b); mean(c)]];

        % mean costs
        a = []; b = []; c = [];
        for i = 1:1:length(served_surgeries)
            if served_surgeries(i) == surgery_types.A
                a = [a surgery_costs(served_surgeries(i), server_hospitals(i))];
            elseif served_surgeries(i) == surgery_types.B
                b = [b surgery_costs(served_surgeries(i), server_hospitals(i))];
            else
                c = [c surgery_costs(served_surgeries(i), server_hospitals(i))];
            end
        end
        mean_costs = [mean_costs [mean(a); mean(b); mean(c)]];

        % max rate
        max_rates = [max_rates max_rate];
    end

    if ~isempty(mean_surgery_times_k)
        mean_surgery_times_k = cat(3, mean_surgery_times_k, mean_surgery_times);
        mean_request_completion_times_k = cat(3, mean_request_completion_times_k, mean_request_completion_times);
        mean_costs_k = cat(3, mean_costs_k, mean_costs);
    else
        mean_surgery_times_k = mean_surgery_times;
        mean_request_completion_times_k = mean_request_completion_times;
        mean_costs_k = mean_costs;
    end
    k
    mean_surgery_times_ = [mean_surgery_times(1, end) mean_surgery_times(2, end) mean_surgery_times(3, end)]
    mean_request_completion_times_ = [mean_request_completion_times(1, end) mean_request_completion_times(2, end) mean_request_completion_times(3, end)]
    mean_costs_ = [mean_costs(1, end) mean_costs(2, end) mean_costs(3, end)]
    max_rates_ = max_rates(end)
end
%% figures
figure
subplot(3, 1, 1)
for i = 1:1:k+1
    plot(simulation_times, mean_surgery_times_k(1,:, i))
    hold on
end
hold off
title('Mean Duration of Surgery A')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')
subplot(3, 1, 2)
for i = 1:1:k+1
    plot(simulation_times, mean_surgery_times_k(2,:, i))
    hold on
end
hold off
title('Mean Duration of Surgery B')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')
subplot(3, 1, 3)
for i = 1:1:k+1
    plot(simulation_times, mean_surgery_times_k(3,:, i))
    hold on
end
hold off
title('Mean Duration of Surgery C')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')

figure
subplot(3, 1, 1)
for i = 1:1:k+1
    plot(simulation_times, mean_request_completion_times_k(1,:, i))
    hold on
end
hold off
title('Mean Request Completion Duration of Surgery A')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')
subplot(3, 1, 2)
for i = 1:1:k+1
    plot(simulation_times, mean_request_completion_times_k(2,:, i))
    hold on
end
hold off
title('Mean Request Completion Duration of Surgery B')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')
subplot(3, 1, 3)
for i = 1:1:k+1
    plot(simulation_times, mean_request_completion_times_k(3,:, i))
    hold on
end
hold off
title('Mean Request Completion Duration of Surgery C')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')

figure
subplot(3, 1, 1)
for i = 1:1:k+1
    plot(simulation_times, mean_costs_k(1,:, i))
    hold on
end
hold off
title('Mean Cost of Surgery A')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')
subplot(3, 1, 2)
for i = 1:1:k+1
    plot(simulation_times, mean_costs_k(2,:, i))
    hold on
end
hold off
title('Mean Cost of Surgery B')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')
subplot(3, 1, 3)
for i = 1:1:k+1
    plot(simulation_times, mean_costs_k(3,:, i))
    hold on
end
hold off
title('Mean Costs of Surgery C')
legend('k = 0', 'k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5')
xlabel('Simulation Time')

figure;
plot(simulation_times, mean_surgery_times(:,:,1))
title('Mean Surgery Durations for k = 0')
legend('A', 'B', 'C')
xlabel('Simulation Time')

figure
plot(simulation_times, mean_request_completion_times(:,:,1))
title('Mean Request Completion Durations for k = 0')
legend('A', 'B', 'C')
xlabel('Simulation Time')

figure
plot(simulation_times, mean_costs(:,:,1))
title('Mean Costs for k = 0')
legend('A', 'B', 'C')
xlabel('Simulation Time')

figure
plot(simulation_times, max_rates)
title('Maximum Rates')
xlabel('Simulation Time')
%% functions
function [queue, queue_k] = wait_in_queue(surgeries_preference_relations, queue_k, queue, surgery_index, best_hospitals, free_hospitals)
    out = [];
    wait = [];
    l = length(queue);
    counter = surgery_index;
    while counter <= l && length(out) <=3
        if queue_k(counter) ~= 0
            [~, I] = min(surgeries_preference_relations(queue(counter), :));
            if free_hospitals(I) ~= best_hospitals(queue(counter))
                wait = [wait queue(counter)];
                queue_k(counter) = queue_k(counter) - 1;
            else
                out = [out queue(counter)];
%                queue_k(surgery_index) = 0; % should it keep the k?
            end
        else
            out = [out queue(counter)];
        end
        counter = counter + 1;
    end
    queue = [queue(1, 1:surgery_index - 1) out wait queue(1, length(out)+length(wait)+surgery_index:length(queue))];
end
%%
function [served_surgeries, server_hospitals, unserved_requests, hospitals_busy_times, serving_times, serving_surgeries_indexes_in_queue, serving_clock] = Gale_Shapley_Alg(requests, requests_indexes, serving_surgeries_indexes_in_queue, free_hospitals, surgery_durations, surgeries_preference_relations, hospitals_preference_relations, served_surgeries, server_hospitals, hospitals_busy_times, serving_times, serving_clock, clock)
    R = 4; hospitals = struct('a', 1, 'b', 2, 'c', 3);
    unserved_requests = requests;
    choices = zeros(1, length(free_hospitals));
    choice_indexes = zeros(1, length(free_hospitals));
    a_requests = []; b_requests = []; c_requests = [];
    a_request_index = []; b_request_index = []; c_request_index = [];

    visited_hospitals = zeros(length(requests), length(free_hospitals));

    while true
        % stage 1
        for i = 1:1:length(unserved_requests)
            if unserved_requests(i) ~= R
                for j = 1:1:length(free_hospitals)
                    if visited_hospitals(i, j) == 1
                        surgeries_preference_relations(unserved_requests(i), j) = R;
                    end
                end

                [~, I] = min(surgeries_preference_relations(unserved_requests(i), :));

                if free_hospitals(1, I) == hospitals.a
                    a_requests = [a_requests unserved_requests(i)];
                    a_request_index = [a_request_index i];
                    visited_hospitals(i, I) = 1;
                elseif free_hospitals(1, I) == hospitals.b
                    b_requests = [b_requests unserved_requests(i)];
                    b_request_index = [b_request_index i];
                    visited_hospitals(i, I) = 1;
                else
                    c_requests = [c_requests unserved_requests(i)];
                    c_request_index = [c_request_index i];
                    visited_hospitals(i, I) = 1;
                end
            end
        end

        % stage 2
        for i = 1:1:length(free_hospitals)
            if free_hospitals(1, i) == hospitals.a
                p = [];
                for j = 1:1:length(a_requests)
                    p = [p hospitals_preference_relations(a_requests(j), hospitals.a)]; 
                end
                if ~isempty(p)
                    [~, I] = min(p);
                    if choices(i) ~= a_requests(1, I) && choices(i) ~= 0
                        unserved_requests(choice_indexes(i)) = choices(i);
                    end
                    choices(i) = a_requests(1, I);
                    choice_indexes(i) = a_request_index(1, I);
                    a_requests = []; a_requests = choices(i);
                    a_request_index = []; a_request_index = choice_indexes(i);
                    unserved_requests(choice_indexes(i)) = R;
                end
            elseif free_hospitals(1, i) == hospitals.b
                p = [];
                for j = 1:1:length(b_requests)  
                    p = [p hospitals_preference_relations(b_requests(j), hospitals.b)]; 
                end
                if ~isempty(p)
                    [~, I] = min(p);
                    if choices(i) ~= b_requests(1, I) && choices(i) ~= 0
                        unserved_requests(choice_indexes(i)) = choices(i);
                    end
                    choices(i) = b_requests(1, I);
                    choice_indexes(i) = b_request_index(1, I);
                    b_requests = []; b_requests = choices(i);
                    b_request_index = []; b_request_index = choice_indexes(i);
                    unserved_requests(choice_indexes(i)) = R;
                end
                 
            else
                p = [];
                for j = 1:1:length(c_requests)
                    p = [p hospitals_preference_relations(c_requests(j), hospitals.c)]; 
                end
                if ~isempty(p)
                    [~, I] = min(p);
                    if choices(i) ~= c_requests(1, I) && choices(i) ~= 0
                        unserved_requests(choice_indexes(i)) = choices(i);
                    end
                    choices(i) = c_requests(1, I);
                    choice_indexes(i) = c_request_index(1, I);
                    c_requests = []; c_requests = choices(i);
                    c_request_index = []; c_request_index = choice_indexes(i);
                    unserved_requests(choice_indexes(i)) = R;
                end
            end
        end
        
        all_served = 0;
        for i = 1:1:length(unserved_requests)
            if unserved_requests(i) == R
                all_served = all_served + 1;
            end
        end

        c = 0;
        for i= 1:1:length(choices)
            if choices(i) == 0
                c = c + 1;
            end
        end

        if c == 0 || all_served == length(unserved_requests)
            for i = 1:1:length(free_hospitals)
                if choices(i) ~= 0
                    served_surgeries = [served_surgeries choices(i)];
                    serving_surgeries_indexes_in_queue = [serving_surgeries_indexes_in_queue requests_indexes(choice_indexes(i))];
                    server_hospitals = [server_hospitals free_hospitals(i)];
                    serving_times = [serving_times surgery_durations(choices(i), i)];
                    hospitals_busy_times(free_hospitals(i)) = surgery_durations(choices(i), i);
                    serving_clock = [serving_clock clock];
                end
            end
            break;
        end
    end
end
%%
function [queue, queue_k, arrival_times, surgery_index] = queue_rearrangement(queue, queue_k, arrival_times, requests, requests_arrival_times, requests_k, unserved_requests, surgery_index)
    q = []; a = []; k = []; R = 4; counter = 0;
    for i = 1:1:length(unserved_requests)
        if unserved_requests(i) ~= R
            counter = counter + 1;
            q = [q unserved_requests(i)];
            a = [a requests_arrival_times(i)];
            k = [k requests_k(i)];
        end
    end

    for i = length(unserved_requests):-1:1
        if unserved_requests(i) == R
            q = [requests(i) q];
            a = [requests_arrival_times(i) a];
            k = [requests_k(i) k];
        end
    end
   
    queue = [queue(1, 1:surgery_index-1) q queue(1, surgery_index + length(unserved_requests):length(queue))];
    arrival_times = [arrival_times(1, 1:surgery_index-1) a arrival_times(1, surgery_index + length(unserved_requests):length(arrival_times))];
    queue_k = [queue_k(1, 1:surgery_index-1) k queue_k(1, surgery_index + length(unserved_requests):length(queue_k))];
    surgery_index = surgery_index + length(unserved_requests) - counter;
end
%%
function [surgeries_preference_relations, free_surgery_durations, best_hospitals] = find_preference(mus, free_hospitals)
    surgery_durations = normrnd(mus, ones(size(mus, 1), size(mus, 2)));
    best_hospitals = [];
    surgeries_preference_relations_mus = [];
    for i = 1:1:3
        [~, ~, ranks] = unique(surgery_durations(i,:));
        surgeries_preference_relations_mus  = [surgeries_preference_relations_mus;ranks'];
    end
    for i = 1:1:length(surgeries_preference_relations_mus)
        [~, I] = min(surgeries_preference_relations_mus(i,:));
        best_hospitals(i) = I;
    end

    free_surgery_durations = [];
    for i = 1:1:length(free_hospitals)
        free_surgery_durations = [free_surgery_durations surgery_durations(:,free_hospitals(i))];
    end
    surgeries_preference_relations = [];
    for i = 1:1:3
        if ~isempty(free_surgery_durations)
            [~, ~, ranks] = unique(free_surgery_durations(i,:));
            surgeries_preference_relations  = [surgeries_preference_relations;ranks'];
        end
    end
end
%%
function [arrivals, oprs] = sort_arrivals(arrivals, oprs)
    for i = 1:1:length(arrivals)
        for j = 1:1:i
            if arrivals(i) < arrivals(j)
                swap = arrivals(i);
                arrivals(i) = arrivals(j);
                arrivals(j) = swap;
                swap = oprs(i);
                oprs(i) = oprs(j);
                oprs(j) = swap;
            end
        end
    end
end
%%
function [arrival_times, queue, queue_k] = create_queue(arrival_times, queue, rate, clock, queue_k, k)
    slot = 8;
    surgery_types = struct('A', 1, 'B', 2, 'C', 3);
    lambda_A = 0.2*rate*slot;
    lambda_B = 0.1*rate*slot;
    lambda_C = 0.05*rate*slot;
    
    oprs = poissrnd([lambda_A lambda_B lambda_C]);

    A_arrivals = cumsum(exprnd(1/lambda_A, 1, oprs(1)));
    B_arrivals = cumsum(exprnd(1/lambda_B, 1, oprs(2)));
    C_arrivals = cumsum(exprnd(1/lambda_C, 1, oprs(3)));
    arrivals = [A_arrivals B_arrivals C_arrivals] + clock;

    surgeries = [];
    for i = 1:1:oprs(1)
        surgeries = [surgeries surgery_types.A];
    end
    for i = 1:1:oprs(2)
        surgeries = [surgeries surgery_types.B];
    end
    for i = 1:1:oprs(3)
        surgeries = [surgeries surgery_types.C];
    end

    [arrivals, surgeries] = sort_arrivals(arrivals, surgeries);

    slot = 8;
    for i = 1:1:length(arrivals)
        if  arrivals(i) <= clock + slot && arrivals(i) >= clock - slot
            queue = [queue surgeries(i)];
            queue_k = [queue_k k];
            arrival_times = [arrival_times arrivals(i)];
        end
    end
end
clear; clc;

sigma = 0.5; % resource dependent MAXIMUM spore yield
mz = 0.8; % loss rate/background mortality rate of spores
a = 0.005; % space clearance rate
epsilon = 0.7156; % host suscpetibilty/ per-parasite infectivity
e = 0.02; % conversion efficiency
ms = 0.05; % host background death rate
v = 0.05; % virulent effects of infection on host survival
mi = ms + v; % host background death rate
r = 2; % maximal growth rate of the algal resource
K = 1500; % algal carrying capacity

y0 = [100 0 400 500]; % initial prey and predator densities
t_max = 300;
tspan = [0 t_max]; % start end times
ode = @(t,y) RSIZ_model(t,y,sigma,mz,a,epsilon,mi,e,ms,r,K); % compile function and call
    [t1,y1] = ode45(ode, tspan, y0); % return time and population density vectors

figure(1); clf(1);
    box on; hold on;
    h1 = plot(t1,y1(:,1),'-k','LineWidth',2);
    h2 = plot(t1,y1(:,2),'-r','LineWidth',2);
    h3 = plot(t1,y1(:,3),'-g','LineWidth',2);
    h4 = plot(t1,y1(:,4),'-b','LineWidth',2);
    legend([h1 h2 h3 h4],'Z','I','S','R');
    xlim([0 t_max]);
    shg;

%% define key evolution values and traits

IS_host_1 = 9; % value of species 1 in immunological space
IS_par_A = 10; % area under the curve
IS_par_p = 1; % height of the peak
IS_par_range_mean = 5; % location of peak on the x axis
epsilon_mean = IS_par_p - (16*IS_par_p^3*(IS_host_1 - IS_par_range_mean).^2)./(9*IS_par_A.^2); % calc epsilon

cv = 0.3; % set coefficient of variation in traits
h_2 = 0.75; % set heritability of traits

%% Gillespie algorithm

num_replicates = 1; % number of simulations
stand_times = 0:1:t_max; % standardized time steps for storing time series
num_time_steps = length(stand_times);

R_stand = nan(num_replicates,num_time_steps); % preallocate matrix for standardized population size
S_stand = nan(num_replicates,num_time_steps);
I_stand = nan(num_replicates,num_time_steps);
Z_stand = nan(num_replicates,num_time_steps);
p_stand = nan(num_replicates,num_time_steps);
rm_stand = nan(num_replicates,num_time_steps);
epsilon_stand = nan(num_replicates,num_time_steps);
x_var_stand = nan(num_replicates,num_time_steps);

for i = 1:num_replicates % start Gillespie algorithm
    i
    % preallocate for whole time series
        R = zeros(1,1e7); % 1e7 is just a large number to ensure the vector is long enough
        S = zeros(1,1e7);
        I = zeros(1,1e7);
        Z = zeros(1,1e7);
        x_mean = nan(1e7,3);
        x_var = nan(1e7,3);
        t = nan(1,1e7);
    % define initial states
        t(1) = 0; % initial time
        R(1) = y0(4); % initial resource (prey) population size
        S(1) = y0(3); % initial susceptible population size
        I(1) = y0(2); % initial infected population size
        Z(1) = y0(1); % initial parasite population size

    % create initial distribution for parameter
        rng('shuffle'); % change random number seed

        % distribution of trait p
        MU = log(IS_par_p^2 / sqrt((cv*IS_par_p)^2+IS_par_p^2)); % mean for lognormal
        SIGMA1 = sqrt(log((cv*IS_par_p)^2/IS_par_p^2 + 1)); % std for lognormal
        x_dist_init(:,1) = lognrnd(MU,SIGMA1,Z(1),1); % specify initial distribution of traits

        % distribution of trait range mean
        MU = log(IS_par_range_mean^2 / sqrt((cv*IS_par_range_mean)^2+IS_par_range_mean^2)); % mean for lognormal
        SIGMA2 = sqrt(log((cv*IS_par_range_mean)^2/IS_par_range_mean^2 + 1)); % std for lognormal
        x_dist_init(:,2) = lognrnd(MU,SIGMA2,Z(1),1); % specify initial distribution of traits

        x_dist_init(:,3) = x_dist_init(:,1) - (16.*x_dist_init(:,1).^3.*(IS_host_1 - x_dist_init(:,2)).^2)./(9*IS_par_A.^2);
            x_dist_init(x_dist_init(:,3)<0,3) = 0; % swap out negatives for zeros

        x_dist = x_dist_init; % reset trait distribution at the start of each simulation
        x_mean(1,1:3) = mean(x_dist(:,1:3)); % initial mean trait
        x_var(1,1:3) = var(x_dist(:,1:3)); % initial variance in trait

        infected_dist = []; % open up empty vector for the epsilons of infected hosts

    count = 1; % start counter to index steps while inside loop
    while t(count) < t_max
        if size(x_dist,1) > 0 % as long as population size is > 0, pick another individual
            whosnext = randi(size(x_dist,1),1); % randomly choose individual from the vector
            p_next = x_dist(whosnext,1); % pick the p for that individual
            range_mean_next = x_dist(whosnext,2); % pick the range mean for that individual
            epsilon_next = max([p_next - (16*p_next^3*(IS_host_1 - range_mean_next).^2)./(9*IS_par_A.^2) 0]); % calc epsilon
            if size(infected_dist,1) > 0 % identify the IS traits of the infecting parasite
                which_infected = randi(size(infected_dist,1)); % randomly choose individual from the vector
                infected_p_next = infected_dist(which_infected,1);
                infected_rm_next = infected_dist(which_infected,2);
            end
        end

        % set up rates of each possible event, given by ODE in RSIZ_model.m
        % birth of prey
            b_R = r*R(count);
        % natural death of prey
            d_R_nat = r*R(count)^2/K;
        % consumption of resource
            d_R_pred = a*R(count)*(S(count) + I(count));
        % birth of a consumer (into susceptible class)
            b_S = e*a*R(count)*(S(count) + I(count));
        % natural death of consumer
            d_S_nat = ms*S(count);
        % parasite induced death of consumer
            d_I_par = (ms + v)*I(count);
        % infection of susceptible
            t_S_to_I = epsilon_next*a*Z(count)*S(count);
        % consumption of parasite
            d_Z_pred = a*Z(count)*(S(count) + I(count));
        % natural death of parasite
            d_Z_nat = mz*Z(count);
        % birth (shedding) rate of parsites
            b_Z = sigma*I(count);

    % sum the events to make wheel of fortune
        CS_vector = cumsum([b_R d_R_nat d_R_pred b_S d_S_nat d_I_par t_S_to_I d_Z_pred d_Z_nat b_Z]);
        Slice_widths = CS_vector./CS_vector(end);
        LI = rand < Slice_widths;
        Event_index = find(LI,1,'first');

    % now choose actual events
        if Event_index == 1 % choose birth of prey
            R(count+1) = R(count)+1; % add a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant
            Z(count+1) = Z(count); % hold parasite population constant

        elseif Event_index == 2 % choose natural death of prey
            R(count+1) = R(count)-1; % remove a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant
            Z(count+1) = Z(count); % hold parasite population constant

        elseif Event_index == 3 % choose death of prey by predator
            R(count+1) = R(count)-1; % remove a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant
            Z(count+1) = Z(count); % hold parasite population constant

        elseif  Event_index == 4 % choose birth of susceptible
            R(count+1) = R(count); % remove a resource individual
            S(count+1) = S(count)+1; % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant
            Z(count+1) = Z(count); % hold parasite population constant

        elseif  Event_index == 5 % choose natural death of susceptible
            R(count+1) = R(count); % remove a resource individual
            S(count+1) = S(count)-1; % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant
            Z(count+1) = Z(count); % hold parasite population constant

        elseif  Event_index == 6 % choose death of infected
            R(count+1) = R(count); % remove a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            Z(count+1) = Z(count); % hold parasite population constant

             infected_dist = infected_dist([1:which_infected-1,which_infected+1:end],:); % reduce dist by lost individual
             I(count+1) = size(infected_dist,1); % recount the infecteds

        elseif  Event_index == 7 % infection of susceptible
            R(count+1) = R(count); % hold resource population constant
            S(count+1) = S(count)-1; % remove a susceptible individual
            Z(count+1) = Z(count); % hold parasite population constant
            I(count+1) = I(count)+1; % hold parasite population constant

            infected_dist(size(infected_dist,1)+1,1) = p_next;
            infected_dist(size(infected_dist,1),2) = range_mean_next;
            I(count+1) = size(infected_dist,1); % recount the infecteds

        elseif  Event_index == 8 % choose death of a parasite by consumption
            R(count+1) = R(count); % remove a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant

            x_dist = x_dist([1:whosnext-1,whosnext+1:end],:); % reduce dist by lost individual
            Z(count+1) = size(x_dist,1); % hold parasite population constant

        elseif  Event_index == 9 % choose natural death of parasite
            R(count+1) = R(count); % remove a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant

            x_dist = x_dist([1:whosnext-1,whosnext+1:end],:); % reduce dist by lost individual
            Z(count+1) = size(x_dist,1); % hold parasite population constant

        elseif  Event_index == 10 % choose birth of parasite
            R(count+1) = R(count); % remove a resource individual
            S(count+1) = S(count); % hold susceptible population constant
            I(count+1) = I(count); % hold infected population constant

            % offspring trait distribution mean
            x_parent(1) = (1-h_2)*mean(x_dist_init(:,1)) + h_2*infected_p_next;
            x_parent(2) = (1-h_2)*mean(x_dist_init(:,2)) + h_2*infected_rm_next;
            % offspring trait distribution std
            off_std(1) = sqrt(1-h_2^2)*((1-h_2)*std(x_dist_init(:,1))+h_2*std(x_dist(:,1)));
            off_std(2) = sqrt(1-h_2^2)*((1-h_2)*std(x_dist_init(:,2))+h_2*std(x_dist(:,2)));

            MU_p = log(x_parent(1)^2 / sqrt(off_std(1)^2+x_parent(1)^2));
            MU_rm = log(x_parent(2)^2 / sqrt(off_std(2)^2+x_parent(2)^2));

            SIGMA_p = sqrt(log(off_std(1)^2/x_parent(1)^2 + 1));
            SIGMA_rm = sqrt(log(off_std(2)^2/x_parent(2)^2 + 1));

            x_dist(size(x_dist,1)+1,1) = lognrnd(MU_p,SIGMA_p,1,1); % pick offspring from lognormal dist
            x_dist(size(x_dist,1),2) = lognrnd(MU_rm,SIGMA_rm,1,1); % pick offspring from lognormal dist
            x_dist(size(x_dist,1),3) = max([x_dist(end,1) - (16*x_dist(end,1)^3*(IS_host_1 - x_dist(end,2))^2)/(9*IS_par_A^2) 0]); % calc epsilon

            Z(count+1) = size(x_dist,1); % recount the number of parasites
        end

            x_mean(count+1,:) = mean(x_dist(:,1:3)); % calculate new mean trait
            x_var(count+1,:) = var(x_dist(:,1:3)); % calculate new variance trait

        t(count+1) = t(count) + exp(-1/CS_vector(end))/CS_vector(end); % updating time
        count = count+1;
    end

    % find standardized times and corresponding densities (need for ci's)
    for q = 1:num_time_steps
        val = stand_times(q); %value to find
        tmp = abs(t-val);
        [idx idx] = min(tmp); %index of closest value
        closest = t(idx); %closest value
        R_stand(i,q) = R(idx); % resource abundance at standard time
        S_stand(i,q) = S(idx); % susceptible abundance at standard time
        I_stand(i,q) = I(idx); % infected abundance at standard time
        Z_stand(i,q) = Z(idx); % parasite abundance at standard time
        p_stand(i,q) = x_mean(idx,1); % trait values at standard time
        rm_stand(i,q) = x_mean(idx,2); % trait values at standard time
        epsilon_stand(i,q) = x_mean(idx,3); % trait values at standard time
    end

% figure(1);
%     box on; hold on;
%     h1 = plot(t,Z,'-k','LineWidth',1);
%     h2 = plot(t,I,'-r','LineWidth',1);
%     h3 = plot(t,S,'-g','LineWidth',1);
%     h4 = plot(t,R,'-b','LineWidth',1);
%     legend([h1 h2 h3 h4],'Z','I','S','R');
end

% calculate ci's for time series
    upper_ci_level = 75; % choose ci levels
    lower_ci_level = 25; % choose ci levels
    % resource abundance
        test(:,:) = R_stand(:,:);
        ci_R_up = prctile(test,lower_ci_level);
        ci_R_down = prctile(test,upper_ci_level);
        median_R = prctile(test,50);
    % susceptible abundance
        test(:,:) = S_stand(:,:);
        ci_S_up = prctile(test,lower_ci_level);
        ci_S_down = prctile(test,upper_ci_level);
        median_S = prctile(test,50);
    % infected abundance
        test(:,:) = I_stand(:,:);
        ci_I_up = prctile(test,lower_ci_level);
        ci_I_down = prctile(test,upper_ci_level);
        median_I = prctile(test,50);
    % parasite abundance
        test(:,:) = Z_stand(:,:);
        ci_Z_up = prctile(test,lower_ci_level);
        ci_Z_down = prctile(test,upper_ci_level);
        median_Z = prctile(test,50);
    % parasite trait p
        test(:,:) = p_stand(:,:);
        ci_p_up = prctile(test,lower_ci_level);
        ci_p_down = prctile(test,upper_ci_level);
        median_p = prctile(test,50);
    % parasite trait range mean
        test(:,:) = rm_stand(:,:);
        ci_rm_up = prctile(test,lower_ci_level);
        ci_rm_down = prctile(test,upper_ci_level);
        median_rm = prctile(test,50);
    % parasite trait range mean
        test(:,:) = epsilon_stand(:,:);
        ci_epsilon_up = prctile(test,lower_ci_level);
        ci_epsilon_down = prctile(test,upper_ci_level);
        median_epsilon = prctile(test,50);

%% Solve ODE

realized_epsilon = mean(epsilon_stand(:,1)); % get the average of the realized starting epsilons

ode = @(t,y) RSIZ_model(t,y,sigma,mz,a,realized_epsilon,mi,e,ms,r,K); % compile function and call
    [t1,y1] = ode45(ode, tspan, y0); % return time and population density vectors

%% plotting

figure(2); clf(2); % plot medians and ci's overtop individual lines
    subplot(4,2,1); box on;
        jbfill(stand_times,ci_Z_up,ci_Z_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        h1 = plot(stand_times,median_Z,'-k','LineWidth',2);
        h2 = plot(t1,y1(:,1),'--k','LineWidth',2);
        ylabel('Parasite abundance','FontSize',12);
        %axis([0 t_max 0 250]);
    subplot(4,2,3); box on;
        jbfill(stand_times,ci_I_up,ci_I_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        plot(stand_times,median_I,'-r','LineWidth',2);
        plot(t1,y1(:,2),'--k','LineWidth',2);
        ylabel('Infected abundance','FontSize',12);
        %axis([0 t_max 0 25]);
    subplot(4,2,5); box on;
        jbfill(stand_times,ci_S_up,ci_S_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        plot(stand_times,median_S,'-g','LineWidth',2);
        plot(t1,y1(:,3),'--k','LineWidth',2);
        %axis([0 t_max 8000 12000]);
        ylabel('Susceptible abundance','FontSize',12);
    subplot(4,2,7); box on;
        jbfill(stand_times,ci_R_up,ci_R_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        plot(stand_times,median_R,'-b','LineWidth',2);
        plot(t1,y1(:,4),'--k','LineWidth',2);
        ylabel('Prey abundance','FontSize',12);
    subplot(4,2,2); box on;
        jbfill(stand_times,ci_p_up,ci_p_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        plot(stand_times,median_p,'-b','LineWidth',2);
        plot([0 t_max],[IS_par_p IS_par_p],'--k');
        ylabel('p','FontSize',12);
    subplot(4,2,4); box on;
        jbfill(stand_times,ci_rm_up,ci_rm_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        plot(stand_times,median_rm,'-b','LineWidth',2);
        plot([0 t_max],[IS_par_range_mean IS_par_range_mean],'--k');
        ylabel('rm','FontSize',12);
    subplot(4,2,6); box on;
        jbfill(stand_times,ci_epsilon_up,ci_epsilon_down,[0.5 0.5 0.5],'w',1,0.2);
        hold on;
        plot(stand_times,median_epsilon,'-b','LineWidth',2);
        plot([0 t_max],[realized_epsilon realized_epsilon],'--k');
        ylabel('epsilon','FontSize',12);


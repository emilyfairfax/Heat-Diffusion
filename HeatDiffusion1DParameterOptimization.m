% 1D Heat Diffusion Model
% written by Emily Fairfax 
% March 16th, 2016

clear all
figure(1)
clf
figure(2)
clf
%% Model Optimization
%Number of model runs in each parameter
number_of_runs = 100;

%Model Parameters for Manipulation
min_years = 0; %least time since temp change to test
max_years = 200; %longest time since temp change to test
year_range = max_years - min_years; %range of years elapsed to test
year_param_step = year_range/number_of_runs; %step between years 

min_Ts = -7; %lowest "hot" surface temperature
max_Ts = 0; %highest "hot" surface temperature
Ts_range = abs(max_Ts - min_Ts); %range of temperatures to test
Ts_param_step = Ts_range/number_of_runs; %step between temps

fprintf('Running model iterations for parameter optimization...')
for j = 2:number_of_runs %loop through time parameter
    years_elapsed(1) = min_years; %initial time is the min tim
    years_elapsed(j) = years_elapsed(j-1) + year_param_step; % years since the warming step added to min time
    for l = 2:number_of_runs %loop through all the temp parameters for each loop of time parameter to get all temp/time combos
        Ts_hot(1) = min_Ts; %initial temp is the low temp
        Ts_hot(l) =Ts_hot(l-1) + Ts_param_step; %new surface temperature after heating temperature step change
        
        %% Model within Optimization
        %all time parameter related terms are tracked with j
        %all temperature parameter related terms are tracked with l
        %% Initialize

        %Data
            %Load the Experimental Data
            load cape_thompson_copy.dat; %this is the data file in text format, 2 columns
            data_depth = cape_thompson_copy(:,1); %make array of depth data from text file (1st column in file)
            data_T = cape_thompson_copy(:,2); %make array of temperature data from text file (2st column in file)

            %Experimental Data Geotherm Determination
            data_bottom_depth = data_depth(end); %deepest depth point in the data file (least affected by temp change at surface)
            data_bottom_T = data_T(end); %deepest temperature point in the data file (least affected by temp change at surface)
            slope = diff(data_T)./diff(data_depth); %calculate all slopes between data points
            dTdz_data = mean(slope(end-5:end)); %take the average of the slopes in the last handful of data points (ones likely unaffected by surface T change)
            Ts_initial = data_bottom_T -(dTdz_data*data_bottom_depth); %solve for intercept at surface, Ts_initial
            data_geotherm = Ts_initial + (dTdz_data*data_depth); %final equation of straight line geotherm for the data set

        %Constants
        k = 2.2; %Thermal Conductivity of the ground in W/mK
        rho = 2700; %Density of the ground in kg/m^3
        Cp = 1000; % Specific Heat Capacity

        %Time and Space
            %Space Array in Z
            N = 200; % number of nodes
            maxdepth = 400; % maximum depth in m, slightly larger than deepest experimental depth
            dz = maxdepth/N; % z spacing in meters, determined based on desired number of nodes
            z = 0:dz:dz*(N-1); %depth array, need the N-1 to make correct length accounting for the 0 at the beginning 

            %Time Array
            tmax(j) = years_elapsed(j)*24*3600*365; % years to run the code, converted to seconds
            years_time_step = 1/24; % time increments to use in loop, measured in fraction of a year
            dt = years_time_step*24*3600*365; %increment of time step, converted to seconds
            t = 0:dt:tmax(j); %time array
            imax(j) = length(t); %for use in time loop max time
            nplots=100; %number of plots generated
            tplot(j) = tmax(j)/nplots; %time at which to generate plots

        %Variable Arrays
        Q = zeros(N,1); %Flux Array, empty
        T = zeros(N,1); %Temperature Array, empty

        %Boundaries and Initial Conditions
        dTdz_initial = dTdz_data; %create a copy of calculated slope of data
        T = Ts_initial+(dTdz_initial*z); %initial condition of temperature is the geotherm extracted from the base of the experimental data
        T_initial_condition=T; %create a copy of the original geotherm
        dTdz(N) = dTdz_data; % keep model heat flux from base same as heat flux determined from data set

        %% Run

        for i = 1:imax(j) 

            %Surface Temperature
            T(1) = Ts_hot(l); %Update the surface temperature to be the step changed temperature, hold constant

            %Heat Diffusion
            dTdz(1:N-1) = diff(T)/dz; %calculate temperature gradient between each cell
            Q = -k*dTdz; %calculate heat flux via diffusion equation (no source or sink terms)
            dqdz = diff(Q)/dz; %calculate rate of temperature change between cells

            %Update Temperatures
            T(2:N) = T(2:N) - (1/(rho*Cp))*dqdz*dt; %Update all nodes below the surface

        end

        %% Finalize
        %Perform Chi Squared Analysis
        T_model = interp1(z,T,data_depth); %interpolate temperatures from the model at depths from the experimental data file
        chi_squared(j,l) = sum((data_T - T_model).^2); %chi squared test, print out chi squared in window
    end
end

%% Evaluate Chi-Squared Values and Find Best Parameters
chi_squared_plot = chi_squared(2:end,2:end); %only plot chi squared in the loop range
figure(2)
image([2 number_of_runs], [2 number_of_runs], chi_squared_plot, 'CDataMapping', 'scaled') %make the heat map
title('Chi Squared Map')

[M,I] = min(chi_squared_plot(:)); %find the lowest chi squared value
[row,col] = ind2sub(size(chi_squared_plot),I); %locate it's indices in the matrix
Best_Chi_Squared = chi_squared_plot(row,col) %locate and print the value of the lowest chi squared
Best_Elapsed_Time = years_elapsed(row+1) %locate and print the best elapsed time parameter
Best_Ts = Ts_hot(col+1) %locate and print the best hot surface temperature parameter



%% Run Model with Best Parameters
%% Initialize
fprintf('Plotting optimized model...')

%Parameters from Optimization
best_years_elapsed = Best_Elapsed_Time;
best_Ts_hot = Best_Ts;

%Data
    %Load the Experimental Data
    best_data_depth = cape_thompson_copy(:,1); %make array of depth data from text file (1st column in file)
    best_data_T = cape_thompson_copy(:,2); %make array of temperature data from text file (2st column in file)

    %Experimental Data Geotherm Determination
    best_data_bottom_depth = best_data_depth(end); %deepest depth point in the data file (least affected by temp change at surface)
    best_data_bottom_T = best_data_T(end); %deepest temperature point in the data file (least affected by temp change at surface)
    best_slope = diff(best_data_T)./diff(best_data_depth); %calculate all slopes between data points
    best_dTdz_data = mean(best_slope(end-5:end)); %take the average of the slopes in the last handful of data points (ones likely unaffected by surface T change)
    best_Ts_initial = best_data_bottom_T -(best_dTdz_data*best_data_bottom_depth); %solve for intercept at surface, Ts_initial
    best_data_geotherm = best_Ts_initial + (best_dTdz_data*best_data_depth); %final equation of straight line geotherm for the data set

%Time and Space
    %Space Array in Z
    best_N = 200; % number of nodes
    best_maxdepth = 400; % maximum depth in m, slightly larger than deepest experimental depth
    best_dz = best_maxdepth/best_N; % z spacing in meters, determined based on desired number of nodes
    best_z = 0:best_dz:best_dz*(best_N-1); %depth array, need the N-1 to make correct length accounting for the 0 at the beginning 

    %Time Array
    best_nplots=100; %number of plots generated
    best_tmax = best_years_elapsed*24*3600*365; % years to run the code, converted to seconds
    best_years_time_step = best_years_elapsed/(best_nplots*50); % time increments to use in loop, measured in fraction of a year
    best_dt = best_years_time_step*24*3600*365; %increment of time step, converted to seconds
    best_t = 0:best_dt:best_tmax; %time array
    best_imax = length(best_t); %for use in time loop max time
    best_tplot = (1/best_nplots)*best_tmax; %time at which to generate plots

%Variable Arrays
best_Q = zeros(best_N,1); %Flux Array, empty
best_T = zeros(best_N,1); %Temperature Array, empty

%Boundaries and Initial Conditions
best_dTdz_initial = best_dTdz_data; %create a copy of calculated slope of data
best_T = best_Ts_initial+(best_dTdz_initial*best_z); %initial condition of temperature is the geotherm extracted from the base of the experimental data
best_T_initial_condition=best_T; %create a copy of the original geotherm
best_dTdz(best_N) = best_dTdz_data; % keep model heat flux from base same as heat flux determined from data set

%% Run

for p = 1:best_imax

    %Surface Temperature
    best_T(1) = best_Ts_hot; %Update the surface temperature to be the step changed temperature, hold constant

    %Heat Diffusion
    best_dTdz(1:best_N-1) = diff(best_T)/best_dz; %calculate temperature gradient between each cell
    best_Q = -k*best_dTdz; %calculate heat flux via diffusion equation (no source or sink terms)
    best_dqdz = diff(best_Q)/best_dz; %calculate rate of temperature change between cells

    %Update Temperatures
    best_T(2:best_N) = best_T(2:best_N) - (1/(rho*Cp))*best_dqdz*best_dt; %Update all nodes below the surface
    %For plotting
    str = {'years ago'};
    %Plot the results as a movie
            if rem(best_t(p),best_tplot)==0 
                figure(1)
                clf
                plot(best_T,best_z,'Color', [0 .447 .333],'linewidth', 1)
                hold all
                plot(best_T_initial_condition,best_z,'--k','linewidth', 1)
                plot(best_data_T,best_data_depth,'o','MarkerEdgeColor',[0 .337 .447],'MarkerFaceColor', [.875 .588 .757],'MarkerSize', 6)

                xlabel('Temperature (°C)','fontname','arial','fontsize', 21)
                ylabel('Depth (m)', 'fontname', 'arial', 'fontsize', 21)
                ht=text(-3,100,['  ',num2str(ceil(best_years_elapsed-round(best_t(p)/(24*3600*365)))), '    years ago'],'fontsize',18); %print time in animation
                set(gca, 'fontsize', 18, 'fontname', 'arial') 
                set(gca, 'YDIR', 'reverse') %change y-axis direction
                drawnow
                pause(0.1)
                hold off
            end

end
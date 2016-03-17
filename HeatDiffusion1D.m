% 1D Heat Diffusion Model
% written by Emily Fairfax 
% March 16th, 2016

clear all
figure(1)
clf
%% Initialize

%Parameters
years_elapsed = 93.2;
Ts_hot = -5.26;

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
    tmax = years_elapsed*24*3600*365; % years to run the code, converted to seconds
    years_time_step = 1/100; % time increments to use in loop, measured in fraction of a year
    dt = years_time_step*24*3600*365; %increment of time step, converted to seconds
    t = 0:dt:tmax; %time array
    imax = length(t); %for use in time loop max time
    nplots=100; %number of plots generated
    tplot = 1/100*tmax; %time at which to generate plots

%Variable Arrays
Q = zeros(N,1); %Flux Array, empty
T = zeros(N,1); %Temperature Array, empty

%Boundaries and Initial Conditions
dTdz_initial = dTdz_data; %create a copy of calculated slope of data
T = Ts_initial+(dTdz_initial*z); %initial condition of temperature is the geotherm extracted from the base of the experimental data
T_initial_condition=T; %create a copy of the original geotherm
dTdz(N) = dTdz_data; % keep model heat flux from base same as heat flux determined from data set

%% Run

for i = 1:imax

    %Surface Temperature
    T(1) = Ts_hot; %Update the surface temperature to be the step changed temperature, hold constant

    %Heat Diffusion
    dTdz(1:N-1) = diff(T)/dz; %calculate temperature gradient between each cell
    Q = -k*dTdz; %calculate heat flux via diffusion equation (no source or sink terms)
    dqdz = diff(Q)/dz; %calculate rate of temperature change between cells

    %Update Temperatures
    T(2:N) = T(2:N) - (1/(rho*Cp))*dqdz*dt; %Update all nodes below the surface

    %Plot the results as a movie
            if(rem(t(i),tplot)==0) 
                figure(1)
                plot(T,z,'Color', [0 .447 .333],'linewidth', 1)
                hold on
                plot(T_initial_condition,z,'--k','linewidth', 1)
                plot(data_T,data_depth,'o','MarkerEdgeColor',[0 .337 .447],'MarkerFaceColor', [.875 .588 .757],'MarkerSize', 6)

                xlabel('Temperature (°C)','fontname','arial','fontsize', 21)
                ylabel('Depth (m)', 'fontname', 'arial', 'fontsize', 21)
                set(gca, 'fontsize', 18, 'fontname', 'arial') 
                set(gca, 'YDIR', 'reverse') %change y-axis direction
                ht=text(-3,100,['  ',num2str(ceil(years_elapsed-round(t(i)/(24*3600*365)))), '    years ago'],'fontsize',18); 
                pause(0.1)
                hold off
            end

end

%% Finalize
%Perform Chi Squared Analysis
T_model = interp1(z,T,data_depth); %interpolate temperatures from the model at depths from the experimental data file
chi_squared = sum((data_T - T_model).^2) %chi squared test, print out chi squared in window
 
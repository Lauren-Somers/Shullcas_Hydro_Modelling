%% Glacier_Melt_Module_bands_2019
% by Lauren Somers, McGill University 
% January 18, 2019
%
% Calculate melt water production from the Huaytapallana Glaciers using a
% modified temperature index approach
% Further improvements on Glacier_Melt_Module_bands.m in order to simplify the overall calibration approach
%
% Divides Huaytapallana into 3 drainage sections: Chuspicocha, Lashuntay
% and Duraznoyoc, and calculates melt and area for each one

%% Import Data

%Import observed discharge for Chuspi Outlet Channel
load('D:\Glacier_Melt\Stream_flow_records\Chuspi_Out_Time.mat');
load('D:\Glacier_Melt\Stream_flow_records\Chuspi_Out_discharge_exp.mat');
Chuspi_Out_Time = datetime(datestr(Chuspi_Out_Time));

% Import Huaytapallana met data
met_data = xlsread('Met_Data_new_withoutspike_Updated_Aug2_2018.xlsx');
Huay_precip = met_data(:,10);
Huay_Tmin = met_data(:,12);
Huay_Tmax = met_data(:,11);
Huay_Tavg = (Huay_Tmin + Huay_Tmax)/2;
Huay_elev = 4684; % Elevation of the huaytapallana station masl

%% Define parameters:

% Melt factors
MF_D = 0.0028; % Melt factor for the dry season 0.0028
MF_R = 0.0030; % Melt factor for the rainy season 0.0030
T_crit = -1; % temperature at which melt starts

% Temperature and precip lapse rates
lapse_T = 0.62; % Average temperature lapse rate per 100m
lapse_P = 0.08; % Precip lapse rate as a proportion

% Sublimation factors
SF_D = 0.0003 ;% Sublimation factor for the dry season ~1/10 of melt factors(More sublimation during dry season)
SF_R = 0.0002; % Sulimation factor for the rainy season

p_ice = 0.8; % Relative density of ice
Thickness_i = 55; %initial ice thickness for all glaciers
n_bands = 8; % Number of 100m elevation bands

%Time variable
Start_time = datetime(1968,1,1);
End_time = datetime(2018,8,1);
Time_mod = (Start_time:End_time)';
num_steps = length(Time_mod);

%% CHUSPICOCHA

% 2015 Hypsometric data
%G_A_Chuspi = 1419607; % Total glacier area in square meters contributing to Chuspicocha Outlet
% Band_A = [234050, 390092, 341315, 348279, 105875]; % Band areas for bands 1:5 in square m
% Cum_A = [1419611, 1185561, 795469, 454154, 105875, 0]; %1419611 Cumulative band areas (area in and above given band)
% A_HRUs = 2262293.912; % Area of the three Glacier HRUs for Chuspi Outlet (HRUs 3,11,12)
% Band_elevs = [5050, 5150, 5250, 5350, 5450];

% 1989 Hypsometric data
% total chuspi glacier area = 1 930 778 square m 
Band_A = [192173, 132926, 146870, 263994, 477792, 346724, 351084, 107276]; %Areas of bands calculated in ArcGIS
Cum_A = [2018841, 1826667, 1693741, 1546871, 1282876, 805084, 458359, 107276, 0]; %Cumulative areas
A_HRUs = 2262293.912 %Area of the HRUs over which glacier melt will be applied
Band_elevs = [4750, 4850, 4950, 5050, 5150, 5250, 5350, 5450]; %Elevation of mittle of bands

% Create the output variables for the melt loop 
Melt = zeros(num_steps,n_bands + 1); % Bands 1-5 plus the total
Sublimation = zeros(num_steps,n_bands + 1); % Bands 1-5 plus the total
Accumulation = zeros (num_steps,n_bands + 1); % Bands 1-5 plus the total
Mass = zeros(num_steps,1); % For the glacier as a whole
Area = zeros(num_steps,1);
Thickness = zeros(num_steps,1);
A_rain = zeros(num_steps,1);
A_snow = zeros(num_steps,1);

% Initial conditions in 1968 for spin-up
Area(1,1) = Cum_A(1); %initial area in square meters (for non glacier version use 1) 
Thickness(1,1) = Thickness_i; %initial thickness in meters
Mass(1,1) = Area(1,1)*Thickness(1,1)*p_ice; % initial mass of ice in tonnes

% Run the loop!
for i = 2:num_steps
      for band = 1:n_bands
          temp = Huay_Tavg(i)-(lapse_T*(Band_elevs(band)-Huay_elev)/100); % Use lapse rate to adjust temp
          precip = Huay_precip(i)*(lapse_P*((Band_elevs(band)-Huay_elev)/100)+1); % Use lapse rate to adjust precip
            if   Area(i-1)> Cum_A(band) % is the glacier front beyond this elev band?
                 area = Band_A(band); % If so, just use the current area
            elseif Area(i) <= Cum_A(band)
                area = Area(i-1)-Cum_A(band+1); % if not, use total area minus area of the other bands, cumulative are of the next highest band
                if area<0
                   area = 0 ;
                end
            end 
            
         % Choose Melt Factor and sublimation factor         
            if  month(Time_mod(i))>= 5 && month(Time_mod(i))<10 % Choose a melt factor depending on the month
                MF = MF_D; % then it's the dry season!
                SF = SF_D;
            else
                MF = MF_R; % Then it's the rainy season!
                SF = SF_R;
            end
            
         % Calc Melt
           if temp >= T_crit
              Melt(i , band) = (temp-(-1)) * MF * area; % Need to add a sumnmer v winter if
           else
              Melt(i,band) = 0;
           end
           
         % Calc Sublimation
           if temp >= -2
              Sublimation(i , band) = (temp-(-2)) * SF * area; % Need to add a summer v winter if
           else
              Sublimation(i,band) = 0;
           end
           
       % Calc Accumulation
           if temp < 0 
                Accumulation (i,band) = (precip/1000) * area; % accumulation in m of SWE to a given elev band.
                A_snow(i) = A_snow(i)+area;
           else
                Accumulation (i,band) = 0;
           end
      end
    A_rain(i) = A_HRUs - A_snow(i); %Rain area cannot be negative
    if A_rain(i) < 0 
       A_rain(i) = 0;
    end
    Melt(i,n_bands+1) = sum(Melt (i,1:n_bands));
    Accumulation (i,n_bands+1) = sum(Accumulation (i,1:n_bands));
    %Sublimation (i,6) = sum(Sublimation(i,1:5));
      
    % Mass balance for the whole galcier!  
    Mass(i,1) = Mass(i-1,1) + Accumulation(i,n_bands+1) - Melt(i,n_bands+1) - Sublimation (i,n_bands+1); % - sub_rat*Melt(i,1); % in cubic meters of water of kg of ice (cubic m of swe)
    Area(i,1) = Area(i-1,1) + (((Mass(i,1) - Mass(i-1))/p_ice)*0.8)/Thickness(i-1,1); % old area plus change in mass converted to volume
    Thickness(i,1)= Thickness(i-1,1) + (((Mass(i,1) - Mass(i-1))/p_ice)*0.2)/Area(i,1);% 20% goes to thinning and 80% goes to edge ablation
end

% Plot the melt water
figure;
plot(Time_mod,Melt(:,n_bands+1));
title('Melt')

% Plot the melt water
figure;
plot(Time_mod,Area);
title('Area')

% Create the precip data that we can input into GSFLOW
Precip_and_melt_Chuspi = ((((Huay_precip/1000).* A_rain*(1+lapse_P))+ Melt(:,n_bands+1))/A_HRUs)/(25.4/1000); % Apply lapse rate to huaytapallan precip assumin the rain area is aroun 4700 masl.
Area_Chuspi = Area;

%% LAZOHUNTAY

% 1989 Hypsometric data
% total chuspi glacier area = 1 930 778 square m 
Band_A = [161498, 120842, 284444, 462919, 525200, 487088, 317647, 107276];
Cum_A = [2466915, 2305417, 2184574, 1900130, 1437211, 912011, 424923, 107276, 0]; %Correct
A_HRUs = 2514449.5; % Area of the three Glacier HRUs for Lazo (HRUs 7,13)
Band_elevs = [4750, 4850, 4950, 5050, 5150, 5250, 5350, 5450];

% Create the output variables for the melt loop 
Melt = zeros(num_steps,n_bands + 1); % Bands 1-5 plus the total
Sublimation = zeros(num_steps,n_bands + 1); % Bands 1-5 plus the total
Accumulation = zeros (num_steps,n_bands + 1); % Bands 1-5 plus the total
Mass = zeros(num_steps,1); % For the glacier as a whole
Area = zeros(num_steps,1);
Thickness = zeros(num_steps,1);
A_rain = zeros(num_steps,1);
A_snow = zeros(num_steps,1);

% Initial conditions in 1968
Area(1,1) = Cum_A(1); %initial area in square meters (for non glacier version use 1) 
Thickness(1,1) = Thickness_i; %initial thickness in meters
Mass(1,1) = Area(1,1)*Thickness(1,1)*p_ice; % initial mass of ice in tonnes

% Run the loop!
for i = 2:num_steps
      for band = 1:n_bands
          temp = Huay_Tavg(i)-(lapse_T*(Band_elevs(band)-Huay_elev)/100); % Use lapse rate to adjust temp
          precip = Huay_precip(i)*(lapse_P*((Band_elevs(band)-Huay_elev)/100)+1); % Use lapse rate to adjust precip
            if   Area(i-1)> Cum_A(band) % is the glacier front beyond this elev band?
                 area = Band_A(band); % If so, just use the current area
            elseif Area(i) <= Cum_A(band)
                area = Area(i-1)-Cum_A(band+1); % if not, use total area minus area of the other bands, cumulative are of the next highest band
                if area<0
                   area = 0 ;
                end
            end 
            
         % Choose Melt Factor and sublimation factor         
            if  month(Time_mod(i))>= 5 && month(Time_mod(i))<10 % Choose a melt factor depending on the month
                MF = MF_D; % then it's the dry season!
                SF = SF_D;
            else
                MF = MF_R; % Then it's the rainy season!
                SF = SF_R;
            end
            
         % Calc Melt
           if temp >= T_crit
              Melt(i , band) = (temp-(-1)) * MF * area; % Need to add a sumnmer v winter if
           else
              Melt(i,band) = 0;
           end
           
         % Calc Sublimation
           if temp >= -2
              Sublimation(i , band) = (temp-(-2)) * SF * area; % Need to add a summer v winter if
           else
              Sublimation(i,band) = 0;
           end
           
       % Calc Accumulation
           if temp < 0 
                Accumulation (i,band) = (precip/1000) * area; % accumulation in m of SWE to a given elev band.
                A_snow(i) = A_snow(i)+area;
           else
                Accumulation (i,band) = 0;
           end
      end
    A_rain(i) = A_HRUs - A_snow(i); %Rain area cannot be negative
    if A_rain(i) < 0 
       A_rain(i) = 0;
    end
    Melt(i,n_bands+1) = sum(Melt (i,1:n_bands));
    Accumulation (i,n_bands+1) = sum(Accumulation (i,1:n_bands));
    %Sublimation (i,6) = sum(Sublimation(i,1:5));
      
    % Mass balance for the whole galcier!  
    Mass(i,1) = Mass(i-1,1) + Accumulation(i,n_bands+1) - Melt(i,n_bands+1) - Sublimation (i,n_bands+1); % - sub_rat*Melt(i,1); % in cubic meters of water of kg of ice (cubic m of swe)
    Area(i,1) = Area(i-1,1) + (((Mass(i,1) - Mass(i-1))/p_ice)*0.8)/Thickness(i-1,1); % old area plus change in mass converted to volume
    Thickness(i,1)= Thickness(i-1,1) + (((Mass(i,1) - Mass(i-1))/p_ice)*0.2)/Area(i,1);% 20% goes to thinning and 80% goes to edge ablation
end

% Plot the melt water
figure;
plot(Time_mod,Melt(:,n_bands+1));
title('Melt')

% Plot the melt water
figure;
plot(Time_mod,Area);
title('Area')

% Create the precip data that we can input into GSFLOW
Precip_and_melt_Lazo = ((((Huay_precip/1000).* A_rain*(1+lapse_P))+ Melt(:,n_bands+1))/A_HRUs)/(25.4/1000); % Apply lapse rate to huaytapallan precip assumin the rain area is aroun 4700 masl.
Area_Lazo = Area;

%% DURAZNOYOC & CHUCO CHUSPI

% 1989 Hypsometric data
% total chuspi glacier area = 1 930 778 square m 
Band_A = [137329, 489877, 298388, 329993, 251910, 103181, 84985, 1393];
Cum_A = [1697055, 1559726, 1069849, 771461, 441469, 189559, 86378, 1393, 0]; %Correct
A_HRUs = 1608992; % Area of the three Glacier HRUs for Lazo (HRUs 7,13)
Band_elevs = [4750, 4850, 4950, 5050, 5150, 5250, 5350, 5450];

% Create the output variables for the melt loop 
Melt = zeros(num_steps,n_bands + 1); % Bands 1-5 plus the total
Sublimation = zeros(num_steps,n_bands + 1); % Bands 1-5 plus the total
Accumulation = zeros (num_steps,n_bands + 1); % Bands 1-5 plus the total
Mass = zeros(num_steps,1); % For the glacier as a whole
Area = zeros(num_steps,1);
Thickness = zeros(num_steps,1);
A_rain = zeros(num_steps,1);
A_snow = zeros(num_steps,1);

% Initial conditions in 1968
Area(1,1) = Cum_A(1); %initial area in square meters (for non glacier version use 1) 
Thickness(1,1) = Thickness_i; %initial thickness in meters
Mass(1,1) = Area(1,1)*Thickness(1,1)*p_ice; % initial mass of ice in tonnes

% Run the loop!
for i = 2:num_steps
      for band = 1:n_bands
          temp = Huay_Tavg(i)-(lapse_T*(Band_elevs(band)-Huay_elev)/100); % Use lapse rate to adjust temp
          precip = Huay_precip(i)*(lapse_P*((Band_elevs(band)-Huay_elev)/100)+1); % Use lapse rate to adjust precip
            if   Area(i-1)> Cum_A(band) % is the glacier front beyond this elev band?
                 area = Band_A(band); % If so, just use the current area
            elseif Area(i) <= Cum_A(band)
                area = Area(i-1)-Cum_A(band+1); % if not, use total area minus area of the other bands, cumulative are of the next highest band
                if area<0
                   area = 0 ;
                end
            end 
            
         % Choose Melt Factor and sublimation factor         
            if  month(Time_mod(i))>= 5 && month(Time_mod(i))<10 % Choose a melt factor depending on the month
                MF = MF_D; % then it's the dry season!
                SF = SF_D;
            else
                MF = MF_R; % Then it's the rainy season!
                SF = SF_R;
            end
            
         % Calc Melt
           if temp >= T_crit
              Melt(i , band) = (temp-(-1)) * MF * area; % Need to add a sumnmer v winter if
           else
              Melt(i,band) = 0;
           end
           
         % Calc Sublimation
           if temp >= -2
              Sublimation(i , band) = (temp-(-2)) * SF * area; % Need to add a summer v winter if
           else
              Sublimation(i,band) = 0;
           end
           
       % Calc Accumulation
           if temp < 0 
                Accumulation (i,band) = (precip/1000) * area; % accumulation in m of SWE to a given elev band.
                A_snow(i) = A_snow(i)+area;
           else
                Accumulation (i,band) = 0;
           end
      end
    A_rain(i) = A_HRUs - A_snow(i); %Rain area cannot be negative
    if A_rain(i) < 0 
       A_rain(i) = 0;
    end
    Melt(i,n_bands+1) = sum(Melt (i,1:n_bands));
    Accumulation (i,n_bands+1) = sum(Accumulation (i,1:n_bands));
    %Sublimation (i,6) = sum(Sublimation(i,1:5));
      
    % Mass balance for the whole galcier!  
    Mass(i,1) = Mass(i-1,1) + Accumulation(i,n_bands+1) - Melt(i,n_bands+1) - Sublimation (i,n_bands+1); % - sub_rat*Melt(i,1); % in cubic meters of water of kg of ice (cubic m of swe)
    Area(i,1) = Area(i-1,1) + (((Mass(i,1) - Mass(i-1))/p_ice)*0.8)/Thickness(i-1,1); % old area plus change in mass converted to volume
    Thickness(i,1)= Thickness(i-1,1) + (((Mass(i,1) - Mass(i-1))/p_ice)*0.2)/Area(i,1);% 20% goes to thinning and 80% goes to edge ablation
end

% Plot the melt water
figure;
plot(Time_mod,Melt(:,n_bands+1));
title('Melt')

% Plot the melt water
figure;
plot(Time_mod,Area);
title('Area')

% Create the precip data that we can input into GSFLOW
Precip_and_melt_Dur = ((((Huay_precip/1000).* A_rain*(1+lapse_P))+ Melt(:,n_bands+1))/A_HRUs)/(25.4/1000); % Apply lapse rate to huaytapallan precip assumin the rain area is aroun 4700 masl.
Area_Dur = Area;

%% Sum the areas all together and plot them agains the observations
Total_area = Area_Chuspi + Area_Lazo + Area_Dur
L_M_areas = readtable('Lopez-Moreno_Areas.xlsx')

figure
plot(Time_mod, Total_area);
title('Total Area')
hold on
scatter(L_M_areas.Date, L_M_areas.AreaInBasin_sq_M_)

%% New input files for GSFLOW

% Import the existing precip file (sans glacier)
Old_precip_file = xlsread('precip_updated_Aug3_2018.xlsx');
% 
%Modify the precip file to include glacier melt.
New_precip_file = Old_precip_file;

%CHUSPI:
New_precip_file(:,(3+6)) = Precip_and_melt_Chuspi ; % for HRU 3
New_precip_file(:,(11+6)) = Precip_and_melt_Chuspi ; % for HRU 11
New_precip_file(:,(12+6)) = Precip_and_melt_Chuspi ; % for HRU 12

%LAZO:
New_precip_file(:,(7+6)) = Precip_and_melt_Lazo ; % for HRU 7
New_precip_file(:,(13+6)) = Precip_and_melt_Lazo ; % for HRU 13

%DUR:
New_precip_file(:,(9+6)) = Precip_and_melt_Dur ; % for HRU 9
New_precip_file(:,(4+6)) = Precip_and_melt_Dur ; % for HRU 4
New_precip_file(:,(6+6)) = Precip_and_melt_Dur ; % for HRU 6

% Export the new precip file which includes glacier melt
% dlmwrite('precip_melt.day', New_precip_file, 'delimiter','\t');

Old = sum(Old_precip_file(:,9))
New = sum(New_precip_file(:,9))

Area_Chuspi(17330)
Total_area(17330)


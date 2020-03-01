% Wave field processing in the time domain
% By: Tiffany Vlaar 

clear all; 
close all; 
clc; 

fid1 = fopen('concrete2.dat','r+');
H = fread(fid1,6,'float64',0,'s'); % machineformat = 's' (since each number in the header is stored in binary format as a 64 bit big-endian real number)

%Define which part of the header (matrix H) corresponds to which data characteristic
Trace_spacing = H(1); % in m
Amount_of_traces = H(2);
Time_sampling_interval = H(3); %in s, 
Number_of_time_samples_per_trace = H(4);
Center_frequency_of_the_transmitted_signal = H(5); % in s^-1
Horizontal_distance_between_source_and_receiver = H(6);

%So the actual amount of data is given by H(2) and H(4).
%I want to put this data into 1 matrix using fread.
%Since I want to plot the data in a certain way later on I use sizedata = [H(4),H(2)]
data = fread(fid1,[Number_of_time_samples_per_trace,Amount_of_traces],'float64',0,'s');

frewind(fid1); %resets file position indicator to beginning of file
fclose(fid1); 

t = (0:1:Number_of_time_samples_per_trace-1)*Time_sampling_interval*(10^9); %time in nanoseconds
%since it now runs from 0, I have to substract 1 from the total amount
x = (0:1:Amount_of_traces-1)*Trace_spacing; %distance in meter
imagesc(x,t,data) %plots the data
colormap(flipud(gray)); %the extra flipud command reverses the colormap such that
%the smallest number is white, and the largest number black instead of vice versa

colorbar;
xlabel('Distance (in meter)', 'Fontsize', 12) 
ylabel('Time (in nanoseconds)', 'Fontsize',12) 
title('Raw data', 'Fontsize',12) 

t_exp = Horizontal_distance_between_source_and_receiver./300000000; %expected time in seconds

%As long as you didn't receive any data, the wave has not arrived yet so there's no onset of the direct wave.
%So using this fact I can find the amount of time samples per trace it takes till the signal is received 
%by counting how many time samples are needed in order for the data to be not equal to zero anymore.
a = 0; %counter parameter
for  k = 1:H(4) %k runs from 1 to the total possible amount of time samples
    if data(k,54) == 0 % as long as the data for a certain randomly chosen trace (it's for all traces the same) is not zero..
        a = a+1; 
    end
end

%Now a is the total amount of time samples till the direct wave arrives at the receiver
Time_samples_till_arrival = a;

t_direct = Time_samples_till_arrival*Time_sampling_interval; %the time it takes the direct wave to arrive at the receiver

if (t_exp - t_direct) ~= 0 %so if these times are not equal
    lag = t_exp - t_direct; %then there is a lag in the data
    lagt = ceil(lag/Time_sampling_interval) %the lag in the amount of time samples rounded to above
end

if lagt > 0 %So if the expected time is larger than the real time
    Adjusted_data = data(lagt:end,:); %Eliminates the lag in the data by 
                                      %removing all values with lower indices than of the time lag
    Adjusted_t = (0:1:Number_of_time_samples_per_trace-lagt-1)*Time_sampling_interval*(10^9); %Adjusted time in nanoseconds
    figure
    imagesc(x,Adjusted_t,Adjusted_data) 
    colormap(flipud(gray)); colorbar;
    xlabel('Distance (in meter)', 'Fontsize', 12) 
    ylabel('Adjusted Time (in nanoseconds)', 'Fontsize',12)
    title('Concrete Data which is timelag corrected', 'Fontsize',12) 
end

hold on 
text(0.02,2,'Pick two time levels between which','FontSize',10)
text(0.02,3,'you want to set all data to zero. ','FontSize',10)
[a,t] = ginput(2); %Gives you to chance to pick 2 time levels 

Adjusted_data2 = Adjusted_data;
%The following loop sets all data between the chosen 2 time levels to zero
if t(2)>t(1)
   Adjusted_data2((Adjusted_t > t(1)) & (Adjusted_t < t(2)),:) = 0;
   else if t(1)>t(2)
       Adjusted_data2((Adjusted_t > t(2)) & (Adjusted_t < t(1)),:) = 0;
   end
end

Median_data = median(Adjusted_data,2); %because I want to find the median value of each row
Expanded_median = Median_data*ones(1,Amount_of_traces); %to make the median the same size as the data set
Data_substracted_by_median = Adjusted_data - Expanded_median; %substract the median from the data

figure
imagesc(x,Adjusted_t,Adjusted_data2)
colormap(flipud(gray)); colorbar;
xlabel('Distance (in meter)', 'Fontsize', 12) 
ylabel('Adjusted Time (in nanoseconds)', 'Fontsize',12)
title('Data which is timelag and groundroll corrected', 'Fontsize',12) 

v_initial = 0.10; %initial guess for the velocity in meter/nanoseconds
hold on 
text(0.02,2,'Pick an apex of a diffraction hyperbola','FontSize',10)
[x_apex,t_apex] = ginput(1); %picking a random apex of the diffraction hyperbola

% d = t_apex*v_initial/2, x = x-xapex and x0 = Horizontal_distance_between_source_and_receiver
t_hyperbola = ((1/v_initial)*(sqrt((((t_apex*v_initial)/2).^2)+(((x-x_apex)-...
Horizontal_distance_between_source_and_receiver/2).^2))+sqrt((((t_apex*v_initial)/2).^2)+...
(((x-x_apex)+Horizontal_distance_between_source_and_receiver/2).^2)))); 

Hyperbola = plot(x,t_hyperbola,'b-');
title('Hyperbola', 'Fontsize',12) 

choice = menu('Adjust velocity to find the best fit of the hyperbola to the hyperbolic shape of the diffraction pattern:',...
    'increase velocity','decrease velocity','end changing velocity');
dv = 0.0033; %velocity step in m/ns 

v = v_initial;
while choice ~=3 %so as long as the user wants to change the velocity
   if choice == 1 %if the user wants to increase the velocity
       v = v + dv; %add 0.33 cm/ns to the former velocity
       delete(Hyperbola); %remove the former hyperbola plot
   else                   %if the user wants to decrease the velocity 
       v = v - dv; %substract 0.33cm/ns from the former velocity
       delete(Hyperbola);
   end
   t_hyperbola = ((1/v)*(sqrt((((t_apex*v)/2).^2)+(((x-x_apex)-...
   Horizontal_distance_between_source_and_receiver/2).^2))+sqrt((((t_apex*v)/2).^2)+...
   (((x-x_apex)+Horizontal_distance_between_source_and_receiver/2).^2)))); 
   Hyperbola = plot(x,t_hyperbola,'c-'); %plots the new hyperbola for the changed velocity
   
   %again gives the user the choice to either adjust the velocity of move on to the next step
   choice = menu('Adjust velocity to find the best fit of the hyperbola to the hyperbolic shape of the diffraction pattern:',...
   'increase velocity','decrease velocity','end changing velocity'); 
end

% Now that the correct velocity has been found, migration is applied:

z = 0:0.01:0.5; %depth vector since I investigate all depths between 
%the surface, z = 0 and 50 cm below the surface with image depth intervals of 1 cm

%call in the function
[Migration_output] = migration(z,x,v,Horizontal_distance_between_source_and_receiver,Time_sampling_interval, Adjusted_data2);

figure
imagesc(x,z,Migration_output);
colormap(flipud(gray)); colorbar;
xlabel('Distance (in meter)', 'Fontsize', 12) 
ylabel('Depth (in meter)', 'Fontsize',12)
title('Migrated Dataset with correct velocity', 'Fontsize',12) 

v_high = 1.1*v; %velocity 10% too high
[Migration_output] = migration(z,x,v_high,Horizontal_distance_between_source_and_receiver,Time_sampling_interval, Adjusted_data2);

figure
imagesc(x,z,Migration_output);
colormap(flipud(gray)); colorbar;
xlabel('Distance (in meter)', 'Fontsize', 12) 
ylabel('Depth (in meter)', 'Fontsize',12)
title('Migrated dataset with velocity 10% too high', 'Fontsize',12) 

%If you compare this figure with the correct velocity migration plot,
%you can see that, if you use a velocity which is 10% too high,
%the hyperbolas become a bit blurred and tend to go upwards.
%Using a different velocity changes both the slope of the reflector and the position of the diffractors.
%An increasing velocity leads to a steeper slope and to lower positions of the diffractors than where they should be.

v_low = 0.9*v; %velocity 10% too low
[Migration_output] = migration(z,x,v_low,Horizontal_distance_between_source_and_receiver,Time_sampling_interval, Adjusted_data2);

figure
imagesc(x,z,Migration_output);
colormap(flipud(gray)); colorbar;
xlabel('Distance (in meter)', 'Fontsize', 12) 
ylabel('Depth (in meter)', 'Fontsize',12)
title('Migrated dataset with velocity 10% too low', 'Fontsize',12) 

%In the last figure it's clearly visible that if you apply migration with a
%velocity which is 10% too low the hyperbolas are also blurred and
%tend to become broader downwards.
%Also there is a less steep slope and the diffractors are shifted upwards.

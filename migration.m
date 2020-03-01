% Function file to support: Wave field processing in the time domain
% By: Tiffany Vlaar 

% Migration is one of the most important processing steps, since its goal is 
% to make the stacked section resemble the true geologic cross section along the seismic line. 
% Migration is necessary since travel times are always plotted vertically down, 
% while rays might come from the right of the section. 
% Migration does therefore not affect horizontal events, 
% but it moves dipping reflectors to their true subsurface position.
% Migration also collapses diffraction patterns, and thereby delineates fault planes and other detailed subsurface features. 
% This increases spatial resolution and makes structural interpretation easier and more accurate

function [Migration_output] = migration(z,x,v,Horizontal_distance_between_source_and_receiver,Time_sampling_interval,Adjusted_data2)
    
    Migration_output = zeros(length(z),length(x)); %makes the output of the function an array filled with zeros
    
    for i = 1:length(z) %loop over all different depths
        
        for j = 1:length(x) %summing for each horizontal distance
        
            %equation (3) with d = z, x = x-x(j) and x0 = Horizontal_distance_between_source_and_receiver
            t_hyperbola = ((1/v)*(sqrt((z(i).^2)+(((x-x(j))-...
            Horizontal_distance_between_source_and_receiver/2).^2))+sqrt((z(i).^2)+...
            (((x-x(j))+Horizontal_distance_between_source_and_receiver/2).^2)))); 
            
            time = (t_hyperbola./Time_sampling_interval)*(10^-9);
            Sum1 = zeros(1,length(t_hyperbola)); %creates an array filled with zeros to store the values for the sum in
            Sum2 = zeros(1,length(t_hyperbola));
            time1 = floor(time) + 1; %lower neighbouring time point
            time2 = ceil(time) + 1; %higher neighbouring time point
            
            for k=1:length(x)          %for each different horizontal distance
                if time2(k) <= size(Adjusted_data2,1)  
                %to make sure that I don't call elements from my data matrix that do not exist
                    
                    Sum1(k) = Adjusted_data2(time1(k),k); 
                    %Fills the empty array's Sum1 and Sum2 with the data
                    %for resp. the lower neighbouring time point and for time2.
                    Sum2(k) = Adjusted_data2(time2(k),k);            
                end
            end
            
        Sum=Sum1+(Sum2-Sum1).*((time-(time1-1))/(time2-time1)); %weighted data values
        
        Migration_output(i,j)=sum(Sum); %sum of all weighted data values that fit the hyperbola 
        
        end
    end
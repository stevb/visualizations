%steve brust 4049195
%simple PIV script to show the theory behind the experimental technique.
%The plotted images during processing show the correlation map in the upper
%left, the two upper right images are frames A and B (with resulting vector
%overlayed) of the image pair separated by the time step, dt, and lastly, 
%the bottom image shows the interrogation window travelling across the
%image. After processing a full vector field over an original image is
%shown. Airfoil used was a NACA 0012.
clear all, close all, clc

%1 will result in processing display (slow)
%0 will go straight to results
code_processing_visualization = 0;

%load image
%good quality recording from experiment
image_data = imread('B00003.tif');

dt = 80e-6; %80 [us]
mag_fac = 0.059168; %magnification factor [-]

%settings
AoA = 15;

%find number of rows for splitting
[n_rows, ~] = size(image_data);

%separate the two images
image_a = image_data(1:n_rows/2,:);
image_b = image_data(n_rows/2+1:end,:);

%find the number of rows and columns
[n_rows, n_cols] = size(image_a);

%corresponds to window size of 32
window_size = 32;

row_index = [0:n_rows/window_size];
col_index = [0:n_cols/window_size];

%storing cen_x cen_y dx dy
vector_vals = zeros(length(col_index)-1,4,length(row_index)-1);

for j = 1:length(row_index)-1
    for i = 1:length(col_index)-1
    
        window_a = image_a(row_index(j)*window_size+1:row_index(j+1)*window_size,...
            col_index(i)*window_size+1:col_index(i+1)*window_size);

        window_b = image_b(row_index(j)*window_size+1:row_index(j+1)*window_size,...
            col_index(i)*window_size+1:col_index(i+1)*window_size); 
 
        window_a = double(window_a);
        window_b = double(window_b);

        window_a = window_a - mean(window_a(:));
        window_b = window_b - mean(window_b(:));

        cross_cor_result = xcorr2(window_b, window_a);
        
        [max_x,max_y] = find(cross_cor_result == max(cross_cor_result(:)));
        dx = max_x - window_size;
        dy = max_y - window_size;

        vector_vals(i,1,j) = window_size + (i-1)*2*window_size;
        vector_vals(i,2,j) = (window_size + (j-1)*2*window_size);

        vector_vals(i,3,j) = dx;
        vector_vals(i,4,j) = dy;
        
        if code_processing_visualization == 1
            %figure display
            %row_index(j), col_index(i)
            subplot(3,3,1); surf(cross_cor_result)
            subplot(3,3,2); imagesc(window_a)
            subplot(3,3,3); imagesc(window_b)
            hold on
            quiver(16,16,dy,dx,'color',[1 .5 0],'LineWidth',2)
            hold off
            subplot(3,3,4:9); image(image_a)
            hold on
            plot([col_index(i)*window_size col_index(i+1)*window_size+1],...
                [row_index(j)*window_size row_index(j)*window_size],...
                'r','LineWidth',1) %top horiz. line
            plot([col_index(i)*window_size col_index(i+1)*window_size+1],...
                [row_index(j+1)*window_size row_index(j+1)*window_size],...
                'r','LineWidth',1) %bottom horiz. line
            plot([col_index(i)*window_size col_index(i)*window_size+1],...
                [row_index(j)*window_size row_index(j+1)*window_size],...
                'r','LineWidth',1) %left vert. line
            plot([col_index(i+1)*window_size col_index(i+1)*window_size+1],...
                [row_index(j)*window_size row_index(j+1)*window_size],...
                'r','LineWidth',1) %right vert. line
            hold off
            pause
            %pause(0.01)
        end
    end
end

for i = 1:size(vector_vals,3)
    for j = 1:size(vector_vals,1)
        magnitude = sqrt((vector_vals(j,3,i)*(6.45e-6)/((dt)*mag_fac))^2 ...
            + (vector_vals(j,4,i)*(6.45e-6)/((dt)*mag_fac))^2);
        if magnitude > 15
            vector_vals(j,3,i) = 0;%vector_vals(j,3,i)*0.01
            vector_vals(j,4,i) = 0;%vector_vals(j,4,i)*0.01
        end
    end
end

%erase the leading values that are messed up
vector_vals(1:6,3:4,:) = 0;

fig = figure;

imagesc(imresize(image_a, 2))
hold on

%dy and dx have been flipped
quiver(vector_vals(:,1,:),vector_vals(:,2,:),...
    vector_vals(:,4,:)*(6.45e-6)/((dt)*mag_fac),...
    vector_vals(:,3,:)*(6.45e-6)/((dt)*mag_fac),'color',[1 .5 0])

[airfoil_x,airfoil_y] = read_airfoil('seligdatfile.dat');

%For 0 AoA
if AoA == 0
    airfoil_x = airfoil_x * 2400 + 10;
    airfoil_y = airfoil_y * 2400 + 1160;
    fill(airfoil_x,airfoil_y,[0.5 0.51 0.5])
end
%For 5 AoA
if AoA == 10
    rotAngle = deg2rad(-5);
    xRot     = airfoil_x*cos(rotAngle) - airfoil_y*sin(rotAngle);
    yRot     = airfoil_x*sin(rotAngle) + airfoil_y*cos(rotAngle);
    xRot = xRot*2400 - 0;
    yRot = yRot*2400 + 920;
    fill(xRot,yRot,[0.5 0.51 0.5]);
end

%For 15 AoA
if AoA == 15
    rotAngle = deg2rad(-12);
    xRot     = airfoil_x*cos(rotAngle) - airfoil_y*sin(rotAngle);
    yRot     = airfoil_x*sin(rotAngle) + airfoil_y*cos(rotAngle);
    xRot = xRot*2400 + 46;
    yRot = yRot*2400 + 1316;
    fill(xRot,yRot,[0.5 0.51 0.5]);
end

ylabel('y location')
xlabel('x location')
title('PIV image with vector and airfoil overlay for \alpha = 15^\circ')

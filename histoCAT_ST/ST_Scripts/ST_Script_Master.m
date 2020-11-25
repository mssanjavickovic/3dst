% Denis Schapiro
% Independent Fellow Harvard Medical School and Broad Institute
% 2019
clc
clear all
%% Load everything needed
% Load Mask
Mask = imread('/Users/denis/Desktop/ST/RA2_9.jpg');
% Load Image
Image = imread('/Users/denis/Desktop/ST/RA2_HE_9.jpg');
% Load spot location
[~, ~,raw] = tsvread('/Users/denis/Desktop/ST/RA2_9_xy.tsv');
Spot_XY = cellfun(@str2double,raw(2:end,2:3));

%% Visualize for QC (Comment if not needed in batch processing)
% Extract X and Y locations with a small margin
Spot_size_radius_with_margin = 200;

Crop_coordinates_X_small = Spot_XY(44,1)-Spot_size_radius_with_margin;
Crop_coordinates_X_large = Spot_XY(44,1)+Spot_size_radius_with_margin;

Crop_coordinates_Y_small = Spot_XY(44,2)-Spot_size_radius_with_margin;
Crop_coordinates_Y_large = Spot_XY(44,2)+Spot_size_radius_with_margin;

Crop_mask = Mask(Crop_coordinates_Y_small:Crop_coordinates_Y_large,...
    Crop_coordinates_X_small:Crop_coordinates_X_large);

Crop_HE = Image(Crop_coordinates_Y_small:Crop_coordinates_Y_large,...
    Crop_coordinates_X_small:Crop_coordinates_X_large,:);
figure()
imshow(Crop_HE);
figure()
imshow(Crop_mask);

%% Create new tiff files using for-loop
% Decide on margin
Spot_size_radius_with_margin = 200;

for i=1:size(Spot_XY,1)
    
    % Extract X and Y locations with a small margin
    Crop_coordinates_X_small = Spot_XY(i,1)-Spot_size_radius_with_margin;
    Crop_coordinates_X_large = Spot_XY(i,1)+Spot_size_radius_with_margin;
    
    Crop_coordinates_Y_small = Spot_XY(i,2)-Spot_size_radius_with_margin;
    Crop_coordinates_Y_large = Spot_XY(i,2)+Spot_size_radius_with_margin;
    
    % Create cropped image
    Crop_HE_loop = Image(Crop_coordinates_Y_small:Crop_coordinates_Y_large,...
        Crop_coordinates_X_small:Crop_coordinates_X_large,:);
    
    % Create cropped mask for QC
    Crop_mask_loop = Mask(Crop_coordinates_Y_small:Crop_coordinates_Y_large,...
        Crop_coordinates_X_small:Crop_coordinates_X_large);
    
    imwrite(Crop_HE_loop,strcat('RA2_9_Spot_',raw{i+1,1},'.jpg'));
    imwrite(Crop_mask_loop,strcat('RA2_9_mask_Spot_',raw{i+1,1},'.jpg'));
    
end

clc
clear all
RA_which = 'RA1';
Sample_which = '1';

% Load spot location
[~, ~,raw] = tsvread(strcat('/Users/denis/Desktop/ST/',RA_which,'_',Sample_which,'_xy.tsv'));
Spot_XY = cellfun(@str2double,raw(2:end,2:3));

for k=2:size(raw,1)
    %% Create histoCAT folder structure
    % If script starts here:
    
    Image_name = strcat(...
        '/Users/denis/Desktop/ST/RA1_1_Spots_for_histoCAT/',RA_which,'_',Sample_which,'_Spot_',raw{k,1},'.jpg');
    New_folder_path = strcat(...
        '/Users/denis/Desktop/ST/', RA_which,'_',Sample_which,'_Spots_for_histoCAT/',RA_which,'_',Sample_which,'_Spot_',raw{k,1});
    
    JPG = imread(Image_name);
    
    % Extract RGB
    R = JPG(:,:,1);
    G = JPG(:,:,2);
    B = JPG(:,:,3);
    
    % Create folder
    mkdir(New_folder_path);
    imwrite(R,fullfile(New_folder_path,strcat('Spot_Red.tif')));
    imwrite(G,fullfile(New_folder_path,strcat('Spot_Green.tif')));
    imwrite(B,fullfile(New_folder_path,strcat('Spot_Blue.tif')));
    
    % Move mask to folder
    Mask_location = fullfile(strcat(...
        '/Users/denis/Desktop/ST/RA1_1_Spots_for_histoCAT/',RA_which,'_',Sample_which,'_Spot_',raw{k,1},'_Probabilities_mask.tiff'));
    Mask_destination = fullfile(New_folder_path,strcat(RA_which,'_',Sample_which,'_Spot_',raw{k,1},'_Probabilities_mask.tiff'));
    copyfile(Mask_location,Mask_destination);
end
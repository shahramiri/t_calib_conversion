classdef Undistortion < handle % Second Version
                                                            
    properties
        LUTList
        GrayThreshold
        SizeThreshold
        SizeOfNeighborhood
        UndistortMask
        UndistortRadius
        NumInterpPoints         % Number of points used for interpolation
        IDWPowerParameter       % Power parameter for inverse distance weighting function.
                                % Increase for greater influence from closer points
        CropImage               % Option to crop image to show only undistorted region.
        ParentFolderName        % Parent Data Folder Name
        InputPoints             % Input Points Digitized        
        BasePoints              % Base Points 
    end
    
    methods
        function Obj = Undistortion(CalibFolderName)
            Obj.LUTList = List('LUT List');
            Obj.GrayThreshold = 100;
            Obj.SizeThreshold = 300;
            Obj.SizeOfNeighborhood = 20;
            Obj.UndistortRadius = 500;
            Obj.NumInterpPoints = 5;
            Obj.IDWPowerParameter = 2;
            Obj.CropImage = 1;
            Obj.ParentFolderName = CalibFolderName; %fullfile('Calibration\',CarmModelName);            
            if ~exist(Obj.ParentFolderName,'dir')
                mkdir(Obj.ParentFolderName);
            end
        end
        
        function [pixLUT] = createLUT(Obj,imageFile)
            % Code from Brown University XROMM Undistorter
            [I,Iinfo] = calimread(imageFile,0);
            rawImax=max(max(I(:,:))); % maximum raw image pixel value
            %do some initial processing (Adjusted strel value to accomodate
            %image size)
            Ibackground = imopen(I,strel('disk',round(Iinfo.Height*(15/1024))));
            I = imsubtract(I,Ibackground);
            I = imadjust(I);
            
            levelMax = (2^8)-1;
            if(isa(I,'uint16'))
                levelMax = (2^16)-1;
            end
            level = Obj.GrayThreshold/levelMax;
            
            % We don't want to include insufficiently sized cells that may alter
            % our calculations, so color the pixels from such cells black.
            sizeLim = Obj.SizeThreshold;
            badPixels = removesmalls(I,level,sizeLim);
            I(badPixels) = 0;
            
            [InputPoints, BasePoints, dY, dX, dOffY, gridOrient] = ...
                getControlPoints(I, level, sizeLim);
            iResolution = [size(I,2) size(I,1)];
            
            undistortTform = cp2tform(InputPoints, BasePoints, 'lwm', Obj.SizeOfNeighborhood);
            pixLUT = xCreateLookupTable(undistortTform,iResolution);
            
            Obj.InputPoints=InputPoints;
            Obj.BasePoints=BasePoints;
            
            % Undistortion preview
%             preview = xFastDistortionCorrection(pixLUT, I);
%             figure('NumberTitle', 'off', ...
%                 'Name', 'Undistorted Calibration Grid Image');
%             imshow(preview);
        end
        
        function undistortGrid(Obj,imageFile,matFile)
            pixLUT = Obj.createLUT(imageFile);
            
            load(matFile);
            tilt = ImgPlane.CarmPose.Tilt;
            rotation = ImgPlane.CarmPose.Rotation;
            lut = LUT(pixLUT,tilt,rotation);
            lut.save(Obj.ParentFolderName);
            Obj.LUTList.addItem(lut);
        end
        
        function undistortGridWoTracking(Obj,imageFile,tilt,rotation)
            
            pixLUT = Obj.createLUT(imageFile);
            
%             load(matFile);
%             tilt = ImgPlane.CarmPose.Tilt;
%             rotation = ImgPlane.CarmPose.Rotation;
            lut = LUT(pixLUT,tilt,rotation);
            lut.save(Obj.ParentFolderName);
            Obj.LUTList.addItem(lut);
        end
        
        function batchUndistortGrids(Obj,imagesFolder)
            if ~strcmp(imagesFolder(end),'\')
                imagesFolder = [imagesFolder '\'];  % Add \ to folder name if not included
            end
            d1=dir([imagesFolder '*.tif']);
            imageFiles = sort_nat({d1.name});   % Sort image files in numerical order
            d2=dir([imagesFolder '*.mat']);
            matFiles = sort_nat({d2.name});
            for i = 1:size(d1,1)
                disp(['Generating Lookup Table for ' imageFiles{i}(1:end-4)]);
                imageFile = [imagesFolder imageFiles{i}];
                matFile = [imagesFolder matFiles{i}];
                Obj.undistortGrid(imageFile,matFile);
            end
        end
        
        % Save the LUTList
        function save(Obj)
            % Remove PixLUT and PixChangeTable before saving.
            % Only Tilt, Rotation, PixChangeMatFile need to be saved.
            for i = 1:Obj.LUTList.Count
                Obj.LUTList.Items(i).PixLUT = [];
                Obj.LUTList.Items(i).PixChangeTable = [];
            end
            LUTList = Obj.LUTList;
            save(fullfile(Obj.ParentFolderName,'UndistortionLUTs.mat'),'LUTList');
        end
        
        % Load the LUTList
        function load(Obj)
            load(fullfile(Obj.ParentFolderName,'UndistortionLUTs.mat'));
            Obj.LUTList = LUTList;
        end
        
        % Find 3 LUTs with tilt/rotation closest to current image tilt/rot
        function varargout = findNClosestOrientations(Obj,numberOfPoints,currTilt,currRot)
            tilts = zeros(Obj.LUTList.Count,1);
            rots = zeros(Obj.LUTList.Count,1);
            for i = 1:Obj.LUTList.Count
                tilts(i) = Obj.LUTList.Items(i).Tilt;
                rots(i) = Obj.LUTList.Items(i).Rotation;
            end
            tiltDiffs = tilts - repmat(currTilt,size(tilts,1),1);   % Difference between current tilt and LUT tilts
            rotDiffs = rots - repmat(currRot,size(rots,1),1);       % Difference between current rot and LUT rots
            orientDiffs = (tiltDiffs.^2 + rotDiffs.^2).^(1/2);      % Overall difference in orientation
            [min,minIndex] = sort(orientDiffs);     % minIndex returns indices of orientation differences sorted from least to greatest
            varargout = cell(numberOfPoints,1);
            for n = 1:numberOfPoints
                varargout{n} = Obj.LUTList.Items(minIndex(n));  % Return N closest LUTs
            end
        end
        
        % Create new LUT interpolated from 3 LUTs with closest orientations
        function pixLUT = createInterpLUT(Obj,imageTilt,imageRot)
            numPoints = Obj.NumInterpPoints;
            [luts{1:numPoints}] = Obj.findNClosestOrientations(numPoints,imageTilt,imageRot);
            [defaultGrid(:,:,2),defaultGrid(:,:,1)]=meshgrid(1:1024,1:1024);
            tilts = zeros(size(luts,2),1);
            rots = zeros(size(luts,2),1);
            for i = 1:numPoints
                tilts(i) = luts{i}.Tilt;        % Create vector of tilts
                rots(i) = luts{i}.Rotation;     % Create vector of rotations                
                luts{i}.load(Obj.ParentFolderName);% Load pixel change table for each LUT
            end
            
            % Create undistort mask if not already created
            if isempty(Obj.UndistortMask)
                Obj.UndistortMask = Obj.createUndistortMask(size(luts{1}.PixChangeTable,1),Obj.UndistortRadius);
            end
            
            % Inverse distance weighting interpolation method
            tiltDist = tilts - repmat(imageTilt,size(tilts,1),1);   % Difference between current tilt and LUT tilts
            rotDist = rots - repmat(imageRot,size(rots,1),1);       % Difference between current rot and LUT rots
            pointDist = (tiltDist.^2 + rotDist.^2).^(1/2);          % Overall difference in orientation
            if ~any(pointDist==0)   
                w = 1./(pointDist).^(Obj.IDWPowerParameter);        % Weight factors
                A = 0;
                for i = 1:numPoints
                    A = A + w(i)*luts{i}.PixChangeTable;    % Sum of weight factors times corresponding values
                end
                B = sum(w);                                 % Sum of weight factors
                interpPixChange = A/B;                      % IDW interpolating function
            else
                % If orientation matches calibration point, use calibration point's pixel change table.
                ind = (pointDist==0);
                interpPixChange = luts{ind}.PixChangeTable;
            end
            mask = cat(3,Obj.UndistortMask,Obj.UndistortMask);
            interpPixChange = interpPixChange.*mask;    % Apply undistort radius mask;
            pixLUT = defaultGrid + interpPixChange;     % Create pixel position lookup table
        end
        
        % Create mask for radius to apply undistortion to
        function mask = createUndistortMask(Obj,imageSize,radius)
            center = imageSize/2;
            [W,H] = meshgrid(1:imageSize,1:imageSize);
            mask = sqrt((W-center).^2 + (H-center).^2) < radius;
        end
        
        function setUndistortRadius(Obj,radius)
            Obj.UndistortMask = [];
            Obj.UndistortRadius = radius;
        end
        
        % Create undistorted image using interpolated LUT
        function im = undistortImage(Obj,imageFile,varargin)
            if nargin == 3                  % mat file was provided
                matFile = varargin{1};
                load(matFile);
                imageTilt = ImgPlane.CarmPose.Tilt;
                imageRot = ImgPlane.CarmPose.Rotation;
            elseif nargin == 4              % tilt and rotation were manually provided
                imageTilt = varargin{1};
                imageRot = varargin{2};
            end
            
            if strcmp(class(imageFile),'char') % if the imageFile is a String                
                I = imread(imageFile);
            else
                I = imageFile; % the provided variable is an image
            end
            
            pixLUT = Obj.createInterpLUT(imageTilt,imageRot);
            
            if Obj.CropImage == 1
                % Crop image to show only undistorted region
                mask = cat(3,Obj.UndistortMask,Obj.UndistortMask,Obj.UndistortMask);                
                if size(I,3)==1
                    I = I.*uint8(mask(:,:,1));
                else
                    I = I.*uint8(mask);
                end
            end
            
            im = xFastDistortionCorrection(pixLUT,I);
        end
        
        
                % Create undistorted image using interpolated LUT
        function im = undistortImageSimple(Obj,I,imrot)

            Rot=[-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90];
            dRot=abs(Rot-imrot);
            [i,j]=min(dRot);
            rotation=Rot(j);
            
            RotStr=num2str(rotation);
            switch numel(RotStr)
                case 1
                    RotStr=['  ', RotStr];
                case 2
                    RotStr=[' ', RotStr];
            end
            MatFName=sprintf('Dist\\S%s_DistCalib.mat',RotStr);
            ImgFName=sprintf('Dist\\S%s_DistCalib.tif',RotStr);            
            LutFName=sprintf('Pixel Change Tables\\PixChangeTable_0tilt_%irot.mat',rotation);
            disp(sprintf('Distortion File Used: %s',LutFName'));
            FFName=fullfile(['Calibration\' Obj.ParentFolderName],LutFName);
            
            % load the LUTfile and correct distortion if the file is found
            if exist(FFName, 'file')
                try
                    load(FFName);
                    PixChangeTable = double(pixChangeTable)./100;
                    [defaultGrid(:,:,2),defaultGrid(:,:,1)]=meshgrid(1:1024,1:1024);
                    % Create undistort mask if not already created
                    if isempty(Obj.UndistortMask)
                        Obj.UndistortMask = Obj.createUndistortMask(size(PixChangeTable,1),Obj.UndistortRadius);
                    end
                    mask = cat(3,Obj.UndistortMask,Obj.UndistortMask);
                    pixLUT = defaultGrid + PixChangeTable; 

                    if Obj.CropImage == 1
                        % Crop image to show only undistorted region
                        mask = cat(3,Obj.UndistortMask,Obj.UndistortMask,Obj.UndistortMask);                
                        if size(I,3)==1
                            I = im2uint8(I).*uint8(mask(:,:,1));
                        else
                            I = im2uint8(I).*uint8(mask);
                        end
                    end
                    im = xFastDistortionCorrection(pixLUT,I);
                catch
                    im=[]; %return an empty image in case there is an error in processing.
                end
            else
                im=[]; %return an empty image in case there is file not found
            end
            
        end
        
        
        function plotClosestPoints(Obj,varargin)
            if nargin == 2                  % mat file was provided
                matFile = varargin{1};
                load(matFile);
                imageTilt = ImgPlane.CarmPose.Tilt;
                imageRot = ImgPlane.CarmPose.Rotation;
            elseif nargin == 3              % tilt and rotation were manually provided
                imageTilt = varargin{1};
                imageRot = varargin{2};
            end
            
            numPoints = Obj.NumInterpPoints;
            [luts{1:numPoints}] = Obj.findNClosestOrientations(numPoints,imageTilt,imageRot);
            tilts = zeros(Obj.LUTList.Count,1);
            rots = zeros(Obj.LUTList.Count,1);
            for i = 1:Obj.LUTList.Count
                tilts(i) = Obj.LUTList.Items(i).Tilt;
                rots(i) = Obj.LUTList.Items(i).Rotation;
            end
            
            figure;
            plot(imageTilt,imageRot,'o','MarkerFaceColor','y');
            hold on
            for i = 1:numPoints
                plot(luts{i}.Tilt,luts{i}.Rotation,'s','MarkerFaceColor','r');
            end
            plot(tilts,rots,'+');
            xlabel('Tilt');
            ylabel('Rotation');
        end
        
        
        function pixDistortion = showImageDistortion(Obj,imageFile)
            % Code from Brown University XROMM Undistorter
            [I,Iinfo] = calimread(imageFile,0);
            rawImax=max(max(I(:,:))); % maximum raw image pixel value
            %do some initial processing (Adjusted strel value to accomodate
            %image size)
            Ibackground = imopen(I,strel('disk',round(Iinfo.Height*(15/1024))));
            I = imsubtract(I,Ibackground);
            I = imadjust(I);
            
            levelMax = (2^8)-1;
            if(isa(I,'uint16'))
                levelMax = (2^16)-1;
            end
            level = Obj.GrayThreshold/levelMax;
            
            % We don't want to include insufficiently sized cells that may alter
            % our calculations, so color the pixels from such cells black.
            sizeLim = Obj.SizeThreshold;
            badPixels = removesmalls(I,level,sizeLim);
            I(badPixels) = 0;
            
            [InputPoints, BasePoints, dY, dX, dOffY, gridOrient] = ...
                getControlPoints(I, level, sizeLim);
            
            udim = size(I,1);
            vdim = size(I,2);
            distFromCenter = rnorm([InputPoints(:,1)-udim/2 InputPoints(:,2)-vdim/2]);
            [distSorted,order] = sort(distFromCenter);
            
            pixDiff = BasePoints - InputPoints;
            pixDistortion = (pixDiff(:,1).^2 + pixDiff(:,2).^2).^(1/2);
            pixDistortion = pixDistortion(order);
            
%             figure;
%             plot(distSorted,pixDistortion);
%             xlabel('Distance from center (pixels)');
%             ylabel('Distortion (pixels)');
        end
        
        function showDistortionMap(Obj,origImageFile,undImageFile)
            [defaultGrid(:,:,2),defaultGrid(:,:,1)]=meshgrid(1:1024,1:1024);
            pixLUT = Obj.createLUT(origImageFile);
            origPixChange = pixLUT - defaultGrid;
            origPixChange(abs(origPixChange)>200) = 0;
            
            % Create undistort mask if not already created
            if isempty(Obj.UndistortMask)
                Obj.UndistortMask = Obj.createUndistortMask(size(origPixChange,1),Obj.UndistortRadius);
            end
            mask = cat(3,Obj.UndistortMask,Obj.UndistortMask);
            
            origPixChange = origPixChange.*mask;    % Apply undistort radius mask;
            origPixDistortion = (origPixChange(:,:,1).^2 + origPixChange(:,:,2).^2).^(1/2);
            figure;
            origMap = imagesc(origPixDistortion);        % draw image and scale colormap to values range
            colorbar;          % show color scale
            title('Distortion (in pixels) for Original Image');
            
            undPixLUT = Obj.createLUT(undImageFile);
            undPixChange = undPixLUT - defaultGrid;
            undPixChange(abs(undPixChange)>200) = 0;
            undPixChange = undPixChange.*mask;    % Apply undistort radius mask;
            undPixDistortion = (undPixChange(:,:,1).^2 + undPixChange(:,:,2).^2).^(1/2);
            figure;
            clim = [min(min(get(origMap,'CData'))), max(max(get(origMap,'CData')))];
            undMap = imagesc(undPixDistortion,clim);        % draw image and scale colormap to values range
            colorbar;          % show color scale
            title('Distortion (in pixels) for Undistorted Image');
        end
        
    end
    
end



% Undistorter subfunctions
% Code from Brown University XROMM Undistorter

function [pixelList] = removesmalls(I,level,sizeLim)

%level = get(hsGray,'Value')/levelMax;
udim = size(I,2); %horizontal
vdim = size(I,1); %vertical

% create black/white image and extract information
bw = im2bw(I,level);
bwlabeled = bwlabel(bw,4);
celldata = regionprops(bwlabeled,'Area','BoundingBox','PixelIdxList');

%sizeLim = get(hsSize,'value');% get the size limit from hsSize

%get areas and bounding boxes
areas = cat(1,celldata.Area);
bb= cat(1,celldata.BoundingBox);

% select pixels with edge touching cells -- i.e., edge of rectangle
% that bounds entire image.
idx = find(bb(:,1) == 0.5 | bb(:,1)+bb(:,3) > udim-0.5 |...
    bb(:,2) == 0.5 | bb(:,2)+bb(:,4) > vdim-0.5 );
edgeTouchingCells = celldata(idx);
pixlist = cat(1,edgeTouchingCells.PixelIdxList);

% select small cell pixels
idx = find(areas < sizeLim);
smallcells = celldata(idx);

%stack the two lists
pixelList = cat(1,pixlist, cat(1,smallcells.PixelIdxList));
%       set(handles.hsSize,'UserData',pixelList); % <--caller should do this

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Non-Staggered Grid Procedure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BasePoints,InputPoints,dY,dX,dOffY] = getPts4(cp,centroids,four)

hwait = waitbar(0, 'Creating/checking Control Points');

% set first point to the point closest to the center of the image
BasePoints = cp;
InputPoints = cp;

% set the cell center to center distance
dY = mean(rnorm([four(:,1)-cp(1) four(:,2)-cp(2)]));
dX = NaN;
dOffY = NaN;

% establish base points to undisort the quadrants.
r = four(2,:); %get the initial right cell for 1st quadrant
f.q1 = four(1,:);
f.q2 = four(3,:);
f.q3 = four(4,:);

cp = four(1,:);
[four] = adj4Cells(cp,centroids);
f.q4 = four(4,:);

cp = four(3,:);


i = 2;
ychk = 1; % used to check for edge on y image axis
xchk = 1; % used to check for edge on y image axis
col = 1; a = 1; b = 2; c = -1; d = 1;

for j = 1:4
    
    while xchk == 1
        
        while ychk == 1
            
            [four] = adj4Cells(cp,centroids);
            
            angles = (180/pi)*(atan2(four(:,2)-cp(2),four(:,1)-cp(1)));
            ychk = ycheck(angles);
            if ychk == 0
                break
            end
            
            cp = four(a,:);
            InputPoints(i,:) = cp;
            BasePoints(i,1) = BasePoints(i-1,1);
            BasePoints(i,2) = BasePoints(i-1,2)+c*dY;
            i = i+1;
        end
        
        cp = r;
        [four] = adj4Cells(cp,centroids);
        
        angles = (180/pi)*(atan2(four(:,2)-cp(2),four(:,1)-cp(1)));
        xchk = ycheck(angles);
        if xchk == 0
            break
        end
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1)+(d*dY)*col;
        if j ==1 || j == 3
            BasePoints(i,2) = BasePoints(1,2);
        elseif j == 2
            BasePoints(i,2) = BasePoints(1,2)+dY;
        elseif j == 4
            BasePoints(i,2) = BasePoints(1,2)-dY;
        end
        
        r = four(b,:);
        
        i = i + 1;
        col = col + 1;
        ychk = 1;
        
    end
    
    if j == 1
        
        cp = f.q2;
        [four] = adj4Cells(cp,centroids);
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1);
        BasePoints(i,2) = BasePoints(1,2)+dY;
        i = i+1;
        
        ychk = 1;
        xchk = 1;
        col = 1; a = 3; b = 2; c = 1; d = 1;
        
        r = four(b,:);
        
    elseif j == 2
        
        cp = f.q3;
        [four] = adj4Cells(cp,centroids);
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1)-dY;
        BasePoints(i,2) = BasePoints(1,2);
        i = i+1;
        
        ychk = 1;
        xchk = 1;
        col = 2; a = 3; b = 4; c = 1; d = -1;
        
        r = four(b,:);
        
    elseif j == 3
        
        cp = f.q4;
        [four] = adj4Cells(cp,centroids);
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1)-dY;
        BasePoints(i,2) = BasePoints(1,2)-dY;
        i = i+1;
        
        ychk = 1;
        xchk = 1;
        col = 2; a = 1; b = 4; c = -1; d = -1;
        
        r = four(b,:);
        
    end
    waitbar(j*0.25,hwait);
end
close(hwait);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Staggered Grid Procedure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% called by getControlPoints.
function [BasePoints,InputPoints,dY,dX,dOffY] = getPts6(cp,centroids,six)

hwait = waitbar(0, 'Creating/checking Control Points');

% get matrix of points for staggered cells
BasePoints = cp;
InputPoints = cp;

% set the cell center to center distance.  Use average of distances
% between centerpoint (cp) and six surrounding points.
dY = mean(rnorm([six(:,1)-cp(1) six(:,2)-cp(2)]));
dX = cosd(30)*dY;
dOffY = sind(30)*dY;

% establish base points to undisort the quadrants.
r = six(2,:); %get the initial right cell for 1st quadrant
f.q1 = six(1,:);
f.q2 = six(4,:);
f.q3 = six(6,:);
f.q4 = six(5,:);

i = 2;
ychk = 1; % used to check for edge on y image axis
xchk = 1; % used to check for edge on y image axis
col = 1; a = 1; b = 2; c = -1; d = 1;

for j = 1:4
    while xchk == 1 % travel down column
        while ychk == 1 %(repeat) until edge (i.e., bottom)reached.
            
            [six,sm] = adjCells(cp,centroids);
            
            angles = (180/pi)*(atan2(six(:,2)-cp(2),six(:,1)-cp(1)));
            ychk = ycheck(angles);
            if ychk == 0
                break  % break out of while ychk==1 loop
            end
            
            cp = six(a,:); % change centerpoint to point closest to
            % bottom of hexagon (a=1) or top of hexagon (a=4)
            InputPoints(i,:) = cp;
            BasePoints(i,1) = BasePoints(i-1,1);
            BasePoints(i,2) = BasePoints(i-1,2)+c*dY;
            i = i+1;
        end   % end while ychk
        
        cp = r; %change centerpoint--move right (b=2) or left (b=5)
        
        [six,sm] = adjCells(cp,centroids);
        
        angles = (180/pi)*(atan2(six(:,2)-cp(2),six(:,1)-cp(1)));
        xchk = ycheck(angles);
        if xchk == 0 %break out of while xchk==1 loop.  This stops
            % us when edge of grid is reached as we traverse in
            % (roughly) the x-direction.
            break
        end
        
        
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1)+(d*dX)*col;
        if j == 3 || j == 4
            BasePoints(i,1) = BasePoints(1,1)+(d*dX)*(col+1);
        end
        %Next column:
        col = col + 1;
        
        if j == 1
            if rem(col,2) == 0 % even column number is offset
                BasePoints(i,2) = BasePoints(1,2)-dOffY; %offset down
            else % odd column number is same as first column.
                BasePoints(i,2) = BasePoints(1,2);
            end
        elseif j == 2
            if rem(col,2) == 0
                BasePoints(i,2) = BasePoints(1,2)+dOffY; %offset up
            else
                BasePoints(i,2) = BasePoints(1,2)+dY;
            end
        elseif j == 3
            if rem(col,2) == 0
                BasePoints(i,2) = BasePoints(1,2);
            else
                BasePoints(i,2) = BasePoints(1,2)+dOffY;
            end
        elseif j == 4
            if rem(col,2) == 0
                BasePoints(i,2) = BasePoints(1,2)-dY;
            else
                BasePoints(i,2) = BasePoints(1,2)-dOffY;
            end
        end
        
        if rem(col,2) == 0
            r = six(b+1,:);
        else
            r = six(b,:); %initialize for next column--to right
            %(b=2) or left (b=5)
        end
        
        
        i = i + 1;
        ychk = 1;
        
    end % end while xchk==1
    
    if j == 1 %initialize for pass j=2 for-loop.
        
        cp = f.q2; % above centerpoint and up
        [six] = adjCells(cp,centroids);
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1);
        BasePoints(i,2) = BasePoints(1,2)+dY;
        i = i+1;
        
        ychk = 1;
        xchk = 1;
        col = 1; a = 4; b = 2; c = 1; d = 1;
        % travel up the column (a=4),
        % next column is to the right (b=2),
        % up (c=1)
        % right (d=1)
        r = six(b,:);
        
    elseif j == 2 %initialize for pass j=3 for-loop
        
        cp = f.q3; % left and up
        [six] = adjCells(cp,centroids);
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1)-dX;
        BasePoints(i,2) = BasePoints(1,2)+dOffY;
        i = i+1;
        
        ychk = 1;
        xchk = 1;
        col = 1; a = 4; b = 5; c = 1; d = -1;
        %         up     left   up    left
        r = six(b,:);
        
    elseif j == 3 % initialize for pass j=4 for-loop
        
        cp = f.q4; % left and down
        [six] = adjCells(cp,centroids);
        
        InputPoints(i,:) = cp;
        BasePoints(i,1) = BasePoints(1,1)-dX;
        BasePoints(i,2) = BasePoints(1,2)-dOffY;
        i = i+1;
        
        ychk = 1;
        xchk = 1;
        col = 1; a = 1; b = 5; c = -1; d = -1;
        %         down   left   down    left
        r = six(b,:);
        
    end % end of if else...
    waitbar(j*0.25,hwait);
end % end of for j=1:4
close(hwait);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% called by getControlPoints
function [gridOrient] = gridOrientation(cp,centroids,gridType)

dist = rnorm([centroids(:,1)-cp(1) centroids(:,2)-cp(2)]);
r = sort(dist);

if strcmp(gridType,'square')==1
    
    r = r(2:5);
    
    idx = find(dist == r(1) | dist == r(2) | dist == r(3) | dist == r(4));
    
    four = centroids(idx,:);
    
    % reorder four
    angles = (180/pi)*(atan2(four(:,2)-cp(2),four(:,1)-cp(1)));
    perfectlyPlacedGrid = [0;0;0;0];
    
    idx = find(abs(angles+90)==min(abs(angles+90)));
    perfectlyPlacedGrid(idx) = 90;
    idx = find(abs(angles)==min(abs(angles)));
    perfectlyPlacedGrid(idx) = 0;
    idx = find(abs(angles-90)==min(abs(angles-90)));
    perfectlyPlacedGrid(idx) = -90;
    idx = find(abs(abs(angles)-180) == min(abs(abs(angles)-180)));
    perfectlyPlacedGrid(idx) = -180;
    
elseif strcmp(gridType,'hexagonal')==1
    
    r = r(2:7);
    
    idx = find(dist == r(1) | dist == r(2) | dist == r(3)...
        | dist == r(4)| dist == r(5)| dist == r(6));
    
    six = centroids(idx,:);
    
    % reorder six
    % angles in radians from the center point to each of the six points
    angles = (180/pi)*(atan2(six(:,2)-cp(2),six(:,1)-cp(1)));
    perfectlyPlacedGrid = [0;0;0;0;0;0];
    
    for n = 1:2
        % First try
        if n == 1
            modifier = 0;
        % Second try - grid must be closer to horizontal orientation
        elseif n == 2
            % Subtract or add 30 to angles to get closer to vertical orientation
            if any((angles+30)>180)
                modifier = -30;
            else
                modifier = 30;
            end
        end
        angles = angles + modifier;
        
        idx = find(abs(angles+90)==min(abs(angles+90)));
        perfectlyPlacedGrid(idx) = 90;
        idx = find(abs(angles+30)==min(abs(angles+30)));
        perfectlyPlacedGrid(idx) = 30;
        idx = find(abs(angles-30)==min(abs(angles-30)));
        perfectlyPlacedGrid(idx) = -30;
        idx = find(abs(angles-90)==min(abs(angles-90)));
        perfectlyPlacedGrid(idx) = -90;
        idx = find(abs(angles-150)==min(abs(angles-150)));
        perfectlyPlacedGrid(idx) = -150;
        idx = find(abs(angles+150)==min(abs(angles+150)));
        perfectlyPlacedGrid(idx) = 150;
        % If calculated angles are similar, continue. 
        % If not, retry and add/subtract 30 from angles
        if diff(angles+perfectlyPlacedGrid) < 2
            break   
        end
    end
    
    if strcmp(gridType,'square')
        gridOrient = mean(angles+perfectlyPlacedGrid);
    elseif strcmp(gridType,'hexagonal')
        gridOrient = mean(angles+perfectlyPlacedGrid) - modifier; % Remove modifier to get true grid orientation
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [four] = adj4Cells(cp,centroids)

% returns the 4 closest cells to cp reorders so that 1st row is
% always up then moves clockwise 2=right 3=down 4=left

dist = rnorm([centroids(:,1)-cp(1) centroids(:,2)-cp(2)]);
r = sort(dist);
r = r(2:5);

idx = find(dist == r(1) | dist == r(2) | dist == r(3) | dist == r(4));

four = centroids(idx,:);

% reorder four
angles = (180/pi)*(atan2(four(:,2)-cp(2),four(:,1)-cp(1)));

idx = find(abs(angles+90)==min(abs(angles+90)));
U = four(idx,:);
idx = find(abs(angles)==min(abs(angles)));
R = four(idx,:);
idx = find(abs(angles-90)==min(abs(angles-90)));
D = four(idx,:);
idx = find(abs(abs(angles)-180) == min(abs(abs(angles)-180)));
L = four(idx,:);

four = [U;R;D;L];



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% called by getControlPoints, getPts6
function [six,sm] = adjCells(cp,centroids)

% returns the 6 closest cells to cp
% [note: dist is for all centroids, not just neighborhood of cp.]
dist = rnorm([centroids(:,1)-cp(1) centroids(:,2)-cp(2)]);
r = sort(dist);
r = r(2:7); % r is distance--i.e., radius--from cp to point.

idx = find(dist == r(1) | dist == r(2) | dist == r(3)...
    | dist == r(4)| dist == r(5)| dist == r(6));

six = centroids(idx,:);

sm = std(r)/mean(r);

% reorder six
% angles from the center point to each of the six points
angles = (180/pi)*(atan2(six(:,2)-cp(2),six(:,1)-cp(1)));

idx = find(abs(angles+90)==min(abs(angles+90)));
U = six(idx,:);
idx = find(abs(angles+30)==min(abs(angles+30)));
UR = six(idx,:);
idx = find(abs(angles-30)==min(abs(angles-30)));
DR = six(idx,:);
idx = find(abs(angles-90)==min(abs(angles-90)));
D = six(idx,:);
idx = find(abs(angles-150)==min(abs(angles-150)));
DL = six(idx,:);
idx = find(abs(angles+150)==min(abs(angles+150)));
UL = six(idx,:);

six = [U;UR;DR;D;UL;DL];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% called by getPts6
function [ychk] = ycheck(angles)

ychk = 1;

if length(angles) == 4 % rectilinear grid
    
    a = [90;0;90;180];
    
    for i = 1:4
        
        hdiff = a(i) - abs(angles(i));
        if  hdiff > 20 || hdiff < -20
            ychk = 0;
        end
    end
    
elseif length(angles) == 6 %hexagonal grid
    
    a = [90;30;30;90;150;150]; % ideal absolute angles
    
    for i = 1:6
        
        hdiff = a(i) - abs(angles(i));
        if  hdiff > 20 || hdiff < -20
            ychk = 0;
        end
    end
    
else
    return  %error? Grid of neither type.
    ychk = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ip, bp, dY, dX, dOffY, gridOrient] = getControlPoints(I, level, sizeLim, dY, dX, dOffY)
% getControlPoints is called when "Set" is clicked (called with 3
% parameters, not six!).
% ip is input points; bp is base points.
% returns ip--coordinates of blue asterisks.
udim = size(I,2);
vdim = size(I,1);
% We don't want to include insufficiently sized cells that may alter
% our calculations, so color the pixels from such cells black.
%  Before calling getControlPoints, caller must do this:
%   badPixels = removesmalls(I,level,sizeLim);
%           set(handles.hsSize,'UserData',badPixels);
%   I(badPixels) = 0;

% convert original image to black and white
bw = im2bw(I,level);

% get centroid data from grid image
bwlabeled = bwlabel(bw,4);
celldata = regionprops(bwlabeled,'Centroid');
centroids = cat(1,celldata.Centroid);

%find the centroid closest to the center of the image
% Note: if bounding box dimensions are even, "centroid" cannot be true
% center of mass.
dist = rnorm([centroids(:,1)-udim/2 centroids(:,2)-vdim/2]);
idx = find(dist == min(dist));

cp = centroids(idx,:); % set current point to

% determine if grid cells are staggered or square
% if sm is below 0.08, then it is the staggered grid
[six,sm] = adjCells(cp,centroids);

% set the grid type and get an orientation
if sm > 0.08
    
    % set grid type
    gridType = 'square';
    gridOrient = gridOrientation(cp,centroids,gridType);
    
else
    
    % set grid type
    gridType = 'hexagonal';
    gridOrient = gridOrientation(cp,centroids,gridType);
    
end



% rotate the centroids for processing (they will be returned to
% original location after "true" grid is determined)
% move rotation point to center of the image
centroids(:,1) = centroids(:,1) - (udim/2); %translate center to cp
centroids(:,2) = centroids(:,2) - (vdim/2);

% rotate centroids to match true orientation
rotm = [cosd(gridOrient) -sind(gridOrient); sind(gridOrient) cosd(gridOrient)];
centroids = centroids*rotm;

% move back to original image coordinate system
centroids(:,1) = centroids(:,1) + (udim/2); %trans. center to original pt.
centroids(:,2) = centroids(:,2) + (vdim/2);

cp = centroids(idx,:);

if sm > 0.08
    
    [four] = adj4Cells(cp,centroids);
    % Get "true grid" points
    if(nargin > 3)
        [bp,ip,dY,dX,dOffY] = getPts4(cp,centroids,...
            four, dY);
    else
        [bp,ip,dY,dX,dOffY] = getPts4(cp,centroids,...
            four);
    end
    
    dX = NaN;
    dOffY = NaN;
    
else
    
    [six,sm] = adjCells(cp,centroids);
    if(nargin < 6)
        [bp,ip,dY,dX,dOffY] = getPts6(cp,centroids,six);
    else
        [bp,ip,dY,dX,dOffY] = getPts6(cp,centroids,...
            six, dY, dX, dOffY);
    end
    
end

% move rotation point to center of the image
bp(:,1) = bp(:,1) - (udim/2);
bp(:,2) = bp(:,2) - (vdim/2);
ip(:,1) = ip(:,1) - (udim/2);
ip(:,2) = ip(:,2) - (vdim/2);

% rotate the data ang degrees (clockwise)
rotm = [cosd(-gridOrient) -sind(-gridOrient); sind(-gridOrient) cosd(-gridOrient)];
bp = bp*rotm;
ip = ip*rotm;

% move the data back to image coordinate system
bp(:,1) = bp(:,1) + (udim/2);
bp(:,2) = bp(:,2) + (vdim/2);
ip(:,1) = ip(:,1) + (udim/2);
ip(:,2) = ip(:,2) + (vdim/2);

end


classdef Exam < handle
    %EXAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ImagePlanes
        FolderName
        ImageFolder
        ImageFiles
        FolderToSave
        FileToShow
        SelectedFileToShow
        GreenFileToShow
        SelectedGreenFileToShow
        %mainly for save/load purposes
        OriginalFolder
        OriginalImageFolder

        LastCarmPose
        IgnoreList
        ExamsSubFolderName
        
        CorruptTag
    end
    
    methods
        %  old Load - Use May 23 one below
        function load(Obj)
            % this is to get the IP data
            FldName=Obj.FolderName;
%             images = dir(fullfile(Obj.FolderName,'TCarm_Data/*.tif'));
%             calibdata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.txt'));
            images = dir(fullfile(Obj.FolderName,'TCarm_Data/*.tif'));
            calibdata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.txt'));
                        
            numIPs = size(images, 1);
            IPs = [];
            
            filenumbers = zeros(size(images));
            for i = 1:numIPs % sort data
                
                switch size(images(i).name,2)
                    
                    case 33 % single digit
                        filenumbers(i) = cast(str2double(images(i).name(29)), 'uint8');
                        
                    case 34 % double digit
                        filenumbers(i) = cast(str2double(images(i).name(29:30)), 'uint8');
                end
                
                
            end
            
            [~,indices]=sort(filenumbers);
            images = images(indices);
            calibdata = calibdata(indices);
            
            for i = 1:numIPs
                
                IPs = [IPs; ImagePlane(num2str(i), '', [FldName '\TCarm_Data\' calibdata(i).name], [FldName '\TCarm_Data\' images(i).name])];
                Obj.ImageFiles{i} = [Obj.ImageFolder, images(i).name];
            end
            
            Obj.ImagePlanes = IPs;
        end
        
        function Obj = Exam(Folder,CorruptTag)
            Obj.FolderName = Folder;
            Obj.ExamsSubFolderName='Exams';
            % change the folder name for HTML coding
            ChangedFolderName = strrep(Folder,'\','/');
            Obj.ImageFolder = [ChangedFolderName,'/TCarm_Data/'];
            Obj.ImagePlanes = [];
            Obj.LastCarmPose=0;
            Obj.IgnoreList=[];
            Obj.CorruptTag=CorruptTag;
            
        end
        
        function resetExam(Obj,Folder)
            % reset the Exam like a new class. 
            
            % first delete all the ImagePlanes
            if numel(Obj.ImagePlanes)>0
                for f=numel(Obj.ImagePlanes):-1:1
                    delete(Obj.ImagePlanes(f));
                    Obj.ImagePlanes(f)=[];
                end
            end
            
            %then reset the rest of the parameters
            Obj.FolderName = Folder;
            % change the folder name for HTML coding
            ChangedFolderName = strrep(Folder,'\','/');
            Obj.ImageFolder = [ChangedFolderName,'/TCarm_Data/'];
            Obj.ImagePlanes = [];
            Obj.LastCarmPose=0;
        end
        
        function eraseExam(Obj)
            % reset the Exam like a new class. 
            
            % first delete all the ImagePlanes
            if numel(Obj.ImagePlanes)>0
                for f=numel(Obj.ImagePlanes):-1:1
                    delete(Obj.ImagePlanes(f));
                    Obj.ImagePlanes(f)=[];
                end
            end            
            % erase the folder and subfolders  
            rmdir(Obj.FolderName, 's')                        
            %then reset the rest of the parameters
            Obj.FolderName = [];
            Obj.ImageFolder =[];
            Obj.ImagePlanes = [];
            Obj.LastCarmPose=0;            
        end
                
        function ConvertJPG(Obj)
            
            numimages = max(size(Obj.ImageFiles));
            Obj.FolderToSave = [Obj.FolderName '\Resources\'];
            %removecirc = imread('BlueTemplate3.png');
            %removecirc1 = removecirc(:,:,1);
            bluefil = imread('bluefilter.png');
            
            % check if the folder exists
            if ne(exist(Obj.FolderToSave,'dir'),0)
                D = dir([Obj.FolderToSave '*.jpg']);
                numfiles = length(D(not([D.isdir])))/4; %divide the number of files by 4, since there are 4 different sets
                % check the number of images in folder, if not the same then
                % delete the current images in the folder
                if ne(numimages,numfiles)
                    delete([Obj.FolderToSave '*.jpg']);
                    
                    % create the new images
                    for i = 1:numimages
                        a = imread(Obj.ImageFiles{i});
                        j = num2str(i);
                        
                        Obj.FileToShow{i} = [Obj.FolderToSave j '.jpg'];
                        imwrite(a(:,:,1), Obj.FileToShow{i}); %because some tif may be multipage
                        
                        % this is to create the selected images views
%                         imgselected = a(:,:,1);
%                         newIMG = cat(3, imgselected,imgselected,imgselected);
                        newIMG=a;
                        %Resize the new image to 1024X1024 if it is undersized
                        newIMG=imresize(newIMG,'OutputSize',[1024 1024]);
                        kkk = cat(3, newIMG(:,:,1),newIMG(:,:,1),newIMG(:,:,1));
                        %RGBselected = newIMG + bluefil;
                        RGBselected = kkk + bluefil;
                        Obj.SelectedFileToShow{i} = [Obj.FolderToSave 'selected_' j '.jpg'];
                        imwrite(RGBselected,Obj.SelectedFileToShow{i});
                        
                        % this is to create the assigned images views
                        RGBassigned = cat(3, zeros(size(newIMG(:,:,1))),newIMG(:,:,1),zeros(size(newIMG(:,:,1))));
                        Obj.GreenFileToShow{i} = [Obj.FolderToSave 'green_' j '.jpg'];
                        imwrite(RGBassigned,Obj.GreenFileToShow{i});
                        
                        % this is to create the selected assigned image view
                        Greenselected = RGBassigned + bluefil;
                        Obj.SelectedGreenFileToShow{i} = [Obj.FolderToSave 'greenselect_' j '.jpg'];
                        imwrite(Greenselected,Obj.SelectedGreenFileToShow{i});
                        
                    end
                else
                    % if the same, then just set the paths
                    for i = 1:numimages
                        j = num2str(i);
                        Obj.FileToShow{i} = [Obj.FolderToSave j '.jpg'];
                        Obj.SelectedFileToShow{i} = [Obj.FolderToSave 'selected_' j '.jpg'];
                        Obj.GreenFileToShow{i} = [Obj.FolderToSave 'green_' j '.jpg'];
                        Obj.SelectedGreenFileToShow{i} = [Obj.FolderToSave 'greenselect_' j '.jpg'];
                        
                    end
                end
            else
                % create the folder
                mkdir(Obj.FolderToSave);
                % creat the new images
                for i = 1:numimages
                    a = imread(Obj.ImageFiles{i});
                    j = num2str(i);
                    
                    Obj.FileToShow{i} = [Obj.FolderToSave j '.jpg'];
                    imwrite(a(:,:,1), Obj.FileToShow{i});
                    
                    % this is to create the selected images views
%                     imgselected = a(:,:,1);
%                     newIMG = cat(3, imgselected,imgselected,imgselected);
                    newIMG=a;
                    %Resize the new image to 1024X1024 if it is undersized
                    newIMG=imresize(newIMG,'OutputSize',[1024 1024]);
                    kkk = cat(3, newIMG(:,:,1),newIMG(:,:,1),newIMG(:,:,1));
                    %RGBselected = newIMG + bluefil;
                    RGBselected = kkk + bluefil;
                    Obj.SelectedFileToShow{i} = [Obj.FolderToSave 'selected_' j '.jpg'];
                    imwrite(RGBselected,Obj.SelectedFileToShow{i});
                    
                    % this is to create the assigned images views
                    RGBassigned = cat(3, zeros(size(newIMG(:,:,1))),newIMG(:,:,1),zeros(size(newIMG(:,:,1))));
                    Obj.GreenFileToShow{i} = [Obj.FolderToSave 'green_' j '.jpg'];
                    imwrite(RGBassigned,Obj.GreenFileToShow{i});
                    
                    % this is to create the selected assigned image view
                    Greenselected = RGBassigned + bluefil;
                    Obj.SelectedGreenFileToShow{i} = [Obj.FolderToSave 'greenselect_' j '.jpg'];
                    imwrite(Greenselected,Obj.SelectedGreenFileToShow{i});
                end 
            end    
        end 
        
        %  old add_IP - Use May 23 one below
        function add_IP(Obj,folderpath,imgcell,calibcell)
            % this is to get the IP data
            %FldName=Obj.FolderName;
            %images = dir(fullfile(Obj.FolderName,'TCarm_Data/*.tif'));
            %calibdata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.txt'));
            FldName = folderpath;
            
            %see number of items in the input
            numCell = max(size(imgcell));
            
            for ii = 1:numCell
                images(ii) = dir(imgcell{ii});
                calibdata(ii) = dir(calibcell{ii});
            end
            
            ImgSize = max(size(Obj.ImagePlanes)); % the number of files previously
            %ImgSize = max(size(dir(fullfile(exam_name,'\TCarm_Data\'))));
            
            % change the folder name for HTML coding
            ChangedFolderName = strrep(folderpath,'\','/');
            Obj.ImageFolder = [ChangedFolderName,'/TCarm_Data/'];
            
            numIPs = size(images, 2);
            IPs = [];
            
            filenumbers = zeros(size(images));
            for i = 1:numIPs % sort data
                
                filenumbers(i) = ImgSize + i;
                
                s = images(i).name;
                rep = strfind(s,'_');
                %s = s(1:end-5);
                s = s(1:rep(end));
                thenumber = num2str(filenumbers(i));
                new_img(i).name = [s thenumber '.tif'];
                new_calib(i).name = [s thenumber '.txt'];
                
            end
            
            %[~,indices]=sort(filenumbers);
            %images = images(indices);
            %calibdata = calibdata(indices);
            
            % copy the file over now with the new name
            for i = 1:numIPs
                copyfile(imgcell{i},[FldName '\TCarm_Data\' new_img(i).name],'f');
                copyfile(calibcell{i},[FldName '\TCarm_Data\' new_calib(i).name],'f');
            end
            
            
            
            for i = 1:numIPs
                
                thenumber = num2str(i + ImgSize);
                
                IPs = [IPs; ImagePlane(thenumber, '', [FldName '\TCarm_Data\' new_calib(i).name], [FldName '\TCarm_Data\' new_img(i).name])];
                Obj.ImageFiles{i+ImgSize} = [Obj.ImageFolder, new_img(i).name];
            end
            
            % add the IPs to the end of the previous IP list
            if ~isempty(Obj.ImagePlanes)
                Obj.ImagePlanes = [Obj.ImagePlanes; IPs];
            else
                Obj.ImagePlanes = IPs;
            end
            
        end
        
        % DO NOT USE - added May 13 to include the patient pose and stage info - Use May 23 one below
        function new_add_IP(Obj,folderpath,imgcell,calibcell,AS,PP)
            % this is to get the IP data
            %FldName=Obj.FolderName;
            %images = dir(fullfile(Obj.FolderName,'TCarm_Data/*.tif'));
            %calibdata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.txt'));
            FldName = folderpath;
            
            %see number of items in the input
            numCell = max(size(imgcell));
            
            for ii = 1:numCell
                images(ii) = dir(imgcell{ii});
                calibdata(ii) = dir(calibcell{ii});
            end
            
            ImgSize = max(size(Obj.ImagePlanes)); % the number of files previously
            %ImgSize = max(size(dir(fullfile(exam_name,'\TCarm_Data\'))));
            
            % change the folder name for HTML coding
            ChangedFolderName = strrep(folderpath,'\','/');
            Obj.ImageFolder = [ChangedFolderName,'/TCarm_Data/'];
            
            numIPs = size(images, 2);
            IPs = [];
            
            filenumbers = zeros(size(images));
            for i = 1:numIPs % sort data
                
                filenumbers(i) = ImgSize + i;
                
                s = images(i).name;
                rep = strfind(s,'_');
                %s = s(1:end-5);
                s = s(1:rep(end));
                thenumber = num2str(filenumbers(i));
                new_img(i).name = [s thenumber '.tif'];
                new_calib(i).name = [s thenumber '.txt'];
                
            end
            
            %[~,indices]=sort(filenumbers);
            %images = images(indices);
            %calibdata = calibdata(indices);
            
            % copy the file over now with the new name
            for i = 1:numIPs
                copyfile(imgcell{i},[FldName '\TCarm_Data\' new_img(i).name],'f');
                copyfile(calibcell{i},[FldName '\TCarm_Data\' new_calib(i).name],'f');
            end
            
            for i = 1:numIPs
                thenumber = num2str(i + ImgSize);
                IPs = [IPs; ImagePlane(thenumber, '', [FldName '\TCarm_Data\' new_calib(i).name], [FldName '\TCarm_Data\' new_img(i).name])];
                Obj.ImageFiles{i+ImgSize} = [Obj.ImageFolder, new_img(i).name];
            end
            
            % add the IPs to the end of the previous IP list
            if ~isempty(Obj.ImagePlanes)
                Obj.ImagePlanes = [Obj.ImagePlanes; IPs];
            else
                Obj.ImagePlanes = IPs;
            end
            
        end
        
        % added May 23 to include the patient pose and stage info and ref
        function new_add_IP_with_ref(Obj,folderpath,imgcell,calibcell,matcell,AS,PP,NumSteps)
            % this is to get the IP data

            FldName = folderpath;
            
            %see number of items in the input
            numCell = max(size(imgcell));
            
            for ii = 1:numCell
                images(ii) = dir(imgcell{ii});
                calibdata(ii) = dir(calibcell{ii});
                if ~isempty(matcell)
                    matfiles(ii) = dir(matcell{ii});
                end
            end
            
            ImgSize = max(size(Obj.ImagePlanes)); % the number of files previously

            % change the folder name for HTML coding
            ChangedFolderName = strrep(folderpath,'\','/');
            Obj.ImageFolder = [ChangedFolderName,'/TCarm_Data/'];
            
            numIPs = size(images, 2);
            IPs = [];
            
            % Check what the largest number is for existing filenames
            largestNum = 0;
            d1=dir(fullfile([FldName '\TCarm_Data\*.tif']));
            if ~isempty(d1)
                imageNames = sort_nat({d1.name});   % Sort image files in numerical order
                ind1 = max(strfind(imageNames{end},'_'))+1;
                largestNum = str2double(imageNames{end}(ind1:end-4));
            end
            
            filenumbers = zeros(size(images));
            for i = 1:numIPs % sort data
                
                % Start file numbering from largest existing file number or number of ImagePlanes
                if largestNum > ImgSize
                    filenumbers(i) = largestNum + i;
                else
                    filenumbers(i) = ImgSize + i;
                end
                
                s = images(i).name;
                rep = strfind(s,'_');
                %s = s(1:end-5);
                s = s(1:rep(end));
                thenumber = num2str(filenumbers(i));
                new_img(i).name = [s thenumber '.tif'];
                new_calib(i).name = [s thenumber '.txt'];
                if ~isempty(matcell)
                    new_mat(i).name = [s thenumber '.mat'];
                end
            end
            

            % copy the file over now with the new name
            for i = 1:numIPs
                
                % Check if image file is corrupt before copying
                while(1)
                    try
                        imread(imgcell{i});
                        break
                    catch
                    end
                end
                copyfile(imgcell{i},[FldName '\TCarm_Data\' new_img(i).name],'f');
                
                % Copy the txt calib file
                copyfile(calibcell{i},[FldName '\TCarm_Data\' new_calib(i).name],'f');
                
                % Check if mat file is corrupt before copying
                if ~isempty(matcell)
                    while(1)
                        try
                            load(matcell{i});
                            break
                        catch
                        end
                    end
                    copyfile(matcell{i},[FldName '\TCarm_Data\' new_mat(i).name],'f');
                end
                
            end
            
            for i = 1:numIPs
                thenumber = num2str(i + ImgSize);
                IPs = [IPs; ImagePlane(thenumber, '', [FldName '\TCarm_Data\' new_calib(i).name], [FldName '\TCarm_Data\' new_img(i).name])];
                Obj.ImageFiles{i+ImgSize} = [Obj.ImageFolder, new_img(i).name];
            end
            
            % add the IPs to the end of the previous IP list
            if ~isempty(Obj.ImagePlanes)
                Obj.ImagePlanes = [Obj.ImagePlanes; IPs];
            else
                Obj.ImagePlanes = IPs;
            end
            
        end
        
        % added May 23 to include ref
        function load_with_ref(Obj,NumSteps)
            % this is to get the IP data
            FldName=Obj.FolderName;
            convertTIFtoJPG(FldName);
            % images = dir(fullfile(Obj.FolderName,'TCarm_Data/*.tif'));            
            % calibdata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.txt'));
            posedata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.mat'));
            % thumbimages = dir(fullfile(Obj.FolderName,'TCarm_Data/*.jpg'));            
              
            %numIPs = numel(images);
            numIPs = numel(posedata);
            IPs = [];
                        
            LstPose=0;
            LstPPose=0;
            
            %[~,idx] = sort([images.datenum]);
            [~,idx] = sort([posedata.datenum]);
            Obj.ImageFiles={};
            for ii = 1:numIPs

                i=idx(ii);
                PlaneName=num2str(i);
                Path='';
                JTPoseFileName=[FldName '\TCarm_Data\' posedata(i).name];                
                clear ImgPlane;
                load(JTPoseFileName,'ImgPlane');
                Im=ImgPlane;
                Im.Name=PlaneName;
                Im.Path=Path;
                Im.JTPoseFileName=JTPoseFileName;                                                             
                Im.CarmPose.PoseIndex=i; % correct the PoseIndex

                % Updating the last patient pose and Carm pose according to
                % the loaded data:
                if Im.CarmPose.PoseIndex> LstPose
                    LstPose=Im.CarmPose.PoseIndex;
                end
                IPs=[IPs; Im];
            end
            Obj.ImagePlanes = IPs;            
            Obj.LastCarmPose=LstPose;
            % Now load the .ref files and update the refgrid fields for each ivdividual ImagePlane:
            Obj.LoadImagePlanesRefGridFiles;            
        end
        
        
        % added May 23 to include ref
        function AddtoFilesGridInfo(Obj,GridOrg,GridVx,GridVy,GridVz)
            % this is to get the IP data
            FldName=Obj.FolderName;
            convertTIFtoJPG(FldName);
            images = dir(fullfile(Obj.FolderName,'TCarm_Data/*.tif'));            
            calibdata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.txt'));
            posedata = dir(fullfile(Obj.FolderName,'TCarm_Data/*.mat'));
            thumbimages = dir(fullfile(Obj.FolderName,'TCarm_Data/*.jpg'));            
              
            numIPs = numel(images);
            IPs = [];
                        
            LstPose=0;
            LstPPose=0;
            
            [~,idx] = sort([images.datenum]);
            Obj.ImageFiles={};
            for ii = 1:numIPs
                i=idx(ii);
                PlaneName=num2str(i);
                Path='';
                JTPoseFileName=[FldName '\TCarm_Data\' posedata(i).name];                
                clear ImgPlane;
                load(JTPoseFileName,'ImgPlane');
                Im=ImgPlane;
                ImgPlane=ImagePlane(Im.Name);
               
                ImgPlane.Name=Im.Name;
                ImgPlane.JTCalibFileName=Im.JTCalibFileName;
                ImgPlane.JTImageFileName=Im.JTImageFileName;
                ImgPlane.JTPoseFileName=Im.JTPoseFileName;
                ImgPlane.PrincipalDist=Im.PrincipalDist;
                ImgPlane.Xoffset=Im.Xoffset;
                ImgPlane.Yoffset=Im.Yoffset;
                ImgPlane.PixelScale=Im.PixelScale;
                ImgPlane.SourcePosition=Im.SourcePosition;
                ImgPlane.ImageNormal=Im.ImageNormal;
                ImgPlane.ImageUp=Im.ImageUp;
                ImgPlane.CentrePoint=Im.CentrePoint;
                ImgPlane.Features=Im.Features;
                ImgPlane.Path=Im.Path;
                ImgPlane.Image=Im.Image;

                ImgPlane.CarmPose=Im.CarmPose;                
                ImgPlane.RefGridOrg=GridOrg;
                ImgPlane.RefGridVx=GridVx;
                ImgPlane.RefGridVy=GridVy;
                ImgPlane.RefGridVz=GridVz;
                
                save(JTPoseFileName,'ImgPlane');
            end            
            
        end
        
   
        % added June 21 2018 to speed up the calibration process
        function RefreshFilesWithNewGridInfo_wo_LoadSave(Obj,SM,Old_GridOrg,Old_GridVx,Old_GridVy,Old_GridVz,GridOrg,GridVx,GridVy,GridVz,ShowGraphics)
            % Update the ImagePlane GridInfo wo Load and Save; Only apply to the Obj.ImagePlanes.
            if (ShowGraphics)
                wb_h=waitbar(0,'Revising ImagePlans with New Grid Info ..');
%                 changeFigIcon(wb_h);
            end
            numIPs=numel(Obj.ImagePlanes);
            for ii = 1:numIPs
                if (ShowGraphics)
                    waitbar(ii/numIPs);
                end
                Obj.ImagePlanes(ii)=adjustImagePlaneGrid(SM,Obj.ImagePlanes(ii),Old_GridOrg,Old_GridVx,Old_GridVy,Old_GridVz,GridOrg,GridVx,GridVy,GridVz);
            end
            if (ShowGraphics)
                close(wb_h);
            end
        end
        
        % added Aug 16 2018: overwriting the ImagePlane files with the new
        % grid info:
        function OverwriteFilesWithNewGridInfo(Obj,ShowGraphics)
            % Update the ImagePlane GridInfo wo Load and Save; Only apply to the Obj.ImagePlanes.
            if (ShowGraphics)
                wb_h=waitbar(0,'Overwriting ImagePlans with New Grid Info ..');
                changeFigIcon(wb_h);
            end
            numIPs=numel(Obj.ImagePlanes);           
            for ii = 1:numIPs
                if (ShowGraphics)
                    waitbar(ii/numIPs);
                end
                Filename=sprintf('%s\\MPJT_Calib_PhantomObj_Frame_%i.mat',Obj.ImageFolder,ii);
                ImgPlane=Obj.ImagePlanes(ii);
                save(Filename,'ImgPlane');                
            end                        
            if (ShowGraphics)
                close(wb_h);
            end
        end        

        % added Sep 28 2018: save the individual imageplane grid info:
        function UpdateImagePlanesRefGridFiles(Obj)
            numIPs=numel(Obj.ImagePlanes);           
            for ii = 1:numIPs
                Filename=sprintf('%s\\MPJT_Calib_PhantomObj_Frame_%i.ref',Obj.ImageFolder,ii);
                % if the filename exists then load the .ref file
                if exist(Filename,'file')
                   load(Filename,'-mat');
                   sz=numel(IPRef);
                else
                    IPRef=[];
                    IPRef.RefGridOrg=[];
                    IPRef.RefGridVx=[];
                    IPRef.RefGridVy=[];
                    IPRef.RefGridVz=[];
                    IPRef.DateTime=[]; % Add the date information for a complete record
                    sz=0;
                end
                % now add the latest PlaneRef Info on the very top                
                IPRef(sz+1).RefGridOrg=Obj.ImagePlanes(ii).RefGridOrg;
                IPRef(sz+1).RefGridVx=Obj.ImagePlanes(ii).RefGridVx;
                IPRef(sz+1).RefGridVy=Obj.ImagePlanes(ii).RefGridVy;
                IPRef(sz+1).RefGridVz=Obj.ImagePlanes(ii).RefGridVz;
                IPRef(sz+1).DateTime=datestr(now);
                save(Filename,'IPRef');                
            end                        
        end                     
                
        % added Sep 28 2018: load the individual imageplane grid info:
        function LoadImagePlanesRefGridFiles(Obj)
            numIPs=numel(Obj.ImagePlanes);           
            for ii = 1:numIPs
                Filename=sprintf('%s\\MPJT_Calib_PhantomObj_Frame_%i.ref',Obj.ImageFolder,ii);
                % if the filename exists then load the .ref file
                if exist(Filename,'file')
                   load(Filename,'-mat');
                   if (~isempty(IPRef(end).RefGridOrg) && ...
                       ~isempty(IPRef(end).RefGridVx) && ...
                       ~isempty(IPRef(end).RefGridVy) && ...
                       ~isempty(IPRef(end).RefGridVz))
                            % load the latest saved reference for each
                            % individual imageplane:
                            Obj.ImagePlanes(ii).RefGridOrg=IPRef(end).RefGridOrg;
                            Obj.ImagePlanes(ii).RefGridVx=IPRef(end).RefGridVx;
                            Obj.ImagePlanes(ii).RefGridVy=IPRef(end).RefGridVy;
                            Obj.ImagePlanes(ii).RefGridVz=IPRef(end).RefGridVz;                            
                   end
                end
            end                        
        end         
        
        
        % added June 21 2018 to speed up the calibration process
        function SaveImagePlanes(Obj)
            % Update the ImagePlane GridInfo wo Load and Save; Only apply to the Obj.ImagePlanes.
            
            numIPs=numel(Obj.ImagePlanes);
            for ii = 1:numIPs
                waitbar(ii/numIPs);
                ImgPlane=Obj.ImagePlanes(ii);
                JTPoseFileName=ImgPlane.JTPoseFileName;
                save(JTPoseFileName,'ImgPlane');
            end
            close(wb_h);
        end        
        
        
      
        
        function add_Convert(Obj)
            
            numimages = max(size(Obj.ImageFiles));
            Obj.FolderToSave = [Obj.FolderName '\Resources\'];
            %removecirc = imread('BlueTemplate3.png');
            %removecirc1 = removecirc(:,:,1);
            bluefil = imread('bluefilter.png');
            
            % check if the folder exists
            if ne(exist(Obj.FolderToSave,'dir'),0)
                D = dir([Obj.FolderToSave '*.jpg']);
                numfiles = length(D(not([D.isdir])))/4; %divide the number of files by 4, since there are 4 different sets
                % check the number of images in folder, if not the same then
                % delete the current images in the folder
                if ne(numimages,numfiles)
                    delete([Obj.FolderToSave '*.jpg']);
                    
                    % create the new images
                    for i = 1:numimages
                        a = imread(Obj.ImageFiles{i});
                        j = num2str(i);
                        
                        Obj.FileToShow{i} = [Obj.FolderToSave j '.jpg'];
                        imwrite(a(:,:,1), Obj.FileToShow{i}); %because some tif may be multipage
                        
                        % this is to create the selected images views
%                         imgselected = a(:,:,1);
%                         newIMG = cat(3, imgselected,imgselected,imgselected);
                        newIMG=a;
                        %Resize the new image to 1024X1024 if it is undersized
                        newIMG=imresize(newIMG,'OutputSize',[1024 1024]);                        
                        kkk = cat(3, newIMG(:,:,1),newIMG(:,:,1),newIMG(:,:,1));
                        %RGBselected = newIMG + bluefil;
                        RGBselected = kkk + bluefil;
                        Obj.SelectedFileToShow{i} = [Obj.FolderToSave 'selected_' j '.jpg'];
                        imwrite(RGBselected,Obj.SelectedFileToShow{i});
                        
                        % this is to create the assigned images views
                        RGBassigned = cat(3, zeros(size(newIMG(:,:,1))),newIMG(:,:,1),zeros(size(newIMG(:,:,1))));
                        Obj.GreenFileToShow{i} = [Obj.FolderToSave 'green_' j '.jpg'];
                        imwrite(RGBassigned,Obj.GreenFileToShow{i});
                        
                        % this is to create the selected assigned image view
                        Greenselected = RGBassigned + bluefil;
                        Obj.SelectedGreenFileToShow{i} = [Obj.FolderToSave 'greenselect_' j '.jpg'];
                        imwrite(Greenselected,Obj.SelectedGreenFileToShow{i});
                        
                    end
                else
                    % if the same, then just set the paths
                    for i = 1:numimages
                        j = num2str(i);
                        Obj.FileToShow{i} = [Obj.FolderToSave j '.jpg'];
                        Obj.SelectedFileToShow{i} = [Obj.FolderToSave 'selected_' j '.jpg'];
                        Obj.GreenFileToShow{i} = [Obj.FolderToSave 'green_' j '.jpg'];
                        Obj.SelectedGreenFileToShow{i} = [Obj.FolderToSave 'greenselect_' j '.jpg'];
                        
                    end
                end
            else
                % create the folder
                mkdir(Obj.FolderToSave);
                % creat the new images
                for i = 1:numimages
                    a = imread(Obj.ImageFiles{i});
                    j = num2str(i);
                    
                    Obj.FileToShow{i} = [Obj.FolderToSave j '.jpg'];
                    imwrite(a(:,:,1), Obj.FileToShow{i});
                    
                    % this is to create the selected images views
%                     imgselected = a(:,:,1);
%                     newIMG = cat(3, imgselected,imgselected,imgselected);
                    newIMG=a;
                    %Resize the new image to 1024X1024 if it is undersized
                    newIMG=imresize(newIMG,'OutputSize',[1024 1024]);                    
                    kkk = cat(3, newIMG(:,:,1),newIMG(:,:,1),newIMG(:,:,1));
                    %RGBselected = newIMG + bluefil;
                    RGBselected = kkk + bluefil;
                    Obj.SelectedFileToShow{i} = [Obj.FolderToSave 'selected_' j '.jpg'];
                    imwrite(RGBselected,Obj.SelectedFileToShow{i});
                    
                    % this is to create the assigned images views
                    RGBassigned = cat(3, zeros(size(newIMG(:,:,1))),newIMG(:,:,1),zeros(size(newIMG(:,:,1))));
                    Obj.GreenFileToShow{i} = [Obj.FolderToSave 'green_' j '.jpg'];
                    imwrite(RGBassigned,Obj.GreenFileToShow{i});
                    
                    % this is to create the selected assigned image view
                    Greenselected = RGBassigned + bluefil;
                    Obj.SelectedGreenFileToShow{i} = [Obj.FolderToSave 'greenselect_' j '.jpg'];
                    imwrite(Greenselected,Obj.SelectedGreenFileToShow{i});
                end 
            end    
        end 
        
        function updateIgnoreList(Obj, FrameNo)
            
            %update the IgnoreList first
            [Found, Indx]=ismember(FrameNo,Obj.IgnoreList);
            if Found % if the item exists then remove from the list
                Obj.IgnoreList(Indx)=[];
            else     % if not found then add it to the list
                Obj.IgnoreList(end+1)=FrameNo;
            end
            
            %Now write it in the file directory:
            %FName=[Obj.FolderName,'\IngnoreList.IgInfo'];
            IgInfoFileName=Obj.returnIgnoreFile(Obj.FolderName);
            FName=[Obj.FolderName,'\',IgInfoFileName];
            
            fid=fopen(FName,'w');                        
            if ~isempty(Obj.IgnoreList)
                fprintf(fid,'%i ', int16(Obj.IgnoreList));
            end
            fclose(fid);                        
        end
        

        function [Found]=ShallIgnore(Obj, FrameNo)
            % Returns the answer if the FrameNo shall be ignored:
            [Found, Indx]=ismember(FrameNo,Obj.IgnoreList);
        end
        
        function loadIgnoreList(Obj)            
            %Load the ingorelist if the file exists            
            % FName=[Obj.FolderName,'\IngnoreList.IgInfo'];
            IgInfoFileName=Obj.returnIgnoreFile(Obj.FolderName);
            FName=[Obj.FolderName,'\',IgInfoFileName];
            if exist(FName,'file');
                fid=fopen(FName,'r');       
                Obj.IgnoreList=fscanf(fid,'%i ');
                fclose(fid);
            else
                Obj.IgnoreList=[];
            end
            
        end

        function eraseIgnoreList(Obj)
            %Load the ingorelist if the file exists
            % FName=[Obj.FolderName,'\IngnoreList.IgInfo'];
            IgInfoFileName=Obj.returnIgnoreFile(Obj.FolderName);
            FName=[Obj.FolderName,'\',IgInfoFileName];            
            if exist(FName,'file');
                delete(FName);
            end
            Obj.IgnoreList=[];            
        end

        function str=shortExamPath(Obj)
            % return the short form of the exam path
            indx=strfind(Obj.FolderName,[Obj.ExamsSubFolderName,'\']);
            if isempty(indx)
                str=[];
            else
                str=Obj.FolderName(indx(end):end);
            end                        
        end
        
        function str=minExamPath(Obj)
            % return the short form of the exam path
            indx=strfind(Obj.FolderName,[Obj.ExamsSubFolderName,'\']);
            if isempty(indx)
                str=[];
            else
                str=Obj.FolderName((indx(end)+numel(Obj.ExamsSubFolderName)+1):end);
            end                        
        end
        
        function FileName=returnIgnoreFile(Obj,FolderName)
            FileName='IgnoreList.IgInfo';
            try
                FName=[FolderName,'\*.IgInfo'];
                d=dir(FName);
                if ~isempty(d)
                   FileName=d(1).name; 
                end
            catch
                fprintf('\n\rIgInfo file not found. Using the defailt IgInfo file.\n\r');
                FileName='IgnoreList.IgInfo';
            end
        end

        
    end
    
end


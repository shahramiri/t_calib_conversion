classdef FrameGrabber < handle
% FrameGrabber Class is responsible to connecting to a framecapture device
% and acquiring, scaling, cropping, masking, and undistorting of the frame
    properties
        ImageFolder;    % the folder the images are saved in
        ImageName;      % the name of the images
        ImageCount;     % the current total of images
        VideoIndex;
        DefaultSize;
        CropSize;
        IsInitialized;
        Video;
        DeviceName;
        DeviceFormat;
        PlaneIndex;         % To determine whether AP or ML must be used for adding images
        CarmModelImagePars; % Carm-model specific rotation and flipping
                            % CarmModelImagePars(1): Rotation (ccw) degrees
                            % CarmModelImagePars(2): 1 or 0 for flipping
                            % about the vertical axis after rotation
        CFGFileName
        MSKFileName
        OriginalCapture
        CroppedCapture
        UnDistortObject
        UnDistortFlag
        %DAQ % retire the DAQ
        Carm
        SimMode             % Simulation mode: -> 1: Used SimExam folder to randomly select a fluoro view instead of receiving an image
        ImageMask           % Masking any area of the image before processing
        
        CaptureRes          % Resultion of the Captured Image (should match CaptureRes is available)
        OriginalRes         % [Opional] Intended Resolution of the FrameGrabber
        ForceOriginalRes    % [Optional] Force Resolution
        
        AutoCheck                   % Autocheck Flag
        %AutoCheckBusy               % Autocheck busy flag for internal use to avoid running interferring commands 
        AutoCheckRefImage           % Autocheck Reference Image
        AutoCheckLastImage          % Autocheck Last Image Check
        AutoCheckStatus             % 0: No Auto Image
                                    % 1: Connection good - New Image Available
                                    % 2: Attempting to Connect to the framegrabber
                                    % 3: Failed in Reinitiating Connection
                                    % 4: Connection good - No new Image Available
        AutoCheckTriggerFolder      % Folder that if contains any file triggers auto imager capture is the SimMode        
        AutoCheckImageThreshold     % Image discrepancy threshold in percentage for picking up a new new image through autocheck
        AutoCheckTimeThreshold      % Time threshold for taking
        AutoCheckLoops              % Number of times snapshots have to be taken one after another in a loop every time a new image is in; this is prevent half images (Default is 2)
        AutoCheckLastTime           % Clock when last image was autocollected
    end
    
    methods
        
        function Obj = FrameGrabber(folder,name)
            Obj.ImageFolder = folder;
            Obj.ImageName = name;
            Obj.ImageCount = 0;
            Obj.VideoIndex=0;
            Obj.DefaultSize=[1024 1024];  % Set as standard for all systems
            Obj.CropSize=[8 0 1023 1024]; % Initial  Crop Size
            Obj.IsInitialized=0;
            Obj.Video=[];
            Obj.DeviceName='';
            Obj.DeviceFormat='';
            Obj.CFGFileName=[];
            Obj.MSKFileName=[];
            Obj.OriginalCapture=[];
            Obj.CroppedCapture=[];
            Obj.UnDistortFlag=0;
            % Obj.DAQ=[];
            Obj.Carm=[];
            Obj.SimMode=0;
            Obj.ForceOriginalRes=0;
            Obj.AutoCheck=0;
            Obj.AutoCheckRefImage=[];
            Obj.AutoCheckLastImage=[];
            Obj.AutoCheckStatus=0;
            Obj.AutoCheckTriggerFolder='SimExam\AutoImageTrigger';
            Obj.AutoCheckImageThreshold=0.1;    % Default is 0.1%
            Obj.AutoCheckTimeThreshold=1;       % Time difference required to register an image as a new image
            Obj.AutoCheckLastTime=clock;        % Clock when the last image was autochecked
            Obj.AutoCheckLoops=2;

        end
                
        % Reset the Autocheck after an image is successfully retrieved
        function resetAuto(Obj)              
            try
                if (Obj.SimMode)
                    if exist(Obj.AutoCheckTriggerFolder,'dir')
                        delete([Obj.AutoCheckTriggerFolder,'\*.*']); % delete the file after triggering the AutoCheckImage;
                    end
%                     Obj.AutoCheckStatus=0;                
                else

                end
                Obj.AutoCheckRefImage=[];
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.resetAuto()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end
        end
        
        % Start the AutoCheck 
        function startAuto(Obj)
            try 
                Obj.AutoCheck=1;
                Obj.resetAuto;
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.startAuto()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end                
        end
        
        % Stop the AutoCheck 
        function stopAuto(Obj)
            try
                Obj.AutoCheck=0;
                Obj.resetAuto;                
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.stopAuto()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end 
        end
        
        % Wait for logging
        function waitForLogging(Obj)   
            try
                if islogging(Obj.Video)
                    while islogging(Obj.Video)  % Pause while video is being logged.
                        pause(0.001);      % This seems to prevent skips in data for very fast reloads.
                    end
                end
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.waitForLogging()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end
        end
        
        function [VidAvailable, Img, VideoTimeout]=checkVideo(Obj,debugtext)  
            % debugtext used for debugging purposes
            tic;
            if ~isempty(debugtext)
                fprintf('<%s ,',debugtext);
            end
            VidAvailable=0;
            VideoTimeout=0;
            Img=[];
            if ~isempty(Obj.Video)
                if isvalid(Obj.Video)
                    if (isrunning(Obj.Video))
                        VidRunning=1;
                    else
                        try
                            start(Obj.Video);
                            VidRunning=isrunning(Obj.Video);
                        catch
                            VidRunning=0;
                        end
                    end
                    if (VidRunning)
                        Obj.waitForLogging;
                        try
                            for f=1:Obj.AutoCheckLoops % Repeat snapshot to make sure image is good; this is to prevent incomplete images
                                Img=getsnapshot(Obj.Video);
                            end
                        catch
                           VideoTimeout=1;
                        end
                        if (~VideoTimeout)
                            Obj.waitForLogging;
                            if max(max(Img))>0
                                VidAvailable=1;
                            end
                            Obj.waitForLogging;
                            flushdata(Obj.Video); % flush data if blank image
                            Obj.waitForLogging;
                        end
                    end
                    
                end
            end
            tm=toc;
            if ~isempty(debugtext)
                fprintf('%s >',num2str(tm));
            end
        end
                
        % Check the availability of the new image to be run using
        % clientBytesReceived
        function checkAuto(Obj,SetRef,handles)

        try
            Obj.AutoCheckStatus=2; % Busy Status
            if ~(Obj.SimMode)
                %% Non-Simulation Mode:
                    [VidAvailable,Img,VideoTimeout]=Obj.checkVideo([]); 
                    sendCmdProgress(handles, 0.2);
                    if (~VidAvailable || isempty(Img) || VideoTimeout) % Try to reinitialize connection
                            fprintf('\r\nHard Connection Re-initialization of the Framegrabber..\n');
                            Obj.Initialize(1,handles);
                            sendCmdProgress(handles, 0.21);
                            [VidAvailable,Img, VideoTimeout]=Obj.checkVideo([]);
                    end
                    
                    ContinueProcessing=0;
                    if (~VidAvailable || isempty(Img) || VideoTimeout)
                        Obj.AutoCheckStatus=3; % Error
                        sendCmdProgress(handles, 0.22);
                    else
                        try
                            % Img is available from the above snapshot
                            Obj.AutoCheckLastImage=imresize(Img(:,:,1),0.1);
                            ContinueProcessing=1;
                        catch 
                            ContinueProcessing=0;
                            Obj.AutoCheckStatus=3; % Error
                        end
                    end
                    
                    if (ContinueProcessing)  % Image is available and processing can continue
                        if SetRef
                            Obj.AutoCheckRefImage=Obj.AutoCheckLastImage;
                            Obj.AutoCheckStatus=4;              % Status: (4) No Image is available now but framegrabber is ready
                        else
                            if max(max(Obj.AutoCheckRefImage))==0 % If reference is blank then set it to empty
                                Obj.AutoCheckRefImage=[];
                            end
                            if isempty(Obj.AutoCheckRefImage)
                                Obj.AutoCheckRefImage=Obj.AutoCheckLastImage;
                                Obj.AutoCheckStatus=4;          % Status: (4) No Image is available
                            else
                                sub=imabsdiff(Obj.AutoCheckRefImage,Obj.AutoCheckLastImage);
                                diff_im=im2bw(sub,.1);
                                [xxi,~]=find(diff_im==1);
                                % if > imagethereshold is changed the trigger && % and if enoughtime is passed since the last image
                                if (numel(xxi)/numel(diff_im)>(Obj.AutoCheckImageThreshold/100)) && (etime(clock,Obj.AutoCheckLastTime)>Obj.AutoCheckTimeThreshold)
                                        Obj.AutoCheckStatus=1;      % Status: (1) New Image is available
                                        etime(clock,Obj.AutoCheckLastTime)
                                        Obj.AutoCheckLastTime=clock;
                                else
                                    Obj.AutoCheckStatus=4;      % Status: (4) No Image is available
                                end
                            end
                        end
                    end

            else
            %% Simulation Mode:
            % Check the content of the folder to figure if the gramegrabber needs to be triggered:
                        display('C!');
                if ~exist(Obj.AutoCheckTriggerFolder,'dir')
                    mkdir(Obj.AutoCheckTriggerFolder);
                end
                d=dir([Obj.AutoCheckTriggerFolder,'\*.txt']);
                if numel(d)>0
                    switch (d(1).name)                        
                        case '1.txt'
                            Obj.AutoCheckStatus=1;      % Status: (1) New Image  is available
                        case '2.txt'
                            Obj.AutoCheckStatus=2;      % Status: (2) Re-initiating the VideoCapture Card
                        case '3.txt'
                            Obj.AutoCheckStatus=3;      % Status: (3) Failed Re-initiating the Video Capture Card
                    end
                else
                    if (Obj.AutoCheck)
                        Obj.AutoCheckStatus=4;      % Status: (4) No New Image available but framegrabber is active
                    else
                        Obj.AutoCheckStatus=0;      % Status: (0) No new Image is available and framegrabber in inactive
                    end
                end
            end
        catch e %e is an MException struct            
            fprintf('\r\n---------------------------------------');
            fprintf('\r\nError in FrameGrabber.checkAuto()');
            fprintf('\r\n(SimMode: %i)',Obj.SimMode);
            fprintf('\r\nThe identifier was:\n%s',e.identifier);
            fprintf('\r\nThere was an error! The message was:\n%s',e.message);
            fprintf('\r\n---------------------------------------');
        end
%         fprintf('>%i',Obj.AutoCheckStatus);
        end
        
        % Connecting to the FrameGrabber hardware:
        function Initialize(Obj,TotalReset,handles)            
           try      
           Obj.IsInitialized=1;
               if Obj.SimMode
                   display('Using a Simulated Framegrabber. No Attempt made to connect to framegrabber device.');
               else
                   Obj.readConfigFile;                  
                   
                   if (TotalReset)
                       imaqreset;
                       sendCmdProgress(handles,0.01);
                   end
                   
                   % get the list of all available image acquisition hardware
                   imin=imaqhwinfo('winvideo');
                   sendCmdProgress(handles,0.015);
                   DevN=numel(imin.DeviceInfo);
                   if DevN>0
                       DevIDs=zeros(DevN,1);
                       DevNames=cell(DevN,1);
                       for f=1:DevN
                           DevIDs(f)=imin.DeviceInfo(f).DeviceID;
                           DevNames{f}=imin.DeviceInfo(f).DeviceName;
                       end
                   end
                   
                   % see if an active device already exists
                   VidObjs=imaqfind; % List of all image acquisition objects
                   DevFound=0;
                   for g=1:numel(VidObjs) % Loop through all the available VideoObjects
                       if isvalid(VidObjs(g))
                           for f=1:numel(DevIDs) % Find the matching device name by index:
                               if VidObjs(g).DeviceID==DevIDs(f);
                                   if strcmpi(DevNames{f},Obj.DeviceName) %If there is a match assign it:
                                       Obj.VideoIndex=DevIDs(f);
                                       Obj.Video=VidObjs(g);
                                       if ~isrunning(Obj.Video)
                                           start(Obj.Video);
                                       end
                                       DevFound=1;
                                       flushdata(Obj.Video);
                                       Obj.waitForLogging;
                                       Obj.IsInitialized=1;
                                   end
                               end
                           end
                       end
                   end
                   
                   if ~DevFound % If the object does not exist then create obj
                       Obj.VideoIndex=0;
                       imin=imaqhwinfo('winvideo');
                       % Search for the right DeviceInfo
                       for f=1:numel(imin.DeviceInfo)
                           if strcmpi(imin.DeviceInfo(f).DeviceName,Obj.DeviceName)
                               Obj.VideoIndex=f;
                           end
                       end
                       
                       if ~Obj.VideoIndex
                           fprintf('No video framegrabber has been detected.\n');
                       else
                           fprintf('Framegrabber detected: %s\n',Obj.DeviceName);
                           % if a video object exists delete it first
                           if ~isempty(Obj.Video)
                               delete(Obj.Video);
                               Obj.Video=[];
                           end
                           sendCmdProgress(handles,0.02);
                           % establish new video capture
                           if ~isempty(Obj.DeviceFormat)
                               Obj.Video=videoinput('winvideo',Obj.VideoIndex,Obj.DeviceFormat);
                           else
                               Obj.Video=videoinput('winvideo',Obj.VideoIndex);
                           end
                           set(Obj.Video,'TriggerRepeat',Inf);
                           set(Obj.Video,'FrameGrabInterval',1);
                           % Wait for a fraction of a second to ensure the capture
                           % frame is not chopped:
                           set(Obj.Video,'TimerPeriod',0.1);
                           set(Obj.Video,'FramesPerTrigger',1);
                           set(Obj.Video,'TriggerRepeat',0);
                           set(Obj.Video,'Timeout',2.0);
                           %                                set(Obj.Video,'Timeout',3);
                           triggerconfig(Obj.Video, 'manual');
                           set(Obj.Video,'ReturnedColorSpace','grayscale');
                           start(Obj.Video);
                           Obj.IsInitialized=1;
                           flushdata(Obj.Video);
                           Obj.waitForLogging;
                       end
                   end
                   display('Completed framegrabber initializing ');
               end
                              
           catch
               display('... Error is connecting to the video framegrabber..');
               Obj.IsInitialized=0;
           end
        end

        function release(Obj)
            try 
                stop(Obj.Video);
                Obj.IsInitialized=0;
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.release()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end    
        end            
            
        function [thefile,ErrorCode,ErrorMessage]=framecap(Obj,MeantForCalibration,SaveToFile,CarmPhysicalAngle, handles)
            
            ErrorCode=1;
            ErrorMessage='';
            StartTime=cputime;
            try
                thefile='';    
                UseAvailableCarmPhysicalAngle=0;
                if exist('CarmPhysicalAngle','var')
                    if ~isempty(CarmPhysicalAngle)                            
                        UseAvailableCarmPhysicalAngle=1;
                    else
                        UseAvailableCarmPhysicalAngle=0;
                    end
                end                
                if SaveToFile
                    if ~exist('TCarmCalib','dir')
                       mkdir('TCarmCalib');
                    end
                    if ~exist('TCarmCalib\LastShot Data','dir')
                       mkdir('TCarmCalib\LastShot Data');
                    end
                    if ~exist('TCarmCalib\LastShot Data\Image','dir')
                        mkdir('TCarmCalib\LastShot Data\Image');
                    end
                    if ~exist('TCarmCalib\LastShot Data\OrigImage','dir')
                        mkdir('TCarmCalib\LastShot Data\OrigImage');
                    end
                    path='TCarmCalib\LastShot Data\Image';
					OrigImageFile='TCarmCalib\LastShot Data\OrigImage\OrigImage.bmp';
                    delcmd=sprintf('del /Q "%s\\*.*',path);
					dos(delcmd);
					% delete(sprintf('%s\\*.*',path));
                    filename = [path '\' Obj.ImageName];
                    Obj.ImageCount = Obj.ImageCount + 1;
                    imgcount = num2str(Obj.ImageCount);
                    filename = [filename '_' imgcount];
                    thefile = [filename '.bmp'];
                end
                
                if isempty(MeantForCalibration)
                    MeantForCalibration=0;
                end
                
                if Obj.SimMode 
                    % Using the FrameGrabber Simulator option
                    ErrorCode=6;
                    dfiles=dir([handles.SimExamDataFolder,'SimExam\*.tif']);
                    
                    RandomOrder=1;
                    % Read the order file if available:
                    try
                        OrderFileName=[handles.SimExamDataFolder,'SimExam\order.txt'];
                        if exist(OrderFileName,'file')
                            fid=fopen(OrderFileName);
                            str=fgetl(fid);
                            fclose(fid);
                            order_index=str2num(str);
                            RandomOrder=0;
                            % rewrite the order file after reading by adding one:
                            if order_index<numel(dfiles)
                                delete OrderFileName;
                                fid=fopen(OrderFileName,'w');                            
                                fprintf(fid,'%i',order_index+1);
                                fclose(fid);
                            end
                        end
                    catch
                        RandomOrder=1;
                    end
                                        
                    if ~RandomOrder
                       RndFileName=num2str(order_index);
                        % Add zeros to the beginning of the filename
                        if order_index<10
                            RndFileName=['00',RndFileName,'.tif'];
                        elseif order_index<100
                            RndFileName=['0',RndFileName,'.tif'];
                        end

                        RndFileName=[handles.SimExamDataFolder, 'SimExam\' ,RndFileName];
                        
                        if ~exist(RndFileName,'file')
                            RandomOrder=1; % Go with a random order if file not found
                        end
                    end
                    
                    if RandomOrder
                        RndFileName=[handles.SimExamDataFolder, 'SimExam\', dfiles(randi([1,numel(dfiles)],1)).name];
                    end
                    
                    Obj.CroppedCapture=imread(RndFileName);
                    Obj.OriginalCapture=Obj.CroppedCapture;
                    if numel(size(Obj.CroppedCapture))>2
                        Obj.CroppedCapture=Obj.CroppedCapture(:,:,1);
                    end
                    if ~isequal([size(Obj.CroppedCapture,1),size(Obj.CroppedCapture,2)],Obj.DefaultSize)
                        Obj.CroppedCapture=imresize(Obj.CroppedCapture,'OutputSize',Obj.DefaultSize);
                    end
                    Obj.CroppedCapture=uint8(Obj.CroppedCapture.*Obj.ImageMask); % apply mask                    
                    % Capture Complete
                    ErrorCode=0;                
                else 
                    % FrameGrabber Files Set up
                    ErrorCode=2;

                    [VidAvailable,im, VideoTimeout]=Obj.checkVideo([]);
                    sendCmdProgress(handles, 0.05);
                    
                    if (~VidAvailable || isempty(im) || VideoTimeout) % Try to reinitialize connection
                        Obj.IsInitialized=0;
                        fprintf('\r\nSoft Connection Re-initialization of the Framegrabber..\n');
                        Obj.Initialize(0,handles);
                        [VidAvailable,im, VideoTimeout]=Obj.checkVideo([]);
                        if (~VidAvailable || isempty(im) || VideoTimeout) % Try to reinitialize connection
                            sendCmdProgress(handles, 0.08);
                            fprintf('\r\nHard Connection Re-initialization of the Framegrabber..\n');
                            Obj.Initialize(1,handles);
                            [VidAvailable,im, VideoTimeout]=Obj.checkVideo([]);
                            sendCmdProgress(handles, 0.10);
                        end
                    end
                    
                    % FrameGrabber initialized
                    ErrorCode=3;
                   
                    if ~(~VidAvailable || isempty(im) || VideoTimeout)
                        
                        fprintf('\ngot here!!! ');
                        
                        validframe=1;
                        if isempty(im)
                            fprintf('\nEmpty Frame captured! ');
                            validframe=0;
                        end
                        
                        if (sum(sum(im(:,:,1)))==0)
                            fprintf('\nZero Frame captured! ');
                            validframe=0;
                        end
                        
                        % Only process the rest if the frame received
                        % is valid:
                        if validframe
                            if SaveToFile
                                imwrite(im, thefile,'bmp');
                                imwrite(im, OrigImageFile,'bmp');
                            end
                            
                            % Scale the image if ForceOriginalRes
                            % option is selected:
                            if Obj.ForceOriginalRes
                                im=imresize(im(:,:,1),[Obj.OriginalRes(2),Obj.OriginalRes(1)]);
                                fprintf('\nFrame is forced to rescale. ');
                            end
                            
                            % set the OriginalRes
                            Obj.CaptureRes=[size(im,2),size(im,1)];
                            
                            % Frame Captured
                            ErrorCode=4;
                            
                            % Crop and pad the image as suggested by the framegrabber parameters and apply
                            % mirroring and rotation:
                            new_im = processScaling(im, [1024,1024,Obj.CropSize],Obj.CarmModelImagePars(1),Obj.CarmModelImagePars(2));
                            Obj.OriginalCapture=im;
                            
                            DistortionCorrectionHasIssues=0;
                            %Apply Undistortion if necessary:
                            if (~MeantForCalibration) && (Obj.UnDistortFlag)
                                if (UseAvailableCarmPhysicalAngle)
                                    Obj.initUnDistortObject([handles.RootDataFolder,Obj.Carm.CalibrationFile]);
                                    new_im=Obj.UnDistortObject.undistortImageSimple(new_im,CarmPhysicalAngle);
                                else
                                    if Obj.PlaneIndex==1 % AP:
                                        Obj.initUnDistortObject([handles.RootDataFolder,Obj.Carm.CalibrationFile]);
                                        new_im=Obj.UnDistortObject.undistortImageSimple(new_im,0);
                                    end
                                    if Obj.PlaneIndex==2 % Lateral:
                                        Obj.initUnDistortObject([handles.RootDataFolder,Obj.Carm.CalibrationFile]);
                                        new_im=Obj.UnDistortObject.undistortImageSimple(new_im,-90);
                                    end
                                end
                                if isempty(new_im)
                                    DistortionCorrectionHasIssues=1;
                                    fprintf('\nError applying distortion correction to the captured frame. ');
                                else
                                    fprintf('Distortion correction applied to the captured frame. ');
                                end
                            end
                            
                            if ~(DistortionCorrectionHasIssues)
                                % Undistortion applied
                                ErrorCode=5;
                                
                                Obj.CroppedCapture=uint8(new_im.*Obj.ImageMask);
                                
                                if (SaveToFile)
                                    imwrite(Obj.CroppedCapture, thefile,'bmp');
                                end
                                
                                % Cropping Applied with no error
                                ErrorCode=0;
                            end
                            
                        end
                        
                    end
                end
            catch err
                st=sprintf('<< FrameCapture Error>> \r\n%s',err.message);
                disp(st);
            end
            
            if ErrorCode
                switch ErrorCode
                    case 1
                        ErrorMessage='Error: Settting up the FrameGrabber!';
                    case 2
                        ErrorMessage='Error: Initializing the FrameGrabber!';
                    case 3
                        ErrorMessage='Error: Capturing Image from FrameGrabber!';
                    case 4
                        ErrorMessage='Error: Undistorting Image Acquired from FrameGrabber!';
                    case 5                        
                        ErrorMessage='Error:: Cropping Image Acquired from FrameGrabber!';
                    case 6                        
                        ErrorMessage='Error:: Preparing Simulated FrameGrabber!';
                end
                Obj.IsInitialized=0;
            else
                if ~isempty(Obj.OriginalRes) && (~Obj.SimMode)
                    if ~((Obj.CaptureRes(1)==Obj.OriginalRes(1)) && (Obj.CaptureRes(2)==Obj.OriginalRes(2)))
                        ErrorCode=7;
                        ErrorMessage='Error: FrameGrabber Resolution does not Match the Intended Resolution!';
                        Obj.IsInitialized=0;
                    end
                end
                
            end 
            StopTime=cputime;
            fprintf('(Capture Time: %.3f seconds)\r\n%s',StopTime-StartTime);
            
        end
        
        % loads the hardware specific parameters of the capturecard to the
        % properties of the class:
        function out=readConfigFile(Obj)
            try
                out=0;
                if exist(Obj.CFGFileName,'file')
                    %Pars=load('FrameGrabber-Config.txt','-ascii');
                    CStr = textread(Obj.CFGFileName, '%s', 'delimiter', ';');

%                     Obj.DeviceName=CStr{1}; 
                    formatindex = strfind(CStr{1},'#');
                    Obj.DeviceFormat='';
                    if isempty(formatindex)
                        Obj.DeviceName=CStr{1};                    
                    else
                        Obj.DeviceName=CStr{1}(1:formatindex-1);
                        if (formatindex(end))<length(CStr{1})
                            Obj.DeviceFormat=CStr{1}(formatindex(end)+1:end);
                        end
                    end
                    Pars(1)=0;
                    % Set optional parameters to default
                    Obj.ForceOriginalRes=0;
                    Obj.AutoCheckImageThreshold=0.1;
                    Obj.AutoCheckTimeThreshold=1;
                    Obj.AutoCheckLoops=2;
                    Obj.AutoCheckLastTime=clock;
                    for f=2:numel(CStr)
                        if isempty(CStr{f})
                            Pars(f)=0;
                        else
                            Pars(f)=str2num(CStr{f});
                        end
                    end
                    Obj.VideoIndex=Pars(1);
                    % Obj.DefaultSize=Pars(2:3);
                    Obj.CropSize=Pars(4:7);
                    Obj.CarmModelImagePars=[Pars(8),Pars(9)];
                    Obj.UnDistortFlag=Pars(10);
                    if numel(Pars)>11
                       Obj.OriginalRes=[Pars(11),Pars(12)];
                    end
                    if numel(Pars)>12                       
                        if (~isempty(Pars(13)))
                            if (Pars(13)~=0)
                               Obj.ForceOriginalRes=1;
                            end
                        end
                    end
                    if numel(Pars)>13
                       if (Pars(14)>0) && (Pars(14)<20)
                           Obj.AutoCheckImageThreshold=Pars(14);
                       end
                    end
                    if numel(Pars)>14
                       if (Pars(15)>0) && (Pars(15)<20)
                           Obj.AutoCheckTimeThreshold=Pars(15);
                       end
                    end     
                    if numel(Pars)>15
                       if (Pars(16)>1) && (Pars(16)<10)
                           Obj.AutoCheckLoops=round(Pars(16));
                       end
                    end                                             

                    out=1;
                end
                % Read the Maskfile and apply it to the ImageMask;
                if exist(Obj.MSKFileName,'file')
                    im=imread(Obj.MSKFileName,'png');
                    Obj.ImageMask=uint8(im(:,:,1)>1);
                else
                    im=ones(Obj.DefaultSize(1),Obj.DefaultSize(2));
                    Obj.ImageMask=uint8(im);
                end
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.readConfigFile()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end    
        end
        
        % Method used to instatiate the Distortion correction class if
        % necessary
        function out=initUnDistortObject(Obj,FldrName)            
            try
                NeedForObject=0;
                if ~isempty(Obj.UnDistortObject)
                    if ~isvalid(Obj.UnDistortObject)
                        NeedForObject=1;
                    end
                else
                    NeedForObject=1;
                end
                if NeedForObject
                    Obj.UnDistortObject=Undistortion(FldrName);
                end
                Obj.UnDistortFlag=1;
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in FrameGrabber.initUnDistortObject()');
                fprintf('\r\n(SimMode: %i)',Obj.SimMode);
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end
            
        end
        
    end
    
end

function sendCmdProgress(handles, progressPerc)

if ~isempty(handles)
    if ~isempty(handles.tc)
        if strcmpi(handles.tc.status,'open')
            handles.CmdP.CmdProgress=progressPerc;
%              display(['>> Send update!! ', num2str(progressPerc)]);
            clientBytesReceived([],[],handles);
        end        
    end
end

end


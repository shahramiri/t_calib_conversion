function [handles,output_args] = TLoadCalibration(handles,Data, DoNOtWriteUsedCalibReg,FromWithinExam)
% XML Command: Load Calibration and Registration files
% Inputs:
%   Data: (struct) from XML
%   DoNOtWriteUsedCalibReg (Optional) if TLoadCalibration is indirectly
%   used within the TLoadExam_Full and TDuplicateExam

SuccessFlag=0;
ErrorCode=[];
TrackingMarkersLocalCoords=[];
ErrorMessage='';
fprintf('Load Calibration..\n');

CalibFolder=Data.CalibFolder;
RegFile=Data.RegFile;

if ~exist('DoNOtWriteUsedCalibReg','var')
    DoNOtWriteUsedCalibReg=0;
end

if ~exist('FromWithinExam','var')
    FromWithinExam=0;
end

% 1. Calibration> b. load the calibration and marker registration files

% try
    [Out]=handles.W1.Carm.loadCalib(handles.CalibsRootDataFolder,CalibFolder);
    if (Out)       
        % Check and Apply the custom Oblique Angles
        ObliqueAngles=checkObliqueAngles([handles.CalibsRootDataFolder,'Calibration\',handles.W1.Carm.CalibrationFile]); %Check if ObliqueAngles are defined
        if isempty(ObliqueAngles)
            handles.W1.StitchMaster.APObliqueAngle=0;
            handles.W1.StitchMaster.MLObliqueAngle=-90;
        else            
            handles.W1.StitchMaster.APObliqueAngle=ObliqueAngles(1);
            handles.W1.StitchMaster.MLObliqueAngle=ObliqueAngles(2);
        end
        
        %LOAD FRAME-GRABBER INFO:
        FldrName=handles.W1.Carm.CalibrationFile;
        FNM=dir([pwd '\' handles.CalibsRootDataFolder 'Calibration\' FldrName '\*.FGB']);
        FGBFileName=[pwd '\' handles.CalibsRootDataFolder 'Calibration\' FldrName '\' FNM(1).name];
        handles.W1.FrameGrabber.CFGFileName=FGBFileName;
        FNM=dir([pwd '\' handles.CalibsRootDataFolder 'Calibration\' FldrName '\*.FGBMSK']);
        FGBMSKFileName=[pwd '\' handles.CalibsRootDataFolder 'Calibration\' FldrName '\' FNM(1).name];
        if ~exist(FGBMSKFileName,'file')
            FGBMSKFileName='';
        end
        handles.W1.FrameGrabber.MSKFileName=FGBMSKFileName;
        if exist(handles.W1.FrameGrabber.CFGFileName,'file')
            handles.W1.FrameGrabber.readConfigFile;
        end

        handles.W1.Carm.RegistrationFile=RegFile;
        [LROut,TrackingRomFileName,TrackingMarkersLocalCoords]=handles.W1.Carm.loadRegistration(handles.CalibsRootDataFolder, handles.W1.Carm.RegistrationFile); % Load Registration:
        if ~(LROut==1)
            ErrorCode=1;
            ErrorMesssage='Unable to Load the Marker Registration';
        end
        if (LROut==1) % only procees if registration is loaded and there is ROM file to load
            %//---------------------------------------------------------------------
            %/Creating or Updating the .USE file to make sure the calibration shows on the top of the list the next time.
            FldrName=handles.W1.Carm.CalibrationFile;
            delete(sprintf('%s\\%sCalibration\\%s\\*.USE',pwd, handles.CalibsRootDataFolder,FldrName)); % delete any *.use file
            UseFileName=sprintf('%s\\%sCalibration\\%s\\%s.USE',pwd, handles.CalibsRootDataFolder, FldrName, urlencode(FldrName));
            use_fid=fopen(UseFileName,'w');
            now_str=datestr(now);
            time_str=now_str(end-7:end);
            ORDate=date;
            ORDateTime=[ORDate, ' ' , time_str];
            fprintf(use_fid,'%s\r\n',ORDateTime);
            fclose(use_fid);
            %//---------------------------------------------------------------------
            % Now search for .jcal file and extract the description from the .jcalinfo file
            RegUseFile=sprintf('%sCalibration\\%s\\%suse',handles.CalibsRootDataFolder,FldrName,handles.W1.Carm.RegistrationFile);
            delete(RegUseFile); % delete reg use
            reg_use_fid=fopen(RegUseFile,'w');
            now_str=datestr(now);
            time_str=now_str(end-7:end);
            ORDateTime=[ORDate, ' ' , time_str];
            fprintf(reg_use_fid,'%s\r\n',ORDateTime);
            fclose(reg_use_fid);
            %//---------------------------------------------------------------------
            % Save the Calibration and Registration files into the exam folder at the time of calibration:
            if ~DoNOtWriteUsedCalibReg
                saveCalRegFileName(handles.exam_name,handles.CalibsRootDataFolder,CalibFolder,RegFile,[],[],[]);
            end
            fprintf(['Calibration: ',strrep(strrep(CalibFolder,'\','/'),'%','_'),' loaded.\n']);
            fprintf(['Registration: ',strrep(strrep(RegFile,'\','/'),'%','_'),' loaded.\n']);
            %//---------------------------------------------------------------------
            % Preparing the attachment for the client to pick-up if the
            % TrackingRomFileName is found:
            if ~isempty(TrackingRomFileName)
                temp_name = tempname;
                [~, FileName_1] = fileparts(temp_name);
                FullFileName_1=['SendingAttachments\',FileName_1,'.dat'];
                copyfile(TrackingRomFileName,FullFileName_1);
                TrackingRomFileBytesSize=checkFileBytes(FullFileName_1);
                fid=fopen(FullFileName_1);
                TrackingRomFileBytesData=fread(fid);
                fclose(fid);
                fprintf('RomFile ready for pickup.\n');
            else
                TrackingRomFileBytesSize=0;
                TrackingRomFileBytesData=[];
                FileName_1='';
                FullFileName_1='';
                fprintf('Warning: No matching tracker markers found.\n');
            end
            
            if (~handles.W1.FrameGrabber.SimMode)
                % try to connect to the framegrabber and take the first image
                try
                    FNM=dir([pwd '\' handles.CalibsRootDataFolder 'Calibration\' handles.W1.Carm.CalibrationFile '\*.FGB']);
                    FGBFileName=[pwd '\' handles.CalibsRootDataFolder 'Calibration\' handles.W1.Carm.CalibrationFile '\' FNM(1).name];
                    if exist(FGBFileName, 'file')
                        handles.W1.FrameGrabber.CFGFileName=FGBFileName;
                    else
                        handles.W1.FrameGrabber.CFGFileName='';
                    end
                    if (~FromWithinExam)
                        sendCmdProgress(handles, 0.01);
                        fprintf('\nInitializing Connection to the framegrabber..');
                        handles.W1.FrameGrabber.Initialize(1,handles);
                        fprintf('\nFramegrabber Initialized. Collecting a sample image..');
                        sendCmdProgress(handles, 0.5);
                        [FG_thefile,FG_ErrorCode,FG_ErrorMessage]=handles.W1.FrameGrabber.framecap(0,0, [],[]);
                        fprintf('\nSample image successfully collected.');
                        sendCmdProgress(handles, 0.75);
                    else
                        FG_ErrorCode=0;
                    end
                    if FG_ErrorCode
                        fprintf('Frame-grabber Error!\r\n');
                        fprintf('Frame-Grabber Error Code: %i; Frame-Grabber Error Message :%s\r\n',FG_ErrorCode,FG_ErrorMessage);
                    else
                        handles.W1.FrameGrabber.OriginalCapture=[];
                        handles.W1.FrameGrabber.CroppedCapture=[];
                    end
                    SuccessFlag=1;
                catch
                    fprintf('Unable to connect to the frame-grabber!');
                    if (~FromWithinExam)
                        sendCmdProgress(handles, 0.75);
                    end
                end
            else
                SuccessFlag=1;
            end
                
        else
            fprintf('Error: Loading registration Failed. Registration file: %s\n',RegFile);
        end
        
    end
% catch
%     if isempty(CalibFolder)
%         fprintf('Error finding the calibration folder: %s\n',CalibFolder);
%     else
%         fprintf('Error loading calibration and registration: %s\n',RegFile);
%     end
% end

if (~FromWithinExam)
    sendCmdProgress(handles, 0.99);
end

if SuccessFlag    
    if ~isempty(FileName_1)
        if handles.SaveOnDiskMode
            handles.AttachmentFileName=FileName_1;
        else
            handles.AttachmentFileName=FileName_1;
            handles.Memory.MemoryFileName=FileName_1; %SavedThumbNailFileName;
            handles.Memory.Files.(FileName_1)=TrackingRomFileBytesData;
        end
    end
    output_args=[];
    output_args.Success=1;
    if ~isempty(FileName_1)
        output_args.HasAttachment=1;
        output_args.AttachmentType='TrackingRomFile';
        output_args.AttachmentFileName=handles.AttachmentFileName;
        output_args.ByteSize=TrackingRomFileBytesSize;
        output_args.MarkersLocalCoords=TrackingMarkersLocalCoords;
    else
        output_args.HasAttachment=0;
        output_args.AttachmentType='';
        output_args.AttachmentFileName='';
        output_args.ByteSize=0;
        output_args.MarkersLocalCoords=[];
        output_args.ErrorCode=1;
        output_args.ErrorMessage='Unable to find and load the Marker .ROM file.';         
    end 
else
    output_args=[];
    output_args.Success=1;
    if isempty(ErrorCode)
        ErrorCode=2;
        ErrorMessage='Unable to load the calibration & registration file.';
    end
    output_args.ErrorCode=ErrorCode;
    output_args.ErrorMessage=ErrorMessage;
    output_args.HasAttachment=0;
    output_args.AttachmentType='TrackingRomFile';    
    output_args.AttachmentFileName='';
    output_args.ByteSize=0;    
end

end 

function sendCmdProgress(handles, progressPerc)

    if ~isempty(handles.tc)
        if strcmpi(handles.tc.status,'open')
            handles.CmdP.CmdProgress=progressPerc;
            clientBytesReceived([],[],handles);
        end
    end

end
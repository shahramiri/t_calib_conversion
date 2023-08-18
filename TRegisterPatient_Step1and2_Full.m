function [handles,output_args] = TRegisterPatient_Step1and2_Full(handles,Data)
% Register Patient Anatomy

fprintf('Registering Grid - Step 1 & 2..\n');

% check if there is ORSetupOverWrite for debugging purposes
if exist('Assets\ORSetupOverWrite.txt','file')
    display('ORSetupOverwrite detected: Assets\ORSetupOverWrite.txt');
    try
        fileID=fopen('Assets\ORSetupOverWrite.txt');
        InputSetup=fgets(fileID);
        fclose(fileID);
        if (numel(InputSetup)==10) || (numel(InputSetup)==8)
            Data.ORSetupString=InputSetup;
            handles.W1.StitchMaster.ORSetupString=InputSetup;
            display(['ORSetup replaced to: ', InputSetup, ' the file: Assets\ORSetupOverWrite.txt']); 
        else
            display('Error in analyzing text in file: Assets\ORSetupOverWrite.txt');    
        end
        
    catch
        display('Error reading the file: Assets\ORSetupOverWrite.txt');
    end
end

% handles.W1.StitchMaster.analyzeORSetup

switch upper(Data.ORSetupString(6))
    case {'S','P'}
        IsStandardViews=1;
    case {'L','R'}
        IsStandardViews=0;
end

% CurrentPlane=handles.W1.StitchMaster.CurrentPlane;

% TCarm Code to Update the DAQ based on the new Tracker info:

% 1. DAQ> a. Update C-arm tracking based on raw TrackerData

% Add tracker Data to the DAQ so they will be processed with the SwitchPlane happens:
% TrackData=Data.TrackerData_AP;
% handles.W1.DAQ.UDPList.Items(1).feed(str2num(TrackData));

if IsStandardViews
%     Data.PlaneID='CARMAP';
    PlaneID='CARMAP';
else
%     Data.PlaneID='CARMML';
    PlaneID='CARMML';
end
% TSwitchPlane(handles,Data);
handles.W1.Carm.switchPlane(PlaneID);
handles.W1.Carm.lookupCalib(Data.TrackerData_AP);
IPInfo=handles.W1.Carm.IPInfo;
handles.W1.StitchMaster.APIPinfo=IPInfo;

% Add tracker Data to the DAQ so they will be processed with the SwitchPlane happens:
% TrackData=Data.TrackerData_ML;
% handles.W1.DAQ.UDPList.Items(1).feed(str2num(TrackData));

% Receive to camera positions and register the patient:
if IsStandardViews
    PlaneID='CARMML';
else
    PlaneID='CARMAP';
end

% TSwitchPlane(handles,Data);
handles.W1.Carm.switchPlane(PlaneID);
handles.W1.Carm.lookupCalib(Data.TrackerData_ML);
IPInfo=handles.W1.Carm.IPInfo;
handles.W1.StitchMaster.MLIPinfo=IPInfo;

Iap=handles.W1.StitchMaster.APIPinfo;
Iml=handles.W1.StitchMaster.MLIPinfo;

% set the currentplaneback
if IsStandardViews
    PlaneID='ML';
else
    PlaneID='AP';
end
% TSwitchPlane(handles,Data);
handles.W1.Carm.switchPlane(PlaneID);
handles.W1.Carm.lookupCalib(Data.TrackerData_ML);

SkipReg=0;
if isfield(Data,'SkipRegistration')
    if str2num(Data.SkipRegistration);
        SkipReg=1;
    end
end

% 2. StitchMaster> k. create grid based on the given beam orientation (Patient Registration)

if (SkipReg)
    % recaluclate the grid based on existing grid
    [Old_APSign,Old_MLSign,Old_SISign]=handles.W1.StitchMaster.analyzeORSetup;
    if isfield(Data,'ORSetupString')
        handles.W1.StitchMaster.ORSetupString=Data.ORSetupString;
    end
    [APSign,MLSign,SISign]=handles.W1.StitchMaster.analyzeORSetup;
    Org=handles.W1.StitchMaster.Grid.VxLineElement.Coords;
    Vx=handles.W1.StitchMaster.Grid.VxLineElement.Vector * Old_APSign * APSign;
    Vy=handles.W1.StitchMaster.Grid.VyLineElement.Vector * Old_MLSign * MLSign;
    handles.W1.StitchMaster.rebuildGrid(Org,Vx,Vy); 
else
    % Check and see if the ORSetupString is a part of input variables and if so assign it:
    if isfield(Data,'ORSetupString')
        handles.W1.StitchMaster.ORSetupString=Data.ORSetupString;
    end
    handles.W1.StitchMaster.initiateGrid(0,[],[],[],turnIPinfo(Iap),turnIPinfo(Iml));
end
    
fprintf('Registering Grid - Step 1 &2 Completed.\n');

% 3. Exam> f. write last exit mode

handles.ReviewOnlyMode=0;
sysWriteExitStatus(handles.LastExitMode,handles.ReviewOnlyMode,handles.TestData.FolderName);

% 4. StitchMaster> n. get Grid And Planes coordinates

% In addition return the plane and grid info:
[handles,output_args] = TGetGridAndPlanesTrans(handles,Data,1);
output_args.ReviewOnlyMode=handles.ReviewOnlyMode;

% Check and see if the ORSetupString is a part of input variables and if so assign it:
saveCalRegFileName(handles.TestData.FolderName,handles.CalibsRootDataFolder,[],[],Data.ORSetupString,Data.TrackerData_AP,Data.TrackerData_ML)

fprintf('Request Completed: Tracker transformations Sent and other parameters sent to the client.\n');
end


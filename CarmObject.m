classdef CarmObject < handle
    % Holds all the main C-arm Calibration and Tracking Data and
    % Calculations detailed explanation goes here
    
    properties
        World
        CGantry         % Class holding the C gantry X-bram lookup table
        T_Imspace_Ref   % 4x4 transformation matrix; result of Joint System Calibration
        T_Ref_Marker    % 4x4 transformation matrix; result of marker registration
        T_Marker_Camera % 4x4 transformation matrix; marker tracking
        IPInfo          % Structure containing the calculated beam info 
        CalibrationFile
        RegistrationFile
        Orbit
        Tilt
    end
    
    methods
        function Obj=CarmObject()           
            Obj.T_Ref_Marker=[];
            Obj.T_Imspace_Ref=[];
            Obj.T_Marker_Camera=[];
            Obj.CGantry=CGantryObject(['CGantry','-CGantry']);            
            Obj.T_Ref_Marker=eye(4);
            Obj.T_Imspace_Ref=eye(4);
            Obj.T_Marker_Camera=eye(4);
            Obj.Orbit=0;
            Obj.Tilt=0;

        end
        
        % Loading the Calibration (Gantry & Joint System Calibration)
        function [Out]=loadCalib(Obj,CalibsRootDataFolder,FolderName)
            Out=0;
            Out1=loadCGantryCalib(Obj,CalibsRootDataFolder,FolderName);
            Out2=loadJointSysCalib(Obj,CalibsRootDataFolder,FolderName);
            if (Out1 && Out2)
                Out=1;
                Obj.CalibrationFile=FolderName;
            end
        end
        
        % Loading the T_Ref_Marker from the marker registration file
        function [Out,TrackingRomFileName,MarkersLocalCoords]=loadRegistration(Obj,CalibsRootDataFolder,FileName)
            Out=0;
            TrackingRomFileName='';
            MarkersLocalCoords=[];
            try
                file=[pwd '\' CalibsRootDataFolder 'Calibration\' Obj.CalibrationFile '\' FileName];
                if path~=0
                    load(file,'MarkerTransform','-mat'); % editted on 2018-Apr-27
                    if exist('MarkerTransform','var') % assign the marker registration if it exists                        
                        Obj.T_Ref_Marker=MarkerTransform;
                    end
                end
                % Finding the corresponding Marker ROM file to be send to the client; reading the info file:
                try
                    MFileName=[file,'info'];
                    fid=fopen(MFileName,'r');
                    MarkerSet=fgetl(fid);
                    RegQ=fgetl(fid);
                    RegQuality=str2num(RegQ);
                    if feof(fid)
                        Mount='';
                    else
                        Mount=fgetl(fid);
                    end
                    fclose(fid);
                    % looking for the file in the Markers folder
                    RomFiles=dir([pwd,'\', CalibsRootDataFolder, 'Calibration\',Obj.CalibrationFile,'\Markers\',Mount,'*.rom']);
                    
                    if ~isempty(RomFiles)
                        TrackingRomFileName=[CalibsRootDataFolder, 'Calibration\',Obj.CalibrationFile,'\Markers\',RomFiles(1).name];
                        TrackingRomExtFileName=[CalibsRootDataFolder, 'Calibration\',Obj.CalibrationFile,'\Markers\',RomFiles(1).name,'ext'];                    
                        
                        Out=1; % Function returns successful output if ROM file exists
                        % Also look for possible extension file: .romext and return the marker coordinates if applicable
                        if exist(TrackingRomExtFileName,'file')
                            % read the file contents:
                            try
                                MrkCoordsRead=dlmread(TrackingRomExtFileName);
                                if numel(size(MrkCoordsRead))==2
                                    if (size(MrkCoordsRead,1)>3) && (size(MrkCoordsRead,2)==3)
                                        MarkersLocalCoords=MrkCoordsRead;
                                    end
                                end
                            catch
                                MarkersLocalCoords=[];
                            end
                        end
                    end
                catch
                    TrackingRomFileName='';                
                end
            catch
               Out=0;
               TrackingRomFileName='';
            end
        end
        
        % Main function for calculating the IPinfo based on all calibration
        % and tracking data:
        function IPInfo=lookupCalib(Obj,TrackingMatrix)
           
           IPInfo=[];
           try  
               res=Obj.applyTracking(str2num(TrackingMatrix));

               R=Obj.CGantry.simpleLookup(Obj.Orbit);
               % T=Obj.JointSys.getTmat(Obj.T_Ref_Marker); %older version
               T=Obj.T_Marker_Camera*Obj.T_Ref_Marker*Obj.T_Imspace_Ref;
               IPInfo=R;
               %....................................................................
               p=T*[R.SourcePosition,1]';
               IPInfo.SourcePosition=p(1:3)';
               
               p=T*[R.CentrePoint,1]';
               IPInfo.CentrePoint=p(1:3)';
               
               p1=T*[R.SourcePosition,1]';
               p1=p1(1:3)';
               p2=T*[R.SourcePosition+(100*R.ImageNormal),1]';
               p2=p2(1:3)';
               IPInfo.ImageNormal=(p2-p1)/norm(p2-p1);
               
               %            p1=T*[R.ImageUp,1]';
               p1=T*[R.SourcePosition,1]';
               p1=p1(1:3)';
               p2=T*[R.SourcePosition+(100*R.ImageUp),1]';
               p2=p2(1:3)';
               IPInfo.ImageUp=(p2-p1)/norm(p2-p1);
               
               % make sure the vectors are perfectly orthogonal:
               if abs(acosd(dot(IPInfo.ImageUp,IPInfo.ImageNormal))-90)>0.0001
                   teV=cross(IPInfo.ImageNormal,IPInfo.ImageUp);
                   IPInfo.ImageUp=cross(teV,IPInfo.ImageNormal);
                   IPInfo.ImageUp=IPInfo.ImageUp/norm(IPInfo.ImageUp);
               end
               Obj.IPInfo=IPInfo;
 
           catch
              disp('Error calculating the IPInfo data!');
              IPInfo=[];
           end

        end        
                     
        % Register markers (T_Ref_Marker) based on the input tracking file
        function [Tcor_ref_mrk, Err, Success]=register(Obj,RecCamFile,SaveTmpRegMarkerData)
            %set the DAQ variable temporary, and remove it to avoid trouble in saving and loading
            Err=[];
            PlotIt=0;
            [Tcor_ref_mrk, Err, Success]=registerMarkerSet(RecCamFile, PlotIt, SaveTmpRegMarkerData);
            Obj.T_Ref_Marker=Tcor_ref_mrk;
        end
        
        % Save the result of marker registration (T_Ref_Marker) in the
        % calibration folder
        function Out=saveRegistration(Obj,CalibsRootDataFolder,CalibFolder,RegFileName) %Saving the registration file into the desired folder:
            Out=0;
            try
                file=[pwd '\' CalibsRootDataFolder 'Calibration\' CalibFolder '\' RegFileName];
                MarkerTransform=Obj.T_Ref_Marker;
                save(file,'MarkerTransform','-mat'); % editted on 2018-Apr-27
            catch
                disp('Error Saving Marker Registration File.');
            end            
        end
        
        
        function [output_args]=switchPlane(Obj,PlaneID)

            output_args=0;
            fprintf('Switch Plane..\n');
            SM=Obj.World.StitchMaster;
            
            switch upper(PlaneID)
                case {'1','AP','APOB'}
                    SelCarmView=SM.APCarmViewStr;
                    SelPlaneIndex=1;
                case {'2','ML','MLOB'}
                    SelCarmView=SM.MLCarmViewStr;
                    SelPlaneIndex=2;
                case 'CARMAP'
                    SelCarmView='AP';
                    if SM.IsStandardViews
                        SelPlaneIndex=1;
                    else
                        SelPlaneIndex=2;
                    end
                case 'CARMML'
                    SelCarmView='ML';
                    if SM.IsStandardViews
                        SelPlaneIndex=2;
                    else
                        SelPlaneIndex=1;
                    end
            end
            
            switch upper(SelCarmView)
                case 'AP'
                    Obj.World.FrameGrabber.PlaneIndex=SelPlaneIndex;
                    Obj.World.Carm.Orbit=0;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Rotation').ConstantValue=0;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Tilt').ConstantValue=0;
                case 'APOB'
                    Obj.World.FrameGrabber.PlaneIndex=SelPlaneIndex;
                    Obj.World.Carm.Orbit=SM.APObliqueAngle;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Rotation').ConstantValue=SM.APObliqueAngle;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Tilt').ConstantValue=0;
                case 'ML'
                    Obj.World.FrameGrabber.PlaneIndex=SelPlaneIndex;
                    Obj.World.Carm.Orbit=-90;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Rotation').ConstantValue=-90;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Tilt').ConstantValue=0;
                case 'MLOB'
                    Obj.World.FrameGrabber.PlaneIndex=SelPlaneIndex;
                    Obj.World.Carm.Orbit=SM.MLObliqueAngle;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Rotation').ConstantValue=SM.MLObliqueAngle;
                    %handles.W1.DAQ.ProcessList.Items(1).InputUnitList.getByName('Tilt').ConstantValue=0;
            end
            
            if SelPlaneIndex==1
                Obj.World.StitchMaster.CurrentPlane='AP';
            else
                Obj.World.StitchMaster.CurrentPlane='ML';
            end
            
            %handles.W1.DAQ.refreshProcessList;
            %IPinfo=handles.W1.DAQ.IPinfo;
            
            fprintf(['Plane Switched to: ',Obj.World.StitchMaster.CurrentPlane,'\n']);
            output_args=1;
        end
        
        
        %% Lower-level Methods
        % Load C-Gantry calibration data
        function [Out]=loadCGantryCalib(Obj,CalibsRootDataFolder,FolderName)
            Out=0;
            FolderName=[CalibsRootDataFolder,'Calibration\',FolderName];            
                        
            %Find the *.ccal files:
            C=dir([FolderName,'\*.ccal']);
            if size(C,1)>0 
                CFileName=[FolderName,'\',C(1).name];
                Obj.CGantry.load(CFileName);
                disp(['Gantry Calibration File loaded: ',CFileName]);
                Out=1;
            else
                disp('No *.ccal file found!');
            end            
        end
        % Load JointSystem calibration data
        function [Out]=loadJointSysCalib(Obj,CalibsRootDataFolder,FolderName)
            Out=0;
            FolderName=[CalibsRootDataFolder,'Calibration\',FolderName];
            %search for the.jscal file
            JS=dir([FolderName,'\*.jscal']);
            if size(JS,1)==0 % If the file was not found th4n 
                % look for the old *.jcal file is *.jscal does not exist
                J=dir([FolderName,'\*.jcal']);
                if size(J,1)>0
                    JFileName=[FolderName,'\',J(1).name];
                    load(JFileName,'JSystem','-mat');
                    disp(['Old Joint Calibration File loaded: ',JFileName]);
                    T_Ref_Imspace=JSystem.JointList(1).Items.Trans;
                    JScalFilename=[FolderName,'\',[J(1).name(1:strfind(J(1).name,'.')-1),'.jscal']];
                    save(JScalFilename,'T_Ref_Imspace','-mat');
                    disp(['Joint Calibration File saved in the new format: ',JFileName]);
                    Obj.T_Imspace_Ref=inv(T_Ref_Imspace);
                    Out=1;
                else
                    disp('No .jscal or old .jcal file found in the folder!');
                    Out=0;
                end
            else
                JScalFilename=[FolderName,'\',JS(1).name];
                load(JScalFilename,'T_Ref_Imspace','-mat');
                Obj.T_Imspace_Ref=inv(T_Ref_Imspace);
                Out=1;
            end
        end   
        % Apply tracking data
        function [Out]=applyTracking(Obj,TrackingMatrix)
            Out=0;
            Obj.T_Marker_Camera=[];
            if ~isempty(TrackingMatrix)
                %convert from quaternion to the matrix form:
                TrackingMatrix(1:3)=TrackingMatrix(1:3)*1000;
                matrix=x2t([TrackingMatrix(1:3),1,TrackingMatrix(4:7)]','qua');
                Obj.T_Marker_Camera=matrix;
                Out=1;
            end
        end
        
        
    end
   
end


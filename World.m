classdef World < handle
    %BODY class to define bodies such as instruments, jig, patient, C-arm,
    %etc
    %   Detailed explanation goes here
    
    properties
        Name
        CarmPoseList
        CarmPatchList
        CurrentPoseIndex
        LogObject
        FrameGrabber
        CalibInProcess        
        CarmPoseListCount
        
        CalibFileName   
        StitchMaster    
        DebugInfo
        Users
        ReportFormat    % Allowable Parameters TIF, DCM; If DCM not is use, then TIF is going used.
        CorruptTag      % Tag to be used for Incomplete (corrupt) stages.
        
        Carm            % Holds the main calibration and tracking calculations
        Handles         % Handle to the parent HandleClass.
    end
    
    properties (Dependent)
    end
    
    methods
        
        function Obj=World(WName)
            Obj.Name=WName;
            Obj.CarmPoseList=List('Saved C-arm Pose List');
            Obj.CurrentPoseIndex=0;
            Obj.CarmPoseListCount=0;
            Obj.CalibFileName = 'MPJT_Calib_PhantomObj_Frame';
            Obj.Users=UserGroup;
            Obj.ReportFormat='TIF'; %'DCM': alternative format
            Obj.CorruptTag='[Incomplete]';
            Obj.Carm=CarmObject; 
        end
        
        function addDisplay(Obj,D)
            D.World=Obj;
            Obj.DisplayList.addItem(D);
        end
        
        function addDisplayMeasure(Obj,D)
            D.World=Obj;
            Obj.DisplayMeasureList.addItem(D);
        end
        
        function addDisplayActivity(Obj,D)
            D.World=Obj;
            Obj.DisplayActivityList.addItem(D);
        end
        
        function addDisplayChannel(Obj,D)
            D.World=Obj;
            Obj.DisplayChannelList.addItem(D);
        end
        
        function OutList=getAllObjects(Obj)
            OList=[];
            % Search through all the branches for the children and grandchildren:
            if Obj.ObjectList.Count>0
                for f=1:Obj.ObjectList.Count
                    OList=[OList,Obj.ObjectList.Items(f)];
                    OList=[OList,Obj.ObjectList.Items(f).getDecendants];
                end
            end
            %Put the results in a list:
            OutList=List('All Objects in the World');
            for f=1:size(OList,2)
                OutList.addItem(OList(f));
            end
        end
        
        function r=takeSnapShot(Obj)
            c=Obj.CarmPoseList.Count;
            st=sprintf('SnapShot-%i',c+1);
            CP=CarmPose(st);
            CP.World=Obj;
            CP.Tilt=0;
            CP.Rotation=Obj.Carm.Orbit;
            CP.First6DOF=Obj.Carm.T_Marker_Camera;
            Obj.CurrentPoseIndex=Obj.CurrentPoseIndex+1;
            CP.PoseIndex=Obj.CurrentPoseIndex;
            Obj.CarmPoseList.addItem(CP);
            r=[];
        end
        
        function cs=getSnapShotString(Obj)
            Cnt=Obj.CarmPoseList.Count;
            CPI=0;
            PCnt=0;
            if (Cnt)
                for f=1:Cnt
                    % add the Patch description to the list:
                    if not(Obj.CarmPoseList.Items(f).PatchIndex==CPI)
                        CPI=Obj.CarmPoseList.Items(f).PatchIndex;
                        cs{f+PCnt}=Obj.CarmPatchList.Items(CPI).getString;
                        PCnt=PCnt+1;
                    end
                    cs{f+PCnt}=sprintf('%i: %s',f,Obj.CarmPoseList.Items(f).getString);
                end
            end
        end
        
        function saveCarmPoseList(Obj,FileName)
            CarmPoseList=Obj.CarmPoseList;
            CarmPatchList=Obj.CarmPatchList;
            save(FileName,'CarmPoseList','CarmPatchList');
        end
        
        function addStatusItem(Obj,S)
            S.World = Obj;
            Obj.StatusItemList.addItem(S);
        end
        
        
        function [ImgPlane,ErrorCode,ErrorMessage]=collectPosition(Obj,FG_thefile)
            
            ImgPlane=[];
            ErrorCode=0;
            ErrorMessage='';
                           
            disp('Collect>>');            
            path='TCarmCalib\LastShot Data\';                        
            
            % Remove snapshot
            Obj.takeSnapShot;
            
            % Revise the Cnt based on the CarmPoseListCount
            if isempty(Obj.Handles.TestData.ImagePlanes)
                LastCount=0;
            else
                LastCount=numel(Obj.Handles.TestData.ImagePlanes);
            end
            
            Obj.CarmPoseListCount=LastCount+1;
            Cnt=Obj.CarmPoseListCount;       
            
            % Create JointTrack Data .tif file
            JT_path = [path 'JointTrack Data'];
            im=Obj.FrameGrabber.CroppedCapture;

            % Create JointTrack Data .mat file
            ImgPlane = ImagePlane(['IP_' num2str(Cnt)]);
            ImgPlane.readIPinfo(Obj.Carm.IPInfo);

            ImgPlane.Image=im;
            % Storing the original captured file as well:
            ImgPlane.OriginalImage=Obj.FrameGrabber.OriginalCapture;
            ImgPlane.CarmPose = Obj.CarmPoseList.Items(Obj.CarmPoseList.Count);
            % ------------------------------------------------------
            ImgPlane.CarmPose.PoseIndex=Obj.CarmPoseList.Count;
            % ------------------------------------------------------
            % Add Grid Info to the stored imageplane:
            GOrg=Obj.StitchMaster.Grid.PlaneOrigin;
            GVx=Obj.StitchMaster.Grid.VxLineElement.Vector;
            GVy=Obj.StitchMaster.Grid.VyLineElement.Vector;
            GVz=Obj.StitchMaster.Grid.VzLineElement.Vector;

            ImgPlane.RefGridOrg=GOrg;
            ImgPlane.RefGridVx=GVx;
            ImgPlane.RefGridVy=GVy;
            ImgPlane.RefGridVz=GVz;
            %---------------------------------------------------------

            filename = [Obj.CalibFileName '_' num2str(Cnt) '.mat'];
            calib_matfile = fullfile(JT_path,filename);
            % Do not save the World properties it creates issues:
            ImgPlane.CarmPose.World=[];
            save(calib_matfile,'ImgPlane','-mat');

            %Turn off the calibration in-process:
            Obj.CalibInProcess=0;
            
        end
        
        % this is used to update the stage indexes
        function updateStages(Obj)            
            % the indexes are native to the application and hence should be
            % obtained from the application
        end
        
    end
    
end
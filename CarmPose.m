classdef CarmPose < handle
    %BODY class to define bodies such as instruments, jig, patient, C-arm,
    %etc
    %   Detailed explanation goes here
    
    properties
        Name
        World
        X
        Y
        Z
        WigWag
        Tilt
        Rotation
        First6DOF
        Display % list of all display values
        Process % list of all processunits
        PoseIndex
        PatchIndex
    end
    
    properties (Dependent)
        String        
    end
    
    methods
        function Obj=CarmPose(Name)%,X,Y,Rotation,Tilt,WigWag)
            Obj.Name=Name;
        end
        
        function takeSnapShot(Obj) %store snapshot of the current World parameters
%             for str = {'X','Y','Z','WigWag','Tilt','Rotation'}
%              measureName = str{1};
%               - Store all processunit data instead of only selected ones
%               as in the older code
                for f=1:Obj.World.DAQ.ProcessList.Count
                    for g=1:Obj.World.DAQ.ProcessList.Items(f).ProcessUnitList.Count
                        % Save uncalibrated values (directly from ProcessUnits)
                        if ~isempty(Obj.World.DAQ.ProcessList.Items(f).ProcessUnitList.Items(g))
                            measureName=Obj.World.DAQ.ProcessList.Items(f).ProcessUnitList.Items(g).Name;
                            Obj.Process.(measureName)=Obj.World.DAQ.ProcessList.Items(f).ProcessUnitList.Items(g).Output;
                        end
                    end
                    
                    for g=1:Obj.World.DAQ.ProcessList.Items(f).InputUnitList.Count
                        % Save uncalibrated values (directly from ProcessUnits)
                        if ~isempty(Obj.World.DAQ.ProcessList.Items(f).InputUnitList.Items(g))
                            measureName=Obj.World.DAQ.ProcessList.Items(f).InputUnitList.Items(g).Name;
                            Obj.Process.(measureName)=Obj.World.DAQ.ProcessList.Items(f).InputUnitList.Items(g).Output;
                        end
                    end
                    
                    for g=1:Obj.World.DAQ.ProcessList.Items(f).FuseUnitList.Count
                        % Save uncalibrated values (directly from ProcessUnits)
                        if ~isempty(Obj.World.DAQ.ProcessList.Items(f).FuseUnitList.Items(g))
                            measureName=Obj.World.DAQ.ProcessList.Items(f).FuseUnitList.Items(g).Name;
                            Obj.Process.(measureName)=Obj.World.DAQ.ProcessList.Items(f).FuseUnitList.Items(g).Output;
                        end
                    end                    
                    
                end
                
                for f=1:Obj.World.DisplayMeasureList.Count
                    % Save display values (calculated from Sensor Tmat)
                    if ~isempty(Obj.World.DisplayMeasureList.Items(f))
                        measureName=Obj.World.DisplayMeasureList.Items(f).Name;
                        Obj.Display.(measureName)=Obj.World.DisplayMeasureList.Items(f).Value;
                    end
                end
                        
        end
        
        function displayString=getString(Obj)
            % Create string to display in Pose List
            displayString = '';
            measureNames = fields(Obj.Display);
            for i = 1:numel(measureNames)
                measureName = measureNames{i};
                displayString = [displayString measureName(1) sprintf('=%2.2f',Obj.Display.(measureName))];
                % Add a comma to separate measures
                if i~=numel(measureNames)
                    displayString = [displayString ', '];
                end
            end
        end        
        

  end
    
end
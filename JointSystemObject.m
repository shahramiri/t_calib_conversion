classdef JointSystemObject < handle
    %VISOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        JointList    
    end
    
    methods        
        function Obj=JointSystemObject(BName)
            Obj.Name=BName;
            Obj.JointList=List('List of Joints');
        end
        
        function T=getTmat(Obj,varargin)
            T=eye(4);
            for f=Obj.JointList.Count:-1:1
                T=T*Obj.JointList.Items(f).getTmat(varargin{1});
            end
        end  
        
        function readPose(Obj,CCarmPose)
            for f=1:Obj.JointList.Count
                Str=Obj.JointList.Items(f).PoseFieldName;
                Obj.JointList.Items(f).SensorRead=CCarmPose.(Str);
            end
        end
        
        function readPoseFile(Obj,FileName,Index)            
            if nargin==2
                load(FileName,'-mat');
                P1=CarmPoses(end);
                Obj.readPose(P1);
            end
            if nargin==3
                load(FileName,'-mat');
                P1=CarmPoses(Index);
                Obj.readPose(P1);                
            end
        end
        
        function save(Obj,FileName)            
            Name=Obj.Name;
            JointList=Obj.JointList;
            save(FileName,'Name','JointList','-mat');
        end
        
        function load(Obj,FileName)            
            load(FileName,'-mat');
            Obj.Name=Name;
            Obj.JointList=JointList;        
        end
        
        function vals = getCalibratedValues(Obj)
            for f=1:Obj.JointList.Count
                name=Obj.JointList.Items(f).PoseFieldName;
                vals.(name) = Obj.JointList.Items(f).getCalibratedValue;
            end
        end
        
        function T=getCorrectedTmat(Obj,varargin)
            T=eye(4);
            for f=Obj.JointList.Count:-1:1
                T=T*Obj.JointList.Items(f).getCorrectedTmat(varargin{1});
            end
        end  
    end
    
end

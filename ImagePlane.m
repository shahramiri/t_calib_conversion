classdef ImagePlane < handle
    % Holds the 3D description of a projection with full description of the
    % source and detector position and orientation.
    properties
        Name
        JTCalibFileName
        JTImageFileName
        JTPoseFileName
        PrincipalDist
        Xoffset
        Yoffset
        PixelScale
        SourcePosition
        ImageNormal
        ImageUp
        ImageRight
        CentrePoint
        Features
        Path
        Image
        OriginalImage
        CarmPose
        RefGridOrg
        RefGridVx
        RefGridVy
        RefGridVz
        InitMag % Initial Magnification calculated at the time of calibration
        PiercingPoint %Calculated based on the other properties
        PiercingPointLocal
    end
    
    properties (Dependent)
        FeatureCount
        Features3D
    end
    
    methods(Static)
%        function obj = loadobj(s)
%           if ~isstruct(s)
%               newObj=ImagePlane('');
%               sp=properties(s);
%               
%               for f=1:numel(sp)                  
%                   if ~(strcmpi(sp{f},'FeatureCount') || strcmpi(sp{f},'Features3D')                  )
%                       if strcmpi(sp{f},'CarmPose')
%                           cp=CarmPose('');
%                           c=properties(s.(sp{f}));
%                           for g=1:numel(c)-1
%                               if ~strcmpi(c{g},'World')
%                                   cp.(c{g})=s.(sp{f}).(c{g});                              
%                               end                          
%                           end
%                           newObj.(sp{f})=cp;
%                       else
%                           newObj.(sp{f})=s.(sp{f});
%                       end
%                       
%                   end
%               end              
%               obj=newObj;
%           else
%               obj=s; % copy otherwise
%           end
%        end
 
       end
   
    methods
        
        function Count=get.FeatureCount(Obj)
            Count=size(Obj.Features,2);
        end
        
        function Out=get.Features3D(Obj)
            
            % Define unit vectors with x and y along the plane of the image frame, z
            % normal to the image
            j = Obj.ImageUp/norm(Obj.ImageUp);
            k = Obj.ImageNormal/norm(Obj.ImageNormal);
            i = cross(j,k);
            
            Xsize = size(Obj.Image,1);
            Ysize = size(Obj.Image,2);
            
            % Transformation matrix to go from 2D to 3D:
            Pnts2D=[0,0,0,1;10,0,0,1;0,20,0,1;0,0,30,1]';
            Pnts3D=zeros(4,4);
            for f=1:4
                Pnts3D(1:3,f) = Obj.CentrePoint + (Pnts2D(3,f))*k + (Pnts2D(2,f))*j + (Pnts2D(1,f))*i;
            end
            Pnts3D(4,:)=1;
            Trans2D3D=Pnts3D/Pnts2D;
            F3D=[];
            for f=1:Obj.FeatureCount
                Q =Obj.Features(f).Points;
%                 P(:,1)=-(Q(:,2)-Xsize/2)*Obj.PixelScale;
%                 P(:,2)=(-Q(:,1)+Ysize/2)*Obj.PixelScale;
                P(:,1)=(Q(:,1)-Xsize/2)*Obj.PixelScale;
                P(:,2)=(Q(:,2)-Ysize/2)*Obj.PixelScale;
                P(:,3)=0;
                P(:,4)=1;
                P=Trans2D3D*P';
                F3D(f).Points=P(1:3,:)';
                F3D(f).FitPoints=P(1:3,:)';
                F3D(f).FitDots=P(1:3,:)';
                F3D(f).Name=Obj.Features(f).Name;                
            end
            Out=F3D;
        end
        
        function Obj=ImagePlane(PlaneName,Path,JTCalibFileName,JTImageFileName,JTPoseFileName)
            switch nargin
                case 1
                    Obj.Name=PlaneName;
                    Obj.ImageRight=[1 0 0];
                    Obj.ImageUp=[0 1 0];
                    Obj.ImageNormal=[0 0 1];
                case 4
                    Obj.Name=PlaneName;
                    Obj.Path=Path;
                    Obj.JTCalibFileName=JTCalibFileName;
                    Obj.JTImageFileName=JTImageFileName;
                    Obj.loadJTCalib(Obj,fullfile(Path,JTCalibFileName));
                    Obj.Features=[];                    
                    
                    Obj.Image = imread(fullfile(Path,JTImageFileName));
                    if max(size(size(Obj.Image)))==2
                        Obj.Image=Obj.Image(:,:);
                    else
                        Obj.Image=Obj.Image(:,:,1:3);
                    end
                    
                    % Check it a -mat file exist and load it if possible
                    FName=strcat(JTCalibFileName(1:size(JTCalibFileName,2)-3),'mat');
                    if exist(FName,'file')
                        load(FName,'-mat');
                        Obj.CarmPose=ImgPlane.CarmPose;
                    end
                    
                case 5
                    Obj.Name=PlaneName;
                    Obj.Path=Path;
                    Obj.JTCalibFileName=JTCalibFileName;
                    Obj.JTImageFileName=JTImageFileName;
                    Obj.loadJTCalib(Obj,fullfile(Path,JTCalibFileName));                    
                    Obj.Features=[];
                    
                    Obj.Image = imread(fullfile(Path,JTImageFileName));
                    if max(size(size(Obj.Image)))==2
                        Obj.Image=Obj.Image(:,:);
                    else
                        Obj.Image=Obj.Image(:,:,1:3);
                    end                    
            
                    if ~isempty(JTPoseFileName)
                        FName=fullfile(Path,JTPoseFileName);
                        load(FName,'-mat');
                        Obj.CarmPose=ImgPlane.CarmPose;
                    end
            end
            Obj.InitMag=0;
        end
        
        function readIPinfo(Obj,IPinfo)
            Obj.PrincipalDist=IPinfo.PrincipalDist;
            Obj.CentrePoint=IPinfo.CentrePoint';
            Obj.Xoffset=IPinfo.Xoffset;
            Obj.Yoffset=IPinfo.Yoffset;
            Obj.PixelScale=IPinfo.PixelScale;
            Obj.SourcePosition=IPinfo.SourcePosition';
            Obj.ImageNormal=IPinfo.ImageNormal';
            Obj.ImageUp=IPinfo.ImageUp';
            Obj.ImageRight=cross(Obj.ImageUp,Obj.ImageNormal);
        end
        
        function addFeature(Obj,Feature)
            Obj.Features=[Obj.Features,Feature];
        end
        
        function [Feat,FeatIndex]=findFeature(Obj,FName)
            FeatIndex=0;
            Feat=[];
            for f=1:Obj.FeatureCount
                if strcmp(FName,Obj.Features(f).Name)
                    Feat=Obj.Features(f);
                    FeatIndex=f;
                end
            end
        end
        
        function deleteFeature(Obj,FName)
            Feats = [];
            for f=1:Obj.FeatureCount
                if f < Obj.FeatureCount
                    if strcmp(FName,Obj.Features(f).Name)
                        Obj.Features(f) = [];
                    else
                        Feats = [Feats,Obj.Features(f)];
                    end
                end
            end
            Obj.Features = Feats;
        end
            
        
        function [Feat,FeatIndex]=findFeatbyType(Obj,FType)
            FeatIndex=[];
            Feat=[];
            cnt=0;
            for f=1:Obj.FeatureCount
                if strcmp(FType,Obj.Features(f).Type)
                    cnt=cnt+1;
                    if cnt==1
                        Feat=Obj.Features(f);
                    else
                        Feat(cnt)=Obj.Features(f);
                    end
                    FeatIndex(cnt)=f;
                end
            end
        end
        
        function loadJTCalib(Obj,JTFileName)
            calibdata = textread(JTFileName,'%f', 13, 'headerlines', 1);
            caldata = struct('principdist', calibdata(1), 'centrepoint',[0,0,0],'xoffset', calibdata(2), 'yoffset', calibdata(3), 'pixconv', calibdata(4), 'sourcepos', calibdata(5:7), 'imgnorm', calibdata(8:10), 'imgup', calibdata(11:13));            
            Obj.PrincipalDist=caldata.principdist;
            ctrpt = Obj.imgfindcentre(caldata);
            Obj.CentrePoint=ctrpt';
            Obj.Xoffset=caldata.xoffset;
            Obj.Yoffset=caldata.yoffset;
            Obj.PixelScale=caldata.pixconv;
            Obj.SourcePosition=caldata.sourcepos';
            Obj.ImageNormal=caldata.imgnorm';
            Obj.ImageUp=caldata.imgup';
            Obj.ImageRight=cross(Obj.ImageUp,Obj.ImageNormal);
        end
        
        % find the image centre
        function ctrpt = imgfindcentre(Obj,frame)
            % Define unit vectors with x and y along the plane of the image frame, z
            % normal to the image/intensifier/detector
            j = frame.imgup/norm(frame.imgup);
            k = frame.imgnorm/norm(frame.imgnorm);
            i = cross(j,k);
            
            % Convert the pixel points on the image into coordinates in 3D space.
            % We assume that the origin of the image coordinate system is at the
            % centre of the image
            
            ctrpt = frame.sourcepos + frame.principdist*-k - frame.xoffset*i - frame.yoffset*j;
        end
        
        % Plot the image in 3D
        function imshow3d(Obj,ShowAxis,ShowFeatures)

            if ~isempty(Obj.Image)
                imageData=Obj.Image;
            else
                imageData=imread(fullfile(Obj.Path,Obj.JTImageFileName));
            end
            
            if size(imageData,3)==4
                imageData=imageData(:,:,1:3);
            end
            CentPoint=Obj.CentrePoint;

            Vy=Obj.ImageUp;
            Vz=Obj.ImageNormal;
            Vx=cross(Vy,Vz);
            Vx=Vx/norm(Vx);
            Vy=Vy/norm(Vy);
            Vz=Vz/norm(Vz);
            
            Dx=(size(imageData,1))*Obj.PixelScale;
            Dy=(size(imageData,2))*Obj.PixelScale;
                       
            P1=CentPoint-Vx*(Dx/2)-Vy*(Dy/2);
            P2=CentPoint+Vx*(Dx/2)-Vy*(Dy/2);
            P3=CentPoint-Vx*(Dx/2)+Vy*(Dy/2);
            P4=CentPoint+Vx*(Dx/2)+Vy*(Dy/2);            
            
            Xdata=[P1(1),P2(1);P3(1),P4(1)];
            Ydata=[P1(2),P2(2);P3(2),P4(2)];
            Zdata=[P1(3),P2(3);P3(3),P4(3)];
            
            colormap(bone);
            
%             surface('XData',Xdata,'YData',Ydata,...
%                 'ZData',Zdata,'CData',rot90(imageData,2),...
%                 'FaceColor','texturemap','EdgeColor','none','FaceAlpha',0.5);           
  
            surface('XData',Xdata,'YData',Ydata,...
                'ZData',Zdata,'CData',imageData,...
                'FaceColor','texturemap','EdgeColor','c','FaceAlpha',0.5);           
            
            if ShowAxis
                hold on;
                zdir=CentPoint+(Vz*100);
                xdir=CentPoint+(Vx*100);
                ydir=CentPoint+(Vy*100);
                plot3([CentPoint(1),zdir(1)],[CentPoint(2),zdir(2)],[CentPoint(3),zdir(3)],'b');
                plot3([CentPoint(1),xdir(1)],[CentPoint(2),xdir(2)],[CentPoint(3),xdir(3)],'r');
                plot3([CentPoint(1),ydir(1)],[CentPoint(2),ydir(2)],[CentPoint(3),ydir(3)],'g');
                plot3(Obj.SourcePosition(1),Obj.SourcePosition(2),Obj.SourcePosition(3),'m*');
                plot3([Obj.SourcePosition(1),CentPoint(1)],[Obj.SourcePosition(2),CentPoint(2)],[Obj.SourcePosition(3),CentPoint(3)],'m');
            end
            
            if ShowFeatures
                hold on;
                for f=1:Obj.FeatureCount
                    Pnts=[];
                    PPnts=Obj.Features3D;
                    Pnts=PPnts(f).FitDots;
                    dv=Vx/norm(Vx)*0;
                    Pnts(:,1)=Pnts(:,1)+dv(1);
                    Pnts(:,2)=Pnts(:,2)+dv(2);
                    Pnts(:,3)=Pnts(:,3)+dv(3);
                    plot3(Pnts(:,1),Pnts(:,2),Pnts(:,3),'go','LineWidth',2);                    
                end
            end
        end
                     
        % Translate an ImagePlane by T matrix:
        function Translate(Obj,T)
            
            if isrow(Obj.ImageNormal)
                Vz=Obj.ImageNormal;
            else
                Vz=Obj.ImageNormal';
            end
            if isrow(Obj.ImageUp)
                Vy=Obj.ImageUp;
            else
                Vy=Obj.ImageUp';
            end                  
            if isrow(Obj.CentrePoint)
                cent=Obj.CentrePoint;
            else
                cent=Obj.CentrePoint';
            end
                        
            Vx=cross(Vy,Vz);                
            NewT=T*([Vx,0;Vy,0;Vz,0;cent,1]'); % New position of the ImagePlane
            
            Obj.CentrePoint=NewT(1:3,4)';
            NVx=NewT(1:3,1)';
            NVy=NewT(1:3,2)';
            NVz=NewT(1:3,3)';
            
            Obj.ImageNormal=NVz;
            Obj.ImageUp=NVy;
            Obj.ImageRight=cross(Obj.ImageUp,Obj.ImageNormal);
            Obj.SourcePosition=Obj.CentrePoint + NVx*Obj.Xoffset + NVy*Obj.Yoffset + NVz*Obj.PrincipalDist;
            Obj.PiercingPoint=Obj.CentrePoint + NVx*Obj.PiercingPointLocal(1) + NVy*Obj.PiercingPointLocal(2);
        end
        
        function [IsValid,AngleError,DistError]=validate(Obj)
            % Validate the ImagePlane for the integrity of the properties.
            % Retun Angular and Distance discrepancies as optional outputs
            % Validate using the tolerance values below.
            TolAngle=0.001;
            TolDist=0.01;
            
            AngleError=abs(acosd(dot(Obj.ImageNormal,Obj.ImageUp))-90);
            if AngleError<TolAngle
                AngleValid=1;
            else
                AngleValid=0;
            end
            
            DistVec=Obj.SourcePosition-Obj.CentrePoint;
            DistError=abs(dot(DistVec,Obj.ImageNormal)-Obj.PrincipalDist);
            if DistError<TolDist
                DistValid=1;
            else
                DistValid=0;
            end
            
            IsValid=DistValid && AngleValid;
        end
    end
end

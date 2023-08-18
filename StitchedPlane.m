classdef StitchedPlane < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        ImagePlanes
        Grid
        Origin
        Vector1
        Vector2
        Vector3
        Image
        ImageSize
        ImageSrcIndx
        ImageDepthIndx % Channel containing the DepthSegment Indices
        
        PlaneType
        PixelPerMM
        PixelScale
        
        GridImageOriginOffset
        GridImageSize
        GridOriginOffetX
        GridOriginOffetY
        GridOriginOffetZ
        CentrePoint
        ImageUp
        ImageNormal
        
        GridFlag
        GridImageFileRead
        BasedOnGrid  % is 1 if the StitchedPlane is already initialized
        UndoInfo
        CroppedImage
        CroppedCoords % Adjusted to a small image if there is no nozero values to crop the image to the size
        OrigCroppedCoords % Orginial Crop coords, which remains largest possible amount if there are no nonzero pixels to calculate the cropping area
        
        DrawingData
        Vis2D
        DownSampleSize
        SaveVis2D_CroppedImage
        SaveVis2D_ImageOutput
        SaveVis2D_Image
        SaveVis2D_CroppedImage_Bare
        DrawingDataOffsetX % if there needs to be an offset after cropping update
        DrawingDataOffsetY
        BlankGrid
        ImageFuseCirRads
        
        AnatomySegments    % 2D points defining the depth of the anatomy on the same plane in 2D
        AnatomySegments3D  % 3D points defining the depth of the anatomy on the same plane in 3D
        
        DepthSegments      % 3D points defining the depth from the paired stitchedplane
        DepthSegmentsNormal% 3D vector defining the depth surface (same as the normal vector to the paried stitchedplane)
        
        CorrectInsideParallax % Correct parallax based on the calculated inside control points
        ParralaxCorrectionRes % Parallax Correction Resolution (smallest divisions in pixels)
        
        
        CNodes          % 2D points defining the nodes of the calibration tool % These are parameters needed to reposition the calibration nodes
        CNodes3D        % 3D points defining the nodes of the calibration tool
        CNodesSource3D  % 3D points defining the nodes of the calibration tool
        CNodes3DReg     % 3D proj points calculated based on depth
        
        DNodes          % 2D points defining the nodes of the depth nodes
        DNodes3D        % 3D points defining the nodes of the depth nodes
        DNodesSource3D  % 3D points defining the nodes of the depth nodes
        DNodes3DReg     % 3D proj points calculated based on depth
        
        
        
        
    end
    
    
    methods
        
        function Obj = StitchedPlane(EName,varargin)
            Obj.Name = EName;
            Obj.GridOriginOffetX=0;
            Obj.GridOriginOffetY=0;
            Obj.GridOriginOffetZ=0;
            Obj.GridImageOriginOffset=[0,0];
            Obj.GridImageSize=[1024*3,1024*10];
            Obj.BasedOnGrid=0; % Initially have this zero, but change it to one as soon as a grid is initialized
            Obj.UndoInfo=[]; %structure to hold the undo add image info
            Obj.CroppedImage=[];
            Obj.CroppedCoords=[];
            Obj.Vis2D=[0,0];
            Obj.SaveVis2D_CroppedImage=[];
            Obj.SaveVis2D_ImageOutput=[];
            Obj.SaveVis2D_Image=[];
            Obj.DrawingDataOffsetX=0;
            Obj.DrawingDataOffsetY=0;
            Obj.ImageFuseCirRads=[100,100];
            Obj.AnatomySegments=[];
            Obj.AnatomySegments3D=[];
            
            Obj.DepthSegments=[];
            Obj.DepthSegmentsNormal=[];
            Obj.CorrectInsideParallax=1;
            Obj.ParralaxCorrectionRes=200;
            
            % Reset the calibration nodes
            Obj.CNodes=[];
            Obj.CNodes3D=[];
            Obj.CNodesSource3D=[];
            Obj.DNodes=[];
            Obj.DNodes3D=[];
            Obj.DNodesSource3D=[];
            
            if ~isempty(varargin)
                Obj.Grid=varargin{1};
                Obj.PlaneType=varargin{2};
                Obj.GridImageFile=varargin{3};
                Obj.GridImageSize=varargin{4};
                Obj.PixelPerMM=varargin{5};
                
                IM=imread(Obj.GridImageFile);
                
                Obj.GridImage=imresize(IM, Obj.GridImageSize);
                Obj.Image=Obj.GridImage*0; % Temporarily set to zero: Jan 18, 2017
                
                Obj.Origin=Obj.Grid.Origin;
                Obj.GridOriginOffetX=0;
                Obj.GridOriginOffetY=0;
                
                switch Obj.PlaneType
                    case 'XY'
                        Obj.Vector1=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VyLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VzLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                    case 'YZ'
                        Obj.Vector1=Obj.Grid.VyLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VxLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                    case 'XZ'
                        Obj.Vector1=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VyLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                end
            end
        end
        
        function readGridImageFile(Obj)
            Obj.GridImageFileRead=imread(Obj.GridImageFile);
            %Obj.GridImageFileRead=double(rgb2gray(IM))/255;
        end
        
        function Obj = reconStitchedPlane(Obj,EName,varargin)
            Obj.Name = EName;
            
            if isempty(varargin)
                Obj.GridOriginOffetX=0;
                Obj.GridOriginOffetY=0;
                Obj.GridImageSize=[1024*3,1024*10];
                Obj.GridImage(1:Obj.GridImageSize(1),1:Obj.GridImageSize(2))=0;
                Obj.GridImageFile='';
            end
            
            if ~isempty(varargin)
                
                Obj.PlaneType=varargin{1};
                Obj.GridImageSize=varargin{2};
                Obj.PixelScale=varargin{3};
                Obj.PixelPerMM=1/Obj.PixelScale;
                Obj.GridImageOriginOffset=varargin{4};
                Obj.BasedOnGrid=varargin{5};
                
                if ~Obj.BasedOnGrid % Run this part only during initialization (since this is time consuming)
                    
                    Sz=ceil(Obj.GridImageSize/Obj.PixelScale);
                    switch Obj.PlaneType
                        case 'XY'
                            Obj.Image=uint8(zeros(Sz(2),Sz(1)));
                            Obj.ImageSrcIndx=uint8(zeros(Sz(2),Sz(1)));
                            Obj.ImageDepthIndx=uint8(zeros(Sz(2),Sz(1)));
                        case 'YZ'
                            Obj.Image=uint8(zeros(Sz(1),Sz(2)));
                            Obj.ImageSrcIndx=uint8(zeros(Sz(1),Sz(2)));
                            Obj.ImageDepthIndx=uint8(zeros(Sz(1),Sz(2)));
                    end
                    
                end
                
                %                 Obj.ImageOutput=cat(3,Obj.Image*0,Obj.Image*0,Obj.Image*0);
                Obj.Origin=Obj.Grid.Origin;
                
                % Old definitions
                switch Obj.PlaneType
                    case 'XY'
                        Obj.Vector1=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VyLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VzLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=Obj.Vector2;
                        Obj.CentrePoint=Obj.Origin;
                        
                    case 'YZ'
                        Obj.Vector1=Obj.Grid.VyLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VxLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=Obj.Vector2;
                        Obj.CentrePoint=Obj.Origin;
                    case 'XZ'
                        Obj.Vector1=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VyLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=1*Obj.Vector1;
                        Obj.CentrePoint=Obj.Origin;
                end
                
                
                
                
                
            end
        end
        
        
        function Obj = updateGrid(Obj,StandardViews)
            
            Obj.Origin=Obj.Grid.Origin;
            
            if StandardViews
                FirstPlaneStr='XY';
                SecondPlaneStr='YZ';
            else
                FirstPlaneStr='YZ';
                SecondPlaneStr='XY';
            end
            
            Obj.Grid.ScoutAPNormal=[];
            Obj.Grid.ScoutMLNormal=[];
            if (~isempty(Obj.Grid.ScoutAPNormal)) && (~isempty(Obj.Grid.ScoutMLNormal))
                switch Obj.PlaneType
                    
                    case (FirstPlaneStr)
                        v3=Obj.Grid.ScoutAPNormal;
                        v1=Obj.Grid.ScoutMLNormal;
                        v2=cross(v3,v1);
                        v1=cross(v2,v3);
                        
                        Obj.Vector1=v1;
                        Obj.Vector2=v2;
                        Obj.Vector3=v3;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=Obj.Vector2;
                        Obj.CentrePoint=Obj.Origin;
                        
                    case (SecondPlaneStr)
                        v3=Obj.Grid.ScoutMLNormal;
                        v2=Obj.Grid.ScoutAPNormal;
                        v1=cross(v2,v3);
                        v2=cross(v3,v1);
                        
                        Obj.Vector1=v1; %Obj.Grid.VyLineElement.Vector;
                        Obj.Vector2=v2; %Obj.Grid.VzLineElement.Vector;
                        Obj.Vector3=v3; %Obj.Grid.VxLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=Obj.Vector2;
                        Obj.CentrePoint=Obj.Origin;
                        
                    case 'XZ'
                        Obj.Vector1=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VyLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=1*Obj.Vector1;
                        Obj.CentrePoint=Obj.Origin;
                end
                
            else
                
                switch Obj.PlaneType
                    case 'XY'
                        Obj.Vector1=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VyLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VzLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=Obj.Vector2;
                        Obj.CentrePoint=Obj.Origin;
                        
                    case 'YZ'
                        Obj.Vector1=Obj.Grid.VyLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VxLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=Obj.Vector2;
                        Obj.CentrePoint=Obj.Origin;
                    case 'XZ'
                        Obj.Vector1=Obj.Grid.VzLineElement.Vector;
                        Obj.Vector2=Obj.Grid.VxLineElement.Vector;
                        Obj.Vector3=Obj.Grid.VyLineElement.Vector;
                        Obj.Origin=Obj.Grid.Origin;
                        Obj.ImageNormal=Obj.Vector3;
                        Obj.ImageUp=1*Obj.Vector1;
                        Obj.CentrePoint=Obj.Origin;
                end
            end
            
            
            
            
        end
        
        
        function [ghandles,Window,ErrorCode,XD,YD]=renderImage_Old(Obj,ImagePlane, ~, ProcessImage, ShowGraphics, ~, ImagePlaneIndex)
            ghandles={};
            ErrorCode=0;
            XD=[];
            YD=[];
            if isempty(ShowGraphics)
                ShowGraphics=0;
            end
            Pp=Obj.Origin;
            Pv=-Obj.Vector3;
            Vz=ImagePlane.ImageNormal;
            Vy=ImagePlane.ImageUp;
            Vx=cross(Vy,Vz);
            Vx=Vx/norm(Vx);
            Vy=Vy/norm(Vy);
            
            O=ImagePlane.CentrePoint;
            PS=ImagePlane.PixelScale;
            Sx=1024;
            Sy=1024;
            ISFLD=0;
            if isfield(ImagePlane,'Image') || isprop(ImagePlane,'Image')
                if ~isempty(ImagePlane.Image)
                    [Sy,Sx]=size(ImagePlane.Image(:,:,1));
                    ISFLD=1;
                end
            end
            if ~ISFLD
                ImagePlane.Image=uint8(zeros(Sx,Sy)); % Set the imageplane
            end
            
            % Finding the intersecting points:
            P1=O+(PS*Sy/2)*(Vy)-(PS*Sx/2)*Vx;
            P2=O+(PS*Sy/2)*(Vy)+(PS*Sx/2)*Vx;
            P3=O-(PS*Sy/2)*(Vy)+(PS*Sx/2)*Vx;
            P4=O-(PS*Sy/2)*(Vy)-(PS*Sx/2)*Vx;
            
            %%
            if ProcessImage
                if (Obj.CorrectInsideParallax && ~isempty(Obj.DepthSegments))
                    % calculate the inbetween points
                    
                    Po0=O+(PS*Sy/2)*(Vy)-(PS*Sx/2)*Vx;
                    Po1=O+(PS*Sy/2)*(Vy)-(PS*Sx/4)*Vx;
                    Po2=O+(PS*Sy/2)*(Vy);
                    Po3=O+(PS*Sy/2)*(Vy)+(PS*Sx/4)*Vx;
                    Po4=O+(PS*Sy/2)*(Vy)+(PS*Sx/2)*Vx;
                    
                    Pa0=O+(PS*Sy/4)*(Vy)-(PS*Sx/2)*Vx;
                    Pa1=O+(PS*Sy/4)*(Vy)-(PS*Sx/4)*Vx;
                    Pa2=O+(PS*Sy/4)*(Vy);
                    Pa3=O+(PS*Sy/4)*(Vy)+(PS*Sx/4)*Vx;
                    Pa4=O+(PS*Sy/4)*(Vy)+(PS*Sx/2)*Vx;
                    
                    Pb0=O-(PS*Sx/2)*Vx;
                    Pb1=O-(PS*Sx/4)*Vx;
                    Pb2=O;
                    Pb3=O+(PS*Sx/4)*Vx;
                    Pb4=O+(PS*Sx/2)*Vx;
                    
                    Pc0=O-(PS*Sy/4)*(Vy)-(PS*Sx/2)*Vx;
                    Pc1=O-(PS*Sy/4)*(Vy)-(PS*Sx/4)*Vx;
                    Pc2=O-(PS*Sy/4)*(Vy);
                    Pc3=O-(PS*Sy/4)*(Vy)+(PS*Sx/4)*Vx;
                    Pc4=O-(PS*Sy/4)*(Vy)+(PS*Sx/2)*Vx;
                    
                    Pd0=O-(PS*Sy/2)*(Vy)-(PS*Sx/2)*Vx;
                    Pd1=O-(PS*Sy/2)*(Vy)-(PS*Sx/4)*Vx;
                    Pd2=O-(PS*Sy/2)*(Vy);
                    Pd3=O-(PS*Sy/2)*(Vy)+(PS*Sx/4)*Vx;
                    Pd4=O-(PS*Sy/2)*(Vy)+(PS*Sx/2)*Vx;
                    
                    InBtwPnts=[Po0,Po1,Po2,Po3,Po4,...
                        Pa0,Pa1,Pa2,Pa3,Pa4,...
                        Pb0,Pb1,Pb2,Pb3,Pb4,...
                        Pc0,Pc1,Pc2,Pc3,Pc4,...
                        Pd0,Pd1,Pd2,Pd3,Pd4];
                    
                    CalcInBtwPnts=zeros(2,25);
                    LocalDepthMesh=zeros(2,25);
                    CorrInBtwPnts=zeros(2,25);
                    Indx=0;
                    %                     for f=[1,5,21,25,2:4,6:20,22:24]  % first caclulate the four corners of the rays and then the rest of the points inside
                    for f=[1,5,21,25,1:25]  % first caclulate the four corners of the rays and then the rest of the points inside
                        Indx=Indx+1;
                        if Indx==5 % calculate the plane main plane after the 4 corners are calculated:
                            TP1=LocalDepthMesh(1:3,1);
                            TP2=LocalDepthMesh(1:3,5);
                            TP3=LocalDepthMesh(1:3,21);
                            TP4=LocalDepthMesh(1:3,25);
                        end
                        
                        InBtwP=InBtwPnts(1:3,f);
                        SRC=ImagePlane.SourcePosition;
                        In_VV=(InBtwP-SRC);
                        In_VV=In_VV/norm(In_VV);
                        [~,In_PI1]=lineAnatomyIntersect(Pv,Pp,In_VV,InBtwP,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                        
                        % project the inner parallax point over the main local box:
                        if Indx<5
                            % if exactly on the corners take the actual intersection point:
                            In_PP1=In_PI1;
                        else
                            % If inside the box, then project the point over the local box:
                            % Check which of the triangles TP1_TP2_TP3 or TP1_TP3_TP4 is the note in:
                            
                            TRefPoint=TP1;
                            %Trangle 1:
                            TPlaneN1=cross(TP2-TRefPoint,TP3-TRefPoint);
                            TPlaneN1=TPlaneN1/norm(TPlaneN1);
                            %Trangle 2:
                            TPlaneN2=cross(TP4-TRefPoint,TP3-TRefPoint);
                            TPlaneN2=TPlaneN2/norm(TPlaneN2);
                            
                            In_PP1_1=linePlaneIntersect(TPlaneN1,TRefPoint,Pv,In_PI1);
                            In_PP1_2=linePlaneIntersect(TPlaneN2,TRefPoint,Pv,In_PI1);
                            
                            % compare the distances to decide which point should be selected:
                            if norm(In_PP1_1-TP2)<norm(In_PP1_2-TP4)
                                In_PP1=In_PP1_1;
                            else
                                In_PP1=In_PP1_2;
                            end
                            
                        end
                        CalcInBtwPnts(1:3,f)=In_PP1;  % calculation of the loction of the projection point in 3D
                        LocalDepthMesh(1:3,f)=In_PI1; % coordinate of the mesh nodes
                        % Corrected point coordinates:
                        CorrInBtwPnts(1:3,f)=linePlaneIntersect(ImagePlane.ImageNormal,ImagePlane.CentrePoint,(SRC-In_PP1)/norm(SRC-In_PP1),In_PP1);
                    end
                    
                end
                
            end
            
            %%
            OP1=P1;
            OP2=P2;
            OP3=P3;
            OP4=P4;
            
            SRC=ImagePlane.SourcePosition;
            VV1=(P1-SRC);
            VV2=(P2-SRC);
            VV3=(P3-SRC);
            VV4=(P4-SRC);
            
            VV1=VV1/norm(VV1);
            VV2=VV2/norm(VV2);
            VV3=VV3/norm(VV3);
            VV4=VV4/norm(VV4);
            
            if isempty(Obj.DepthSegments)
                % Simple cone/plane intersection without variable depth segments:
                PP1=linePlaneIntersect(Pv,Pp,VV1,P1);
                PP2=linePlaneIntersect(Pv,Pp,VV2,P2);
                PP3=linePlaneIntersect(Pv,Pp,VV3,P3);
                PP4=linePlaneIntersect(Pv,Pp,VV4,P4);
                
                IPs=[];
            else
                % Account for the depth of the anatomy
                [PP1,PI1]=lineAnatomyIntersect(Pv,Pp,VV1,P1,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                [PP2,PI2]=lineAnatomyIntersect(Pv,Pp,VV2,P2,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                [PP3,PI3]=lineAnatomyIntersect(Pv,Pp,VV3,P3,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                [PP4,PI4]=lineAnatomyIntersect(Pv,Pp,VV4,P4,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                IPs=[PI1';PI2';PI3';PI4';PI1'];
            end
            
            %Points of Intersections with the Grid plane:
            PPs=[PP1';PP2';PP3';PP4'];
            
            [OverallCheck,~]=Obj.checkPointsOut(PPs);
            
            if OverallCheck~=0 % if points stay out of the bounds:
                Window=[];     % return an empty window
                ErrorCode=1;   % Image Out of Bounds
                
            else % carry on with normal processing:
                V2=Obj.ImageUp;
                V3=Obj.ImageNormal;
                V1=cross(V2,V3);
                V1=V1/norm(V1);
                V2=V2/norm(V2);
                V3=V3/norm(V3);
                
                Pnt1=PP1-Obj.Origin;
                Pnt2=PP2-Obj.Origin;
                Pnt3=PP3-Obj.Origin;
                Pnt4=PP4-Obj.Origin;
                
                %Points with respect to the Grid [2D]:
                Pr1=[dot(Pnt1,V1),dot(Pnt1,V2),dot(Pnt1,V3)];
                Pr2=[dot(Pnt2,V1),dot(Pnt2,V2),dot(Pnt2,V3)];
                Pr3=[dot(Pnt3,V1),dot(Pnt3,V2),dot(Pnt3,V3)];
                Pr4=[dot(Pnt4,V1),dot(Pnt4,V2),dot(Pnt4,V3)];
                
                Prs=[Pr1;Pr2;Pr3;Pr4];
                Prs_Grid=Prs; %All corners wrt the Grid
                
                % Found rectangular box around the image wrt the Grid CS:
                Cr1_Grid=[min(Prs_Grid(:,1)),min(Prs_Grid(:,2))];
                Cr2_Grid=[max(Prs_Grid(:,1)),min(Prs_Grid(:,2))];
                Cr3_Grid=[min(Prs_Grid(:,1)),max(Prs_Grid(:,2))];
                Cr4_Grid=[max(Prs_Grid(:,1)),max(Prs_Grid(:,2))];
                Crs_Grid=[Cr1_Grid;Cr2_Grid;Cr3_Grid;Cr4_Grid];
                
                % Mapping the corners of the box to the global coordinate:
                Cr1=Obj.Origin + Cr1_Grid(1)*V1 + Cr1_Grid(2)*V2;
                Cr2=Obj.Origin + Cr2_Grid(1)*V1 + Cr2_Grid(2)*V2;
                Cr3=Obj.Origin + Cr3_Grid(1)*V1 + Cr3_Grid(2)*V2;
                Cr4=Obj.Origin + Cr4_Grid(1)*V1 + Cr4_Grid(2)*V2;
                Crs=[Cr1';Cr2';Cr3';Cr4'];
                
                % Sorting the corners:
                Indx=ones(4,1);
                for f=1:size(PPs,1)
                    MinDst=norm(Crs(1,:)-PPs(f,:));
                    for g=1:size(Crs,1)
                        Dst=norm(Crs(g,:)-PPs(f,:));
                        if Dst<MinDst
                            MinDst=Dst;
                            Indx(f)=g;
                        end
                    end
                end
                Crs=Crs(Indx,:); %Sorting the Crs
                Crs_Grid=Crs_Grid(Indx,:);
                
                Sz=round(round(Cr4_Grid-Cr1_Grid)*Obj.PixelPerMM);
                
                SCFR=Sz(2)/Sx; % Scale factor calculated based on the size of the image.
                MF=1/SCFR; % *2; %*10; % Magnification factor:
                
                Crs_Grid_p=Crs_Grid;
                Crs_Grid_p(:,1)=Crs_Grid_p(:,1)-min(Crs_Grid_p(:,1));
                Crs_Grid_p(:,2)=Crs_Grid_p(:,2)-min(Crs_Grid_p(:,2));
                
                Prs_Grid_p=Prs_Grid(:,1:2);
                Prs_Grid_p(:,1)=Prs_Grid_p(:,1)-min(Prs_Grid_p(:,1));
                Prs_Grid_p(:,2)=Prs_Grid_p(:,2)-min(Prs_Grid_p(:,2));
                
                
                GrdPnts=round(Crs_Grid_p*MF);
                OrgPnts=round(Prs_Grid_p*MF);
                %Check if the points calculated are good for the rest of the calculation, if not abort:
                IsCollinear=checkCollinearity(OrgPnts);
                
                if IsCollinear
                    Window=[];
                    ErrorCode=2;
                    
                else
                    
                    Window.Origin=[min(Crs_Grid_p(:,1)),min(Crs_Grid_p(:,2))];
                    Window.ScaleFactor=SCFR;
                    %Calculating the pixel Scaling
                    Window.PixelPerMM=1/ImagePlane.PixelScale;
                    
                    if ShowGraphics
                        h1=ImagePlane.showGraphically(1,0);
                        PPS=[OP1,PP1,PP2,OP2,PP2,PP3,OP3,PP3,PP4,OP4,PP4,PP1];
                        h2=[];
                        h2=plot3(PPS(1,:),PPS(2,:),PPS(3,:),'b');
                        PPS=[PP1,SRC,PP2,SRC,PP3,SRC,PP4];
                        h3=[];
                        h3=plot3(PPS(1,:),PPS(2,:),PPS(3,:),'r');
                        ii=[4,3,1,2];
                        
                        % Correct image boundaries
                        Xdata=[PP4(1),PP3(1);PP1(1),PP2(1)];
                        Ydata=[PP4(2),PP3(2);PP1(2),PP2(2)];
                        Zdata=[PP4(3),PP3(3);PP1(3),PP2(3)];
                        
                        imageData=ImagePlane.Image;
                        surf = surface('XData',Xdata,'YData',Ydata,...
                            'ZData',Zdata,'CData',imageData,...
                            'FaceColor','texturemap','EdgeColor','none','FaceAlpha',0.9); %0.5);
                        colormap(bone);
                        CRNS=[Crs;Crs(1,:)]';
                        h4=plot3(CRNS(1,:), CRNS(2,:), CRNS(3,:), 'm');
                        TXTCR=Crs';
                        TXT={'Crs1','Crs2','Crs3','Crs4'};
                        h5=text(TXTCR(1,:),TXTCR(2,:),TXTCR(3,:),TXT,'color','g');
                        TXT2={'PP1','PP2','PP3','PP4'};
                        h6=text(PPs(1:4,1),PPs(1:4,2),PPs(1:4,3),TXT2,'color','m');
                        if ~isempty(IPs)
                            plot3(IPs(:,1),IPs(:,2),IPs(:,3),'c');
                        end
                        axis equal;
                        ghandles={h1{:},h2,h3,surf,h4,h5,h6};
                    end
                    
                    
                    %%
                    % Calculate the local coordinates of the points:
                    
                    if ProcessImage
                        if (Obj.CorrectInsideParallax && ~isempty(Obj.DepthSegments))
                            
                            % Fix the local parallax:
                            % calculate the coordinates with respect to the initial reference frame:
                            P_1=Crs(1,:);
                            P_2=Crs(2,:);
                            P_3=Crs(3,:);
                            P_4=Crs(4,:);
                            
                            % points calculated wrt the corners
                            IP_Vx=(P_2-P_1)/norm(P_2-P_1);
                            IP_Vy=(P_4-P_1)/norm(P_4-P_1);
                            OutImg_Sx=norm(P_2-P_1);
                            OutImg_Sy=norm(P_4-P_1);
                            
                            % Point calculated wrt the original imageplane
                            IP2_Vx=(P2-P1)/norm(P2-P1);
                            IP2_Vy=(P4-P1)/norm(P4-P1);
                            
                            OutImg2_Sx=norm(CorrInBtwPnts(:,5)-CorrInBtwPnts(:,1));
                            OutImg2_Sy=norm(CorrInBtwPnts(:,21)-CorrInBtwPnts(:,1));
                            
                            CorrInBtwPnts2D=zeros(2,25);
                            % Calculate the 9 inner points CorrInBtwPnts2D for image processing:
                            for f=1:25
                                
                                CorrInBtwPnts2D(2,f)=dot(CorrInBtwPnts(:,f)-P_1',IP_Vx')/OutImg_Sx;
                                CorrInBtwPnts2D(1,f)=dot(CorrInBtwPnts(:,f)-P_1',IP_Vy')/OutImg_Sy;
                            end
                            
                            %normalize the coordinates:
                            Sx_ref=min(CorrInBtwPnts2D(1,:));
                            Sy_ref=min(CorrInBtwPnts2D(2,:));
                            Sx_Len=max(CorrInBtwPnts2D(1,:))-Sx_ref;
                            Sy_Len=max(CorrInBtwPnts2D(2,:))-Sy_ref;
                            CorrInBtwPnts2D(1,:)=CorrInBtwPnts2D(1,:)-Sx_ref;
                            CorrInBtwPnts2D(2,:)=CorrInBtwPnts2D(2,:)-Sy_ref;
                            CorrInBtwPnts2D(1,:)=CorrInBtwPnts2D(1,:)/Sx_Len;
                            CorrInBtwPnts2D(2,:)=CorrInBtwPnts2D(2,:)/Sy_Len;
                            
                            if ShowGraphics
                                hold on;
                                plot3(InBtwPnts(1,:),InBtwPnts(2,:),InBtwPnts(3,:),'+');
                                plot3(LocalDepthMesh(1,:),LocalDepthMesh(2,:),LocalDepthMesh(3,:),'gx');
                                plot3(CalcInBtwPnts(1,1:25),CalcInBtwPnts(2,1:25),CalcInBtwPnts(3,1:25),'ro');
                                plot3(CorrInBtwPnts(1,1:25),CorrInBtwPnts(2,1:25),CorrInBtwPnts(3,1:25),'mo');
                                % Labels={'Pa1','Pa2','Pa3','Pb1','Pb2','Pb3','Pc1','Pc2','Pc3','P1','P2','P3','P4'};
                                Labels={'Po0','Po1','Po2','Po3','Po4',...
                                    'Pa0','Pa1','Pa2','Pa3','Pa4',...
                                    'Pb0','Pb1','Pb2','Pb3','Pb4',...
                                    'Pc0','Pc1','Pc2','Pc3','Pc4',...
                                    'Pd0','Pd1','Pd2','Pd3','Pd4'};
                                for f=1:numel(Labels)
                                    text(InBtwPnts(1,f)+10,InBtwPnts(2,f)+10,InBtwPnts(3,f)+10,Labels{f});
                                end
                            end
                            
                        end
                    end
                    %%
                    
                    CRN1_2D=Crs_Grid(1,1:2);
                    CRN2_2D=Crs_Grid(2,1:2);
                    CRN3_2D=Crs_Grid(3,1:2);
                    CRN4_2D=Crs_Grid(4,1:2);
                    
                    % Corners of the image wrt the grid plane (2D):
                    CRNS_2D=[CRN1_2D; CRN2_2D; CRN3_2D; CRN4_2D];
                    Window.Coords=CRNS_2D;
                    
                    if ~ProcessImage % don't do anything with the image: normally used only for updating the window information.
                        Ip=ImagePlane.Image/255;
                    else
                        
                        %% Calculate the sizr of the image:
                        if isempty(Window) % out of bounds detected
                            ErrorCode=1;
                        else
                            Pnts=zeros(size(Window.Coords,1),3);
                            for f=1:size(Window.Coords,1)
                                P=Window.Coords(f,:);
                                P3D=Obj.CentrePoint+Obj.Vector1*P(1)+Obj.Vector2*P(2);
                                Pnts(f,:)=P3D';
                            end
                        end
                        
                        if ErrorCode==0
                            %GridOrigin coordinates in 3d and in mm:
                            GridOrigin3D=Obj.CentrePoint + Obj.GridImageOriginOffset(1)*Obj.Vector2*Obj.PixelScale+ Obj.GridImageOriginOffset(2)*Obj.Vector1*Obj.PixelScale;
                            %Coordinates of the topleft corner with respect to the origin:
                            CR1x=-(size(Obj.Image,2)/2 + Obj.GridImageOriginOffset(1))*Obj.PixelScale;
                            CR1y=-(size(Obj.Image,1)/2 + Obj.GridImageOriginOffset(2))*Obj.PixelScale;
                            CR13D=GridOrigin3D+Obj.Vector1*CR1x+Obj.Vector2*CR1y;
                            Win_Crns2D=zeros(4,2);
                            Win_Crns3D=zeros(4,3);
                            for f=1:4
                                P=Window.Coords(f,:);
                                Pnt=Obj.CentrePoint+Obj.Vector1*P(1)+Obj.Vector2*P(2); % corner of the window
                                Vec=Pnt-CR13D;
                                CRX=round(dot(Vec,Obj.Vector1)/Obj.PixelScale);
                                CRY=round(dot(Vec,Obj.Vector2)/Obj.PixelScale);
                                Win_Crns2D(f,:)=[CRX,CRY];
                                Win3D=CR13D + Obj.Vector1*Win_Crns2D(f,1)*Obj.PixelScale + Obj.Vector2*Win_Crns2D(f,2)*Obj.PixelScale;
                                Win_Crns3D(f,:)=Win3D';
                            end
                            
                            InPnt=[min(Win_Crns2D(:,1)),min(Win_Crns2D(:,2))];
                            InSize=abs([max(Win_Crns2D(:,2))-min(Win_Crns2D(:,2)-1), max(Win_Crns2D(:,1))-min(Win_Crns2D(:,1)-1)]);
                            
                            % Keep track of the XD and YD elements:
                            SXD1=InPnt(2);
                            SXD2=(InPnt(2)+InSize(1)-1);
                            
                            SYD1=InPnt(1);
                            SYD2=(InPnt(1)+InSize(2)-1);
                            
                            XD=SXD1:SXD2;
                            YD=SYD1:SYD2;
                            
                            % Make sure the image is inside the bounds
                            XD1=min(XD);
                            XD2=max(XD);
                            YD1=min(YD);
                            YD2=max(YD);
                            Adjusted=0;
                            if min(XD)<1
                                XD1=1;
                                Adjusted=1;
                            end
                            if max(XD)>size(Obj.Image,1)
                                XD2=size(Obj.Image,1);
                                Adjusted=1;
                            end
                            
                            if min(YD)<1
                                YD1=1;
                                Adjusted=1;
                            end
                            if max(YD)>size(Obj.Image,2)
                                YD2=size(Obj.Image,2);
                                Adjusted=1;
                            end
                            
                            if Adjusted
                                XD=YD1:YD2; YD=XD1:XD2;
                            end
                            %%
                            DownSampleScale=1;
                            udata = [0 1];  vdata = [0 1];
                            OrgPnts(:,1)=OrgPnts(:,1)/max(OrgPnts(:,1))*(YD2-YD1+1);
                            OrgPnts(:,2)=OrgPnts(:,2)/max(OrgPnts(:,2))*(XD2-XD1+1);
                            tform = maketform('projective',[ 0 1;  1  1; 1 0; 0 0],round(OrgPnts/DownSampleScale));
                            tform=projective2d(tform.tdata.T);
                            
                            xdata=[min(GrdPnts(:,1)),max(GrdPnts(:,1))]; %'XData' and 'YData'
                            ydata=[min(GrdPnts(:,2)),max(GrdPnts(:,2))];
                            UVref=imref2d(size(ImagePlane.Image),udata,vdata);
                            outref=imref2d(round([XD2-XD1+1,YD2-YD1+1]/DownSampleScale),1,1);
                            
                            if (Obj.CorrectInsideParallax && ~isempty(Obj.DepthSegments))
                                % correct the inner parallax based on the
                                % previously calculated control points before applying the imwarp
%                                 Ip(:,:,1) = imwarp(uint8(imCorrect(ImagePlane.Image(:,:,1), CorrInBtwPnts2D', 3,0)*255),UVref,tform,'nearest','OutputView',outref);
%                                 Ip(:,:,1) = imwarp(uint8(imCorrect(ImagePlane.Image(:,:,1), CorrInBtwPnts2D', 3,0)*255),UVref,tform,'cubic','OutputView',outref);
                                  Ip(:,:,1) = imwarp(uint8(imCorrect(ImagePlane.Image(:,:,1), CorrInBtwPnts2D', 3,0)*255),UVref,tform,'linear','OutputView',outref);
                                
                            else
                                Ip(:,:,1) = imwarp(ImagePlane.Image(:,:,1),UVref,tform,'nearest','OutputView',outref);
                            end
                            Ip(isnan(Ip))=0; %replace NaNs with zeros:
                        end
                        
                        % Add the ImagePlaneIndex to the ImageSrcIndx:
                        if ~isempty(ImagePlaneIndex)
                            Window.ImageSrcIndx=uint8((Ip>0)*ImagePlaneIndex);
                        else
                            Window.ImageSrcIndx=[];
                        end
                        
                        Window.Image=Ip;
                        
                    end
                    
                end % end of check linearity of the points
                
            end
            
        end
        
        function [ghandles,Window,ErrorCode,XD,YD]=renderImage(Obj,ImagePlane, ~, ProcessImage, ShowGraphics, ~, ImagePlaneIndex)
            ghandles={};
            ErrorCode=0;
            XD=[];
            YD=[];
            ParallaxCorrected=0;
            if isempty(ShowGraphics)
                ShowGraphics=0;
            end
            
            % Establish whether a different imageprocessing needs to be
            % used in case there is parallax correction
            % In ProcessImage mode the Window.Coords are not used.
            if ProcessImage                
                if (Obj.CorrectInsideParallax && ~isempty(Obj.DepthSegments))
                    ParallaxCorrected=1;
                end
            end
                     
            if ParallaxCorrected
                
                [WinBounds,WinImOut]=Obj.intersectIP(ImagePlane,Obj.ParralaxCorrectionRes,ShowGraphics);
                
                if isempty(WinBounds)
                    Window=[];     % return an empty window if image is out of bounds
                    ErrorCode=1;   
                else
                    XD=WinBounds(3):WinBounds(4);
                    YD=WinBounds(1):WinBounds(2);
                    Ip(:,:,1) =WinImOut;
                    Ip(isnan(Ip))=0; %replace NaNs with zeros:
                    Window.Image=Ip;
                    if ~isempty(ImagePlaneIndex) % Add the ImagePlaneIndex to the ImageSrcIndx:
                        Window.ImageSrcIndx=uint8((Ip>0)*ImagePlaneIndex);
                    else
                        Window.ImageSrcIndx=[];
                    end                    
                end
                
            else
                
                Pp=Obj.Origin;
                Pv=-Obj.Vector3;
                Vz=ImagePlane.ImageNormal;
                Vy=ImagePlane.ImageUp;
                Vx=cross(Vy,Vz);
                Vx=Vx/norm(Vx);
                Vy=Vy/norm(Vy);
                
                O=ImagePlane.CentrePoint;
                PS=ImagePlane.PixelScale;
                Sx=1024;
                Sy=1024;
                ISFLD=0;
                if isfield(ImagePlane,'Image') || isprop(ImagePlane,'Image')
                    if ~isempty(ImagePlane.Image)
                        [Sy,Sx]=size(ImagePlane.Image(:,:,1));
                        ISFLD=1;
                    end
                end
                if ~ISFLD
                    ImagePlane.Image=uint8(zeros(Sx,Sy)); % Set the imageplane
                end
                
                % Finding the intersecting points:
                P1=O+(PS*Sy/2)*(Vy)-(PS*Sx/2)*Vx;
                P2=O+(PS*Sy/2)*(Vy)+(PS*Sx/2)*Vx;
                P3=O-(PS*Sy/2)*(Vy)+(PS*Sx/2)*Vx;
                P4=O-(PS*Sy/2)*(Vy)-(PS*Sx/2)*Vx;
                
                
                SRC=ImagePlane.SourcePosition;
                VV1=(P1-SRC);
                VV2=(P2-SRC);
                VV3=(P3-SRC);
                VV4=(P4-SRC);
                
                VV1=VV1/norm(VV1);
                VV2=VV2/norm(VV2);
                VV3=VV3/norm(VV3);
                VV4=VV4/norm(VV4);
                
                if isempty(Obj.DepthSegments)
                    % Simple cone/plane intersection without variable depth segments:
                    PP1=linePlaneIntersect(Pv,Pp,VV1,P1);
                    PP2=linePlaneIntersect(Pv,Pp,VV2,P2);
                    PP3=linePlaneIntersect(Pv,Pp,VV3,P3);
                    PP4=linePlaneIntersect(Pv,Pp,VV4,P4);                    
                else
                    % Account for the depth of the anatomy
                    [PP1,PI1]=lineAnatomyIntersect(Pv,Pp,VV1,P1,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                    [PP2,PI2]=lineAnatomyIntersect(Pv,Pp,VV2,P2,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                    [PP3,PI3]=lineAnatomyIntersect(Pv,Pp,VV3,P3,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                    [PP4,PI4]=lineAnatomyIntersect(Pv,Pp,VV4,P4,Obj.DepthSegments,Obj.DepthSegmentsNormal,ShowGraphics);
                end
                
                %Points of Intersections with the Grid plane:
                PPs=[PP1';PP2';PP3';PP4'];
                
                [OverallCheck,~]=Obj.checkPointsOut(PPs);
                
                if OverallCheck~=0 % if points stay out of the bounds:
                    Window=[];     % return an empty window
                    ErrorCode=1;   % Image Out of Bounds
                    
                else % carry on with normal processing:
                    V2=Obj.ImageUp;
                    V3=Obj.ImageNormal;
                    V1=cross(V2,V3);
                    V1=V1/norm(V1);
                    V2=V2/norm(V2);
                    V3=V3/norm(V3);
                    
                    Pnt1=PP1-Obj.Origin;
                    Pnt2=PP2-Obj.Origin;
                    Pnt3=PP3-Obj.Origin;
                    Pnt4=PP4-Obj.Origin;
                    
                    %Points with respect to the Grid [2D]:
                    Pr1=[dot(Pnt1,V1),dot(Pnt1,V2),dot(Pnt1,V3)];
                    Pr2=[dot(Pnt2,V1),dot(Pnt2,V2),dot(Pnt2,V3)];
                    Pr3=[dot(Pnt3,V1),dot(Pnt3,V2),dot(Pnt3,V3)];
                    Pr4=[dot(Pnt4,V1),dot(Pnt4,V2),dot(Pnt4,V3)];
                    
                    Prs=[Pr1;Pr2;Pr3;Pr4];
                    Prs_Grid=Prs; %All corners wrt the Grid
                    
                    % Found rectangular box around the image wrt the Grid CS:
                    Cr1_Grid=[min(Prs_Grid(:,1)),min(Prs_Grid(:,2))];
                    Cr2_Grid=[max(Prs_Grid(:,1)),min(Prs_Grid(:,2))];
                    Cr3_Grid=[min(Prs_Grid(:,1)),max(Prs_Grid(:,2))];
                    Cr4_Grid=[max(Prs_Grid(:,1)),max(Prs_Grid(:,2))];
                    Crs_Grid=[Cr1_Grid;Cr2_Grid;Cr3_Grid;Cr4_Grid];
                    
                    % Mapping the corners of the box to the global coordinate:
                    Cr1=Obj.Origin + Cr1_Grid(1)*V1 + Cr1_Grid(2)*V2;
                    Cr2=Obj.Origin + Cr2_Grid(1)*V1 + Cr2_Grid(2)*V2;
                    Cr3=Obj.Origin + Cr3_Grid(1)*V1 + Cr3_Grid(2)*V2;
                    Cr4=Obj.Origin + Cr4_Grid(1)*V1 + Cr4_Grid(2)*V2;
                    Crs=[Cr1';Cr2';Cr3';Cr4'];
                    
                    % Sorting the corners:
                    Indx=ones(4,1);
                    for f=1:size(PPs,1)
                        MinDst=norm(Crs(1,:)-PPs(f,:));
                        for g=1:size(Crs,1)
                            Dst=norm(Crs(g,:)-PPs(f,:));
                            if Dst<MinDst
                                MinDst=Dst;
                                Indx(f)=g;
                            end
                        end
                    end
                    
                    Crs_Grid=Crs_Grid(Indx,:);
                    
                    Sz=round(round(Cr4_Grid-Cr1_Grid)*Obj.PixelPerMM);
                    
                    SCFR=Sz(2)/Sx; % Scale factor calculated based on the size of the image.
                    MF=1/SCFR; % *2; %*10; % Magnification factor:
                    
                    Crs_Grid_p=Crs_Grid;
                    Crs_Grid_p(:,1)=Crs_Grid_p(:,1)-min(Crs_Grid_p(:,1));
                    Crs_Grid_p(:,2)=Crs_Grid_p(:,2)-min(Crs_Grid_p(:,2));
                    
                    Prs_Grid_p=Prs_Grid(:,1:2);
                    Prs_Grid_p(:,1)=Prs_Grid_p(:,1)-min(Prs_Grid_p(:,1));
                    Prs_Grid_p(:,2)=Prs_Grid_p(:,2)-min(Prs_Grid_p(:,2));
                    
                    
%                     GrdPnts=round(Crs_Grid_p*MF);
%                     OrgPnts=round(Prs_Grid_p*MF);

                    GrdPnts=round(Crs_Grid_p*MF);
                    OrgPnts=round(Prs_Grid_p*MF);
                    
                    
                    %Check if the points calculated are good for the rest of the calculation, if not abort:
                    IsCollinear=checkCollinearity(OrgPnts);
                    
                    if IsCollinear
                        Window=[];
                        ErrorCode=2;
                        
                    else
                        
                        Window.Origin=[min(Crs_Grid_p(:,1)),min(Crs_Grid_p(:,2))];
                        Window.ScaleFactor=SCFR;
                        %Calculating the pixel Scaling
                        Window.PixelPerMM=1/ImagePlane.PixelScale;
                        
                        
                        CRN1_2D=Crs_Grid(1,1:2);
                        CRN2_2D=Crs_Grid(2,1:2);
                        CRN3_2D=Crs_Grid(3,1:2);
                        CRN4_2D=Crs_Grid(4,1:2);
                        
                        % Corners of the image wrt the grid plane (2D):
                        CRNS_2D=[CRN1_2D; CRN2_2D; CRN3_2D; CRN4_2D];
                        Window.Coords=CRNS_2D;
                        
                        if ~ProcessImage % don't do anything with the image: normally used only for updating the window information.
                            Ip=ImagePlane.Image/255;
                        else
                            
                            %% Calculate the size of the image:
                            if isempty(Window) % out of bounds detected
                                ErrorCode=1;
                            else
                                
                            end
                            
                            if ErrorCode==0
                                
                                
                                %GridOrigin coordinates in 3d and in mm:
                                GridOrigin3D=Obj.CentrePoint + Obj.GridImageOriginOffset(1)*Obj.Vector2*Obj.PixelScale+ Obj.GridImageOriginOffset(2)*Obj.Vector1*Obj.PixelScale;
                                %Coordinates of the topleft corner with respect to the origin:
                                CR1x=-(size(Obj.Image,2)/2 + Obj.GridImageOriginOffset(1))*Obj.PixelScale;
                                CR1y=-(size(Obj.Image,1)/2 + Obj.GridImageOriginOffset(2))*Obj.PixelScale;
                                CR13D=GridOrigin3D+Obj.Vector1*CR1x+Obj.Vector2*CR1y;
                                Win_Crns2D=zeros(4,2);
                                Win_Crns3D=zeros(4,3);
                                for f=1:4
                                    P=Window.Coords(f,:);
                                    Pnt=Obj.CentrePoint+Obj.Vector1*P(1)+Obj.Vector2*P(2); % corner of the window
                                    Vec=Pnt-CR13D;
                                    CRX=round(dot(Vec,Obj.Vector1)/Obj.PixelScale);
                                    CRY=round(dot(Vec,Obj.Vector2)/Obj.PixelScale);
                                    Win_Crns2D(f,:)=[CRX,CRY]; % corners in pixel (in the image coordinate system)
                                    Win3D=CR13D + Obj.Vector1*Win_Crns2D(f,1)*Obj.PixelScale + Obj.Vector2*Win_Crns2D(f,2)*Obj.PixelScale;
                                    Win_Crns3D(f,:)=Win3D';
                                end
                                
                                InPnt=[min(Win_Crns2D(:,1)),min(Win_Crns2D(:,2))];
                                InSize=abs([max(Win_Crns2D(:,2))-min(Win_Crns2D(:,2)-1), max(Win_Crns2D(:,1))-min(Win_Crns2D(:,1)-1)]);
                                
                                % Keep track of the XD and YD elements:
                                SXD1=InPnt(2);
                                SXD2=(InPnt(2)+InSize(1)-1);
                                
                                SYD1=InPnt(1);
                                SYD2=(InPnt(1)+InSize(2)-1);
                                
                                XD=SXD1:SXD2;
                                YD=SYD1:SYD2;
                                
                                % Make sure the image is inside the bounds
                                XD1=min(XD);
                                XD2=max(XD);
                                YD1=min(YD);
                                YD2=max(YD);
                                Adjusted=0;
                                if min(XD)<1
                                    XD1=1;
                                    Adjusted=1;
                                end
                                if max(XD)>size(Obj.Image,1)
                                    XD2=size(Obj.Image,1);
                                    Adjusted=1;
                                end
                                
                                if min(YD)<1
                                    YD1=1;
                                    Adjusted=1;
                                end
                                if max(YD)>size(Obj.Image,2)
                                    YD2=size(Obj.Image,2);
                                    Adjusted=1;
                                end
                                
                                if Adjusted
                                    XD=YD1:YD2; YD=XD1:XD2;
                                end
                                
                                DownSampleScale=1;
                                udata = [0 1];  vdata = [0 1];
                                OrgPnts(:,1)=OrgPnts(:,1)/max(OrgPnts(:,1))*(YD2-YD1+1);
                                OrgPnts(:,2)=OrgPnts(:,2)/max(OrgPnts(:,2))*(XD2-XD1+1);
                                tform = maketform('projective',[ 0 1;  1  1; 1 0; 0 0],round(OrgPnts/DownSampleScale));
                                tform=projective2d(tform.tdata.T);
                                
                                xdata=[min(GrdPnts(:,1)),max(GrdPnts(:,1))]; %'XData' and 'YData'
                                ydata=[min(GrdPnts(:,2)),max(GrdPnts(:,2))];
                                UVref=imref2d(size(ImagePlane.Image),udata,vdata);
                                outref=imref2d(round([XD2-XD1+1,YD2-YD1+1]/DownSampleScale));
                                tic;
                                Ip(:,:,1) = imwarp(ImagePlane.Image(:,:,1),UVref,tform,'linear','OutputView',outref);
                                % Ip(:,:,1) = imwarp(ImagePlane.Image(:,:,1),UVref,tform,'nearest','OutputView',outref);
                                toc;
                                Ip(isnan(Ip))=0; %replace NaNs with zeros:
                            end
                            
                            % Add the ImagePlaneIndex to the ImageSrcIndx:
                            if ~isempty(ImagePlaneIndex)
                                Window.ImageSrcIndx=uint8((Ip>0)*ImagePlaneIndex);
                            else
                                Window.ImageSrcIndx=[];
                            end
                            
                        end
                        Window.Image=Ip;
                        
                    end % end of check linearity of the points
                    
                end
                
            end
            
        end
        
        function [Image,Origin,ImageSize]=stitchWindows(Obj,varargin)
            P=[];
            if ~isempty(varargin)
                for f=1:numel(varargin{1})
                    S=10; % 5 pixels per mm:
                    if f==1
                        P=varargin{1}{f}.Coords(1,:)*S;
                        i=2;
                        j=1;
                        P(2,:)=varargin{1}{f}.Coords(2,:)*S;
                        P(3,:)=varargin{1}{f}.Coords(3,:)*S;
                        P(4,:)=varargin{1}{f}.Coords(4,:)*S;
                    else
                        P(f*4-3,:)=varargin{1}{f}.Coords(1,:)*S;
                        P(f*4-2,:)=varargin{1}{f}.Coords(2,:)*S;
                        P(f*4-1,:)=varargin{1}{f}.Coords(3,:)*S;
                        P(f*4,:)=varargin{1}{f}.Coords(4,:)*S;
                    end
                end
                P(2,:)=1*P(2,:);
                Cr1=int32(round(max(P(:,1))-min(P(:,1))));
                Cr2=int32(round(max(P(:,2))-min(P(:,2))));
                CrD=[Cr1,Cr2];
                
                Cr1=int32(round(min(P(:,1))));
                Cr2=int32(round(min(P(:,2))));
                Offset=[Cr1,Cr2];
                
                Im(CrD(1)+10,CrD(2)+10)=0;
                CNT=[10,10];
                
                REF=[0,0];
                REF_XA=[50,0];
                REF_YA=[0,50];
                
                for f=1:numel(varargin{1})
                    W=varargin{1}{f};
                    
                    SX=size(W.Image,1);
                    SY=size(W.Image,2);
                    
                    DX=int32(max(varargin{1}{f}.Coords(:,1))*S) - int32(min(varargin{1}{f}.Coords(:,1))*S);
                    DY=int32(max(varargin{1}{f}.Coords(:,2))*S) - int32(min(varargin{1}{f}.Coords(:,2))*S);
                    
                    OX=int32(min(varargin{1}{f}.Coords(:,1))*S);
                    OY=int32(min(varargin{1}{f}.Coords(:,2))*S);
                    
                    ReImage=imresize(flipud((W.Image)),[DX,DY]);
                    
                    XD=CNT(1)-Offset(1)+OX+1:CNT(1)-Offset(1)+OX+DX;
                    YD=CNT(2)-Offset(2)+OY+1:CNT(2)-Offset(2)+OY+DY;
                    
                    REFX=CNT(1)-Offset(1)+1+int32(REF(2)*S);
                    REFY=CNT(2)-Offset(2)+1+int32(REF(1)*S);
                    
                    REFX_XA=CNT(1)-Offset(1)+1+int32(REF_XA(2)*S);
                    REFY_XA=CNT(2)-Offset(2)+1+int32(REF_XA(1)*S);
                    
                    REFX_YA=CNT(1)-Offset(1)+1+int32(REF_YA(2)*S);
                    REFY_YA=CNT(2)-Offset(2)+1+int32(REF_YA(1)*S);
                    
                    clear MM;
                    ROI=Im(XD,YD);
                    
                    Msk1=ROI*0;
                    Msk2=Msk1;
                    OVL=ROI*0+1;
                    Msk1(ROI>0)=1;
                    Msk2(ReImage>0)=1;
                    Wght=Msk1+Msk2;
                    OVL(Wght==2)=1;
                    NOVR(Wght<2)=1;
                    
                    ROI2=max(Msk1.*ROI, Msk2.*ReImage);
                    OVERLAP=min(ROI,ReImage).*OVL;
                    OVERLAP(OVERLAP==1)=0;
                    ROI3=max(ROI2,OVERLAP);
                    Im(XD, YD)=(ROI3);
                end
                %Store the .Image property
                Image=Im;
                Origin=[REFX,REFY;REFX_XA,REFY_XA;REFX_YA,REFY_YA];
                ImageSize=size(Im);
            end
            
        end
        
        
        function [Image,Origin,ImageSize]=stitchWindows_Old(Obj,varargin)
            P=[];
            if ~isempty(varargin)
                for f=1:numel(varargin{1})
                    S=varargin{1}{f}.PixelPerMM / varargin{1}{f}.ScaleFactor;
                    if f==1
                        P=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset;
                        i=2;
                        j=1;
                        P(2,:)=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset+[size(varargin{1}{f}.Image,i),size(varargin{1}{f}.Image,j)]-[1,1];
                        P(3,:)=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset+[size(varargin{1}{f}.Image,i),1]-[1,1];
                        P(4,:)=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset+[1,size(varargin{1}{f}.Image,j)]-[1,1];
                    else
                        P(f*4-3,:)=varargin{1}{f}.Origin*S++varargin{1}{f}.Offset;
                        P(f*4-2,:)=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset+[size(varargin{1}{f}.Image,i),size(varargin{1}{f}.Image,j)]-[1,1];
                        P(f*4-1,:)=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset+[size(varargin{1}{f}.Image,i),1]-[1,1];
                        P(f*4,:)=varargin{1}{f}.Origin*S+varargin{1}{f}.Offset+[1,size(varargin{1}{f}.Image,j)]-[1,1];
                    end
                    figure;
                    imshow(varargin{1}{f}.Image);
                end
                Cr1=int32([min(P(:,1)),min(P(:,2))]);
                Cr2=int32([max(P(:,1)),max(P(:,2))]);
                CrD=[Cr1,Cr2];
                
                Im(CrD(1)+1000,CrD(2)+1000)=0;
                
                for f=1:numel(varargin{1})
                    W=varargin{1}{f};
                    Org=int32(W.Origin*S)+int32(W.Offset)+500;
                    SX=size(W.Image,2);
                    SY=size(W.Image,1);
                    XD=Org(1)-Cr1(1)+1:Org(1)-Cr1(1)+size(W.Image,1);
                    YD=Org(2)-Cr1(2)+1:Org(2)-Cr1(2)+size(W.Image,2);
                    
                    clear MM;
                    ROI=Im(XD,YD);
                    
                    Msk1=ROI*0;
                    Msk2=Msk1;
                    OVL=ROI*0+1;
                    Msk1(ROI>0)=1;
                    Msk2(W.Image>0)=1;
                    Wght=Msk1+Msk2;
                    OVL(Wght==2)=1;
                    NOVR(Wght<2)=1;
                    
                    ROI2=max(Msk1.*ROI, Msk2.*W.Image);
                    OVERLAP=min(ROI,W.Image).*OVL;
                    OVERLAP(OVERLAP==1)=0;
                    ROI3=max(ROI2,OVERLAP);
                    Im(XD, YD)=(ROI3);
                end
                %Store the .Image property
                Image=Im;
                Origin=Org;
                ImageSize=size(Im);
            end
            
        end
        
        % plot the image in 3D
        function surf = imshow3d(Obj,ShowAxis,ShowGrid,HighQuality,ShowCropped)
            
            if ShowCropped && isempty(Obj.CroppedCoords)
                ShowCropped=0;
                disp('No crooped info found!');
            end
            
            if ShowGrid
                if ShowCropped
                    imageData= Obj.CroppedImage; % Add the Grid later if required:
                else
                    imageData = Obj.ImageOutput;
                end
            else
                if ShowCropped
                    imageData= Obj.CroppedImage;
                else
                    imageData = Obj.Image;
                end
            end
            
            if size(imageData,3)==4
                imageData=imageData(:,:,1:3);
            end
            
            if ShowCropped
                Cnt=Obj.returnCroppedCentre;
                if ~isempty(Cnt)
                    CentPoint=Cnt;
                else
                    disp('No Cropped Info found');
                end
            else
                CentPoint=Obj.CentrePoint;
            end
            
            
            PixelScale=Obj.PixelScale;
            
            ImageNormal=Obj.ImageNormal;
            ImageUp=Obj.ImageUp;
            
            Vy=ImageUp;
            Vz=ImageNormal;
            Vx=cross(Vy,Vz);
            Vx=Vx/norm(Vx);
            Vy=Vy/norm(Vy);
            Vz=Vz/norm(Vz);
            
            
            Dx=size(imageData,2)*PixelScale;
            Dy=size(imageData,1)*PixelScale;
            
            P1=CentPoint-Vx*(Dx/2)-Vy*(Dy/2);
            P2=CentPoint+Vx*(Dx/2)-Vy*(Dy/2);
            P3=CentPoint-Vx*(Dx/2)+Vy*(Dy/2);
            P4=CentPoint+Vx*(Dx/2)+Vy*(Dy/2);
            
            Xdata=[P1(1),P2(1);P3(1),P4(1)];
            Ydata=[P1(2),P2(2);P3(2),P4(2)];
            Zdata=[P1(3),P2(3);P3(3),P4(3)];
            
            colormap(bone);
            
            % Reduce the image size to 1/10 for faster visulization
            if HighQuality
                RimageData=imageData;
            else
                Sz2=[size(imageData,1)/10,size(imageData,2)/10];
                RimageData=imresize(imageData,Sz2);
            end
            
            surf = surface('XData',Xdata,'YData',Ydata,...
                'ZData',Zdata,'CData',RimageData,...
                'FaceColor','texturemap','EdgeColor','none','FaceAlpha',0.95);
            
            if ShowAxis
                hold on;
                
                Vy=Obj.ImageUp;
                Vz=Obj.ImageNormal;
                Vx=cross(Vy,Vz);
                Vx=Vx/norm(Vx);
                Vy=Vy/norm(Vy);
                Vz=Vz/norm(Vz);
                
                zdir=CentPoint+(Vz*100);
                xdir=CentPoint+(Vx*100);
                ydir=CentPoint+(Vy*100);
                plot3([CentPoint(1),zdir(1)],[CentPoint(2),zdir(2)],[CentPoint(3),zdir(3)],'b','LineWidth',3);
                plot3([CentPoint(1),xdir(1)],[CentPoint(2),xdir(2)],[CentPoint(3),xdir(3)],'r','LineWidth',3);
                plot3([CentPoint(1),ydir(1)],[CentPoint(2),ydir(2)],[CentPoint(3),ydir(3)],'g','LineWidth',3);
                
                if ~isempty(Obj.AnatomySegments3D)
                    plot3(Obj.AnatomySegments3D(:,1),Obj.AnatomySegments3D(:,2),Obj.AnatomySegments3D(:,3),'r');
                    plot3(Obj.AnatomySegments3D(:,1),Obj.AnatomySegments3D(:,2),Obj.AnatomySegments3D(:,3),'y+');
                end
            end
            
        end
        
        
        % Plot the image in 3D
        function ghandles = showGraphically(Obj,ShowAxis,ShowGrid,HighQuality,ShowCropped,Layer)
            
            ghandles={};
            
            if ShowCropped && isempty(Obj.CroppedCoords)
                ShowCropped=0;
                disp('No crooped info found!');
            end
            
            
            switch Layer
                case 1 % .Image
                    imageData=Obj.Image;
                case 2 % .PreopImage
                    %                     imageData=Obj.PreopImage;
                    % to Make it work with Stthe TCarmScript we consider
                    % empty image:
                    imageData=Obj.Image*0;
                case 3 % .ImageOutput
                    imageData=Obj.Image; %Output;
            end
            
            if ShowCropped
                CC=Obj.CroppedCoords;
                % img=Obj.Image(CC(1):CC(2),CC(3):CC(4));
                switch size(imageData,3)
                    case 1
                        imageData=imageData(CC(1):CC(2),CC(3):CC(4));
                    case 2
                        imageData=imageData(CC(1):CC(2),CC(3):CC(4),1:3);
                    case 3
                        imageData=imageData(CC(1):CC(2),CC(3):CC(4),1:3);
                end
            end
            
            if size(imageData,3)==4
                imageData=imageData(:,:,1:3);
            end
            
            if ShowCropped
                Cnt=Obj.returnCroppedCentre;
                if ~isempty(Cnt)
                    CentPoint=Cnt;
                else
                    disp('No Cropped Info found');
                end
            else
                CentPoint=Obj.CentrePoint;
            end
            
            
            PixelScale=Obj.PixelScale;
            ImageNormal=Obj.ImageNormal;
            ImageUp=Obj.ImageUp;
            Vy=ImageUp;
            Vz=ImageNormal;
            Vx=cross(Vy,Vz);
            Vx=Vx/norm(Vx);
            Vy=Vy/norm(Vy);
            Vz=Vz/norm(Vz);
            
            Dx=size(imageData,2)*PixelScale;
            Dy=size(imageData,1)*PixelScale;
            
            P1=CentPoint-Vx*(Dx/2)-Vy*(Dy/2);
            P2=CentPoint+Vx*(Dx/2)-Vy*(Dy/2);
            P3=CentPoint-Vx*(Dx/2)+Vy*(Dy/2);
            P4=CentPoint+Vx*(Dx/2)+Vy*(Dy/2);
            
            Xdata=[P1(1),P2(1);P3(1),P4(1)];
            Ydata=[P1(2),P2(2);P3(2),P4(2)];
            Zdata=[P1(3),P2(3);P3(3),P4(3)];
            
            colormap(bone);
            
            % Reduce the image size to 1/10 for faster visulization
            if HighQuality
                RimageData=imageData;
            else
                Sz2=[size(imageData,1)/10,size(imageData,2)/10];
                RimageData=imresize(imageData,Sz2);
            end
            surf = surface('XData',Xdata,'YData',Ydata,...
                'ZData',Zdata,'CData',RimageData,...
                'FaceColor','texturemap','EdgeColor','none','FaceAlpha',0.95);
            
            if ShowAxis
                hold on;
                
                Vy=Obj.ImageUp;
                Vz=Obj.ImageNormal;
                Vx=cross(Vy,Vz);
                Vx=Vx/norm(Vx);
                Vy=Vy/norm(Vy);
                Vz=Vz/norm(Vz);
                
                zdir=CentPoint+(Vz*100);
                xdir=CentPoint+(Vx*100);
                ydir=CentPoint+(Vy*100);
                h1=plot3([CentPoint(1),zdir(1)],[CentPoint(2),zdir(2)],[CentPoint(3),zdir(3)],'b','LineWidth',3);
                h2=plot3([CentPoint(1),xdir(1)],[CentPoint(2),xdir(2)],[CentPoint(3),xdir(3)],'r','LineWidth',3);
                h3=plot3([CentPoint(1),ydir(1)],[CentPoint(2),ydir(2)],[CentPoint(3),ydir(3)],'g','LineWidth',3);
            end
            
            ghandles={surf,h1,h2,h3};
        end
        
        function [Img]=returnCroppedImage(Obj)
            CC=Obj.CroppedCoords;
            Img=Obj.Image(CC(1):CC(2),CC(3):CC(4));
        end
        
        function [Centre]=returnCroppedCentre(Obj)
            
            Centre=[];
            IC=[size(Obj.Image,1),size(Obj.Image,2)]/2;
            %SZC=Obj.CroppedImage;
            C=Obj.CroppedCoords;
            UP=Obj.ImageUp;
            NR=Obj.ImageNormal;
            RT=cross(UP,NR);
            RT=RT/norm(RT);
            PixelScale=Obj.PixelScale;
            if isempty(C)
                disp('No Crop Info.');
            else
                CC=[(C(1)+C(2))/2,(C(3)+C(4))/2]; % CroppedCentre
                Dist=CC-IC;
                Dist=Dist*PixelScale;
                switch Obj.PlaneType
                    case 'XY'
                        Dist3D=Dist(1)*UP+Dist(2)*RT;
                    case 'YZ'
                        Dist3D=+Dist(1)*UP+Dist(2)*RT;
                end
                Centre=Obj.CentrePoint+Dist3D;
            end
        end
        
        function [img] = returnVisImage(Obj,options,ProcessImage)
            
            if (isempty(Obj.SaveVis2D_CroppedImage) || ...
                    isempty(Obj.SaveVis2D_CroppedImage) || ...
                    isempty(Obj.SaveVis2D_CroppedImage))
                ProcessImage=1;
            end
            
            if ProcessImage
                
                switch options
                    case 1
                        % Croppedimage
                        img = Obj.CroppedImage;
                    case 2
                        % ImageOutput
                        img = Obj.ImageOutput;
                        %img = cat(3,Obj.Image,Obj.Image,Obj.Image);
                    case 3
                        n=1;
                        img=Obj.Image;
                    case 4
                        % Image
                        CC=Obj.CroppedCoords;
                        II=Obj.Image(CC(1):CC(2),CC(3):CC(4));
                        img(:,:,1)=II;
                        img(:,:,2)=II;
                        img(:,:,3)=II;
                    case 5
                        % ImageOutput
                        % removed:
                        img=Obj.Image;
                        
                end
                
                % DownSample
                
                n=Obj.DownSampleSize;
                
                % Decide whether we will flip the image or not
                switch Obj.Vis2D(2)
                    case 1
                        % Flip the normal
                        img = flipdim(imrotate(img(1:n:end,1:n:end,:),Obj.Vis2D(1)),2);
                        %                         img = imrotate(img(1:n:end,1:n:end,:),Obj.Vis2D(1));
                    case 0
                        % Don't flip
                        img = imrotate(img(1:n:end,1:n:end,:),Obj.Vis2D(1));
                end
                
                
                %store the processed images for later use:
                
                switch options
                    case 1
                        % Croppedimage
                        Obj.SaveVis2D_CroppedImage=img;
                    case 2
                        % ImageOutput
                        Obj.SaveVis2D_ImageOutput=img;
                    case 3
                        % Image
                        Obj.SaveVis2D_Image=img;
                    case 4
                        Obj.SaveVis2D_CroppedImage_Bare=img;
                end
                
                %                 pack;
                
            else % if not ProcessImage:
                
                switch options
                    case 1
                        % Croppedimage
                        img = Obj.SaveVis2D_CroppedImage;
                    case 2
                        % ImageOutput
                        img = Obj.SaveVis2D_ImageOutput;
                    case 3
                        % Image
                        img=Obj.SaveVis2D_Image;
                    case 4
                        % Croppedimage wo the Grid
                        img=Obj.SaveVis2D_CroppedImage_Bare;
                end
                
            end
            
        end
        
        function [CropCoords,ImageCoords]=VisToImageCoords(Obj,VisCoords,hSurfHandle)
            XD=get(hSurfHandle,'XData');
            YD=get(hSurfHandle,'YData');
            ImgXdim=size(Obj.Image,2)/2; % assuming that the Img lines up with the surface patch
            ImgYdim=size(Obj.Image,1)/2;
            Xdim=abs(XD(1,1));
            Ydim=abs(YD(1,1));
            XC=VisCoords(1);
            YC=-VisCoords(2);
            ImageCoords=[XC*ImgXdim/Xdim;YC*ImgYdim/Ydim]+[ImgXdim;ImgYdim];
            CropCoords=ImageCoords-[Obj.CroppedCoords(3);Obj.CroppedCoords(1)];
        end
        
        
        function [CropCoords,ImageCoords]=VisToImageCoords_OLD(Obj,VisCoords)
            
            VisCoords=[-VisCoords(2),VisCoords(1)]; %*Obj.DownSampleSize;
            Img=Obj.returnVisImage(3,1);
            CrpImage = Obj.returnVisImage(1,1);
            Rot=Obj.Vis2D(1);
            Mir=Obj.Vis2D(2);
            
            Patch_Cent=VisCoords;
            
            % Centre of returnVisImage
            VisImg_Cent=[size(Img,1)/2,size(Img,2)/2]*Obj.DownSampleSize;
            CrpImg_Cent=[size(CrpImage,1)/2,size(CrpImage,2)/2];
            
            %Patch Centre wrt CroppedImageCentre
            PC_wrt_VisImg=Patch_Cent-VisImg_Cent;
            
            %flip patch centre vertically (wrt the original orientation)
            if Mir
                VisCoords(2) = -VisCoords(2);
            end
            
            % Now we need to rotate the centre point by the same angle but in the
            % opposite direction as returnVisImage
            
            alpha = -Rot*pi/180;
            R = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
            PC_transformed = R*VisCoords';
            PC_wrt_corner = PC_transformed + VisImg_Cent';
            CropCoords=PC_wrt_corner;
            ImageCoords=CropCoords+[Obj.CroppedCoords(1);Obj.CroppedCoords(3)];
            
        end
        
        
        
        function [CropCoords,ImageCoords]=VisToImageCoords_OBSOLETE(Obj,VisCoords)
            
            VisCoords=[-VisCoords(2),VisCoords(1)]*Obj.DownSampleSize;
            Img=Obj.returnVisImage(1,0);
            % Include preop in cropping
            CrpImage = Obj.CroppedImage;
            Rot=Obj.Vis2D(1);
            Mir=Obj.Vis2D(2);
            Patch_Cent=VisCoords;
            % Centre of returnVisImage
            VisImg_Cent=[size(Img,1)/2,size(Img,2)/2]*Obj.DownSampleSize;
            CrpImg_Cent=[size(CrpImage,1)/2,size(CrpImage,2)/2];
            
            %Patch Centre wrt CroppedImageCentre
            PC_wrt_VisImg=Patch_Cent-VisImg_Cent;
            
            %flip patch centre vertically (wrt the original orientation)
            if Mir
                VisCoords(2) = -VisCoords(2);
            end
            
            % Now we need to rotate the centre point by the same angle but in the opposite direction as returnVisImage
            alpha = -Rot*pi/180;
            R = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
            PC_transformed = R*VisCoords';
            PC_wrt_corner = PC_transformed + CrpImg_Cent';
            CropCoords=PC_wrt_corner;
            ImageCoords=CropCoords+[Obj.CroppedCoords(1);Obj.CroppedCoords(3)];
            
        end
        
        
        function [VisCoords]=ImageToVisCoords(Obj,ImageCoords,hSurfHandle)
            % ImageCoords defined wrt the centre of the image
            Pnts=[0,0;10,0;0,10];
            PIs=[];
            PJs=[];
            
            for f=1:3
                Pi=Pnts(f,:);
                [PjC,Pj]=Obj.VisToImageCoords(Pi,hSurfHandle); %+CrpImg_Cent);
                PIs=[PIs;Pi,1];
                PJs=[PJs;Pj',1];
            end
            T=PIs'*inv(PJs');
            
            Out=T*[size(Obj.Image,2)/2+ImageCoords(1),size(Obj.Image,1)/2+ImageCoords(2),1]';
            VisCoords=[Out(1),Out(2)];
        end
        
        function [Img]= returnVisImageFullSize(Obj)
            
            CC=Obj.CroppedCoords;
            Img=Obj.Image(CC(1):CC(2),CC(3):CC(4));
            
        end
        
        function refreshImageOutput(Obj)
            %             Obj.fuseLayers;
            %             Obj.ImageOutput=cat(3,Obj.Image,Obj.Image,Obj.Image);
        end
        
        function [OverallCheck,PointsStatus]=checkPointsOut(Obj,Points)
            %check if the points listed are inside the imageplane
            %First translate points to a 2D version
            V1=Obj.Vector1;
            V2=Obj.Vector2;
            
            OverallCheck=0;
            PointsStatus=zeros(size(Points,1),1);
            for f=1:size(Points,1)
                xconst=0;
                yconst=0;
                Vec=Points(f,:)'-Obj.Origin;
                Px=dot(Vec,V1);
                Py=dot(Vec,V2);
                if abs(Px)<(size(Obj.Image,2)/2*Obj.PixelScale) % check the X boundary
                    xconst=1;
                end
                if abs(Py)<(size(Obj.Image,1)/2*Obj.PixelScale) % check the Y boundary
                    yconst=1;
                end
                if xconst && yconst
                    PointsStatus(f)=1;
                else
                    PointsStatus(f)=0;
                    OverallCheck=1;
                end
            end
            
            if OverallCheck==0 % if points exists within the bounds check for colliearity of the points
                IsCollinear=checkCollinearity(Points);
                if IsCollinear
                    OverallCheck=2;
                end
            end
            
        end
        
        % apply the anatomy segnments to the stitchedplane
        function applyAnatomySegments(Obj,OrigPoints)
            
            Obj.AnatomySegments=[];
            if isempty(OrigPoints)
                Points=[];
            else
                [Points]=interpDevPoints(OrigPoints,50);
%                 [Points]=interpDevPoints2(OrigPoints,Obj.PlaneType,50);
            end
            
            if strcmpi(Obj.PlaneType,'XY')
                if isempty(Points)
                    Obj.AnatomySegments=[];
                else
                    % AP Plane: sort the points based on the Y coordinates:
                    SPnts=sortrows(Points,2);
                    % add two extensions to the end points to cover the entire length:
                    Y1=SPnts(1,1);
                    X1=SPnts(1,2);
                    Y2=SPnts(end,1);
                    X2=SPnts(end,2);
                    dX=Obj.GridImageSize(2)*Obj.PixelScale;
                    if -dX<X2
                        Pnt1=[Y1,-dX];
                    else
                        Pnt1=[];
                    end
                    if dX>X1
                        Pnt2=[Y2,dX];
                    else
                        Pnt2=[];
                    end
                    Obj.AnatomySegments=[Pnt1;SPnts;Pnt2];
                end
            else
                if isempty(Points)
                    Obj.AnatomySegments=[];
                else
                    % ML Plane: sort the opiints based on the X coordinates:
                    SPnts=sortrows(Points,2);
                    % add two extensions to the end points to cover the entire length:
                    Y1=SPnts(1,1);
                    X1=SPnts(1,2);
                    Y2=SPnts(end,1);
                    X2=SPnts(end,2);
                    dX=Obj.GridImageSize(2)*Obj.PixelScale;
                    if -dX<X2
                        Pnt1=[Y1,-dX];
                    else
                        Pnt1=[];
                    end
                    if dX>X1
                        Pnt2=[Y2,dX];
                    else
                        Pnt2=[];
                    end
                    Obj.AnatomySegments=[Pnt1;SPnts;Pnt2];
                end
            end
            
        end
        
        % Calculate the DepthSegment layer
        function calcImageDepthIndex(Obj,ShowGraphics)
            
            tic;
            if (isempty(Obj.DepthSegments) || isempty(Obj.DepthSegmentsNormal)) %check for valid depth information
                Obj.ImageDepthIndx=Obj.Image*0;
                fprintf('\r\nDepth Segments not included in the image.');
            else
                Obj.ImageDepthIndx=Obj.Image*0+1; % Just consider a constant number (The client does not need this info anymore)
            end
            ti=toc;
            fprintf('\r\n Time for processing Depth: %3.3g\r\n',ti);
            
            
        end
        
        % Calculate the DepthSegment layer
        function calcImageDepthIndex_Obsolete(Obj,ShowGraphics)
            
            tic;
            if (isempty(Obj.DepthSegments) || isempty(Obj.DepthSegmentsNormal)) %check for valid depth information
                Obj.ImageDepthIndx=Obj.Image*0;
                fprintf('\r\nDepth Segments not included in the image.');
            else
                P2Ds=calcDepthPrj2D(Obj);
                n=2; % slightly rougher increments for faster processing
                [X,Y] = meshgrid(1:round(size(Obj.Image,1)/n),1:round(size(Obj.Image,2)/n));
                F=X*0;
                Imo=size((Obj.Image)')/2;
                NImo=Imo/n;
                
                for f=1:numel(P2Ds)
                    P2D=P2Ds{f};
                    
                    Pco=[P2D(3,1),P2D(3,2)];
                    Pao=[P2D(4,1),P2D(4,2)];
                    Pbo=[P2D(1,1),P2D(1,2)];
                    Pdo=[P2D(2,1),P2D(2,2)];
                    
                    Pa=(Pao-Imo)/n+NImo;
                    Pb=(Pbo-Imo)/n+NImo;
                    Pc=(Pco-Imo)/n+NImo;
                    Pd=(Pdo-Imo)/n+NImo;
                    
                    m1=(Pb(1)-Pa(1))/(Pb(2)-Pa(2));
                    if m1<0
                        m1=-1*m1;
                    end
                    F=F+((sign((Y-Pb(1))+m1*(X-Pb(2)))+1)/2);
                    
                end
                if n>1
                    Obj.ImageDepthIndx=size(Obj.DepthSegments,1)-(uint8(imresize(F,[size(Obj.Image,2),size(Obj.Image,1)],'nearest'))+1)';
                else
                    Obj.ImageDepthIndx=size(Obj.DepthSegments,1)-(uint8(F)+1)';
                end
                
                if (ShowGraphics)
                    CatImg=cat(3,Obj.Image,Obj.ImageSrcIndx,Obj.ImageDepthIndx);
                    
                    figure; hold on;
                    imshow(CatImg);hold on;
                    
                    for f=2:numel(P2Ds)-1
                        P2D=P2Ds{f};
                        Pco=[P2D(3,1),P2D(3,2)];
                        Pao=[P2D(4,1),P2D(4,2)];
                        Pbo=[P2D(1,1),P2D(1,2)];
                        Pdo=[P2D(2,1),P2D(2,2)];
                        Pl=[Pao;Pbo];
                        plot(Pl(:,1),Pl(:,2),'w');
                        Pl=[Pco;Pdo];
                        plot(Pl(:,1),Pl(:,2),'w');
                    end
                end
                
            end
            ti=toc;
            fprintf('\r\n Time for processing Depth: %3.3g\r\n',ti);
            
        end
        
        % convert any 3D point from global to local coordinates in mm
        function localPointMM=localMM(Obj, Pnt3D, ProjectBool)
            % Pnt3D: 3XN points in global coordinates (in mm)
            % localPointMM: 3XN points in local coordinates (in mm)
            % ProjectBool: Project the points on the plane if set to 1
            if (size(Pnt3D,2)==3) && ~(size(Pnt3D,1)==3)
                Pnt3D=Pnt3D';
            end
            localPointMM=Pnt3D*0;
            for f=1:size(Pnt3D,2)
                V=Pnt3D(:,f)-Obj.Grid.Origin;
                dx=dot(V,Obj.Vector1);
                dy=dot(V,Obj.Vector2);
                dz=dot(V,Obj.Vector3);
                localPointMM(:,f)=[dx,dy,dz]';
            end
            if (ProjectBool)
                localPointMM(3,:)=0;
            end
        end
        
        % convert any 3D point from global to local coordinates in pixels
        function localPointPXL=localPXL(Obj,Pnt3D)
            % Pnt3D: 3XN points in global coordinates (in mm)
            % localPointPXL: 2XN points in local Image (in Pixels)
            localMM=Obj.localMM(Pnt3D,1)/Obj.PixelScale;
            Pnt=localMM(1:2,:);
            Org=round([size(Obj.Image,2),size(Obj.Image,1)]/2);
            Px=Org(1)+Pnt(1,:);
            Py=Org(2)+Pnt(2,:);
            localPointPXL=[Px;Py];
            
            
        end
        
        % convert any pixel coordinate to 3D local ImagePlane and global coordinates
        function [localPnts3D,globalPnts3D]=PXLtoPnt3D(Obj,PntsPXL)
            % Pnt3D: 3XN points in global coordinates (in mm)
            % localPointPXL: 2XN points in local Image (in Pixels)
            localPnts3D=zeros(3,size(PntsPXL,2));
            globalPnts3D=zeros(3,size(PntsPXL,2));
            Org=round([size(Obj.Image,2),size(Obj.Image,1)]/2);
            localPnts3D(1,:)=PntsPXL(1,:)-Org(1);
            localPnts3D(2,:)=PntsPXL(2,:)-Org(2);
            for f=1:size(localPnts3D,2)
                globalPnts3D(:,f)=Obj.Origin + localPnts3D(1,f)*Obj.Vector1 + localPnts3D(2,f)*Obj.Vector2;
            end
        end
        
        % Return the pixel bounding box from four points on the plane
        function [Bounds,CPoints,NCPoints,CPoints3D]=returnPXLBounds(Obj,Pnts3D,UnitSize)
            % Inputs:
            % ---------
            % Pnts3D: 3xN or Nx3 3D points in 3D
            % UnitSize: Number of pixels to be used as the UnitSize. Example: 10
            % Outputs:
            % ---------
            % Bounds: Bounding Box on ImagePlane that contains the projection of the Pnts3D
            % CPoints: Control Points on the ImagePlane (in Pixel Coordinates)
            % NCPoints: Normal Control Points on the ImagePlane (between 0 and 1 in each direction)
            % CPoints3D: 3D coordinate of the control points
            Org2D=round([size(Obj.Image,2),size(Obj.Image,1)]/2);
            
            Pnts2D=Obj.localPXL(Pnts3D);
            MinX=(floor(min(Pnts2D(1,:))/UnitSize)-1)*UnitSize;
            MaxX=(ceil(max(Pnts2D(1,:))/UnitSize)+1)*UnitSize;
            MinY=(floor(min(Pnts2D(2,:))/UnitSize)-1)*UnitSize;
            MaxY=(ceil(max(Pnts2D(2,:))/UnitSize)+1)*UnitSize;
           
            OutofBounds=0;

            if (MinX<1)
                if (MaxX>1)
                    MinX=1;
                else
                    OutofBounds=1;
                end
            end
            if MaxX>size(Obj.Image,2)
                if MinX<=size(Obj.Image,2)
                    MaxX=size(Obj.Image,2);
                else
                    OutofBounds=1;
                end
            end
            
            if MinY<1
                if(MaxY>1)
                    MinY=1;
                else
                    OutofBounds=1;
                end
            end
            if MaxY>size(Obj.Image,1)
                if MinY<=size(Obj.Image,1)
                    MaxY=size(Obj.Image,1);
                else
                    OutofBounds=1;
                end
            end
            
            if OutofBounds   % return empty if out of bounds found:             
                Bounds=[];
                CPoints=[];
                NCPoints=[];
                CPoints3D=[];
            else
                Bounds=[MinX,MaxX,MinY,MaxY];
                % Calculate Control Points:
                XMat=1:UnitSize:(MaxX-MinX+1);
                YMat=1:UnitSize:(MaxY-MinY+1);
                % Add the last element if necessary:
                if ~(XMat(end)==(MaxX-MinX+1))
                    XMat(end+1)=(MaxX-MinX+1);
                end
                if ~(YMat(end)==(MaxY-MinY+1))
                    YMat(end+1)=(MaxY-MinY+1);
                end
                CPoints=zeros(numel(XMat)*numel(YMat),2);
                CPoints3D=zeros(numel(XMat)*numel(YMat),3);
                Cnt=0;
                for f=XMat
                    for g=YMat
                        Cnt=Cnt+1;
                        CPoints(Cnt,:)=[f,g];
                        CPoints3D(Cnt,:)=Obj.Origin + Obj.Vector1 *(f+MinX-Org2D(1))*Obj.PixelScale + Obj.Vector2 *(g+MinY-Org2D(2))*Obj.PixelScale;
                    end
                end
                NCPoints=CPoints;
                MinX=min(NCPoints(:,1));
                MaxX=max(NCPoints(:,1));
                MinY=min(NCPoints(:,2));
                MaxY=max(NCPoints(:,2));
                NCPoints(:,1)=(NCPoints(:,1)-MinX)/MaxX;
                NCPoints(:,2)=(NCPoints(:,2)-MinY)/MaxY;
            end
        end
        
        
        function [Bounds,ImOut]=intersectIP(Obj,ImgPlane,UnitSize,PlotIt)
            
            %Calculate the 4 corners of the imageplane
            ImVz=ImgPlane.ImageNormal;
            ImVy=ImgPlane.ImageUp;
            ImVx=cross(ImVy,ImVz);
            Sx=size(ImgPlane.Image,1)/2*ImgPlane.PixelScale;
            Sy=size(ImgPlane.Image,2)/2*ImgPlane.PixelScale;
            C1=ImgPlane.CentrePoint - Sx*ImVx - Sy*ImVy;
            C2=ImgPlane.CentrePoint + Sx*ImVx - Sy*ImVy;
            C3=ImgPlane.CentrePoint + Sx*ImVx + Sy*ImVy;
            C4=ImgPlane.CentrePoint - Sx*ImVx + Sy*ImVy;
            Crns3D=[C1,C2,C3,C4];
            
            % Calculate the intersection of the corner rays with the Stitched Plane
            IntPnts3D=Obj.intersectRays(Crns3D,ImgPlane.SourcePosition,Obj.Origin,Obj.Vector3,PlotIt,'c','mx');
            if (PlotIt)
                ImgPlane.imshow3d(1,0);
            end
            % Return the pixel bounds and control points within the bounds:
            [Bounds,CPoints,NCPoints,CPoints3D]=Obj.returnPXLBounds(IntPnts3D,UnitSize);
            
            if isempty(Bounds) %if Image out of Bounds skips the rest of the calculations
                ImOut=[];
            else
                % Project the Control points over to the depth surface
                [SurfPnts3D]=projectOnSurf(Obj,CPoints3D,0,'m','ro');

                %             % Plot the 3D points:
                %             if (PlotIt)
                %                 plot3(SurfPnts3D(1,:),SurfPnts3D(2,:),SurfPnts3D(3,:),'ro');
                %                 for f=1:size(SurfPnts3D,2)
                %                     text(SurfPnts3D(1,f)+3,SurfPnts3D(2,f)+3,SurfPnts3D(3,f)+3,num2str(f),'Color',[0, 1 ,0]);
                %                 end
                %             end

                % Plot the 3D points:
                if (PlotIt)
                    plot3(CPoints3D(:,1),CPoints3D(:,2),CPoints3D(:,3),'r+');
                    for f=1:size(CPoints3D,1)
                        text(CPoints3D(f,1)+3,CPoints3D(f,2)+3,CPoints3D(f,3)+3,num2str(f),'Color',[0, 1 ,0]);
                    end
                end

                % Ray project the points over the imageplane
                AdjCPoints3D=Obj.intersectRays(SurfPnts3D,ImgPlane.SourcePosition,ImgPlane.CentrePoint,ImgPlane.ImageNormal,0,'y','g+')';
                if (PlotIt)            
                    for f=1:size(AdjCPoints3D,1)
                        text(AdjCPoints3D(f,1)+5,AdjCPoints3D(f,2)+5,AdjCPoints3D(f,3)+5,num2str(f),'Color',[1, 0 ,0]);
                    end
                end
                %Convert the adjusted control points to the local coordinates
                AdjCPoints2D=AdjCPoints3D(:,1:2)*0;
                for f=1:size(AdjCPoints3D,1)
                    dxyz=AdjCPoints3D(f,:)'-ImgPlane.CentrePoint;
                    AdjCPoints2D(f,:)=[dot(dxyz,ImVx),dot(dxyz,ImVy)]/ImgPlane.PixelScale+[size(ImgPlane.Image,1)/2;size(ImgPlane.Image,2)/2]';
                end

                ImOutputSize=[Bounds(2)-Bounds(1)+1,Bounds(4)-Bounds(3)+1];
                %             tform=cp2tform(AdjCPoints2D, CPoints,'lwm');
                %             ImOut=imtransform(ImgPlane.Image,tform,'bicubic','XData',[1,ImOutputSize(1)],'YData',[1,ImOutputSize(2)]);
                % Faster Alternative
                Img=ImgPlane.Image;
                % Add a thin black edge around the image:
                Img(1,:)=0;
                Img(end,:)=0;            
                Img(:,1)=0;
                Img(:,end)=0;            
                ImOut = warp_triangle(Img,AdjCPoints2D,CPoints,[ImOutputSize(2),ImOutputSize(1)]);
                %             Obj.Image(Bounds(3):Bounds(4),Bounds(1):Bounds(2))=Obj.Image(Bounds(3):Bounds(4),Bounds(1):Bounds(2))+ImOut;
                if (PlotIt)
                    figure;imshow(ImgPlane.Image);hold on;
                    plot(AdjCPoints2D(:,1),AdjCPoints2D(:,2),'r+');
                    for f=1:size(AdjCPoints2D,1)
                        text(AdjCPoints2D(f,1)+3,AdjCPoints2D(f,2)+3,num2str(f),'Color',[0, 1 ,0]);
                    end
                    axis tight;
                    figure;imshow(ImOut);hold on;
                    plot(CPoints(:,1),CPoints(:,2),'r+');
                    for f=1:size(CPoints,1)
                        text(CPoints(f,1)+3,CPoints(f,2)+3,num2str(f),'Color',[0, 1 ,0]);
                    end
                    axis tight;
                end
            end
        end
        
        % Ray plane intersections
        function IntPnts3D=intersectRays(Obj,Pnts3D,Source3D,PlaneOrg,PlaneNorm,PlotIt,LnCol,PntCol)
            if isempty(LnCol)
                LnCol='c';
            end
            if isempty(PntCol)
                PntCol='mx';
            end
            if (size(Pnts3D,2)==3) && ~(size(Pnts3D,1)==3)
                Pnts3D=Pnts3D';
            end
            IntPnts3D=Pnts3D*0;
            for f=1:size(Pnts3D,2)
                IntPnt=linePlaneIntersect(PlaneNorm,PlaneOrg,Pnts3D(:,f)-Source3D,Source3D);
                IntPnts3D(:,f)=IntPnt;
            end
            
            if PlotIt
                for f=1:size(Pnts3D,2)
                    % Ln=[IntPnts3D(:,f),Source3D];
                    Ln=[Pnts3D(:,f),Source3D];
                    plot3(Ln(1,:),Ln(2,:),Ln(3,:),LnCol);
                end
                plot3(IntPnts3D(1,:),IntPnts3D(2,:),IntPnts3D(3,:),PntCol);
                plot3(Pnts3D(1,:),Pnts3D(2,:),Pnts3D(3,:),'r*');
            end
        end
        
        
        % Ray Surface Plane Intersections
        function [PlanePnts3D, SurfPnts3D]=intersectRaysSurf(Obj,Pnts3D,Source3D,PlotIt,LnCol,PntCol)
            if isempty(LnCol)
                LnCol='c';
            end
            if isempty(PntCol)
                PntCol='mx';
            end
            if (size(Pnts3D,2)==3) && ~(size(Pnts3D,1)==3)
                Pnts3D=Pnts3D';
            end
            SurfPnts3D=Pnts3D*0;
            PlanePnts3D=Pnts3D*0;
            for f=1:size(Pnts3D,2)
                [PPnt,SurfPnt]=lineAnatomyIntersect(Obj.Vector3,Obj.Origin,Pnts3D(:,f)-Source3D,Source3D,Obj.DepthSegments,Obj.DepthSegmentsNormal,0);
                PlanePnts3D(:,f)=PPnt';
                SurfPnts3D(:,f)=SurfPnt';
            end
            
            if PlotIt
                for f=1:size(Pnts3D,2)
                    % Ln=[IntPnts3D(:,f),Source3D];
                    Ln=[Pnts3D(:,f),Source3D];
                    plot3(Ln(1,:),Ln(2,:),Ln(3,:),LnCol);
                end
                plot3(SurfPnts3D(1,:),SurfPnts3D(2,:),SurfPnts3D(3,:),PntCol);
                plot3(PlanePnts3D(1,:),PlanePnts3D(2,:),PlanePnts3D(3,:),'y+');
            end
        end
        
        % Project to Surface
        function [SurfPnts3D]=projectOnSurf(Obj,Pnts3D,PlotIt,LnCol,PntCol)
            if isempty(LnCol)
                LnCol='c';
            end
            if isempty(PntCol)
                PntCol='mx';
            end
            if (size(Pnts3D,2)==3) && ~(size(Pnts3D,1)==3)
                Pnts3D=Pnts3D';
            end
            SurfPnts3D=Pnts3D*0;
            for f=1:size(Pnts3D,2)
                [~,SurfPnt]=lineAnatomyIntersect(Obj.Vector3,Obj.Origin,Obj.Vector3,Pnts3D(:,f),Obj.DepthSegments,Obj.DepthSegmentsNormal,0);
                SurfPnts3D(:,f)=SurfPnt';
            end
            
            if PlotIt
                for f=1:size(Pnts3D,2)
                    % Ln=[IntPnts3D(:,f),Source3D];
                    Ln=[Pnts3D(:,f),SurfPnts3D(:,f)];
                    plot3(Ln(1,:),Ln(2,:),Ln(3,:),LnCol);
                end
                plot3(SurfPnts3D(1,:),SurfPnts3D(2,:),SurfPnts3D(3,:),PntCol);
            end
        end
        
        % Calculate the local 2D coordinates
        function [Pnts2D]=returnLocal2D(Obj,Pnts3D)
           Pnts2D=zeros(size(Pnts3D,1),2);
            for f=1:size(Pnts3D,1)
              Vec=Pnts3D(f,:)-Obj.Origin';
              Pnts2D(f,1:2)=[dot(Vec,Obj.Vector1'),dot(Vec,Obj.Vector2')];
           end
        end
            
    end
    
    
    
end
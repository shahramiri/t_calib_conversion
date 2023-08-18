classdef StitchMaster <handle
    % StitchMaster holds main information about the planes, grid, and
    % references to the World
    % 2018 Dec 12   Last edits
    
    properties
        Name
        ImagePlanes
        Grid
        APImagePlane
        MLImagePlane
        APIPinfo
        MLIPinfo
        APIndex
        MLIndex
        APLabel
        MLLabel
        APLongGrid
        MLLongGrid
        APLongView
        MLLongView
        APPlane
        MLPlane
        BiPlanarCropped
        ORSetupString
%         
        CurrentPlane
        CameraSurfHandle
        CameraCircleHandle
        CameraZoomFactor        
        APViewImageHandle
        APViewCircleHandle
        MLViewImageHandle
        MLViewCircleHandle
        ShowCropped
        DownSampleSize       % Biplanar Down Sample Size
        CameraDownSampleSize % Camra Down Sample Size
        Interface
        ReconVoxels
        IsVoxelBasedCrop
        World
        ImageFuseMode

        TrackBody
        AutDetectTrackBody % Flag for Auto detect track body

        APCarmViewStr % Carm View String corresponding to the APPlane
        MLCarmViewStr % Carm View String corresponding to the MLPlane

        IsStandardViews % determine whether the patient is lateral or not

        % Supporting the AP and Lat Angles:
        APObliqueAngle
        MLObliqueAngle

    end
    
    methods
        
        function Obj = StitchMaster(EName,varargin)                        
            Obj.Name = EName;
            Obj.DownSampleSize=2;
%             Obj.ImageFuseMode=1; % do simple image composition.
%             Obj.ImageFuseMode=3; % Simple Averaging
%             Obj.ImageFuseMode=5; % multi band blending            
            Obj.ImageFuseMode=6; % multi band blending + with gaussian gradient            
            Obj.APLabel='XY';
            Obj.MLLabel='YZ';
            Obj.Grid=GridElement('Empty Element');
            Obj.APPlane=StitchedPlane('AP Plane');
            Obj.APPlane.PlaneType='XY';
            Obj.APPlane.DownSampleSize=Obj.DownSampleSize;            
            Obj.APPlane.Grid=Obj.Grid;

            Obj.MLPlane=StitchedPlane('ML Plane');
            Obj.MLPlane.PlaneType='YZ';
            Obj.MLPlane.DownSampleSize=Obj.DownSampleSize;
            Obj.MLPlane.Grid=Obj.Grid;

            Obj.CurrentPlane='AP';
            Obj.CameraZoomFactor=1;
            Obj.ShowCropped=0;
            Obj.CameraDownSampleSize=2; % Camera Down Sample Size for making GUI interface response faster
            Obj.IsVoxelBasedCrop=0; % Start with non-voxel based cropping first

            % Default Carm View Strings 'AP','ML','OB':
            Obj.APCarmViewStr='AP';
            Obj.MLCarmViewStr='ML';
            Obj.IsStandardViews=0;

            Obj.ORSetupString='CRHDPSTB';
            % String Format: C?H?P?T?
            %   (C)-arm:      (R)ight/(L)eft
            %   (H)ead:       (D)own/(U)p
            %   (P)atient:    (S)upine/(P)rone
            %   (T)orus:      (B)ottom/(M)iddle/(T)op
            % Default Value:  'CRHDPSTB'
            
            % Supporting the oblique angles
            % Default is the same as AP and Lateral views if they are
            % specified
            Obj.APObliqueAngle=0;
            Obj.MLObliqueAngle=-90;
            
%             % fiducial based tracked body to establish a common image-based
%             % reference:
%             Obj.TrackBody=ImageTrackBody('Track Body');
%             Obj.TrackBody.defEngDesign('');
%             
%             Obj.AutDetectTrackBody=1;
            
            if ~isempty(varargin)
             
            end      
            
        end
        
        function showDebug(Obj, hFig)
                figure(hFig);
                hold on;
                Obj.APPlane.imshow3d(1,1,1,1);
                Obj.MLPlane.imshow3d(1,1,1,1);            
        end
        
        function turnCorrectInsideParallax(Obj,Val)    
            Obj.APPlane.CorrectInsideParallax=Val;
            Obj.MLPlane.CorrectInsideParallax=Val;
        end
        
        function switchPlane(Obj)
           switch Obj.CurrentPlane
               case 'AP'
                   Obj.CurrentPlane='ML';
               case 'ML'
                   Obj.CurrentPlane='AP';
           end
                      
        end
        
        function Out=buildGrid(Obj,Iap,Iml)
            
            Out=0;
            Vml=Iml.ImageNormal;
            Vap=Iap.ImageNormal;

            Vsi=cross(Vap,Vml);
            Vap=cross(Vml,Vsi);
            Vap=Vap/norm(Vap);
            Vml=Vml/norm(Vml);
            Vsi=Vsi/norm(Vsi);                        
            
            Cap=Iap.CentrePoint;
            Cml=Iml.CentrePoint;
            
            PAA=[Cap';Cap'+Vap'*500];
            PBB=[Cml';Cml'+Vml'*500];
            
            PA=[Cap';Cml'];
            PB=[Cap'+Vap'*500;Cml'+Vml'*500];
            
            [P_intersect,distances] = lineIntersect3D(PA,PB);
            % now place the point on the Vap:
            V=P_intersect'-Cap;
            dst=dot(V,Vap);
            Org=dst*Vap+Cap; % Origin of the grid.
            
            Obj.rebuildGrid(Org,Vml,Vsi);
            
            Out=1;
        end

        function Out=buildGrid_FullIPInfo(Obj,Iap,Iml)
            
            Out=0;
            Cap=retrurnImagePlanePP(Iap);
            Cml=retrurnImagePlanePP(Iml);
            
            CtrAP=Iap.CentrePoint;
            CtrML=Iml.CentrePoint;
            
            Vap=Cap-Iap.SourcePosition;
            Vap=Vap/norm(Vap);
            
            Vml=Cml-Iml.SourcePosition;
            Vml=Vml/norm(Vml);            
            
      
            PAA=[Cap';Cap'+Vap'*500];
            PBB=[Cml';Cml'+Vml'*500];
            
            % Use the centre of the II as the reference point
            PA=[CtrAP';CtrML'];
            PB=[CtrAP'+Vap'*500;CtrML'+Vml'*500];
            
            [P_intersect,distances] = lineIntersect3D(PA,PB);
            % now place the point on the Vap:
            V=P_intersect'-Cap;
            dst=dot(V,Vap);
            Org=dst*Vap+Cap; % Origin of the grid.
                        
%%          % Keeping the original Scout Image Normals
            Obj.Grid.ScoutAPNormal=-Iap.ImageNormal;
            Obj.Grid.ScoutMLNormal=-Iml.ImageNormal;

%% Modified code:            
            Vap=-Iap.ImageNormal;
            Vml=-Iml.ImageNormal;
            
            Vsi=cross(Vap,Vml);
            Vml=cross(Vsi,Vap);
            
            Vap=Vap/norm(Vap);
            Vml=Vml/norm(Vml);
            Vsi=Vsi/norm(Vsi);
            
            %%
            % Adjust the Grid Signs according to the ORSetup (May 17 2018)
            [APSign,MLSign,SISign]=Obj.analyzeORSetup;
            Vap=Vap*APSign;
            Vml=Vml*MLSign;
            Vsi=Vsi*SISign;
            %%
            Obj.rebuildGrid(Org,Vml,Vsi);
            
            Out=1;
        end
        
        function [APSign,MLSign,SISign]=analyzeORSetup(Obj)
            switch upper(Obj.ORSetupString(1:6))
                case {'CRHDPS','CLHUPS','CRHUPL','CLHDPL'}
                    APSign=1;
                    MLSign=1;
                    SISign=1;                    
                case {'CRHDPP','CLHUPP','CRHUPR','CLHDPR'}
                    APSign=-1;
                    MLSign=-1;
                    SISign=+1;                    
                case {'CRHUPS','CLHDPS','CRHDPR','CLHUPR'}
                    APSign=1;
                    MLSign=-1;
                    SISign=-1;       
                case {'CRHUPP','CLHDPP','CRHDPL','CLHUPL'}
                    APSign=-1;
                    MLSign=+1;
                    SISign=-1;                    
            end
            
            % Setting the Carm View Strings (to support lateral patient position)            
            switch upper(Obj.ORSetupString(6))
                case {'S','P'}
                    Obj.APCarmViewStr='AP';
                    Obj.MLCarmViewStr='ML';
                    Obj.IsStandardViews=1;
                case {'L','R'}
                    Obj.APCarmViewStr='ML';
                    Obj.MLCarmViewStr='AP';
                    Obj.IsStandardViews=0;
            end
            
            % Address the oblique position of the planes:
            if numel(Obj.ORSetupString)>9
                if (upper(Obj.ORSetupString(9))=='B')
                    switch (upper(Obj.ORSetupString(10)))
                        case '1' % AP & Lateral [Default]
                            % Leave as default
                        case '2' % AP & Oblique-Lateral 
                            if strcmpi(Obj.APCarmViewStr,'AP')
                                Obj.MLCarmViewStr='MLOB';
                            else
                                Obj.MLCarmViewStr='APOB';
                            end                                                                                                                 
                        case '3' % Oblique-AP & Lateral
                            if strcmpi(Obj.MLCarmViewStr,'ML')
                                Obj.APCarmViewStr='APOB';
                            else
                                Obj.APCarmViewStr='MLOB';
                            end                              
                    end
                end
            end
            
        end
        
        function Window=rebuildGrid(Obj,Org,Vx,Vy)              
            % create the gridelement:
            Sp1=SphereElement('PO',Org,10);
            Sp2=SphereElement('Px',Org+Vx*100,10);
            Sp3=SphereElement('Py',Org+Vy*100,10);
            Obj.Grid.TypeName='simple';            
            Obj.Grid.recreateElement('simple',Sp1,Sp2,Sp3);
            Obj.Grid.Origin=Obj.Grid.PlaneOrigin;
        end
 
        function initiateGrid(Obj, BuildFromScratch, Org, Vx, Vy, Iap, Iml)            

            PixScale=0.21; %<<<<<<<<<
            
            BuildBasedonCompleteImagePlaneInfo=0;
            if BuildFromScratch
                
                Org=[0,0,0]';
                Vx=[1,0,0]';
                Vy=[0,1,0]';
                Vz=[0,0,1]';
                Iap=ImagePlane('AP ImagePlane');
                Iml=ImagePlane('ML ImagePlane');
                Iap.CentrePoint=Org;
                Iap.ImageUp=Vy;
                Iap.ImageNormal=Vz;
                Iap.Image=uint8(zeros(1024,1024));
                Iap.PixelScale=0.21;
                Iap.SourcePosition=Org-1000*Iap.ImageNormal;
                Iml.CentrePoint=Org;
                Iml.ImageUp=Vz;
                Iml.ImageNormal=Vx;
                Iml.Image=uint8(zeros(1024,1024));
                Iml.PixelScale=0.21;
                Iml.SourcePosition=Org-1000*Iml.ImageNormal;
                Obj.APImagePlane=Iap;
                Obj.MLImagePlane=Iml;
                
            else
                
                if ~isempty(Org)  % if vector info is available
                    % Check to see if the normal vector to the AP plane is going to be flipped
                    % and if so prevent that:
                    Vz=cross(Vx,Vy);
                    Nz=Obj.APPlane.Vector3;
                    Zang=real(acosd(dot(Vz,Nz)));                    
                    if Zang>=90
                          Vz=Vz;
                    end
                    Vy=cross(Vz,Vx);
                    Nx=Obj.MLPlane.Vector3;
                    Xang=real(acosd(dot(Vx,Nz)));                    
                    if Xang>=90

                    end        
                    Vx=Vx/norm(Vx);
                    Vy=Vy/norm(Vy);
                    Vz=Vz/norm(Vz);                    
                    Obj.APImagePlane.CentrePoint=Org;
                    Obj.APImagePlane.ImageUp=Vy;
                    Obj.APImagePlane.ImageNormal=Vz;
                    Obj.APImagePlane.SourcePosition=Org-1000*Vz;
                    Obj.MLImagePlane.CentrePoint=Org;
                    Obj.MLImagePlane.ImageUp=Vz;
                    Obj.MLImagePlane.ImageNormal=Vx;
                    Obj.MLImagePlane.SourcePosition=Org-1000*Vx;
                         
                else % if imageplane info is available:
                    
                    BuildBasedonCompleteImagePlaneInfo=1;
                    
                    Obj.APImagePlane.PrincipalDist=Iap.PrincipalDist;
                    Obj.APImagePlane.Xoffset=Iap.Xoffset;
                    Obj.APImagePlane.Yoffset=Iap.Yoffset;
                    Obj.APImagePlane.PixelScale=Iap.PixelScale;
                    Obj.APImagePlane.CentrePoint=Iap.CentrePoint;
                    Obj.APImagePlane.SourcePosition=Iap.SourcePosition;
                    Obj.APImagePlane.ImageNormal=Iap.ImageNormal;
                    Obj.APImagePlane.ImageUp=Iap.ImageUp;
                    
                    Obj.MLImagePlane.PrincipalDist=Iml.PrincipalDist;
                    Obj.MLImagePlane.Xoffset=Iml.Xoffset;
                    Obj.MLImagePlane.Yoffset=Iml.Yoffset;
                    Obj.MLImagePlane.PixelScale=Iml.PixelScale;
                    Obj.MLImagePlane.CentrePoint=Iml.CentrePoint;
                    Obj.MLImagePlane.SourcePosition=Iml.SourcePosition;
                    Obj.MLImagePlane.ImageNormal=Iml.ImageNormal;
                    Obj.MLImagePlane.ImageUp=Iml.ImageUp;
                end
                
            end
            
            if BuildBasedonCompleteImagePlaneInfo
                Obj.buildGrid_FullIPInfo(Obj.APImagePlane,Obj.MLImagePlane); % for AP & ML Registration
            else
                Obj.buildGrid(Obj.APImagePlane,Obj.MLImagePlane);
            end
            
            if BuildFromScratch                        
                Obj.genStitchPlane(1,'AP long View','XY',[700,1500],PixScale,[0,0],0);
                Obj.genStitchPlane(2,'ML long View','YZ',[700,1500],PixScale,[0,0],0);
            end
%             Obj.cropViews(1);
            Obj.updateGrids;
        end                      
      
        function refreshIndex(Obj)            
            [APIndex,LatIndex]=FindAPandLats( Obj.Grid, Obj.ImagePlanes, Obj.APLabel, Obj.MLLabel);
            Obj.APIndex=APIndex;
            Obj.MLIndex=LatIndex;
        end
        
        function genStitchPlane(Obj,varargin)
            StitchIndex=varargin{1};
            PlaneName=varargin{2};
            PlaneType=varargin{3};
            %GridImageFile=varargin{4};
            GridImageSize=varargin{4};
            PixelPerMM=varargin{5};
            OriginOffset=varargin{6};
            BasedOnGrid=varargin{7}; % Go through the length of generating and scaling long grid views
                                     % in the beginning. Otherwise, just
                                     % revised the grid info in the 'reconstitchedplane below
            if StitchIndex==1  
                Obj.APPlane.Grid=Obj.Grid;
                Obj.APPlane.reconStitchedPlane(PlaneName,PlaneType,GridImageSize,PixelPerMM,OriginOffset,BasedOnGrid);
                Obj.APPlane.Vis2D=[-90,0];
                %Obj.APLongView=Obj.APLongGrid*0;
            else
                Obj.MLPlane.Grid=Obj.Grid;
                Obj.MLPlane.reconStitchedPlane(PlaneName,PlaneType,GridImageSize,PixelPerMM,OriginOffset,BasedOnGrid);
                %Obj.MLPlane.Vis2D=[180,1];
                Obj.MLPlane.Vis2D=[180,0];
                %Obj.MLLongView=Obj.MLLongGrid*0;            
            end
        end
  
        function [Res,UndoImage,UndoImageCoords,UndoImageSrcIndx,UndoImageDepthIndx]=undoAddImage(Obj,SelPlane)            
            
            Res=0; %Undo add incomplete
            SP=[];
            
            switch SelPlane
                case {'AP','1'}
                    SP=Obj.APPlane;
                case {'ML','2'}
                    SP=Obj.MLPlane;
            end
            
            if ~isempty(SP.UndoInfo)
                if ~isempty(SP.UndoInfo.Image)
                    UndoImage=SP.UndoInfo.Image; 
                    UndoImageSrcIndx=SP.UndoInfo.ImageSrcIndx;
                    UndoImageDepthIndx=SP.UndoInfo.ImageDepthIndx;
                    
                    %Swap the Image and UndoImage
                    XD=SP.UndoInfo.CoorX;
                    YD=SP.UndoInfo.CoorY;

                    SP.UndoInfo.Image=SP.Image(XD,YD);
                    SP.Image(XD,YD)=UndoImage;
                    UndoImageCoords=[min(XD),min(YD)];
                    
                    % Replace the ImageSrcIndx as the part of undo
                    SP.ImageSrcIndx(XD,YD)=UndoImageSrcIndx;
                    SP.ImageDepthIndx(XD,YD)=UndoImageDepthIndx;

                    Res=1; %Undo Add complete
                end
            else                
                UndoImage=[];
                UndoImageCoords=[];
                Res=0;
            end

        end
        
        function [AddedImage,AddedImageCoords,AddedImageSrcIndx,AddedImageDepthIndx,ErrCode]=addImage(Obj,varargin)

            AddedImage=[];
            AddedImageCoords=[];
            AddedImageSrcIndx=[];
            AddedImageDepthIndx=[];
            ErrCode=0;

            StitchIndex=varargin{1};
            ImPlane=varargin{2};
            RefreshOutput=varargin{3};
            if nargin>4
                Voxel=varargin{4}; % voxel information determining the depth of the anatomy                
            else
                Voxel=[];           % empty voxel information if nothing is provided
            end

            %%
            ShowGraphics=varargin{5};
            ImagePlaneIndex=varargin{7};
            
            % Detect TrackBody
            if (Obj.AutDetectTrackBody)                
                [Detected,RecT,RecErr,MarkedImage]=Obj.TrackBody.detectTrackBody(StitchIndex,ImPlane,0,0,1,[]);
%                 ImPlane.Image=MarkedImage;
                if (Detected)
                   fprintf('\nTracker Body Detected on the collected Image.\n');                    
                end
            end

            if StitchIndex==1
                SP=Obj.APPlane;
            end            
            if StitchIndex==2
                SP=Obj.MLPlane;
            end

            %%
%             % Debug Help Section to display two imageplanes and the
%             % intersection calculations:
%             % Temporarily added to visualaize Anatomy Depth Calculations
%             if (ImagePlaneIndex>=4) && (ImagePlaneIndex<=5)
%                 ShowGraphics=1;
%             else
%                 ShowGraphics=0;
%             end
            %%
            
            % Main function to calculate the added ImagePlane to the StitchedPlane:
            [~,W,ErrCode,XD,YD]=SP.renderImage(ImPlane,[],1,ShowGraphics,Voxel,ImagePlaneIndex);
                        
            if ~isempty(W)
                Im=W.Image;
                
                ROI=SP.Image(XD,YD);
                % Save the previous Image portion in the Undo Portion so it can be reversed:
                SP.UndoInfo=[];
                SP.UndoInfo.Image=ROI;
                SP.UndoInfo.CoorX=XD;
                SP.UndoInfo.CoorY=YD;
                % -------------------------------------------------------
                % New setup for combining the ROI with With Im (simple overlay):
              
                switch Obj.ImageFuseMode
                    case 1 % simple overwrite
                        Msk1=ROI*0;
                        Msk2=Msk1+1;
                        Msk1(Im>0)=1;
                        Msk2(Im>0)=0;
                        ROI3=(Msk1.*Im)+(Msk2.*ROI);
                    case 2 % simple averaging
                        Msk1=ROI*0;
                        Msk2=Msk1+1;
                        Msk1(Im>0)=1;
                        Msk2(Im>0)=0;
                        ROI3=((Msk1.*Im)+(Msk1.*ROI))/2+(Msk2.*ROI);
                        
                        % Old setup for combining the ROI with With Im:
                        Msk1=ROI*0;
                        Msk2=Msk1;
                        OVL=ROI*0+1;
                        Msk1(ROI>0)=1;
                        Msk2(Im>0)=1;
                        Wght=Msk1+Msk2;
                        OVL(Wght==2)=1;
                        OVERLAP=min(ROI,Im).*OVL;
                        
                        OP=(OVERLAP>0);
                        OOP=~OP;
                        OOOP=uint8(OOP);
                        ROI3=ROI3+OOOP.*Im;
                    case 3
                        Msk1=ROI*0;
                        Msk1(Im>0)=1;
                        % ROI3=((Msk1.*Im)+(Msk1.*ROI))/2;
                        ROI3=(Im+ROI)/2;
                    case 4
                        ROI3=imfuse(ROI,Im,'blend'); 
                    case 5 %MultiBand Blending          
                        Msk1=ROI*0;
                        Msk2=Msk1;
                        OVL=ROI*0+1;
                        Msk1(ROI>0)=1;
                        Msk2(Im>0)=1;
                        Wght=Msk1+Msk2;
                        OVL(Wght==2)=1;
                        OVERLAP=min(ROI,Im).*OVL;
                        OP=(OVERLAP>0);
                        if max(max(OVERLAP))==0
                            % if there is not an overlap:
                            Msk1=ROI*0;
                            Msk2=Msk1+1;
                            Msk1(Im>0)=1;
                            Msk2(Im>0)=0;
                            ROI3=(Msk1.*Im)+(Msk2.*ROI);                             
                        else
                            % if there is an overlap
                            S=regionprops(OP,'centroid'); % centroid of the overlap area
                            pivot=S.Centroid;
                            CP=[floor(size(OP,2)/2),floor(size(OP,1)/2)];
                            
                            dP=CP-pivot;
                            Ang=mod(atan2d(dP(2),dP(1))+270,360);
                            FilterPar1=100;
                            FilterPar2=115;
                            OPD=uint8(OP);                            
                            % SourceIndex
                            SRC=SP.ImageSrcIndx(XD,YD);
                            imgo=uint8(pyrOblqBlend2(Im,ROI,pivot,Ang,FilterPar1,FilterPar2));
                            ROI3=(Msk1-OPD).*ROI + (Msk2-OPD).*Im + OPD.*imgo;
                        end
                    case 6 %Guasian filter edge blurring + MultiBand blending
                        Msk1=ROI*0;
                        Msk2=Msk1;
                        OVL=ROI*0+1;
                        Msk1(ROI>0)=1;
                        Msk2(Im>0)=1;
                        Wght=Msk1+Msk2;
                        OVL(Wght==2)=1;
                        OVERLAP=min(ROI,Im).*OVL;
                        OP=(OVERLAP>0);
                        if max(max(OVERLAP))==0
                            % if there is not an overlap:
                            Msk1=ROI*0;
                            Msk2=Msk1+1;
                            Msk1(Im>0)=1;
                            Msk2(Im>0)=0;
                            ROI3=(Msk1.*Im)+(Msk2.*ROI);                             
                        else
                            OPD=uint8(OP);                            
                            imgo=uint8(gausBlend(Im,ROI,0));
                            ROI3=(Msk1-OPD).*ROI + (Msk2-OPD).*Im + OPD.*imgo;
                        end                        
                end
                % -------------------------------------------------------
                SP.Image(XD,YD)=ROI3;
                
                % Plug in the ImagePlaneIndex
                if ~isempty(W.ImageSrcIndx)
                    %update the image source index
                    SP.ImageSrcIndx(XD,YD)=SP.ImageSrcIndx(XD,YD).*(uint8(W.ImageSrcIndx==0)) + W.ImageSrcIndx.*(uint8(W.ImageSrcIndx>0));                    
                end
                
                %Output the corresponsing parts of the Image
                AddedImage=ROI3;
                AddedImageCoords=[min(XD),min(YD)];
                AddedImageSrcIndx=SP.ImageSrcIndx(XD,YD);
                AddedImageDepthIndx=SP.ImageDepthIndx(XD,YD);
                
                %Store the source index information for the undo.
                SP.UndoInfo.ImageSrcIndx=SP.ImageSrcIndx(XD,YD);
                SP.UndoInfo.ImageDepthIndx=SP.ImageDepthIndx(XD,YD);
                
                if RefreshOutput
                    SP.refreshImageOutput; % Refresh ImageOutput
                end
                
            end
            
        end        
        
        function res= cropViews(Obj,ProcessImage)
            
            % New Version 2018-Apr-09 for fast calculation of cropping
            % viewing wo need for imrotate or mirrors
            % This assumes that the AP image needs to be 90 degrees rotated
            % CCW inorder to make the proper cropping of the views.
            % The result is reflected in the .CroppedCoords of the views
            
            res1=0;
            res2=0;           
            
            scl=10; 
            AP=Obj.APPlane.Image(1:scl:end,1:scl:end,1);
            AP(AP>0)=1;            
            Sz=size(Obj.APPlane.Image(:,:,1));
            W=Sz(2); H=Sz(1);
            r1 = regionprops(AP, 'BoundingBox' );            
            if isempty(r1)
                rec1=[];
            else
                rec1p=r1.BoundingBox+[-2,-2,+4,+4]; % making sure there is enough margin around the cropped area
                rec1p=rec1p*scl;
                rec1=[W-rec1p(1)-rec1p(3)+1,rec1p(2),rec1p(3),rec1p(4)];
            end
                        
            if ~isempty(rec1)
                box1=rec1;%.BoundingBox;
                A1=[box1(1),box1(2)];
                B1=[box1(1)+box1(3)-1,box1(2)+box1(4)-1];
                res1=1;
            end

            ML=Obj.MLPlane.Image(1:scl:end,1:scl:end,1); 
            ML(ML>0)=1;
            r2 = regionprops(ML, 'BoundingBox' );
            if isempty(r2)
                rec2=[];
            else
                rec2=r2.BoundingBox+[-2,-2,+4,+4];
                rec2=rec2*scl;
            end
                        
            if ~isempty(rec2)
                box2=rec2; %.BoundingBox;
                A2=[box2(2),box2(1)];
                B2=[box2(2)+box2(4)-1,box2(1)+box2(3)-1];
                res2=1;
            end

            if (res1 || res2) % if there is bounding box for either views:

                if res2 && res1
                    %coordinate of the cropped image for both views:
                    %which means AP rotated by 90 degreed CCW and ML not rotated:
                    C1=([min([A1(1),A2(1)]),min([A1(2),A2(2)])]);
                    C2=([max([B1(1),B2(1)]),max([B1(2),B2(2)])]);
                end            
            
                % if only image exists on the AP view:
                if res1 && (~res2);
                    C1=[A1(1),A1(2)];
                    C2=[B1(1),B1(2)];                
                end

                % if only image exists on the MLview
                if res2 && (~res1);
                    C1=[A2(1),A2(2)];
                    C2=[B2(1),B2(2)];                 
                end
                
                % Apply coordinates of the cropping coords
                D1=[C1(2),W-C2(1)];
                D2=[C2(2),W-C1(1)];
                Obj.APPlane.CroppedCoords=[D1(1),D2(1),D1(2),D2(2)];
                Obj.MLPlane.CroppedCoords=[C1(1),C2(1),C1(2),C2(2)];
            
            end
            
            res=res1 || res2;
            
            [AP_S1,AP_S2]=size(Obj.APPlane.Image);
            [ML_S1,ML_S2]=size(Obj.MLPlane.Image);

            if ~(res) % if no cropping happened then return different size
                % don't return the full size images
                % Obj.APPlane.CroppedCoords=[1,AP_S1,1,AP_S2];
                % Obj.MLPlane.CroppedCoords=[1,ML_S1,1,ML_S2];                
                
                Obj.APPlane.OrigCroppedCoords=[1,AP_S1,1,AP_S2];
                Obj.MLPlane.OrigCroppedCoords=[1,ML_S1,1,ML_S2];                

                % return a much smaller tiny blank image (50 times smaller)
                RszFactor=50;
                Obj.APPlane.CroppedCoords=[1,round(AP_S1/RszFactor),1,round(AP_S2/RszFactor)];
                Obj.MLPlane.CroppedCoords=[1,round(ML_S1/RszFactor),1,round(ML_S2/RszFactor)];
            else 
                % make sure the APcropping coords stay within bounds
                if Obj.APPlane.CroppedCoords(1)<1
                    Obj.APPlane.CroppedCoords(1)=1;
                end
                if Obj.APPlane.CroppedCoords(3)<1
                    Obj.APPlane.CroppedCoords(3)=1;
                end
                if Obj.APPlane.CroppedCoords(2)>AP_S1
                    Obj.APPlane.CroppedCoords(2)=AP_S1-1;
                end
                if Obj.APPlane.CroppedCoords(4)>AP_S2
                    Obj.APPlane.CroppedCoords(4)=AP_S2-1;
                end
                % make sure the MLcropping coords stay within bounds
                if Obj.MLPlane.CroppedCoords(1)<1
                    Obj.MLPlane.CroppedCoords(1)=1;
                end
                if Obj.MLPlane.CroppedCoords(3)<1
                    Obj.MLPlane.CroppedCoords(3)=1;
                end
                if Obj.MLPlane.CroppedCoords(2)>ML_S1
                    Obj.MLPlane.CroppedCoords(2)=ML_S1-1;
                end
                if Obj.MLPlane.CroppedCoords(4)>ML_S2
                    Obj.MLPlane.CroppedCoords(4)=ML_S2-1;
                end
                
                % OrigCroppedCoords is the same as CroppedCoords if there are nozero pixel values on both views
                Obj.APPlane.OrigCroppedCoords=Obj.APPlane.CroppedCoords;
                Obj.MLPlane.OrigCroppedCoords=Obj.MLPlane.CroppedCoords;
                
            end     
            
            

        end
        
        function res=cropViewsBasedOnVoxels(Obj) %further crop the views after refine to zoom on the target landmarks
                                                 %use the .cropViews to return to the normal box            
            
            % If therewas a previously calculated offset adds it back
            % first:
            MinXX=Obj.APPlane.DrawingDataOffsetX;
            MinYY=Obj.APPlane.DrawingDataOffsetY;
            Obj.APPlane.DrawingData(:,1)=Obj.APPlane.DrawingData(:,1)+MinXX;
            Obj.APPlane.DrawingData(:,2)=Obj.APPlane.DrawingData(:,2)+MinYY;
            MinXX=Obj.MLPlane.DrawingDataOffsetX;
            MinYY=Obj.MLPlane.DrawingDataOffsetY;            
            Obj.MLPlane.DrawingData(:,1)=Obj.MLPlane.DrawingData(:,1)+MinXX;
            Obj.MLPlane.DrawingData(:,2)=Obj.MLPlane.DrawingData(:,2)+MinYY;
            
            DD_AP=Obj.APPlane.DrawingData;
            DD_ML=Obj.MLPlane.DrawingData;
            MinX=min(min(DD_AP(:,1)),min(DD_ML(:,1)));
            MaxX=max(max(DD_AP(:,1)),max(DD_ML(:,1)));            
            MinY=min(min(DD_AP(:,2)),min(DD_ML(:,2)));
            MaxY=max(max(DD_AP(:,2)),max(DD_ML(:,2)));
            
            CenX=(MinX+MaxX)/2;
            CenY=(MinY+MaxY)/2;
            LenX=MaxX-MinX;
            LenY=MaxY-MinY;
            
            % enlarge the area by %10:
            LX=LenX/2;%*1.05;
            LY=LenY/2;%*1.05;
            
            CalRatio=LX/LY; %Calculated Ration
            AxesRatio=98/47; %Axes Aspect Ratio Deriven from the design of the GUI
            
            %Decide on the cropped Size:
            if CalRatio>AxesRatio
                LY=LX/AxesRatio;
            else
                LX=LY*AxesRatio;
            end
            
            MinXX=round(CenX-LX)+1;
            MaxXX=round(CenX+LX)+1;
            MinYY=round(CenY-LY)+1;
            MaxYY=round(CenY+LY)+1;
            

%             P=[MinXX,MinYY;MaxXX,MinYY;MaxXX,MaxYY;MinXX,MaxYY;MinXX,MinYY];
            
            PatchCoords=[MinYY,MaxYY,MinXX,MaxXX];
            PatchImageAP=Obj.APPlane.returnVisImage(1,1);
            
            SZY=size(PatchImageAP,1);
            SZX=size(PatchImageAP,2);
            if MinXX<0;
                MinXX=1;
            end
            if MinYY<0
                MinYY=1;
            end
            
            if MaxXX>SZX;
                MaxXX=SZX;
            end
            if MaxYY>SZY
                MaxYY=SZY;
            end
            
            
            PatchImageAP=imresize(PatchImageAP(MinYY:MaxYY,MinXX:MaxXX,:),Obj.DownSampleSize);
            [PatchImage,OPatchCoords]=Obj.APPlane.returnPatch(PatchCoords*Obj.DownSampleSize,PatchImageAP);
            OP=OPatchCoords;
            SZ=size(PatchImageAP);
            C=[];
            C(1)=OP(1);
            C(2)=OP(1)+SZ(2)-1;
            C(3)=OP(2);
            C(4)=OP(2)+SZ(1)-1;
            Obj.APPlane.CroppedCoords=C;
            Obj.APPlane.CroppedImage=Obj.APPlane.ImageOutput(C(1):C(2),C(3):C(4),:);
            
            PatchCoords=[MinYY,MaxYY,MinXX,MaxXX];
            PatchImageML=Obj.MLPlane.returnVisImage(1,1);
            PatchImageML=imresize(PatchImageML(MinYY:MaxYY,MinXX:MaxXX,:),Obj.DownSampleSize);
            [PatchImage,OPatchCoords]=Obj.MLPlane.returnPatch(PatchCoords*Obj.DownSampleSize,PatchImageAP);
            OP=OPatchCoords;
            SZ=size(PatchImageML);
            C=[];
            C(1)=OP(1);
            C(2)=OP(1)+SZ(1)-1;
            C(3)=OP(2);
            C(4)=OP(2)+SZ(2)-1;
            Obj.MLPlane.CroppedCoords=C;
            Obj.MLPlane.CroppedImage=Obj.MLPlane.ImageOutput(C(1):C(2),C(3):C(4),:);
            res=1;
            %update the Drawing_Data
            Obj.APPlane.DrawingData(:,1)=Obj.APPlane.DrawingData(:,1)-MinXX;
            Obj.APPlane.DrawingData(:,2)=Obj.APPlane.DrawingData(:,2)-MinYY;            
            Obj.MLPlane.DrawingData(:,1)=Obj.MLPlane.DrawingData(:,1)-MinXX;
            Obj.MLPlane.DrawingData(:,2)=Obj.MLPlane.DrawingData(:,2)-MinYY;
            Obj.APPlane.DrawingDataOffsetX=MinXX;
            Obj.APPlane.DrawingDataOffsetY=MinYY;      
            Obj.MLPlane.DrawingDataOffsetX=MinXX;
            Obj.MLPlane.DrawingDataOffsetY=MinYY;              
        end
        
        function showAPPlaneCropped(Obj)
            figure;
            SM.APPlane.imshow3d(0,1,1);
            
        end
        
        function Spot=returnSpot(Obj,ImagePlane,ViewIndex)

            switch ViewIndex
                case 1
                    Plane=Obj.APPlane;
                case 2
                    Plane=Obj.MLPlane;
            end
            Pnts=[];
            Spot=[];
            
            switch ViewIndex
                case 1                    
                    Iap=ImagePlane;
                    [ghandles,W,ErrCode,~,~]=Plane.renderImage(Iap,Iap,0,0,[],[]);                    
                    P=W.Coords/Plane.PixelScale;
                    Cnt=mean(P); % Centre of the view wrt to the grid
                    P(:,1)=P(:,1)-P(1,1);
                    P(:,2)=P(:,2)-P(2,2);
                    P(:,3)=P(:,1).^2+P(:,2).^2;
                    Rad=(median(P(:,3))).^0.5/2; %Radius of the view in Pixels
                    
                    %Rotation of the view:
%                     V=Iap.ImageUp;
%                     Vpr=dot(V,Plane.Vector1)*Plane.Vector1+dot(V,Plane.Vector2)*Plane.Vector2;
%                     Roll=acosd(dot(Vpr,Plane.Vector2))-90;
%                     % Roll=acosd(dot(Iap.ImageUp,Plane.Vector2))-90;
                    
                     V=Iap.ImageUp;
                     Vpr=dot(V,Plane.Vector1)*Plane.Vector1+dot(V,Plane.Vector2)*Plane.Vector2;
                     Roll=acosd(dot(Vpr,Plane.Vector1))-90;                    
                    
                    
                    
                    Tar=[Cnt(1),-Cnt(2),0];
                    Pos=[Cnt(1),-Cnt(2),2000];
                    
                    V=Iap.ImageUp;
                    V1=Plane.Vector2;
                    V2=Plane.Vector1;
                    V3=(Pos-Tar)/norm(Pos-Tar);
                    V1=V1/norm(V1);
                    V2=V2/norm(V2);
                    V3=V3/norm(V3);
                    Vup=dot(V,V1)*V1+dot(V,V2)*V2;
%                     Up=Vup/norm(Vup);
                    Up=cross(V3,Vup);
                    Up=Up/norm(Up);
                    % Always point up:
                    Up=[-1,0,0];

                    
%                     %set the Orientation according to the Grid Refs:
%                     Up=Obj.Grid.VyLineElement.Vector;
%                     Up=Up/norm(Up);

                    
%                     Up=-cross(Iap.ImageUp,Iap.ImageNormal);
%                     Up=Up/norm(Up);
                    Dist=norm(Pos-Tar);
                    VA=2*atand(Rad/Dist);       
                    
                    
                    %Calculate Position2D:
                    Vn=(Pos-Tar)/Dist;
                    Vu=Up;
                    Vh=cross(Up,Vn);
                    Px=dot(Vh,Pos);
                    Py=dot(Vu,Pos);
                    Position2D=round([Px,Py]);%Obj.APPlane.PixelScale);
                    
                    % Position2D=[-Tar(1),-Tar(2)]; %*Plane.PixelScale; %*Plane.PixelScale;
                    Rad2D=Rad; %*Plane.PixelScale; %*Plane.PixelScale;
                    
                case 2
                     Iap=ImagePlane;
                     [ghandles,W,~,~]=Plane.renderImage(Iap,Iap,0,0,[],[]);
                     P=W.Coords/Plane.PixelScale;
                     Cnt=mean(P); % Centre of the view wrt to the grid
                     P(:,1)=P(:,1)-P(1,1);
                     P(:,2)=P(:,2)-P(2,2);
                     P(:,3)=P(:,1).^2+P(:,2).^2;
                     Rad=(median(P(:,3))).^0.5/2*1.1; %Radius of the view in Pixels
                                                               
                     
                     %Rotation of the view:
%                      Roll=acosd(dot(Iap.ImageNormal,Plane.ImageNormal));
                     V=Iap.ImageNormal;
                     Vpr=dot(V,Plane.Vector1)*Plane.Vector1+dot(V,Plane.Vector2)*Plane.Vector2;
                     Roll=acosd(dot(Vpr,Plane.Vector1))-90;
                     
%                     V=Iap.ImageUp;
%                     Vpr=dot(V,Plane.Vector1)*Plane.Vector1+dot(V,Plane.Vector2)*Plane.Vector2;
%                     Roll=acosd(dot(Vpr,Plane.Vector1))-90;                     
                     
                     
                     % Roll=acosd(dot(Iap.ImageNormal,Plane.Vector1))-90;
                     Tar=[Cnt(1),-Cnt(2),0];
                     Pos=[Cnt(1),-Cnt(2),-2000];
                     % Up=Iap.ImageUp;
                     V=Iap.ImageUp;
                     V1=Plane.Vector2;
                     V2=Plane.Vector1;
                     V3=(Pos-Tar)/norm(Pos-Tar);
                     V1=V1/norm(V1);
                     V2=V2/norm(V2);
                     V3=V3/norm(V3);
                     Vup=dot(V,V1)*V1+dot(V,V2)*V2;
                     Up=cross(V3,Vup);
                     Up=Up/norm(Up);
                     % Always point up:
                     Up=[0,1,0];

                     Dist=norm(Pos-Tar);
                     VA=2*atand(Rad/Dist);
                     
                    %Calculate Position2D:
                    Vn=(Pos-Tar)/Dist;
                    Vu=Up;
                    Vh=cross(Up,Vn);
                    Px=dot(Vh,Pos);
                    Py=dot(Vu,Pos);
                    Position2D=round([Px,Py]); %/Obj.MLPlane.PixelScale);
                     
                    %Position2D=[-Tar(1),-Tar(2)]; %*Plane.PixelScale; %*Plane.PixelScale;
                    Rad2D=Rad; %*Plane.PixelScale; %*Plane.PixelScale;
            end
            Spot.CameraTarget=Tar;
            Spot.CameraPosition=Pos;
            Spot.CameraUp=Up;

            % Correct the camera orientation for the ML view
            if ViewIndex==2
                Spot.CameraPosition(3)=-Spot.CameraPosition(3);
                Spot.CameraUp=-Spot.CameraUp;
            end            

            Spot.CameraViewAngle=VA*Obj.CameraZoomFactor;
            Spot.Radius=Rad;
            Spot.Position2D=Position2D;
            Spot.Rad2D=Rad2D;
            Spot.Roll=Roll;

        end
        
        function showImageOutput(Obj,axesHandle)
           axes(axesHandle);
           n=Obj.DownSampleSize;
           switch Obj.CurrentPlane
               case 'AP'
                   imshow(imrotate(Obj.APPlane.ImageOutput(1:n:end,1:n:end,:),90));
               case 'ML'
                   imshow(imrotate(Obj.MLPlane.ImageOutput(1:n:end,1:n:end,:),90));
           end
           java.lang.System.gc();
        end
        
        function [APViewImageHandle,MLViewImageHandle,APCircleHandle,MLCircleHandle]=showBiplanarViews(Obj,ImgPlane, APaxesHandle,MLaxesHandle,WithGraphics)
            
            APViewImageHandle=[];
            MLViewImageHandle=[];
            APCircleHandle=[];
            MLCircleHandle=[];
            
            n=Obj.DownSampleSize;
            
            ShowCropped=Obj.ShowCropped;
            
            if WithGraphics
                axes(APaxesHandle);
                hold(APaxesHandle,'off');
            end
            if isempty(Obj.APPlane.CroppedImage) || isempty(Obj.MLPlane.CroppedImage)
                ShowCropped=0; %Don't show cropped if either views are empty
            end
                
            if ~ShowCropped                
                if WithGraphics
                    APViewImageHandle=imshow(imrotate(Obj.APPlane.ImageOutput(1:n:end,1:n:end,:),90));
                end
            else
                Obj.cropViews(1);
                if WithGraphics
                    APViewImageHandle=imshow(imrotate(Obj.APPlane.CroppedImage(1:n:end,1:n:end,:),90));
                end
            end
            
            if WithGraphics
                axis equal; axis tight;
                set(APaxesHandle,'Units', 'normalized');
                APpos=get(APaxesHandle,'Position');
                hold(APaxesHandle,'on');

                axes(MLaxesHandle);
                hold(MLaxesHandle,'off');
                if ~ShowCropped
                    MLViewImageHandle=imshow(imrotate(Obj.MLPlane.ImageOutput(1:n:end,1:n:end,:),0));
                else
                    %Obj.cropViews;
                    MLViewImageHandle=imshow(imrotate(Obj.MLPlane.CroppedImage(1:n:end,1:n:end,:),0));
                end
                axis equal; axis tight;
                set(MLaxesHandle,'Units', 'normalized');
                MLpos=get(MLaxesHandle,'Position');
                MLpos(1)=APpos(1);
                MLpos(3)=APpos(3);
                set(MLaxesHandle,'Units', 'normalized', 'Position', MLpos);
                hold(MLaxesHandle,'on');
            end

%             APOffset=[size(Obj.APPlane.ImageOutput,1)/n/2,size(Obj.APPlane.ImageOutput,2)/n/2];
%             MLOffset=[size(Obj.MLPlane.ImageOutput,1)/n/2,size(Obj.MLPlane.ImageOutput,2)/n/2];
            
            APOffset=[size(Obj.APPlane.Image,1)/n/2,size(Obj.APPlane.Image,2)/n/2];
            MLOffset=[size(Obj.MLPlane.Image,1)/n/2,size(Obj.MLPlane.Image,2)/n/2];
            
            
            if ShowCropped
                C=Obj.APPlane.CroppedCoords;
                APOffset=APoffset+[(C(1)+C(2))/2,(C(3)+C(4))/2];
                C=Obj.MLPlane.CroppedCoords;
                MLOffset=MLoffset+[(C(1)+C(2))/2,(C(3)+C(4))/2];                
            end
            
            switch Obj.CurrentPlane
                case 'AP'
                    Spot=Obj.returnSpot(ImgPlane,1);
                case 'ML'
                    Spot=Obj.returnSpot(ImgPlane,2);
            end
            
            Tar=Spot.CameraTarget;
            Pos=Spot.CameraPosition;
            Up=Spot.CameraUp;
            VA=Spot.CameraViewAngle;
            Rad=Spot.Radius;
            Pos2D=Spot.Position2D/n;
            Rad2D=Spot.Rad2D/n;
            
            %Circular point:
            a=[0:5:360]';
            c=0.9;
            Pnts(:,1)=c*Rad2D*cosd(a);
            Pnts(:,2)=c*Rad2D*sind(a);
            Pnts(:,3)=Pnts(:,1)*0;
            % axes(AxesHandle);
            
            if WithGraphics
                switch Obj.CurrentPlane
                    case 'AP'
                        APViewCircHandle=plot(APaxesHandle, Pnts(:,1)+Pos2D(1)+APOffset(1),Pnts(:,2)+Pos2D(2)+APOffset(2),'y--','LineWidth',3);
                        MLViewCircHandle=plot(MLaxesHandle, Pnts(:,1)*0,Pnts(:,2)*0,'y--','LineWidth',3);
    %                     axes(APaxesHandle);
    %                     axis equal;
                    case 'ML'
                        APViewCircHandle=plot(APaxesHandle, Pnts(:,1)*0,Pnts(:,2)*0,'y--','LineWidth',3);
                        MLViewCircHandle=plot(MLaxesHandle, Pnts(:,1)-Pos2D(2)+MLOffset(2),Pnts(:,2)+Pos2D(1)+MLOffset(1),'y--','LineWidth',3);
    %                     axes(MLaxesHandle);
    %                     axis equal;
                end   
            end
            
%             Obj.APViewImageHandle=APViewImageHandle;
%             Obj.APViewCircleHandle=APViewCircHandle;
%             Obj.MLViewImageHandle=MLViewImageHandle;
%             Obj.MLViewCircleHandle=MLViewCircHandle;
            
%             java.lang.System.gc();
        end

        function updateBiplanarViews(Obj,ImgPlane, APaxesHandle,MLaxesHandle,UpdateData,UpdateCircle,UpdateTemplate,WithGraphics)
%             pause (0.03);
            n=Obj.DownSampleSize;
            temp_n=1;
            ShowCropped=Obj.ShowCropped;
            % CropResult=0;
            
            NeverCropped=0;
            if isempty(Obj.APPlane.CroppedCoords) || isempty(Obj.MLPlane.CroppedCoords)
                NeverCropped=1;
            end
            
            
            if UpdateData
                OkToCrop=1;
                if ShowCropped && NeverCropped
                    CropRes=Obj.cropViews(0);
                    if isempty(Obj.APPlane.CroppedImage) || isempty(Obj.MLPlane.CroppedImage)
                        ShowCropped=0; %Don't show cropped if either views are empty
                        Obj.ShowCropped=0;
                        %OkToCrop=0;
                    end
                end
                
                if ShowCropped
                    Obj.cropViews(1);
                    if Obj.IsVoxelBasedCrop
                        Obj.cropViewsBasedOnVoxels; %Zoom to the selected portion of the image only
                    end
                    
                    % To include preop in the cropping
                    if UpdateTemplate
                        DATA=Obj.APPlane.returnVisImage(1,1);
                        Im=Obj.APPlane.Image(1:end,1:end);
%                         DATAtemp=cat(3,Im,Im*0.95+10,Im*0.95+10);
                    end
                    
                else
                    
                    if UpdateTemplate
                        DATA=Obj.APPlane.returnVisImage(2,1);
                        Im=Obj.APPlane.Image(1:end,1:end);
%                         DATAtemp=cat(3,Im,Im*0.95+10,Im*0.95+10);
                    end
                end
                
                
                if WithGraphics
                    set(Obj.APViewImageHandle,'XData',[1,size(DATA,2)],'YData',[1,size(DATA,1)]);
                    set(Obj.APViewImageHandle,'CData',DATA);
                
                    axis(APaxesHandle,'equal'); %axis tight;
                    xlim(APaxesHandle,'auto');
                    ylim(APaxesHandle,'auto');
                end
                
%                 %=========================================================================================
%                 if UpdateTemplate
%                     % Update the corresponding AP tempplating view:
%                     Sz=size(DATAtemp);
%                     XD=round(Sz(2)/2);
%                     YD=round(Sz(1)/2);
%                     Xdata=[-XD, XD;-XD,XD];
%                     Ydata=[YD,YD;-YD,-YD];
%                     Zdata=[0,0;0,0];
%                     if WithGraphics
%                         set(Obj.APSurfHandle,'XData',Xdata,'YData',Ydata,'ZData',Zdata);
%                     end
%                     Dat=DATAtemp;
%                     %                         DATAtemp(:,:,2)=Dat;
%                     %                         DATAtemp(:,:,3)=Dat;
%                     if WithGraphics
%                         set(Obj.APSurfHandle,'CData',DATAtemp(1:temp_n:end,1:temp_n:end,:));
%                         xlim(Obj.APTemplatingAxes,[-XD,XD]);
%                         ylim(Obj.APTemplatingAxes,[-YD,YD]);
%                         axis(Obj.APTemplatingAxes,'vis3d'); %axis(Obj.APTemplatingAxes,'tight');
%                         axis(Obj.APTemplatingAxes,'equal');
%                     end
%                 end
%                 %=========================================================================================
                
                
                %                     case 'ML'
                if ShowCropped
                    
                    DATA=Obj.MLPlane.returnVisImage(1,1);
                    
                    if UpdateTemplate
                        Im=Obj.MLPlane.Image(1:end,1:end);
%                         DATAtemp=cat(3,Im,Im*0.95+10,Im*0.95+10);
                    end
                    
                else
                    
                    DATA=Obj.MLPlane.returnVisImage(2,1);
                    
                    if UpdateTemplate
                        Im=Obj.MLPlane.Image(1:end,1:end);
%                         DATAtemp=cat(3,Im,Im*0.95+10,Im*0.95+10);
                    end
                end
                
                
                %                         end
                if WithGraphics
                    set(Obj.MLViewImageHandle,'XData',[1,size(DATA,2)],'YData',[1,size(DATA,1)]);
                    set(Obj.MLViewImageHandle,'CData',DATA);
                end
                
                
%                 %=========================================================================================
%                 if UpdateTemplate
%                     % Update the corresponding ML tempplating view:
%                     Sz=size(DATAtemp);
%                     XD=round(Sz(2)/2);
%                     YD=round(Sz(1)/2);
%                     Xdata=[-XD, XD;-XD,XD];
%                     Ydata=[YD,YD;-YD,-YD];
%                     Zdata=[0,0;0,0];
%                     if WithGraphics
%                         set(Obj.MLSurfHandle,'XData',Xdata,'YData',Ydata,'ZData',Zdata);
%                     end
%                     Dat=DATAtemp;
%                     %                         DATAtemp(:,:,2)=Dat;
%                     %                         DATAtemp(:,:,3)=Dat;
%                     if WithGraphics
%                         set(Obj.MLSurfHandle,'CData',DATAtemp(1:temp_n:end,1:temp_n:end,:));
%                         xlim(Obj.MLTemplatingAxes,[-XD,XD]);
%                         ylim(Obj.MLTemplatingAxes,[-YD,YD]);
%                         axis(Obj.MLTemplatingAxes,'vis3d'); %axis(Obj.MLTemplatingAxes,'tight');
%                         axis(Obj.MLTemplatingAxes,'equal');
%                     end
%                 end
%                 %=========================================================================================
                
                if ShowCropped
                    %                             axes(MLaxesHandle);
                    if WithGraphics
                        axis(MLaxesHandle,'equal'); %axis tight;
                        xlim(MLaxesHandle,'auto');
                        ylim(MLaxesHandle,'auto');
                    end
                end
                
            end

            if UpdateCircle
                
                if UpdateCircle==2 % flag to remove the circles from both views
                    if WithGraphics
                        Pnts=get(Obj.APViewCircleHandle, 'XData');
                        set(Obj.APViewCircleHandle, 'XData', Pnts*0, 'YData', Pnts*0);
                        Pnts=get(Obj.MLViewCircleHandle, 'XData');
                        set(Obj.MLViewCircleHandle, 'XData', Pnts*0, 'YData', Pnts*0);
                    end
                else

%                     if ShowCropped
%                         APOffset=[size(Obj.APPlane.ImageOutput,2)/n/2,size(Obj.APPlane.ImageOutput,1)/n/2];
%                         MLOffset=[size(Obj.MLPlane.ImageOutput,1)/n/2,size(Obj.MLPlane.ImageOutput,2)/n/2];
%                     else                
%                         APOffset=[size(Obj.APPlane.ImageOutput,2)/n/2,size(Obj.APPlane.ImageOutput,1)/n/2];
%                         MLOffset=[size(Obj.MLPlane.ImageOutput,1)/n/2,size(Obj.MLPlane.ImageOutput,2)/n/2];
%                     end
                                        
                    if ShowCropped
                        APOffset=[size(Obj.APPlane.Image,2)/n/2,size(Obj.APPlane.Image,1)/n/2];
                        MLOffset=[size(Obj.MLPlane.Image,1)/n/2,size(Obj.MLPlane.Image,2)/n/2];
                    else                
                        APOffset=[size(Obj.APPlane.Image,2)/n/2,size(Obj.APPlane.Image,1)/n/2];
                        MLOffset=[size(Obj.MLPlane.Image,1)/n/2,size(Obj.MLPlane.Image,2)/n/2];
                    end                    


                    switch Obj.CurrentPlane
                        case 'AP'
                            Spot=Obj.returnSpot(ImgPlane,1);
                            C=Obj.APPlane.CroppedCoords;
                            if isempty(C)
                                C=[0,0,0,0];
                            end
                            CropCent=[(-C(3)+C(4))/2,(-C(1)+C(2))/2]/n;
                            CropCentroid=[(C(3)+C(4))/2,(C(1)+C(2))/2]/n;

                            DD=CropCentroid-APOffset;
                            delta_Offset=CropCent+ [-DD(1),DD(2)];

                            Tar=Spot.CameraTarget;
                            Pos=Spot.CameraPosition;
                            Up=Spot.CameraUp;
                            VA=Spot.CameraViewAngle;
                            Rad=Spot.Radius;

                            APPos2D=delta_Offset+[-Spot.Position2D(2),Spot.Position2D(1)]/n;
                            Rad2D=Spot.Rad2D/n;
                            %Circular point:
                            a=[0:5:360]';
                            c=0.9;
                            Pnts(:,1)=c*Rad2D*cosd(a);
                            Pnts(:,2)=c*Rad2D*sind(a);
                            Pnts(:,3)=Pnts(:,1)*0;
                            
                            if WithGraphics
                                set(Obj.APViewCircleHandle, 'XData',Pnts(:,1)+APPos2D(2), 'YData', Pnts(:,2)+APPos2D(1));
                                set(Obj.MLViewCircleHandle, 'XData', Pnts(:,1)*0, 'YData', Pnts(:,2)*0);
                            end
                        case 'ML'
                            Spot=Obj.returnSpot(ImgPlane,2);
                            C=Obj.MLPlane.CroppedCoords;
                            if isempty(C)
                                C=[0,0,0,0];
                            end
                            CropCent=[(-C(1)+C(2))/2,(-C(3)+C(4))/2]/n;
                            CropCentroid=[(C(1)+C(2))/2,(C(3)+C(4))/2]/n;
                            DD=CropCentroid-MLOffset;
                            delta_Offset=CropCent+ [DD(1),DD(2)];                        

                            Tar=Spot.CameraTarget;
                            Pos=Spot.CameraPosition;
                            Up=Spot.CameraUp;
                            VA=Spot.CameraViewAngle;
                            Rad=Spot.Radius;
                            
                            MLPos2D=delta_Offset+[Spot.Position2D(2),Spot.Position2D(1)]/n;                            
                            Rad2D=Spot.Rad2D/n;
                            
                            %Circular point:
                            a=[0:5:360]';
                            c=0.9;
                            Pnts(:,1)=c*Rad2D*cosd(a);
                            Pnts(:,2)=c*Rad2D*sind(a);
                            Pnts(:,3)=Pnts(:,1)*0;
                            if WithGraphics
                                set(Obj.APViewCircleHandle, 'XData', Pnts(:,1)*0, 'YData', Pnts(:,2)*0);
                                set(Obj.MLViewCircleHandle, 'XData', Pnts(:,1)+MLPos2D(2), 'YData', Pnts(:,2)+MLPos2D(1));
                            end
                    end
                    
                end
                
                drawnow;

            end
            
        end
         
        function [SurfHandle,CircleHandle]= showCameraView(Obj, ImgPlane, AxesHandle, WithGraphics)

            SurfHandle=[];
            CircleHandle=[];
            switch Obj.CurrentPlane
                case 'AP'
                    Sz=size(Obj.APPlane.Image); %Output);
                    XD=round(Sz(2)/2);
                    YD=round(Sz(1)/2);
                    Xdata=[-XD, XD;-XD,XD];
                    Ydata=[YD,YD;-YD,-YD];
                    Zdata=[0,0;0,0];
                case 'ML'
                    Sz=size(Obj.MLPlane.Image); %Output);
                    XD=round(Sz(2)/2);
                    YD=round(Sz(1)/2);
                    Xdata=[-XD, XD;-XD,XD];
                    Ydata=[YD,YD;-YD,-YD];
                    Zdata=[0,0;0,0];
            end            
            %Sz=size(imageData);

            if WithGraphics
                axes(AxesHandle);
                shading FLAT
                axis vis3d;
                axis equal;
                hold on;
            end
            
            nn=Obj.CameraDownSampleSize;
            switch Obj.CurrentPlane
                case 'AP'
                    if WithGraphics
                        SurfHandle = surface('XData',Xdata,'YData',Ydata,...
                            'ZData',Zdata,'CData',Obj.APPlane.ImageOutput(1:nn:end,1:nn:end,:),...
                            'FaceColor','texturemap','EdgeColor','none','FaceAlpha',1,'CDataMapping','direct'); %0.5);
                    end
                    Spot=Obj.returnSpot(ImgPlane,1);
                case 'ML'
                    if WithGraphics
                        SurfHandle = surface('XData',Xdata,'YData',Ydata,...
                            'ZData',Zdata,'CData',Obj.MLPlane.ImageOutput(1:nn:end,1:nn:end,:),...
                            'FaceColor','texturemap','EdgeColor','none','FaceAlpha',1,'CDataMapping','direct'); %0.5);
                        Spot=Obj.returnSpot(ImgPlane,2);
                    end
            end

            Tar=Spot.CameraTarget;
            Pos=Spot.CameraPosition;
            Up=Spot.CameraUp;
            VA=Spot.CameraViewAngle;
            Rad=Spot.Radius;
            
            %Circular point:
            a=[0:5:360]';
            c=0.9;
            Pnts(:,1)=c*Rad*cosd(a);
            Pnts(:,2)=c*Rad*sind(a);
            Pnts(:,3)=Pnts(:,1)*0;
            % axes(AxesHandle);
            if WithGraphics
                CircleHandle=plot3(AxesHandle, Pnts(:,1)+Tar(1),Pnts(:,2)+Tar(2),Pnts(:,3)+Tar(3),'y--','LineWidth',3);
            end
            
            %axes(AxesHandle);
            if WithGraphics
                campos(AxesHandle,Pos);
                camtarget(AxesHandle,Tar);
                camup(AxesHandle,Up);
                camva(AxesHandle,VA);
            end
            
            % viewfinder            
            pixPerMM = Obj.APPlane.PixelPerMM;
            Sz=size(Obj.APPlane.Image); %Output);
%             Obj.APViewFinder.init(pixPerMM,'A',Sz);
            
            Sz=size(Obj.MLPlane.Image); %Output);
%             Obj.MLViewFinder.init(pixPerMM,'M',Sz);
            
            
           % set(handles.ViewFinderAxes,'units','normalized');
           % set(handles.ViewFinderAxes,'OuterPosition',[0.05,0.05,0.9,0.9]);
           % set(AxesHandle,'units','normalized');
           % set(AxesHandle,'OuterPosition',[0.05,0.05,0.9,0.9]);
            
%             % Save Camera Handle properties here:
%             Obj.CameraSurfHandle=SurfHandle;
%             Obj.CameraCircleHandle=CircleHandle;

        end        
               
        function updateCameraView(Obj, ImgPlane, AxesHandle, UpdateData, UpdateCircle, WithGraphics)

            switch Obj.CurrentPlane
                case 'AP'
                    Sz=size(Obj.APPlane.Image); %Output);
                    XD=round(Sz(2)/2);
                    YD=round(Sz(1)/2);
                    Xdata=[-XD, XD;-XD,XD];
                    Ydata=[YD,YD;-YD,-YD];
                    Zdata=[0,0;0,0];
                    ViewIndex=1;
%                     Obj.APViewFinder.enableView('moving',1);
%                     Obj.MLViewFinder.enableView('moving',0);
%                     ViewFinder = Obj.APViewFinder;
                case 'ML'
                    Sz=size(Obj.MLPlane.Image); %Output);
                    XD=round(Sz(2)/2);
                    YD=round(Sz(1)/2);
                    Xdata=[-XD, XD;-XD,XD];
                    Ydata=[YD,YD;-YD,-YD];
                    Zdata=[0,0;0,0];                
                    ViewIndex=2;
%                     Obj.APViewFinder.enableView('moving',0);
%                     Obj.MLViewFinder.enableView('moving',1);
%                     ViewFinder = Obj.MLViewFinder;                    
            end            

            nn=Obj.CameraDownSampleSize;
            
            switch Obj.CurrentPlane
                case 'AP'
%                     SurfHandle = surface('XData',Xdata,'YData',Ydata,...
%                         'ZData',Zdata,'CData',Obj.APPlane.ImageOutput,...
%                         'FaceColor','texturemap','EdgeColor','none','FaceAlpha',1,'CDataMapping','direct'); %0.5);
                    if UpdateData
                        if WithGraphics
                            set(Obj.CameraSurfHandle,'XData',Xdata,'YData',Ydata,'CData',Obj.APPlane.ImageOutput(1:nn:end,1:nn:end,:));
                        end
%                         set(Obj.CameraSurfHandle,'XData',Xdata,'YData',Ydata,'CData',Obj.APPlane.returnVisImage(2));
                    end
                    %set(Obj.CameraSurfHandle,'CData',Obj.APPlane.ImageOutput);
                    if UpdateCircle
                        Spot=Obj.returnSpot(ImgPlane,ViewIndex);
                    end
                    
                    WW_Mark_OffsetAngle=0;
                case 'ML'
%                     SurfHandle = surface('XData',Xdata,'YData',Ydata,...
%                         'ZData',Zdata,'CData',Obj.MLPlane.ImageOutput,...
%                         'FaceColor','texturemap','EdgeColor','none','FaceAlpha',1,'CDataMapping','direct'); %0.5);
                    if UpdateData
                        if WithGraphics
                            set(Obj.CameraSurfHandle,'XData',Xdata,'YData',Ydata,'CData',Obj.MLPlane.ImageOutput(1:nn:end,1:nn:end,:));
                        end
%                         set(Obj.CameraSurfHandle,'XData',Xdata,'YData',Ydata,'CData',Obj.MLPlane.returnVisImage(2));
                    end
                    %set(Obj.CameraSurfHandle,'CData',Obj.MLPlane.ImageOutput);
                    if UpdateCircle                        
%                         Spot=Obj.returnSpot(ImgPlane,2);
                        Spot=Obj.returnSpot(ImgPlane,ViewIndex);
                    end
                    WW_Mark_OffsetAngle=90;
            end

            if UpdateCircle
                
                
                Tar=Spot.CameraTarget;
                Pos=Spot.CameraPosition;
                Up=Spot.CameraUp;
                VA=Spot.CameraViewAngle;
                Rad=Spot.Radius;
                Roll=Spot.Roll;
                % shoudl fix zooming in too far with the marker, if issues
                % remove
                if VA >= 140 
                    VA = 140;
                elseif VA <= 1
                    VA = 1;
                end
                %%%%%%%
                surfSize = abs(Tar(3)+Pos(3)) * 4 * tand(VA/2);
                errFlag = false;
                if Tar(1)> XD
                    Tar(1) = XD;
                    Pos(1) = XD;
                    if WithGraphics
                        ViewFinder.updateViewFinder([Pos(1),Pos(2),Pos(3),Roll],surfSize);
                        ViewFinder.updateWarning(1,'On');
                    end
                    errFlag = true;
                elseif Tar(1)< -XD
                    Tar(1) = -XD;
                    Pos(1) = -XD;
                    if WithGraphics
                        ViewFinder.updateViewFinder([Pos(1),Pos(2),Pos(3),Roll],surfSize);
                        ViewFinder.updateWarning(1,'On');
                    end
                    errFlag = true;
                end
                if Tar(2)> YD
                    Tar(2) = YD;
                    Pos(2) = YD;
                    if WithGraphics
                        ViewFinder.updateViewFinder([Pos(1),Pos(2),Pos(3),Roll],surfSize);
                        ViewFinder.updateWarning(1,'On');
                    end
                    errFlag = true;
                elseif Tar(2)< -YD
                    Tar(2) = -YD;
                    Pos(2) = -YD;
                    if WithGraphics
                        ViewFinder.updateViewFinder([Pos(1),Pos(2),Pos(3),Roll],surfSize);
                        ViewFinder.updateWarning(1,'On');
                    end
                    errFlag = true;
                end
                if ~errFlag
                    if WithGraphics
                        ViewFinder.updateViewFinder([Pos(1),Pos(2),Pos(3),Roll],surfSize);
                        ViewFinder.updateWarning(1,'Off');
                    end
                end
                %Circular point:
                ang=3;
%                 a=[ang:5:360-ang]';
                a=[0:5:360]';
                c=0.9;
                Pnts(:,1)=c*Rad*cosd(a+Roll);
                Pnts(:,2)=c*Rad*sind(a+Roll);
                Pnts(:,3)=Pnts(:,1)*0;
                
%                 ang=5;
                c=1;
                N1=[c*Rad*cosd(ang+Roll+WW_Mark_OffsetAngle),c*Rad*sind(ang+Roll+WW_Mark_OffsetAngle),0];
                N3=[c*Rad*cosd(-ang+Roll+WW_Mark_OffsetAngle),c*Rad*sind(-ang+Roll+WW_Mark_OffsetAngle),0];
                c=1.2;
                N2=[c*Rad*cosd(ang+Roll+WW_Mark_OffsetAngle),c*Rad*sind(ang+Roll+WW_Mark_OffsetAngle),0];
                N4=[c*Rad*cosd(-ang+Roll+WW_Mark_OffsetAngle),c*Rad*sind(-ang+Roll+WW_Mark_OffsetAngle),0];
                
                c=1;
                N1p=[c*Rad*cosd(ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(ang+Roll+WW_Mark_OffsetAngle+180),0];
                N3p=[c*Rad*cosd(-ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(-ang+Roll+WW_Mark_OffsetAngle+180),0];
                c=1.2;
                N2p=[c*Rad*cosd(ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(ang+Roll+WW_Mark_OffsetAngle+180),0];
                N4p=[c*Rad*cosd(-ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(-ang+Roll+WW_Mark_OffsetAngle+180),0];
                
                sth=0.9;          
                thk=1.3;
                N1c=[sth*Rad*cosd(WW_Mark_OffsetAngle),sth*Rad*sind(WW_Mark_OffsetAngle),0];
                N2c=[thk*Rad*cosd(WW_Mark_OffsetAngle),thk*Rad*sind(WW_Mark_OffsetAngle),0];
                
                N3c=[-sth*Rad*cosd(WW_Mark_OffsetAngle),-sth*Rad*sind(WW_Mark_OffsetAngle),0];
                N4c=[-thk*Rad*cosd(WW_Mark_OffsetAngle),-thk*Rad*sind(WW_Mark_OffsetAngle),0];

                
%                 c=1;
%                 N1p=[c*Rad*cosd(ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(ang+Roll+WW_Mark_OffsetAngle+180),0];
%                 N3p=[c*Rad*cosd(-ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(-ang+Roll+WW_Mark_OffsetAngle+180),0];
%                 c=1.3;
%                 N2p=[c*Rad*cosd(ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(ang+Roll+WW_Mark_OffsetAngle+180),0];
%                 N4p=[c*Rad*cosd(-ang+Roll+WW_Mark_OffsetAngle+180),c*Rad*sind(-ang+Roll+WW_Mark_OffsetAngle+180),0];                
%                 
                %             axes(AxesHandle);
                % CircleHandle=plot3(Pnts(:,1)+Tar(1),Pnts(:,2)+Tar(2),Pnts(:,3)+Tar(3),'y--','LineWidth',3);
                
                Pnts=[Pnts;[nan,nan,nan];N1;N2;[nan,nan,nan];N3;N4;[nan nan nan];N2p;N1p;[nan,nan,nan];N3p;N4p;[nan nan nan];N2c;N1c;[nan,nan,nan];N3c;N4c];
                if UpdateCircle==2 % flag to remove the circles from both views
                    Pnts=get(Obj.CameraCircleHandle, 'XData');
                    if WithGraphics
                        set(Obj.CameraCircleHandle, 'XData', Pnts*0, 'YData', Pnts*0);
                    end
                else
                    if WithGraphics
                        set(Obj.CameraCircleHandle,'XData',Pnts(:,1)+Tar(1),'YData',Pnts(:,2)+Tar(2),'ZData',Pnts(:,3)+Tar(3));
                    end
                end
                
                %                 axes(AxesHandle);
                
                if WithGraphics
                    campos(AxesHandle,Pos);
                    camtarget(AxesHandle,Tar);
                    camup(AxesHandle,Up);
                    camva(AxesHandle,VA);                                


                    drawnow;
                end
                
                %             % Save Camera Handle properties here:
                %             Obj.CameraSurfHandle=SurfHandle;
                %             Obj.CameraCircleHandle=CircleHandle;
            end
            
        end     
        
        function limit = updateCameraZoom(Obj,ImgPlane, AxesHandle)
            
            switch Obj.CurrentPlane
                case 'AP'
                    ViewIndex=1;
%                     ViewFinder = Obj.APViewFinder;
                case 'ML'
                    ViewIndex=2;
%                     ViewFinder = Obj.MLViewFinder;
            end  
            Spot=Obj.returnSpot(ImgPlane,ViewIndex);
            VA=Spot.CameraViewAngle;
            if VA >= 140 || VA <= 1
                limit = false;
            else
                limit = true;
                camva(AxesHandle,VA);
                Tar=Spot.CameraTarget;
                Pos=Spot.CameraPosition;
                surfSize = abs(Tar(3)+Pos(3)) * 4 * tand(VA/2);
%                 ViewFinder.updateViewFinder(ViewFinder.position,surfSize);
            end

        end
        
        function emptyViews(Obj)
%             % Reset the image plane info:
%             Obj.APPlane.Image=Obj.APPlane.Image*0;
%             Obj.APPlane.ImageOutput=Obj.APPlane.GridImage*0;
%             Obj.MLPlane.Image=Obj.MLPlane.Image*0;
%             Obj.MLPlane.ImageOutput=Obj.MLPlane.GridImage*0;   
        end
            
        function IPs=breakupImagePlane(Obj,ImPlane)
            
            IPs={};
            %StitchIndex=varargin{1};
            % ImPlane=Obj;
            
            ImPlane.Image=CropCircle(ImPlane.Image,95,85); % get rid of the border first
            
            HDivs=6;
            WDivs=6;
            OVR=0; %5; % pixel overlap between the tiles
            
            ImgW=size(ImPlane.Image,2);
            ImgH=size(ImPlane.Image,1);
            
            %             nDivs=WDivs/2;
            %             mDivs=HDivs/2;
            
            dW=(ImgW/WDivs);
            dH=(ImgH/HDivs);
            dW2=dW/2;
            dH2=dH/2;
            Vz=ImPlane.ImageNormal;
            Vy=ImPlane.ImageUp;
            Vx=cross(Vy,Vz);
            Vx=Vx/norm(Vx);
            PS=ImPlane.PixelScale;
            cc=0;
            %for xx=-nW2:nW2
            for xx=(-ImgW/2+dW2):dW:(ImgW/2-dW2)
                for yy=(ImgH/2-dH2):-dH:(-ImgH/2+dH2)
                    Cnt=ImPlane.CentrePoint + yy*Vy*PS + xx* Vx*PS;
                    % Cnt=ImPlane.CentrePoint + xx*Vy*PS - yy*Vx*PS;
                    cc=cc+1;
                    Str=sprintf('ImagePlane:%.2g - %.2g',xx,yy);
                    IPs{cc}=ImagePlane(Str);
                    IPs{cc}.ImageNormal=ImPlane.ImageNormal;
                    IPs{cc}.ImageUp=ImPlane.ImageUp;
                    IPs{cc}.PixelScale=ImPlane.PixelScale;
                    IPs{cc}.SourcePosition=ImPlane.SourcePosition;
                    IPs{cc}.CentrePoint=Cnt;
                    IPs{cc}.PrincipalDist=ImPlane.PrincipalDist;
                    IPs{cc}.Xoffset=0;
                    IPs{cc}.Yoffset=0;
%                     X1=xx-dW2+ImgW/2+1;
%                     X2=X1+2*dW2-1;
                    X1=xx-dW2+ImgW/2+1;
                    X2=X1+2*dW2-1;
%                     Y1=yy-dH2+ImgH/2+1; 
%                     Y2=Y1+2*dH2-1;
%                     jj=ImgH/2+yy;  
                    Y1=(yy+dH2)+ImgH/2; 
                    Y2=Y1-2*dH2+1;    
                    J2=ImgH-Y2+1;
                    J1=ImgH-Y1+1;
%                     IPs{cc}.Image=ImPlane.Image(X1:X2,Y1:Y2);
%                     [X1,X2,J2,J1]
                    % Add the overlapping pixels before and after the bounds:
                    X1=X1-OVR;
                    X2=X2+OVR;
                    J1=J1-OVR;
                    J2=J2+OVR;
                    if X1<1
                        X1=1;
                    end
                    if J1<1;
                        J1=1;
                    end
                    if X2>ImgW
                        X2=round(ImgW);
                    end
                    if J2>ImgH
                        J2=round(ImgH);
                    end
                    
                    IPs{cc}.Image=ImPlane.Image(J1:J2 ,X1:X2);
%                     IPs{cc}.Image=ImPlane.Image(X1:X2, J1:J2);
                    
%                     IPs{cc}.Image=ImPlane.Image(X1:X2,Y2:Y1);
                end
            end
            
        end        
        

        function [Img]=returnBiPlanarCroppedImage(Obj)
            ML=Obj.MLPlane.returnCroppedImage;
            APR=imrotate(Obj.APPlane.returnCroppedImage,90);
            Spacer=zeros(20,size(ML,2))+50;
            Img=cat(1,ML,Spacer,APR);
        end
        
        function getFolders(Obj,Exam)
            Obj.GridFolder=Exam.GridFolder;
            Obj.TemplateFolder=Exam.TemplateFolder;
            Obj.OutputFolder=Exam.OutputFolder;
            Obj.PreopFolder=Exam.PreopFolder;
        end
        
        function updateGrids(Obj)
            Obj.APPlane.updateGrid(Obj.IsStandardViews);
            Obj.MLPlane.updateGrid(Obj.IsStandardViews);            
        end
        
        function updateGridCoords(Obj,GridOrigin,GridVx, GridVy, GridVz) 
            Obj.Grid.Origin=GridOrigin;
            Obj.Grid.PlaneOrigin=GridOrigin;
            Obj.Grid.VxLineElement.Vector=GridVx;            
            Obj.Grid.VyLineElement.Vector=GridVy;
            Obj.Grid.VzLineElement.Vector=GridVz;            
            Obj.Grid.VxLineElement.Coords=GridOrigin;
            Obj.Grid.VyLineElement.Coords=GridOrigin;
            Obj.Grid.VzLineElement.Coords=GridOrigin;
        end
              
        function ImgPlane=alignImagePlane(Obj,Im)                     
            Org=Obj.Grid.PlaneOrigin;
            Vx=Obj.Grid.VxLineElement.Vector;
            Vy=Obj.Grid.VyLineElement.Vector;
            Vz=Obj.Grid.VzLineElement.Vector;
            ImgPlane=Obj.alignImagePlaneWithGridInfo(Im,Org,Vx,Vy,Vz);
        end
        
        function ImgPlane=alignImagePlaneWithGridInfo(Obj,Im,Org,Vx,Vy,Vz)
            
            % return a version of the ImagePlane that is aligned with
            % current Grid directions. This function uses the Grid*
            % Properties stored previously in the imageplane and uses that
            % to realign with current .Grid properties of the StitchMaster
            ImgPlane=ImagePlane(Im.Name);
            
            ImgPlane.Name=Im.Name;         
            ImgPlane.JTCalibFileName=Im.JTCalibFileName;
            ImgPlane.JTImageFileName=Im.JTImageFileName;
            ImgPlane.JTPoseFileName=Im.JTPoseFileName;
            ImgPlane.PrincipalDist=Im.PrincipalDist;
            ImgPlane.Xoffset=Im.Xoffset;
            ImgPlane.Yoffset=Im.Yoffset;
            ImgPlane.PixelScale=Im.PixelScale;
            ImgPlane.Features=Im.Features;
            ImgPlane.Path=Im.Path;
            ImgPlane.Image=Im.Image;

            ImgPlane.CarmPose=Im.CarmPose;
            
            ImgPlane.RefGridOrg=Im.RefGridOrg;
            ImgPlane.RefGridVx=Im.RefGridVx;
            ImgPlane.RefGridVy=Im.RefGridVy;
            ImgPlane.RefGridVz=Im.RefGridVz;

            P0=Im.RefGridOrg';
            P1=P0+(Im.RefGridVx'/norm(Im.RefGridVx'))*10;
            P2=P0+(Im.RefGridVy'/norm(Im.RefGridVy'))*20;
            P3=P0+(Im.RefGridVz'/norm(Im.RefGridVz'))*30;                        
            Pi=[P0,1;P1,1;P2,1;P3,1]';
            
%             Org=Obj.Grid.PlaneOrigin;
%             Vx=Obj.Grid.VxLineElement.Vector;
%             Vy=Obj.Grid.VyLineElement.Vector;
%             Vz=Obj.Grid.VzLineElement.Vector;
            
            Q0=Org';
            Q1=Q0+(Vx'/norm(Vx'))*10;
            Q2=Q0+(Vy'/norm(Vy'))*20;
            Q3=Q0+(Vz'/norm(Vz'))*30;                        
            Qi=[Q0,1;Q1,1;Q2,1;Q3,1]';
            
            T=Qi*inv(Pi); %Transformation Matrix to convert the coordinates
            
            Src=T*[Im.SourcePosition;1];
            ImgPlane.SourcePosition=Src(1:3);
            
            Cnt=T*[Im.CentrePoint;1];
            ImgPlane.CentrePoint=Cnt(1:3);

            Nrmi=Im.CentrePoint+10*Im.ImageNormal;
            Nrmj=T*[Nrmi;1];
            Nrm=Nrmj(1:3)-Cnt(1:3);
            Nrm=Nrm/norm(Nrm);            
            ImgPlane.ImageNormal=Nrm;
            
            Upi=Im.CentrePoint+10*Im.ImageUp;
            Upj=T*[Upi;1];
            Up=Upj(1:3)-Cnt(1:3);
            Up=Up/norm(Up);
            ImgPlane.ImageUp=Up;
            
        end
                 
        function setCircleCropRads(Obj,Rads) % Set the appearance of cropping circles.
           Obj.APPLane.ImageFuseCirRads=Rads;
           Obj.MLPLane.ImageFuseCirRads=Rads;
        end
         
        function Exam=realignGridBasedOnUI_woLoading(Obj,handles,APTransform,MLTransform,ShowGraphics,SyncOrigin)
            
            % SyncOrigin is the origin determined through the
            % 'synchronized' points on the AP and ML view with a cross
            % calibration tool is in use. Empty SyncOrigin means the
            % function should work as normal. If SyncOrigin has a value
            % then the value should be used as the new coordinate for the
            % grid
            
            if nargin<6
                SyncOrigin=[];  
            end
            

                   
            Exam=handles.TestData;
            
            Old_GOrg=Obj.Grid.Origin;
            Old_GVx=Obj.Grid.VxLineElement.Vector;
            Old_GVy=Obj.Grid.VyLineElement.Vector;
            Old_GVz=Obj.Grid.VzLineElement.Vector;

            AP_V1=Obj.APPlane.Vector1;
            AP_V2=Obj.APPlane.Vector2;
            AP_V3=Obj.APPlane.Vector3;
            AP_Org=Obj.APPlane.Origin;
            
            ML_V1=Obj.MLPlane.Vector1;
            ML_V2=Obj.MLPlane.Vector2;
            ML_V3=Obj.MLPlane.Vector3;
            ML_Org=Obj.MLPlane.Origin;
            
            % if the transformes are empty it could be assumed that the
            % method is being run during load/duplicate exams and therefore
            % the DepthNormals does not need to be calculated.
            if (isempty(APTransform) || isempty(MLTransform))
                IsInCalibnMode=0;
            else
                IsInCalibnMode=1;
                Obj.APPlane.DepthSegments=[];
                Obj.APPlane.DepthSegmentsNormal=[];
                Obj.MLPlane.DepthSegments=[];
                Obj.MLPlane.DepthSegmentsNormal=[];
            end
            % handle empty values for the parmeters:
            if isempty(APTransform)
                APTransform=[];
                APTransform.position=[0,0,0];
                APTransform.rotation=[0,0,0];
            end
            
            if isempty(MLTransform)
                MLTransform=[];
                MLTransform.position=[0,0,0];
                MLTransform.rotation=[0,0,0];
            end
            
            %Magnitude:
            APTransform.position(1)=APTransform.position(1)*-21;
            APTransform.position(2)=APTransform.position(2)*-21;
            MLTransform.position(1)=MLTransform.position(1)*21;
            MLTransform.position(2)=MLTransform.position(2)*-21;
            
            % Calculate the 3D Origin of the grid
            AP_XYZ=AP_Org + APTransform.position(1)*AP_V1 + APTransform.position(2)*AP_V2;
            ML_XYZ=ML_Org + MLTransform.position(2)*ML_V1 + MLTransform.position(1)*ML_V2;            
            
            A0=AP_XYZ;
            A1=AP_XYZ+AP_V3*100;
            B0=ML_XYZ;
            B1=ML_XYZ+ML_V3*100;
            
%             LA=[A0';A1'];
%             LB=[B0';B1'];
%             plot3(LA(:,1),LA(:,2),LA(:,3),'m');
%             plot3(LB(:,1),LB(:,2),LB(:,3),'c');
            
            
            [P1]=intersectLines3D(A0,A1,B0,B1);
            
            
            % Overridde the Grid_Org if the SyncOrigin value is provided.
            if isempty(SyncOrigin)
                Grid_Org=P1(1:3,1);
            else
                if ~iscolumn(SyncOrigin) % ensure a row matrix
                    SyncOrigin=SyncOrigin'; 
                end
                Grid_Org=SyncOrigin;
            end
            
                                    
            %% Refining the rotation angle in Spine Sync mode:  
            % This option was experimented with but later removed. The
            % idea is to keep the hips berfectly kept horizontal accounting
            % for correcting for parallax knowing the fluoro parallax
            % information.
            
%             % Check if CNodes and CNodesSource3D have index 1 and 4 then
%             % refine the Rotation Angle
%             RotAngleNeedsRefinement=0;
%             if ~isempty(handles.W1.StitchMaster.APPlane.CNodes) && ~isempty(handles.W1.StitchMaster.APPlane.CNodesSource3D)
%                 if (size(handles.W1.StitchMaster.APPlane.CNodes,1)>3) && (size(handles.W1.StitchMaster.APPlane.CNodesSource3D,1)>3)
%                     if (~isequal(handles.W1.StitchMaster.APPlane.CNodesSource3D(1,:),[-1,-1,-1])) && (~isequal(handles.W1.StitchMaster.APPlane.CNodesSource3D(4,:),[-1,-1,-1]))
%                         RotAngleNeedsRefinement=1;
%                     end
%                 end
%             end
%              
%             if (RotAngleNeedsRefinement)
%                 P2D_L=handles.W1.StitchMaster.APPlane.CNodes(1,:);
%                 T_LH3D=AP_Org-AP_V1*P2D_L(1)*21-AP_V2*P2D_L(2)*21;
%                 P2D_R=handles.W1.StitchMaster.APPlane.CNodes(4,:);
%                 T_RH3D=AP_Org-AP_V1*P2D_R(1)*21-AP_V2*P2D_R(2)*21;
%                 
%                 SrcPnt_LH=handles.W1.StitchMaster.APPlane.CNodesSource3D(1,:);
%                 SrcPnt_RH=handles.W1.StitchMaster.APPlane.CNodesSource3D(4,:);
%                 
%                 L_Lv=T_LH3D'-SrcPnt_LH;
%                 L_Lv=L_Lv/norm(L_Lv);
%                 N_LH=linePlaneIntersect(AP_V3',Grid_Org',L_Lv,SrcPnt_LH);
%                 
%                 R_Lv=T_RH3D'-SrcPnt_RH;
%                 R_Lv=R_Lv/norm(R_Lv);
%                 N_RH=linePlaneIntersect(AP_V3',Grid_Org',R_Lv,SrcPnt_RH);
%                 
%                 RelV=N_RH-N_LH;
%                 RelV=RelV/norm(RelV);
%                 RotAngle=acosd(dot(RelV,AP_V1'));
%                 AP_Rot=RotAngle*pi/180; %<<< Revised Rotation Angle
%             else
%                 AP_Rot=APTransform.rotation(3)*pi/180;
%             end
            %%                       
            AP_Rot=APTransform.rotation(3)*pi/180;            
            AP_V1n=AP_V1*cos(AP_Rot)+AP_V2*sin(AP_Rot);
            AP_V1n=AP_V1n/norm(AP_V1n);

            % Locate the Y Axis:
            AP_V2n=-AP_V1*sin(AP_Rot)+AP_V2*cos(AP_Rot);
            AP_V2n=AP_V2n/norm(AP_V2n);           
                       
            Grid_Vx=AP_V1n;
            
            %Locate the X Axis:
            ML_Rot=MLTransform.rotation(3)*pi/180;
            ML_V1n=ML_V1*cos(ML_Rot)+ML_V2*sin(ML_Rot);
            ML_V1n=ML_V1n/norm(ML_V1n);

            % Locate the Y Axis:
            ML_V2n=-ML_V1*sin(ML_Rot)+ML_V2*cos(ML_Rot);
            ML_V2n=ML_V2n/norm(ML_V2n);
                   
            
            [~,Grid_Vy,~]=plane_intersect(AP_V1n,AP_XYZ,ML_V2n,ML_XYZ);
            %correcting the sign of the Grid_Vy:
            VySign=1;
            if dot(Grid_Vy,AP_V2n)<0
                VySign=-1;
            end
            Grid_Vy=VySign*Grid_Vy;            
            Grid_Vy=Grid_Vy/norm(Grid_Vy);                   
            
            Grid_Vz=cross(Grid_Vx,Grid_Vy);
            Grid_Vz=Grid_Vz/norm(Grid_Vz);
            
            Obj.APImagePlane.CentrePoint=Grid_Org;
            Obj.APImagePlane.ImageUp=Grid_Vy;
            Obj.APImagePlane.ImageNormal=Grid_Vz;
            Obj.APImagePlane.SourcePosition=Grid_Org-1000*Grid_Vz;

            ML_Nrm=ML_V3-dot(ML_V3,Grid_Vy)*Grid_Vy;
            ML_Nrm=ML_Nrm/norm(ML_Nrm);

            Rz=cross(ML_Nrm,Grid_Vy);
            Rz=Rz/norm(Rz);
            Obj.MLImagePlane.CentrePoint=Grid_Org;
            Obj.MLImagePlane.ImageUp=Rz;
            Obj.MLImagePlane.ImageNormal=ML_Nrm;
            Obj.MLImagePlane.SourcePosition=Grid_Org-1000*ML_Nrm;             
            
            New_GOrg=Grid_Org;
            New_GVx=Grid_Vx;
            New_GVy=Grid_Vy;
            New_GVz=Grid_Vz;   


            %% Calculate Anatomy Segments (if this is in CalibMode)
            if (IsInCalibnMode)
                % AP Plane:
                if ~isempty(Obj.APPlane.AnatomySegments)
                   Obj.MLPlane.DepthSegments=Obj.APPlane.AnatomySegments3D;
                   if strcmpi(Obj.APCarmViewStr,'APOB') % depth segment normal should account for the angle of the images
                       OblqNorm=Obj.returnOblqLatNorm;
                       Obj.MLPlane.DepthSegmentsNormal=OblqNorm(1,:)';
                   else
                       Obj.MLPlane.DepthSegmentsNormal=Obj.APPlane.ImageNormal;
                   end
                else
                   Obj.MLPlane.DepthSegments=[];
                   Obj.MLPlane.DepthSegmentsNormal=[];
                end
                % ML Plane:
                if ~isempty(Obj.MLPlane.AnatomySegments)              
                   Obj.APPlane.DepthSegments=Obj.MLPlane.AnatomySegments3D;
                   if strcmpi(Obj.MLCarmViewStr,'MLOB')  % depth segment normal should account for the angle of the images
                       OblqNorm=Obj.returnOblqLatNorm;
                       Obj.APPlane.DepthSegmentsNormal=OblqNorm(2,:)';
                   else
                       Obj.APPlane.DepthSegmentsNormal=Obj.MLPlane.ImageNormal;
                   end
                else
                   Obj.APPlane.DepthSegments=[];
                   Obj.APPlane.DepthSegmentsNormal=[];
                end
            end
            %%
            
            if (ShowGraphics)
                figure;
                hold on;
                Obj.APPlane.imshow3d(1,0,0,0);
                Px=[Grid_Org';Grid_Org'+1000*Grid_Vx'];
                Py=[Grid_Org';Grid_Org'+1000*Grid_Vy'];
                Pz=[Grid_Org';Grid_Org'+1000*Grid_Vz'];
                plot3(Px(:,1),Px(:,2),Px(:,3),'r');
                plot3(Py(:,1),Py(:,2),Py(:,3),'g');
                plot3(Pz(:,1),Pz(:,2),Pz(:,3),'b');                
                axis equal;
                if ~isempty(Obj.APPlane.AnatomySegments)
                    plot3(Obj.APPlane.AnatomySegments3D(:,1),Obj.APPlane.AnatomySegments3D(:,2),Obj.APPlane.AnatomySegments3D(:,3),'r*');                    
                end               
            end       
            
            
            if (ShowGraphics)
                hold on;
                Obj.MLPlane.imshow3d(1,0,0,0);
                Px=[Grid_Org';Grid_Org'+1000*Grid_Vy'];
                Py=[Grid_Org';Grid_Org'+1000*Rz'];
                Pz=[Grid_Org';Grid_Org'+1000*ML_Nrm'];
                plot3(Px(:,1),Px(:,2),Px(:,3),'r');
                plot3(Py(:,1),Py(:,2),Py(:,3),'g');
                plot3(Pz(:,1),Pz(:,2),Pz(:,3),'b');
                axis equal;
                if ~isempty(Obj.MLPlane.AnatomySegments)
                    plot3(Obj.MLPlane.AnatomySegments3D(:,1),Obj.MLPlane.AnatomySegments3D(:,2),Obj.MLPlane.AnatomySegments3D(:,3),'r*');                    
                end               
            end              
            
            %clear the image planes before restitching the images
            Obj.APPlane.Image=Obj.APPlane.Image*0;
            Obj.MLPlane.Image=Obj.MLPlane.Image*0;
            Obj.APPlane.ImageSrcIndx=Obj.APPlane.ImageSrcIndx*0;
            Obj.MLPlane.ImageSrcIndx=Obj.MLPlane.ImageSrcIndx*0;
            Obj.APPlane.ImageDepthIndx=Obj.APPlane.ImageDepthIndx*0;
            Obj.MLPlane.ImageDepthIndx=Obj.MLPlane.ImageDepthIndx*0;            
                        
            %updating the grid location based on the new Grid coords
            Obj.updateGridCoords(New_GOrg,New_GVx,New_GVy,New_GVz);            
            Obj.updateGrids;
            
            % Faster revision of the imageplanes wo need for writting to the memory:
            if ~isempty(Exam)
                Exam.RefreshFilesWithNewGridInfo_wo_LoadSave(Obj,Old_GOrg,Old_GVx,Old_GVy,Old_GVz,New_GOrg,New_GVx,New_GVy,New_GVz,ShowGraphics);
            end
                      
            
        end
        
        function [SrcPosition3D, CentrePoint3D, ImageNormal3D]=realignGridBasedOnUI(Obj,handles,APTransform,MLTransform,Data,ShowGraphics)

            SrcPosition3D=[];
            CentrePoint3D=[];
            ImageNormal3D=[];
            
            [IsCrossSync,LocalizedOrigin3D,OriginIndices]=Obj.checkCrossSyncCalib(Data);
            
            % Realign the Grid and reprocess Images based on UI inputs
            % (Skip re-aligning for a special testmode where the
            % calibration should not move the Grid around.

            if IsCrossSync && ~isempty(OriginIndices)
                % Correct the cross sync and anatomy segments based on the calculated coordinate of the sync'ed origin:
                if iscolumn(LocalizedOrigin3D)
                    LocalizedOrigin3D=LocalizedOrigin3D';
                end
                % Adjust the grid origin now considering the correct localizeed point
                Exam=Obj.realignGridBasedOnUI_woLoading(handles,APTransform,MLTransform,ShowGraphics,LocalizedOrigin3D);
            else
                Exam=Obj.realignGridBasedOnUI_woLoading(handles,APTransform,MLTransform,ShowGraphics);
            end
                             
            %% Reprocess all the image planes based on the new Grid Info:
            if(ShowGraphics)
                wb_h=waitbar(0,'Reformating the images..');
%                 changeFigIcon(wb_h);
            end
            %Loop through Images and stitch them up:
            for ii= 1:numel(Exam.ImagePlanes)
                if(ShowGraphics)
                    waitbar(ii/numel(Exam.ImagePlanes));
                end
                
                if isempty(Exam.ImagePlanes(ii).CarmPose.Process)
                    Rot=Exam.ImagePlanes(ii).CarmPose.Rotation;
                else  %supporting the old Exams where .Process held the Rotation
                    Rot=Exam.ImagePlanes(ii).CarmPose.Process.Rotation;
                end                
                DR1=abs(Rot);
                DR2=abs(90+Rot);
                if DR1<DR2
                    if Obj.IsStandardViews
                        CP=1; %AP Plane
                    else
                        CP=2;
                    end
                else
                    if Obj.IsStandardViews
                        CP=2; %ML Plane
                    else
                        CP=1;
                    end
                end
                
                ErrCode=-1;
                % stitch the image if it should not be ignored, otherwise skip
                ShallIgnoreFrame=Exam.ShallIgnore(ii);
                % Selective show of the graphics below for debugging
                

                if ~ShallIgnoreFrame
                    
                    SelImagePlane=Obj.alignImagePlane(Exam.ImagePlanes(ii));
                    if isempty(SrcPosition3D)
                        SrcPosition3D=SelImagePlane.SourcePosition';
                    else
                        SrcPosition3D=[SrcPosition3D; SelImagePlane.SourcePosition'];
                    end
                    
                    if isempty(CentrePoint3D)
                        CentrePoint3D=handles.W1.StitchMaster.alignImagePlane(handles.TestData.ImagePlanes(ii)).CentrePoint';
                    else
                        CentrePoint3D=[CentrePoint3D; SelImagePlane.CentrePoint'];
                    end
                    
                    % create the Normal Vector Matrix
                    if isempty(ImageNormal3D)
                        ImageNormal3D=handles.W1.StitchMaster.alignImagePlane(handles.TestData.ImagePlanes(ii)).ImageNormal';
                    else
                        ImageNormal3D=[ImageNormal3D; SelImagePlane.ImageNormal'];
                    end
                    
%                     profile on;
                    [~,~,~,~,ErrCode]=Obj.addImage(CP,Obj.alignImagePlane(Exam.ImagePlanes(ii)),0,[],ShowGraphics,[],ii);
%                     profile viewer;
                    
                    % Update the progress indicator and send a respond back to the client
                    if ~isempty(handles.tc)
                        if strcmpi(handles.tc.status,'open')
                            handles.CmdP.CmdProgress=ii/numel(Exam.ImagePlanes);
                            if mod(ii,handles.CmdP.ClientResponseCycle)==0
                                clientBytesReceived([],[],handles); %,obj)
                            end
                        end
                    end
                                        
                    switch ErrCode
                        case 0
                            st=sprintf('Image # %i added..',ii);
                        case 1
                            st=sprintf('Image # %i NOT added: Out of bounds detected..',ii);
                        case 2
                            st=sprintf('Image # %i NOT added: Collinear corners of the box..',ii);
                    end

                    disp(st);
                end
            
            end
            if (ShowGraphics)
                close(wb_h);
            end
           
            % Update Graphics:
            Obj.cropViews(1);
            
            % calculate the depth index layer to be send to the frontend
            handles.W1.StitchMaster.APPlane.calcImageDepthIndex(0);
            handles.W1.StitchMaster.MLPlane.calcImageDepthIndex(0);

        end
        
        %experimental:
        function realignGridBasedOnUI_new_experimental(Obj,handles,APTransform,MLTransform,ShowGraphics)
            % Realign the Grid and reprocess Images based on UI inputs
            
            Exam=handles.TestData;
            
            Old_GOrg=Obj.Grid.Origin;
            Old_GVx=Obj.Grid.VxLineElement.Vector;
            Old_GVy=Obj.Grid.VyLineElement.Vector;
            Old_GVz=Obj.Grid.VzLineElement.Vector;

            AP_V1=Obj.APPlane.Vector1;
            AP_V2=Obj.APPlane.Vector2;
            AP_V3=Obj.APPlane.Vector3;
            AP_Org=Obj.APPlane.Origin;
            
            ML_V1=Obj.MLPlane.Vector1;
            ML_V2=Obj.MLPlane.Vector2;
            ML_V3=Obj.MLPlane.Vector3;
            ML_Org=Obj.MLPlane.Origin;
            
            %Magnitude:
            APTransform.position(1)=APTransform.position(1)*-21;
            APTransform.position(2)=APTransform.position(2)*-21;
            MLTransform.position(1)=MLTransform.position(1)*21;
            MLTransform.position(2)=MLTransform.position(2)*-21;
            % Calculate the 3D Origin of the grid:
            % AP_XYZ=AP_Org + APTransform.position(1)*AP_V1 + APTransform.position(2)*AP_V2;
            % ML_XYZ=ML_Org + MLTransform.position(2)*ML_V1 + MLTransform.position(1)*ML_V2;
            % dML_XYZ=ML_XYZ-AP_XYZ;
            % Location of the origin
            % Grid_Org=AP_XYZ + dot(dML_XYZ,AP_V3)*AP_V3; 
            
            
            % Calculate the 3D Origin of the grid
            AP_XYZ=AP_Org + APTransform.position(1)*AP_V1 + APTransform.position(2)*AP_V2;
            ML_XYZ=ML_Org + MLTransform.position(2)*ML_V1 + MLTransform.position(1)*ML_V2;            
            
            A0=AP_XYZ;
            A1=AP_XYZ+AP_V3*100;
            B0=ML_XYZ;
            B1=ML_XYZ+ML_V3*100;
            [P1]=intersectLines3D(A0,A1,B0,B1);
            Grid_Org=P1(1:3,1);
                        
            %Locate the X Axis:
            AP_Rot=APTransform.rotation(3)*pi/180;
            AP_V1n=AP_V1*cos(AP_Rot)+AP_V2*sin(AP_Rot);
            AP_V1n=AP_V1n/norm(AP_V1n);

            % Locate the Y Axis:
            AP_V2n=-AP_V1*sin(AP_Rot)+AP_V2*cos(AP_Rot);
            AP_V2n=AP_V2n/norm(AP_V2n);
            
            %AP_V3n=cross(AP_V1n,AP_V2n);
            %AP_V3n=AP_V3n/norm(AP_V3n);
                       
            Grid_Vx=AP_V1n;
            
            %Locate the X Axis:
            ML_Rot=MLTransform.rotation(3)*pi/180;
            ML_V1n=ML_V1*cos(ML_Rot)+ML_V2*sin(ML_Rot);
            ML_V1n=ML_V1n/norm(ML_V1n);

            % Locate the Y Axis:
            ML_V2n=-ML_V1*sin(ML_Rot)+ML_V2*cos(ML_Rot);
            ML_V2n=ML_V2n/norm(ML_V2n);
            
            %ML_V3n=cross(ML_V1n,ML_V2n);
            %ML_V3n=ML_V3n/norm(ML_V3n);            
            
            [~,Grid_Vy,~]=plane_intersect(AP_V1n,AP_XYZ,ML_V2n,ML_XYZ);
            %correcting the sign of the Grid_Vy:
            VySign=1;
            if dot(Grid_Vy,AP_V2n)<0
                VySign=-1;
            end
            Grid_Vy=VySign*Grid_Vy;            
            Grid_Vy=Grid_Vy/norm(Grid_Vy);
                     
            
            Grid_Vz=cross(Grid_Vx,Grid_Vy);
            Grid_Vz=Grid_Vz/norm(Grid_Vz);
            
            Obj.APImagePlane.CentrePoint=Grid_Org;
            Obj.APImagePlane.ImageUp=Grid_Vy;
            Obj.APImagePlane.ImageNormal=Grid_Vz;
            Obj.APImagePlane.SourcePosition=Grid_Org-1000*Grid_Vz;

            ML_Nrm=ML_V3-dot(ML_V3,Grid_Vy)*Grid_Vy;
            ML_Nrm=ML_Nrm/norm(ML_Nrm);

            Rz=cross(ML_Nrm,Grid_Vy);
            Rz=Rz/norm(Rz);
            Obj.MLImagePlane.CentrePoint=Grid_Org;
%             Obj.MLImagePlane.ImageUp=Rz;
%             Obj.MLImagePlane.ImageNormal=ML_Nrm;
%             Obj.MLImagePlane.SourcePosition=Grid_Org-1000*ML_Nrm;

%%
% New experimental equaltion:
            Obj.MLImagePlane.ImageUp=Grid_Vz;
            Obj.MLImagePlane.ImageNormal=Grid_Vx;
            Obj.MLImagePlane.SourcePosition=Grid_Org-1000*Grid_Vx;            
%%            
            
            New_GOrg=Grid_Org;
            New_GVx=Grid_Vx;
            New_GVy=Grid_Vy;
            New_GVz=Grid_Vz;   

            if (ShowGraphics)
                figure;
                subplot(1,2,1);
                hold on;
%                 Obj.APPlane.imshow3d(1,0,0,0);
                Obj.APPlane.imshow3d(1,0,0,0);
                Px=[Grid_Org';Grid_Org'+1000*Grid_Vx'];
                Py=[Grid_Org';Grid_Org'+1000*Grid_Vy'];
                Pz=[Grid_Org';Grid_Org'+1000*Grid_Vz'];
                plot3(Px(:,1),Px(:,2),Px(:,3),'r');
                plot3(Py(:,1),Py(:,2),Py(:,3),'g');
                plot3(Pz(:,1),Pz(:,2),Pz(:,3),'b');
                axis equal;
            end       
            
            
            if (ShowGraphics)
                subplot(1,2,2);
                hold on;
%                 Obj.APPlane.imshow3d(1,0,0,0);
                Obj.MLPlane.imshow3d(1,0,0,0);
                Px=[Grid_Org';Grid_Org'+1000*Grid_Vy'];
                Py=[Grid_Org';Grid_Org'+1000*Rz'];
                Pz=[Grid_Org';Grid_Org'+1000*ML_Nrm'];
                plot3(Px(:,1),Px(:,2),Px(:,3),'r');
                plot3(Py(:,1),Py(:,2),Py(:,3),'g');
                plot3(Pz(:,1),Pz(:,2),Pz(:,3),'b');
                axis equal;
            end              
            
            %clear the image planes before restitching the images
            Obj.APPlane.Image=Obj.APPlane.Image*0;
            Obj.MLPlane.Image=Obj.MLPlane.Image*0;
            Obj.APPlane.ImageSrcIndx=Obj.APPlane.ImageSrcIndx*0;
            Obj.MLPlane.ImageSrcIndx=Obj.MLPlane.ImageSrcIndx*0;
            Obj.APPlane.ImageDepthIndx=Obj.APPlane.ImageDepthIndx*0;
            Obj.MLPlane.ImageDepthIndx=Obj.MLPlane.ImageDepthIndx*0;            
                        
            
            %% Revise all the ImagePlane files one by one:
            % Exam.RefreshFilesWithNewGridInfo(Obj,Old_GOrg,Old_GVx,Old_GVy,Old_GVz,New_GOrg,New_GVx,New_GVy,New_GVz);
            
            %>>>>>>>>>>> The grid needs to be updated before updating the
            %texture? Sep 24,2018:
            %updating the grid location based on the new Grid coords
            Obj.updateGridCoords(New_GOrg,New_GVx,New_GVy,New_GVz);            
            Obj.updateGrids;
            
            % Faster revision of the imageplanes wo need for writting to
            % the memory:
            Exam.RefreshFilesWithNewGridInfo_wo_LoadSave(Obj,Old_GOrg,Old_GVx,Old_GVy,Old_GVz,New_GOrg,New_GVx,New_GVy,New_GVz,ShowGraphics);
            
                 
            %% Reprocess all the image planes based on the new Grid Info:
            if(ShowGraphics)
                wb_h=waitbar(0,'Reformating the images..');
                changeFigIcon(wb_h);
            end
            %Loop through Images and stitch them up:
            for ii= 1:numel(Exam.ImagePlanes)
                if(ShowGraphics)
                    waitbar(ii/numel(Exam.ImagePlanes));
                end
                
                if isempty(Exam.ImagePlanes(ii).CarmPose.Process)
                    Rot=Exam.ImagePlanes(ii).CarmPose.Rotation;
                else  %supporting the old Exams where .Process held the Rotation
                    Rot=Exam.ImagePlanes(ii).CarmPose.Process.Rotation;
                end  
                
                DR1=abs(Rot);
                DR2=abs(90+Rot);
                if DR1<DR2
                    if Obj.IsStandardViews
                        CP=1; %AP Plane
                    else
                        CP=2;
                    end
                else
                    if Obj.IsStandardViews
                        CP=2; %ML Plane
                    else
                        CP=1;
                    end
                end
                
                ErrCode=-1;
                % stitch the image if it should not be ignored, otherwise skip
                ShallIgnoreFrame=Exam.ShallIgnore(ii);
                if ~ShallIgnoreFrame
                    [~,~,~,~,ErrCode]=Obj.addImage(CP,Obj.alignImagePlane(Exam.ImagePlanes(ii)),0,[],0,[],ii);
                    
                    % Update the progress indicator and send a respond back to the client
                    if strcmpi(handles.tc.status,'open')                                            
                        handles.CmdP.CmdProgress=ii/numel(Exam.ImagePlanes);
                        if mod(ii,handles.CmdP.ClientResponseCycle)==0
                            clientBytesReceived([],[],handles); %,obj)
                        end                       
                    end
                                        
                    switch ErrCode
                        case 0
                            st=sprintf('Image # %i added..',ii);
                        case 1
                            st=sprintf('Image # %i NOT added: Out of bounds detected..',ii);
                        case 2
                            st=sprintf('Image # %i NOT added: Collinear corners of the box..',ii);
                    end

                    disp(st);
                end
            
            end
            if (ShowGraphics)
                close(wb_h);
            end
           
            % Update Graphics:
            Obj.cropViews(1);

  

            
        end
         
        function ploInputsFromUI(Obj,Exam,APTransform,MLTransform,ShowGraphics)
            
            Old_GOrg=Obj.Grid.Origin;
            Old_GVx=Obj.Grid.VxLineElement.Vector;
            Old_GVy=Obj.Grid.VyLineElement.Vector;
            Old_GVz=Obj.Grid.VzLineElement.Vector;

            AP_V1=Obj.APPlane.Vector1;
            AP_V2=Obj.APPlane.Vector2;
            AP_V3=Obj.APPlane.Vector3;
            AP_Org=Obj.APPlane.Origin;
            
            ML_V1=Obj.MLPlane.Vector1;
            ML_V2=Obj.MLPlane.Vector2;
            ML_V3=Obj.MLPlane.Vector3;
            ML_Org=Obj.MLPlane.Origin;
            
            %Magnitude:
            APTransform.position(1)=APTransform.position(1)*-21;
            APTransform.position(2)=APTransform.position(2)*-21;
            MLTransform.position(1)=MLTransform.position(1)*21;
            MLTransform.position(2)=MLTransform.position(2)*-21;
            % Calculate the 3D Origin of the grid
            AP_XYZ=AP_Org + APTransform.position(1)*AP_V1 + APTransform.position(2)*AP_V2;
            ML_XYZ=ML_Org + MLTransform.position(2)*ML_V1 + MLTransform.position(1)*ML_V2;

            Grid_Org=AP_XYZ;
            
            %Locate the X Axis:
            AP_Rot=APTransform.rotation(3)*pi/180;
            AP_V1n=AP_V1*cos(AP_Rot)+AP_V2*sin(AP_Rot);
            AP_V1n=AP_V1n/norm(AP_V1n);

            % Locate the Y Axis:
            AP_V2n=-AP_V1*sin(AP_Rot)+AP_V2*cos(AP_Rot);
            AP_V2n=AP_V2n/norm(AP_V2n);
            
            % Locate the Z Axis:
            Grid_Vx=AP_V1n;
            Grid_Vy=AP_V2n;            
            Grid_Vz=cross(Grid_Vx,Grid_Vy);
              
            
            if (ShowGraphics)
                figure;
                subplot(1,2,1);
                hold on;
%                 Obj.APPlane.imshow3d(1,0,0,0);
                Obj.APPlane.imshow3d(1,0,0,0);
                Px=[Grid_Org';Grid_Org'+1000*Grid_Vx'];
                Py=[Grid_Org';Grid_Org'+1000*Grid_Vy'];
                Pz=[Grid_Org';Grid_Org'+1000*Grid_Vz'];
                plot3(Px(:,1),Px(:,2),Px(:,3),'r');
                plot3(Py(:,1),Py(:,2),Py(:,3),'g');
                plot3(Pz(:,1),Pz(:,2),Pz(:,3),'b');
                axis equal;
            end            
            
            ML_XYZ=ML_Org + MLTransform.position(2)*ML_V1 + MLTransform.position(1)*ML_V2;

            Grid_Org=ML_XYZ;
            
            %Locate the X Axis:
            ML_Rot=MLTransform.rotation(3)*pi/180;
            ML_V1n=ML_V1*cos(ML_Rot)+ML_V2*sin(ML_Rot);
            ML_V1n=ML_V1n/norm(ML_V1n);

            % Locate the Y Axis:
            ML_V2n=-ML_V1*sin(ML_Rot)+ML_V2*cos(ML_Rot);
            ML_V2n=ML_V2n/norm(ML_V2n);
            
            % Locate the Z Axis:
            Grid_Vx=ML_V1n;
            Grid_Vy=ML_V2n;            
            Grid_Vz=cross(Grid_Vx,Grid_Vy);
              
            
            if (ShowGraphics)
%                 figure;
                subplot(1,2,2);
                hold on;
%                 Obj.APPlane.imshow3d(1,0,0,0);
                Obj.MLPlane.imshow3d(1,0,0,0);
                Px=[Grid_Org';Grid_Org'+1000*Grid_Vx'];
                Py=[Grid_Org';Grid_Org'+1000*Grid_Vy'];
                Pz=[Grid_Org';Grid_Org'+1000*Grid_Vz'];
                plot3(Px(:,1),Px(:,2),Px(:,3),'r');
                plot3(Py(:,1),Py(:,2),Py(:,3),'g');
                plot3(Pz(:,1),Pz(:,2),Pz(:,3),'b');
                axis equal;
            end               

            
        end
        
        function GnP=exportGridAndPlanes(Obj)
        % Exporing most essential data defining plane and grids
            % copy the structure for saving
            GnP=[];

            GnP.GridOrigin=Obj.Grid.Origin;
            GnP.GridPlaneOrigin=Obj.Grid.PlaneOrigin;
            GnP.GridVx=Obj.Grid.VxLineElement.Vector;
            GnP.GridVy=Obj.Grid.VyLineElement.Vector;
            GnP.GridVz=Obj.Grid.VzLineElement.Vector; 
            
            GnP.APCarmViewStr=Obj.APCarmViewStr;
            GnP.MLCarmViewStr=Obj.MLCarmViewStr;
            GnP.IsStandardViews=Obj.IsStandardViews;
            GnP.ORSetupString=Obj.ORSetupString;
            GnP.APObliqueAngle=Obj.APObliqueAngle;
            GnP.MLObliqueAngle=Obj.MLObliqueAngle;
                        
            %save(Filename,'-mat','GnP');
        end
        
        function Out=importGridAndPlanes(Obj,GnP)
        % Importing most essential data defining plane and grids
            
            Out=0;
            %load the parameters back:
            Obj.Grid.Origin=GnP.GridOrigin;
            Obj.Grid.PlaneOrigin=GnP.GridPlaneOrigin;
            Obj.Grid.VxLineElement.Vector=GnP.GridVx;            
            Obj.Grid.VyLineElement.Vector=GnP.GridVy;
            Obj.Grid.VzLineElement.Vector=GnP.GridVz;            
            Obj.Grid.VxLineElement.Coords=GnP.GridOrigin;
            Obj.Grid.VyLineElement.Coords=GnP.GridOrigin;
            Obj.Grid.VzLineElement.Coords=GnP.GridOrigin;
            
            % Applying the additional GnP if they exist:
            Obj.APCarmViewStr=retunFieldIfExist(GnP,'APCarmViewStr','AP');
            Obj.MLCarmViewStr=retunFieldIfExist(GnP,'MLCarmViewStr','ML');
            Obj.IsStandardViews=retunFieldIfExist(GnP,'IsStandardViews',1);
            Obj.ORSetupString=retunFieldIfExist(GnP,'ORSetupString','CRHDPSTBB1');
            Obj.APObliqueAngle=retunFieldIfExist(GnP,'APObliqueAngle',0);
            Obj.MLObliqueAngle=retunFieldIfExist(GnP,'MLObliqueAngle',-90);
            
            % apply the parameters to both AP and ML planes:
            Obj.updateGrids;
            Out=1;
        end        
        
        function Out=clearImages(Obj)
            Out=[];
            Obj.APPlane.Image=Obj.APPlane.Image*0;
            Obj.APPlane.ImageSrcIndx=Obj.APPlane.ImageSrcIndx*0;
            Obj.APPlane.ImageDepthIndx=Obj.APPlane.ImageDepthIndx*0;
            Obj.APPlane.AnatomySegments=[];
            Obj.APPlane.AnatomySegments3D=[];
            Obj.APPlane.DepthSegments=[];
            Obj.APPlane.DepthSegmentsNormal=[];
            
            Obj.MLPlane.Image=Obj.MLPlane.Image*0;
            Obj.MLPlane.ImageSrcIndx=Obj.MLPlane.ImageSrcIndx*0;
            Obj.MLPlane.ImageDepthIndx=Obj.MLPlane.ImageDepthIndx*0;
            Obj.MLPlane.AnatomySegments=[];
            Obj.MLPlane.AnatomySegments3D=[];
            Obj.MLPlane.DepthSegments=[];
            Obj.MLPlane.DepthSegmentsNormal=[];                   
                       
            Out=1;
        end
        
        function Ang=returnCarmAngle(Obj, PlaneIndex)
                      
            switch (PlaneIndex)
                case 1
                    switch(Obj.APCarmViewStr)
                        case 'AP'
                            Ang=0;
                        case 'APOB'
                            Ang=Obj.APObliqueAngle;
                        case 'ML'
                            Ang=-90;
                        case 'MLOB'
                            Ang=Obj.MLObliqueAngle;
                    end
        
                case 2
                    switch(Obj.MLCarmViewStr)
                        case 'AP'
                            Ang=0;
                        case 'APOB'
                            Ang=Obj.APObliqueAngle;
                        case 'ML'
                            Ang=-90;
                        case 'MLOB'
                            Ang=Obj.MLObliqueAngle;
                    end                  
            end
        end
                
        function OblqNorm=returnOblqLatNorm(Obj)
            % Return Normal vectors to the oblique planes
            Vx=Obj.APPlane.Vector1';
            Vy=Obj.APPlane.Vector2';
            Vz=Obj.APPlane.Vector3';
            APAngle=-Obj.APObliqueAngle;
            MLAngle=-Obj.MLObliqueAngle;
            % Rotate the Z vector of the AP plane to get to the oblique views            
            Vz_AP=Vz*cosd(APAngle) + Vx*sind(APAngle);
            Vz_AP=Vz_AP/norm(Vz_AP);
            Vz_ML=Vz*cosd(MLAngle) + Vx*sind(MLAngle);
            Vz_ML=Vz_ML/norm(Vz_ML);
            OblqNorm=[Vz_AP;Vz_ML];
        end
            
        function [CorrectedAPTransform, CorrectedMLTransform]=calcNodes3D(Obj,APTransform,MLTransform)            
            
            
            if isempty(APTransform)
                APTransform=[];
                APTransform.position=[0,0,0];
                APTransform.rotation=[0,0,0];
            end
            if isempty(MLTransform)
                MLTransform=[];
                MLTransform.position=[0,0,0];
                MLTransform.rotation=[0,0,0];
            end
            
            Obj.APPlane.CNodes3D=calculateNodes3DFlat(Obj.APPlane,APTransform,Obj.APPlane.CNodes);            
            Obj.APPlane.DNodes3D=calculateNodes3DFlat(Obj.APPlane,APTransform,Obj.APPlane.DNodes);            
            Obj.MLPlane.CNodes3D=calculateNodes3DFlat(Obj.MLPlane,MLTransform,Obj.MLPlane.CNodes);
            Obj.MLPlane.DNodes3D=calculateNodes3DFlat(Obj.MLPlane,MLTransform,Obj.MLPlane.DNodes);    

            Obj.APPlane.AnatomySegments3D=calculateNodes3DFlat(Obj.APPlane,APTransform,Obj.APPlane.AnatomySegments);
            Obj.MLPlane.AnatomySegments3D=calculateNodes3DFlat(Obj.MLPlane,MLTransform,Obj.MLPlane.AnatomySegments);            
            
            %% Calculate Anatomy Segments:
            % AP Plane:
            if ~isempty(Obj.APPlane.AnatomySegments)
               Obj.MLPlane.DepthSegments=Obj.APPlane.AnatomySegments3D;
               if strcmpi(Obj.MLCarmViewStr,'APOB') % depth segment normal should account for the angle of the images
                   OblqNorm=Obj.returnOblqLatNorm;
                   Obj.MLPlane.DepthSegmentsNormal=OblqNorm(1,:)';
               else
                   Obj.MLPlane.DepthSegmentsNormal=Obj.APPlane.ImageNormal;
               end
            else
               Obj.MLPlane.DepthSegments=[];
               Obj.MLPlane.DepthSegmentsNormal=[];
            end
            
            % ML Plane:
            if ~isempty(Obj.MLPlane.AnatomySegments)              
               Obj.APPlane.DepthSegments=Obj.MLPlane.AnatomySegments3D;
               if strcmpi(Obj.APCarmViewStr,'MLOB')  % depth segment normal should account for the angle of the images
                   OblqNorm=Obj.returnOblqLatNorm;
                   Obj.APPlane.DepthSegmentsNormal=OblqNorm(2,:)';
               else               
                   Obj.APPlane.DepthSegmentsNormal=Obj.MLPlane.ImageNormal;
               end
            else
               Obj.APPlane.DepthSegments=[];
               Obj.APPlane.DepthSegmentsNormal=[];
            end                                       
            
            % Now calculate the corresponding Registered points onto the  anatomical planes
            Obj.APPlane.CNodes3DReg=calculateNodes3DReg(Obj.APPlane,Obj.APPlane.CNodes3D,Obj.APPlane.CNodesSource3D);
            Obj.APPlane.DNodes3DReg=calculateNodes3DReg(Obj.APPlane,Obj.APPlane.DNodes3D,Obj.APPlane.DNodesSource3D);
            Obj.MLPlane.CNodes3DReg=calculateNodes3DReg(Obj.MLPlane,Obj.MLPlane.CNodes3D,Obj.MLPlane.CNodesSource3D);
            Obj.MLPlane.DNodes3DReg=calculateNodes3DReg(Obj.MLPlane,Obj.MLPlane.DNodes3D,Obj.MLPlane.DNodesSource3D);
            
            
            %% Correct the AP and ML Transformation
            if (isempty(Obj.APPlane.CNodes3D)) && (isempty(Obj.MLPlane.CNodes3D)) 
                 % Cross Calibration Tool in Use
                 if sum(isnan(Obj.APPlane.DNodes3DReg(3,:)))
                     PrjN3Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.DNodes3D(3,:));
                 else                     
                     PrjN3Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.DNodes3DReg(3,:));
                 end
                 if sum(isnan(Obj.APPlane.DNodes3DReg(2,:)))
                     PrjN2Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.DNodes3D(2,:));
                 else
                     PrjN2Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.DNodes3DReg(2,:));
                 end
                 AP_RefAngl_1=acosd(dot(vectnorm(Obj.APPlane.DNodes3D(3,:)-Obj.APPlane.DNodes3D(2,:)),vectnorm(Obj.APPlane.Vector2')));
                 AP_RefAngl_2=acosd(dot(vectnorm(PrjN3Reg-PrjN2Reg),vectnorm(Obj.APPlane.Vector2')));
                 Vpos=PrjN3Reg-Obj.APPlane.DNodes3D(3,:);
                 Vpos_Rel=[dot(Vpos,Obj.APPlane.Vector1'),dot(Vpos,Obj.APPlane.Vector2'),0]/21;
                 rotCorr=AP_RefAngl_1-AP_RefAngl_2;
                 posCorr=Vpos_Rel;
                 CorrectedAPTransform=APTransform;
                 CorrectedAPTransform.position=CorrectedAPTransform.position-posCorr;
                 CorrectedAPTransform.rotation(3)=CorrectedAPTransform.rotation(3)-rotCorr;
                 
                 if (sum(isnan(Obj.MLPlane.DNodes3DReg(3,:)))==0)
                     MLPrjN3Reg=linePlaneIntersect(Obj.MLPlane.Vector3',Obj.MLPlane.Origin',Obj.MLPlane.Vector3',Obj.MLPlane.DNodes3DReg(3,:));
                 else
                     MLPrjN3Reg=linePlaneIntersect(Obj.MLPlane.Vector3',Obj.MLPlane.Origin',Obj.MLPlane.Vector3',Obj.MLPlane.DNodes3D(3,:));
                 end
                 Vpos=MLPrjN3Reg-Obj.MLPlane.DNodes3D(3,:);
                 Vpos_Rel=[dot(Vpos,Obj.MLPlane.Vector1'),dot(Vpos,Obj.MLPlane.Vector2'),0]/21;
                 MLposCorr=Vpos_Rel;
                 CorrectedMLTransform=MLTransform;
                 CorrectedMLTransform.position=CorrectedMLTransform.position-MLposCorr;
            else
                % Spine Calibration in Use:
                if sum(isnan(Obj.APPlane.CNodes3DReg(1,:)))
                    PrjN1Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.CNodes3D(1,:));
                else
                    PrjN1Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.CNodes3DReg(1,:));
                end
                if sum(isnan(Obj.APPlane.CNodes3DReg(4,:)))
                    PrjN4Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.CNodes3D(4,:));
                else
                    PrjN4Reg=linePlaneIntersect(Obj.APPlane.Vector3',Obj.APPlane.Origin',Obj.APPlane.Vector3',Obj.APPlane.CNodes3DReg(4,:));
                end
                AP_RefAngl_1=acosd(dot(vectnorm(Obj.APPlane.CNodes3D(4,:)-Obj.APPlane.CNodes3D(1,:)),vectnorm(Obj.APPlane.Vector1')));                                
                AP_RefAngl_2=acosd(dot(vectnorm(PrjN4Reg-PrjN1Reg),vectnorm(Obj.APPlane.Vector1')));
                rotCorr=AP_RefAngl_1-AP_RefAngl_2;
                
                % Not calculate the correct origin location:
                Pr4=Obj.APPlane.returnLocal2D(PrjN4Reg);
                Pr1=Obj.APPlane.returnLocal2D(PrjN1Reg);
                Pr5=Obj.APPlane.returnLocal2D(Obj.APPlane.DNodes3D(5,:));
                Pr6=Obj.APPlane.returnLocal2D(Obj.APPlane.DNodes3D(6,:));

                Vvert=vectnorm(Pr6-Pr5);
                NVvert=rotateVector(Vvert, rotCorr);
                
                NOrigin=linlinintersect([Pr4; Pr1; Pr5; Pr5+ NVvert*10]);                                
%                 posCorr=([NOrigin,0]-[Pr6,0]+[1,1,0])/21;
                posCorr=([NOrigin,0]-[Pr6,0]+[0,0,0])/21;

                CorrectedAPTransform=APTransform;
%                 CorrectedAPTransform.position=CorrectedAPTransform.position; % + [+posCorr(1),-posCorr(2),0];
%                 CorrectedAPTransform.rotation(3)=CorrectedAPTransform.rotation(3)+rotCorr;
                
                                
                % No change on the lateral view
                CorrectedMLTransform=MLTransform;
%                 CorrectedMLTransform.position=CorrectedMLTransform.position +[0,+posCorr(2),0];
            end
            %%
            

        end
        
        function OutputXY=CalCoordConv(Obj,TypeStr, InputXY, PosTrans, RotTran)
            OutputXY=[];
            switch upper(TypeStr)
                case 'M2U.AP'
                    OutputXY=InputXY*0;
                    Org=size(Obj.APPlane.Image)/2;
                    for f=1:size(InputXY,1)
                        OutputXY(f,1:2)=[(Org(2)-InputXY(f,1)), (Org(1)-InputXY(f,2))]/100;
                    end
                case 'M2U.ML'
                    OutputXY=InputXY*0;                    
                    Org=size(Obj.MLPlane.Image)/2;
                    for f=1:size(InputXY,1)
                        OutputXY(f,1:2)=[-(Org(1)-InputXY(f,2)), Org(2)-InputXY(f,1)]/100;
                    end
                case 'U2M.AP'
                    OutputXY=InputXY*0;  
                    Org=size(Obj.APPlane.Image)/2;                    
                    for f=1:size(InputXY,1)
                        OutputXY(f,1:2)=[Org(2)-InputXY(f,1)*100,Org(1)-InputXY(f,2)*100];
                    end
                case 'U2M.ML'
                    OutputXY=InputXY*0;                  
                    Org=size(Obj.MLPlane.Image)/2;                    
                    for f=1:size(InputXY,1)                    
                        OutputXY(f,1:2)=[Org(2)-InputXY(f,2)*100, Org(1)+InputXY(f,1)*100];
                    end
                case 'M2U.AP.REL'                
                    PO=CalCoordConv('M2U.AP',InputXY);
                    Org=PosTrans(1:2);
                    Vx=[cosd(RotTran(3)),sind(RotTran(3))];
                    Vy=[-sind(RotTran(3)),cosd(RotTran(3))];
                    OutputXY=InputXY*0;
                    for f=1:size(PO,1)
                        V=PO(f,:)-Org;
                        OutputXY(f,1:2)=[dot(V,Vx), dot(V,Vy)];
                    end
                case 'M2U.ML.REL'
                    PO=CalCoordConv('M2U.ML',InputXY);
                    Org=PosTrans(1:2);
                    Vx=[cosd(RotTran(3)),sind(RotTran(3))];
                    Vy=[-sind(RotTran(3)),cosd(RotTran(3))];
                    OutputXY=InputXY*0;
                    for f=1:size(PO,1)
                        V=PO-Org;
                        OutputXY(f,1:2)=[dot(V,Vx), dot(V,Vy)];
                    end
            end
        end
        
        function Out=updateNodes3D(Obj)
            Out=[];
            %Recalculate the Node3D position (on the plane) based on the 3D registered point           
            Obj.APPlane.CNodes3D=calcUpdateNodes3D(Obj.APPlane,Obj.APPlane.CNodes3DReg,Obj.APPlane.CNodes3D);
            Obj.APPlane.DNodes3D=calcUpdateNodes3D(Obj.APPlane,Obj.APPlane.DNodes3DReg,Obj.APPlane.DNodes3D);
            Obj.MLPlane.CNodes3D=calcUpdateNodes3D(Obj.MLPlane,Obj.MLPlane.CNodes3DReg,Obj.MLPlane.CNodes3D);
            Obj.MLPlane.DNodes3D=calcUpdateNodes3D(Obj.MLPlane,Obj.MLPlane.DNodes3DReg,Obj.MLPlane.DNodes3D);
            Out=1;
        end
        
        function [AP_CNodes,AP_DNodes,ML_CNodes,ML_DNodes]=calcNodes2D(Obj)
                        
            % Calculate the 2D positions of the bodes in the unity coordinate system:            
            AP_CNodes=calculateNodes2D(Obj.APPlane,Obj.APPlane.CNodes3D,Obj.APPlane.CNodes);
            AP_DNodes=calculateNodes2D(Obj.APPlane,Obj.APPlane.DNodes3D,Obj.APPlane.DNodes);

            ML_CNodes=calculateNodes2D(Obj.MLPlane,Obj.MLPlane.CNodes3D,Obj.MLPlane.CNodes);
            ML_DNodes=calculateNodes2D(Obj.MLPlane,Obj.MLPlane.DNodes3D,Obj.MLPlane.DNodes);
            
        end

        % Update the calibration dependant nodes:
        function [AP_CNodes, AP_DNodes, ML_CNodes, ML_DNodes]=resetCalibTools(Obj,AP_CNodes, AP_DNodes, ML_CNodes, ML_DNodes)
            AP_CN1=AP_CNodes(1,1:2);
            AP_CN2=AP_CNodes(2,1:2);
            AP_CN3=AP_CNodes(3,1:2);
            AP_CN4=AP_CNodes(4,1:2);
            AP_CN5=AP_CNodes(5,1:2);
            AP_CN6=AP_CNodes(6,1:2);
            AP_CN7=AP_CNodes(7,1:2);
            AP_DN1=AP_DNodes(1,1:2);
            AP_DN2=AP_DNodes(2,1:2);
            AP_DN3=AP_DNodes(3,1:2);
            AP_DN4=AP_DNodes(4,1:2);
            AP_DN5=AP_DNodes(5,1:2);
            AP_DN6=AP_DNodes(6,1:2);
             
            ML_CN1=ML_CNodes(1,1:2);
            ML_CN2=ML_CNodes(2,1:2);
            ML_CN3=ML_CNodes(3,1:2);
            ML_CN4=ML_CNodes(4,1:2);
            ML_CN5=ML_CNodes(5,1:2);            
            ML_DN1=ML_DNodes(1,1:2); 
            ML_DN2=ML_DNodes(2,1:2);   
            ML_DN3=ML_DNodes(3,1:2); 
            ML_DN4=ML_DNodes(4,1:2);               
            ML_DN5=ML_DNodes(5,1:2); 
            ML_DN6=ML_DNodes(6,1:2);   
            
            % Update the AP Nodes
            AP_CN3=(AP_CN1 + AP_CN4)/2;
            V56=[0,-1];
            AP_DN6=linlinintersect([AP_CN1;AP_CN4;AP_DN5;AP_DN5+V56*100]);
            AP_CN6=(AP_DN5+AP_DN6)/2;
            AP_CN7=AP_DN6+(AP_CN6-AP_DN5);
            % Update the ML Nodes
            ML_DN6=(ML_CN1+ML_CN3)/2;
            ML_CN5=(ML_DN5+ML_DN6)/2;  
            
            % Assemble the AP_CNodes:
            AP_CNodes=[AP_CN1,0;AP_CN2,0;AP_CN3,0;AP_CN4,0;AP_CN5,0;AP_CN6,0;AP_CN7,0];
            AP_DNodes=[AP_DN1,0;AP_DN2,0;AP_DN3,0;AP_DN4,0;AP_DN5,0;AP_DN6,0];
            % Assemble the AP_CNodes:
            ML_CNodes=[ML_CN1,0;ML_CN2,0;ML_CN3,0;ML_CN4,0;ML_CN5,0];
            ML_DNodes=[ML_DN1,0;ML_DN2,0;ML_DN3,0;ML_DN4,0;ML_DN5,0;ML_DN6,0];
        end        
        
        function [IsCrossSync,LocalizedOrigin3D,OriginIndices]=checkCrossSyncCalib(Obj,Data)                        
            % Detect if the calibration is cross and the origins are in 'sync'
            IsCrossSync=0;
            LocalizedOrigin3D=[];            
            OriginIndices=[];
            OriginDetectErrorThreshold=0.05; % Detecting the point is on the origin
            SyncErrorThresholdDist=1;
            APFoundIndex=0;
            MLFoundIndex=0;
            
            if (strcmpi(Data.APCalToolName,'cross') && strcmpi(Data.MLCalToolName,'cross')) && (~isempty(Obj.APPlane.DNodesSource3D)) && (~isempty(Obj.MLPlane.DNodesSource3D))
                % determine which index corresponds with the origin of the tool on each plane
                if (~isempty(Obj.APPlane.DNodesSource3D)) && (~isempty(Obj.MLPlane.DNodesSource3D))
                    for f=1:size(Obj.APPlane.DNodes,1)
                        if norm(Obj.APPlane.DNodes(f,:))<OriginDetectErrorThreshold
                            if ~(Obj.APPlane.DNodesSource3D(f,:)==[-1 -1 -1]) % Check for signature of the lack of source (point not on pixel)
                                APFoundIndex=f;
                            end
                        end
                    end
                    for f=1:size(Obj.MLPlane.DNodes,1)
                        if norm(Obj.MLPlane.DNodes(f,:))<OriginDetectErrorThreshold
                            if ~(Obj.APPlane.DNodesSource3D(f,:)==[-1 -1 -1]) % Check for signature of the lack of source (point not on pixel)
                                MLFoundIndex=f;
                            end
                        end
                    end
                end
            end
            
            IsCrossSync=0;
            if (APFoundIndex && MLFoundIndex)
                IsCrossSync=1;
                OriginIndices=[APFoundIndex,MLFoundIndex];
                % Calculate the localized 3D position of the origin:
                AP_S=Obj.APPlane.DNodesSource3D(3,:);
                AP_P=Obj.APPlane.DNodes3D(3,:);
                ML_S=Obj.MLPlane.DNodesSource3D(3,:);
                ML_P=Obj.MLPlane.DNodes3D(3,:);
                [LocalizedOrigin3D,distances] = lineIntersect3D([AP_S;ML_S],[AP_P;ML_P]);
                % If the intersection is not accurate then ignore CrossSync Calibration
                if distances > SyncErrorThresholdDist
                    IsCrossSync=0; % do not use the sync since the input data is inaccurate
                end
            end 
            %%
        end
        
        function Out=plotAllNodes(Obj,h1,h2)
            
            Out=[];
            
            figure(h1);
            hold on;
            Obj.APPlane.imshow3d(1,0,1,0);
            axis equal;
            if ~isempty(Obj.APPlane.CNodes)      
                Pnts=Obj.APPlane.CNodes3D;
                plot3(Pnts(:,1),Pnts(:,2),Pnts(:,3),'mo');
            end
            if ~isempty(Obj.APPlane.DNodes3D)
                Pnts=Obj.APPlane.DNodes3D;
                plot3(Pnts(:,1),Pnts(:,2),Pnts(:,3),'yo');
            end
            axis equal;
            view(Obj.APPlane.Vector3)
            
            figure(h2);
            hold on;
            Obj.MLPlane.imshow3d(1,0,1,0);
            axis equal;
            if ~isempty(Obj.MLPlane.CNodes)
                Pnts=Obj.MLPlane.CNodes3D;
                plot3(Pnts(:,1),Pnts(:,2),Pnts(:,3),'mo');
            end
            if ~isempty(Obj.MLPlane.DNodes3D)
                Pnts=Obj.MLPlane.DNodes3D;
                plot3(Pnts(:,1),Pnts(:,2),Pnts(:,3),'yo');
            end
            axis equal;
            view(Obj.MLPlane.Vector3)
        end
        
        function Out=showNodeDetails(Obj,Planes,indx)
            if exist('indx','var')
                SpecificIndex=1;
            else
                SpecificIndex=0;
            end

            Obj.APPlane.imshow3d(0,0,0,0);
            Obj.MLPlane.imshow3d(0,0,0,0);
            NodeRad=5;
            if ismember(1,Planes)
                if SpecificIndex
                    Pnts=Obj.APPlane.DepthSegments([indx,indx+1],:);
                else
                    Pnts=Obj.APPlane.DepthSegments;
                end
                Nrm=Obj.APPlane.DepthSegmentsNormal;
                plotSegments3D(Pnts,Nrm','ro','c','g',0.4);
                plotSphere(Obj.APPlane.DepthSegments,NodeRad,'none','c',1)               
                plotSphere(Obj.APPlane.CNodes3D,NodeRad,'none','c',1);
                plotSphere(Obj.APPlane.DNodes3D,NodeRad,'none','g',1);
                plotSphere(Obj.APPlane.CNodes3DReg,NodeRad,'none','c',0.5);
                plotSphere(Obj.APPlane.DNodes3DReg,NodeRad,'none','g',0.5);
            end
            
            if ismember(2,Planes)
                if SpecificIndex
                    Pnts=Obj.MLPlane.DepthSegments([indx,indx+1],:);
                else
                    Pnts=Obj.MLPlane.DepthSegments;
                end                
                Nrm=Obj.MLPlane.DepthSegmentsNormal;
                plotSegments3D(Pnts,Nrm','ro','y','r',0.4);
                plotSphere(Obj.MLPlane.DepthSegments,NodeRad,'none','y',1)
                plotSphere(Obj.MLPlane.CNodes3D,NodeRad,'none','y',1);
                plotSphere(Obj.MLPlane.DNodes3D,NodeRad,'none','r',1);
                plotSphere(Obj.MLPlane.CNodes3DReg,NodeRad,'none','y',0.5);
                plotSphere(Obj.MLPlane.DNodes3DReg,NodeRad,'none','r',0.5);                
            end
            
            axis equal;

        end
            
        function SetImageFuseMode(Obj,ReadImageFuseMode)
            AcceptableModes=[1,2,3,4,5,6];            
            try
                if ismember(ReadImageFuseMode,AcceptableModes);
                    Obj.ImageFuseMode=ReadImageFuseMode;
                end
            catch
                fprintf('Error Setting the FuseMode');
            end
        end
        
    end          
end



classdef GridElement < handle
    
    properties
        Name
        Type % not sure why this is here, but don't delete
        Colour % color of the feature
        LayerName % for name of layer that it is on
        PlotType % the type to plot
        SelLayer % cell that keeps track of the selected layer it is on
        InitialAnchorNumber % Number of anchor points to start with
        AnchorPointSpheres = []; % The SphereElements that act as the anchor points
        StepNumber % Number of step that contains this Element
        SelElements % an object that contains all selected elements from other steps
        CentreRadius % the radius of the target circle in the centre of the IP
        TargetRadius % the radius of the target circle that moves dynamically 
                      % depending on the location and orientation of the landmarks
        ProximitySensitivity % the higher this is, the sooner the TargetElement spheres will turn
                            % green as they approach each other
        Mode % 1 - Orientation only (default)
             % 2 - Orientation and position
        ImagesPresent % a flag to indicate whether the target and arrow images are present
        CenterScale % determines the size of the non-moving center circle
        TargetScale % determines the size of the moving, outer target circle
        targetClosePath % holds a string with the path to a png image of the target for the close case
        targetFarPath % holds path to png image for the target in the far case
        IPs % a property for holding all the IPs available in the step
        LastIPInfo % a structure for keeping track of the last colour for each axes
        arrowImgPath % the file path to the arrow png. The arrow must be facing with to point towards the right
        arrowImgData % the image pixel and alpha data in an 1x2 cell array. For easy referencing during update
        arrowScale % determines the size of the arrow png, if imported
        landmarkPath % file path to the landmark indicator png file
        landmarkScale % determines the size of the landmark png, if imported
        lastLineElementVector % required for the initializeLastIPInfo() function to tell whether the coordinates 
                                % have changed since the last run. This
                                % would be the case if the user goes back
                                % a step and changes the position of the
                                % landmarks
        % GridElement Specific Properties                        
        
        % Origin of the Grid Element
        Origin
        OriginSphereElement %Sphere Element marking the origin of the grid
        XSphereElement   %Sphere Element marking the Xdirection of the grid
        YSphereElement   %Sphere Element marking the Ydirection of the grid
        
        VxLineElement
        VyLineElement
        VzLineElement
        PlaneType
        PlaneNormal
        PlaneOrigin
        MajorTickSize
        MinorTickSize
        MajorTickColor
        MinorTickColor
        Axis1Color
        Axis2Color
        Guided
        AutoPlaneSelect
        MeasureAnchorPoint
        RefElements
        TypeName
        ExtraPlots %to hold extra anatomical references that might be created along with the grid during construction
        ExtraElement %Ploygon stickmodel of the pelvis generated during building of the GridElement
        Constants %Added Feb 16:
        StitchedViews
        APIndex
        LatIndex
        
        ScoutAPNormal % To hold the original normal vector from the Scout Views; Added 2018-May-09
        ScoutMLNormal
    end
    
    methods
        function Obj = GridElement(EName,varargin)
            Obj.Name = EName;
            Obj.ExtraPlots=[];
            Obj.Constants=[];
            Obj.StitchedViews=[];
            Obj.ScoutAPNormal=[];
            Obj.ScoutMLNormal=[];
            if ~isempty(varargin)
                Obj.TypeName=varargin{1};              
                TypeName=Obj.TypeName;           

                switch lower(TypeName)
                    case 'simple'
                        SphereElement1 = varargin{2};
                        SphereElement2= varargin{3};
                        SphereElement3 = varargin{4};
                        RefO=SphereElement1.Centre;
                        RefX=SphereElement2.Centre;
                        RefY=SphereElement3.Centre;
                        vx=RefX-RefO;
                        vy=RefY-RefO;
                        vz=cross(vx,vy);
                        vz=vz/norm(vz);
                        RefZ=RefO+50*vz;
                    case 'spinopelvic_old'
                        SphereElement1 = varargin{2};
                        SphereElement2= varargin{3};
                        SplineElement3 = varargin{4};
                        Pel_Nums = varargin{5}; 
                        % [Pel_Preop,Pel_PI,Pel_PT,Pel_S1Size]; 
                        Pel_Preop=Pel_Nums(1);
                        Pel_PI=Pel_Nums(2);
                        Pel_PT=Pel_Nums(3);
                        Pel_S1Size=Pel_Nums(4);

                        H1=SphereElement1.Centre;
                        H2=SphereElement2.Centre;
                        S1_SP=SplineElement3.AnchorPointSpheres(1).Centre;
                        S1_SO=SplineElement3.AnchorPointSpheres(2).Centre;
                        HO=(H1+H2)/2;
    %                     S1_O=(S1_SP+S1_SA)/2;


                        % Add the SplineElement 
                        SPM=SplineElement('Pelvic Stick Model');
                        SPM.Type='Line';
                        SPM.InitialAnchorNumber=8; % Number of anchor points to start with
                        LBS={'1','2','3','4','5','6','7','8'};                        
                        SS1=SphereElement('Anchor 1');
                        SS1.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS1.Radius=1e-5;
                        SS1.Labels=LBS;
                        SS2=SphereElement('Anchor 2');
                        SS2.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS2.Radius=1e-5;
                        SS2.Labels=LBS;                    
                        SS3=SphereElement('Anchor 3');
                        SS3.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS3.Radius=1e-5;
                        SS3.Labels=LBS;                    
                        SS4=SphereElement('Anchor 4');
                        SS4.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS4.Radius=1e-5;                    
                        SS4.Labels=LBS;                    
                        SS5=SphereElement('Anchor 5');
                        SS5.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS5.Radius=1e-5;
                        SS5.Labels=LBS;                    
                        SS6=SphereElement('Anchor 6');
                        SS6.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS6.Radius=1e-5;
                        SS6.Labels=LBS;                    
                        SS7=SphereElement('Anchor 7');
                        SS7.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS7.Radius=1e-5;
                        SS7.Labels=LBS;                    
                        SS8=SphereElement('Anchor 8');
                        SS8.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS8.Radius=1e-5;         
                        SS8.Labels=LBS;
                        %
    %                     SPM.AnchorPointSpheres=[SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8];
                        SPM.AnchorPointSpheres=[SS1;SS2;SS3;SS4;SS5;SS6;SS7;SS8];
                        Obj.ExtraElement=SPM;

                        %sort hips so H1 will be Right and H2 will be left
                        AP=S1_SO-S1_SP;
                        S1HO=S1_SO-HO;
                        AP=AP/norm(AP);
                        S1HO=S1HO/norm(S1HO);
                        ML=cross(S1HO,AP);
                        V=H2-HO;
                        if dot(V,ML)>0 % then the H1 is Right & H2 is Left

                        else  %swap H1 and H2
                            temp=H1;
                            H1=H2;
                            H2=temp;
                        end

                        H1H2=H2-H1;
                        H1H2=H1H2/norm(H1H2);
                        APt=cross(H1H2,S1HO);
                        MLt=cross(S1HO,APt);
                        [APc,S1HOc]=RotZ(HO,APt,S1HO,-Pel_PT);

                        Vx=MLt;
                        Vy=S1HOc;
                        Vz=APc;
                        if Pel_Preop
                            % calculating the Sacral Slope vector from angular
                            % input (PI):
                            [Vsn,Vss]=RotZ(S1_SO,Vy,Vz,(90-Pel_PI));
                            S1_SPP=S1_SO-Vss*Pel_S1Size/2;
                            S1_SAP=S1_SO+Vss*Pel_S1Size/2;                        
                        else
                            %calculating the Sacral Slope vector from input points:
                            VS=S1_SP-S1_SO;
                            VSy=dot(VS,Vy);
                            VSz=dot(VS,Vz);
                            S1_SPP=S1_SO+VSz*Vz+VSy*Vy;
                            Vss=S1_SO-S1_SPP;
                            Vss=Vss/norm(Vss); % Sacral Slope vector (corrected)
                            S1_SAP=S1_SO+Vss*norm(S1_SO-S1_SPP);
                            Vsn=cross(Vss,Vx);  %normal to Sacral slope, directed upwards
                            Pel_S1Size=2*norm(S1_SO-S1_SPP);
                        end

                        RefO=S1_SPP; % Place the origin at the posterior side of the Superior-posterior corner of the S1
                        RefX=RefO+50*Vx;
                        RefY=RefO+50*Vy;
                        RefZ=RefO+50*Vz;

                        % create the ExtraPlots data, here the anatomical
                        % features of pelvis on the sagittal view:
                        % Points 1 to 8 to form an polygon connecting the
                        % references of the spinopelvic anatomy:
                        SZ=norm(S1_SPP-S1_SAP);
                        P2=S1_SPP;
                        P1=P2+Vz*SZ;
                        P3=S1_SO;
                        P5=P3;
                        P7=P3;
                        P4=P3-Vsn*2*SZ;
                        P6=HO;
                        P8=P2+Vss*Pel_S1Size;

                        %Obj.ExtraPlots=[P1,P2,P3,P4,P5,P6,P7,P8];
                        %SPM.AnchorPointSpheres=[S1,S2,S3,S4,S5,S6,S7,S8];
                        Obj.ExtraElement.AnchorPointSpheres(1).Centre=P1;
                        Obj.ExtraElement.AnchorPointSpheres(2).Centre=P2;
                        Obj.ExtraElement.AnchorPointSpheres(3).Centre=P3;
                        Obj.ExtraElement.AnchorPointSpheres(4).Centre=P4;
                        Obj.ExtraElement.AnchorPointSpheres(5).Centre=P5;
                        Obj.ExtraElement.AnchorPointSpheres(6).Centre=P6;
                        Obj.ExtraElement.AnchorPointSpheres(7).Centre=P7;
                        Obj.ExtraElement.AnchorPointSpheres(8).Centre=P8;

                        % Add PI and S1_AP_Width to the constants of the
                        % element, so they can be displayed for the user:
                        Res_APW=SZ;
                        if Pel_Preop
                            Res_PI=Pel_PI;
                        else
                            V1=S1_SAP-S1_SPP;
                            V2=HO-S1_SO;
                            V1=V1/norm(V1);
                            V2=V2/norm(V2);
                            Res_PI=90-180/pi*acos(dot(V1,V2));
                        end
                        Obj.addConstant('S1APW',Res_APW);
                        Obj.addConstant('PI',Res_PI);
                        % =====================================================
                        % plot the spinopelvic stick-model
    %                     figure;
    %                     tt=Obj.ExtraPlots;
    %                     hold on;
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'g');
    %                     TXT={'P1','P2','P3','P4','P5','P6','P7','P8','P9'};
    %                     for f=1:8
    %                         text(tt(1,f)',tt(2,f)',tt(3,f)',TXT{f});                    
    %                     end
    %                     tt=[P1,P2,P8];
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'r');
    % 
    %                     tt=[H1,HO,H2];
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'m');
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'ro');
    %                     axis equal;
    % %                     tt=S1_SPP;
    % %                     plot3(tt(1,:),tt(2,:),tt(3,:),'b+');
    % %                     text(tt(1,:),tt(2,:),tt(3,:),'S1P');
    %                                   
    %                     tt=[RefO,RefX];
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'r');
    %                     tt=[RefO,RefY];
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'g');
    %                     tt=[RefO,RefZ];
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'b');
    %                     
    %                     tt=[RefO,RefO+Vss*100];
    %                     plot3(tt(1,:),tt(2,:),tt(3,:),'c');
    %                     axis equal;        
                        % =====================================================



                    case 'spinopelvic'
                        SphereElement1 = varargin{2};
                        SphereElement2= varargin{3};
                        SplineElement3 = varargin{4};
                        Pel_Nums = varargin{5}; 
                        % [Pel_Preop,Pel_PI,Pel_PT,Pel_S1Size]; 
                        Pel_Preop=Pel_Nums(1);
                        Pel_PI=Pel_Nums(2);
                        Pel_PT=Pel_Nums(3);
                        Pel_S1Size=Pel_Nums(4);

                        H1=SphereElement1.Centre;
                        H2=SphereElement2.Centre;
                        S1_SP=SplineElement3.AnchorPointSpheres(1).Centre;
                        S1_SO=SplineElement3.AnchorPointSpheres(2).Centre;
                        HO=(H1+H2)/2;
    %                     S1_O=(S1_SP+S1_SA)/2;

                        % Update the coordinates of the S1_SP, S1_SO:
                        % by projecting over the sagittal plane:
                        V_ML=(H2-H1)/norm(H2-H1);
                        S1_SP=projPointOnPlane(S1_SP,HO,V_ML);
                        S1_SO=projPointOnPlane(S1_SO,HO,V_ML);

                        % Add the SplineElement 
                        SPM=SplineElement('Pelvic Stick Model');
                        SPM.Type='Line';
                        SPM.InitialAnchorNumber=8; % Number of anchor points to start with
                        LBS={'1','2','3','4','5','6','7','8'};                        
                        SS1=SphereElement('Anchor 1');
                        SS1.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS1.Radius=1e-5;
                        SS1.Labels=LBS;
                        SS2=SphereElement('Anchor 2');
                        SS2.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS2.Radius=1e-5;
                        SS2.Labels=LBS;                    
                        SS3=SphereElement('Anchor 3');
                        SS3.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS3.Radius=1e-5;
                        SS3.Labels=LBS;                    
                        SS4=SphereElement('Anchor 4');
                        SS4.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS4.Radius=1e-5;                    
                        SS4.Labels=LBS;                    
                        SS5=SphereElement('Anchor 5');
                        SS5.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS5.Radius=1e-5;
                        SS5.Labels=LBS;                    
                        SS6=SphereElement('Anchor 6');
                        SS6.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS6.Radius=1e-5;
                        SS6.Labels=LBS;                    
                        SS7=SphereElement('Anchor 7');
                        SS7.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS7.Radius=1e-5;
                        SS7.Labels=LBS;                    
                        SS8=SphereElement('Anchor 8');
                        SS8.Centre=HO; %Obj.OriginSphereElement.Centre;
                        SS8.Radius=1e-5;         
                        SS8.Labels=LBS;
                        %
    %                     SPM.AnchorPointSpheres=[SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8];
                        SPM.AnchorPointSpheres=[SS1;SS2;SS3;SS4;SS5;SS6;SS7;SS8];
                        Obj.ExtraElement=SPM;

                        %sort hips so H1 will be Right and H2 will be left
                        AP=S1_SO-S1_SP;
                        S1HO=S1_SO-HO;
                        AP=AP/norm(AP);
                        S1HO=S1HO/norm(S1HO);
                        ML=cross(S1HO,AP);
                        V=H2-HO;
                        if dot(V,ML)>0 % then the H1 is Right & H2 is Left

                        else  %swap H1 and H2
                            temp=H1;
                            H1=H2;
                            H2=temp;
                        end

                        H1H2=H2-H1;
                        H1H2=H1H2/norm(H1H2);
                        APt=cross(H1H2,S1HO);
                        MLt=cross(S1HO,APt);
                        [APc,S1HOc]=RotZ(HO,APt,S1HO,-Pel_PT);

                        Vx=MLt;
                        Vy=S1HOc;
                        Vz=APc;
                        if Pel_Preop
                            % calculating the Sacral Slope vector from angular
                            % input (PI):
                            [Vsn,Vss]=RotZ(S1_SO,Vy,Vz,(90-Pel_PI));
                            S1_SPP=S1_SO-Vss*Pel_S1Size/2;
                            S1_SAP=S1_SO+Vss*Pel_S1Size/2;                        
                        else
                            %calculating the Sacral Slope vector from input points:
                            VS=S1_SP-S1_SO;
                            VSy=dot(VS,Vy);
                            VSz=dot(VS,Vz);
                            S1_SPP=S1_SO+VSz*Vz+VSy*Vy;
                            Vss=S1_SO-S1_SPP;
                            Vss=Vss/norm(Vss); % Sacral Slope vector (corrected)
                            S1_SAP=S1_SO+Vss*norm(S1_SO-S1_SPP);
                            Vsn=cross(Vss,Vx);  %normal to Sacral slope, directed upwards
                            Pel_S1Size=2*norm(S1_SO-S1_SPP);
                        end

                        RefO=S1_SPP; % Place the origin at the posterior side of the Superior-posterior corner of the S1
                        RefX=RefO+50*Vx;
                        RefY=RefO+50*Vy;
                        RefZ=RefO+50*Vz;

                        % create the ExtraPlots data, here the anatomical
                        % features of pelvis on the sagittal view:
                        % Points 1 to 8 to form an polygon connecting the
                        % references of the spinopelvic anatomy:
                        SZ=norm(S1_SPP-S1_SAP);
                        P2=S1_SPP;
                        P1=P2+Vz*SZ;
                        P3=S1_SO;
                        P5=P3;
                        P7=P3;
                        P4=P3-Vsn*2*SZ;
                        P6=HO;
                        P8=P2+Vss*Pel_S1Size;

                        %Obj.ExtraPlots=[P1,P2,P3,P4,P5,P6,P7,P8];
                        %SPM.AnchorPointSpheres=[S1,S2,S3,S4,S5,S6,S7,S8];
                        Obj.ExtraElement.AnchorPointSpheres(1).Centre=P1;
                        Obj.ExtraElement.AnchorPointSpheres(2).Centre=P2;
                        Obj.ExtraElement.AnchorPointSpheres(3).Centre=P3;
                        Obj.ExtraElement.AnchorPointSpheres(4).Centre=P4;
                        Obj.ExtraElement.AnchorPointSpheres(5).Centre=P5;
                        Obj.ExtraElement.AnchorPointSpheres(6).Centre=P6;
                        Obj.ExtraElement.AnchorPointSpheres(7).Centre=P7;
                        Obj.ExtraElement.AnchorPointSpheres(8).Centre=P8;

                        % Add PI and S1_AP_Width to the constants of the
                        % element, so they can be displayed for the user:
                        Res_APW=SZ;
                        if Pel_Preop
                            Res_PI=Pel_PI;
                        else
                            V1=S1_SAP-S1_SPP;
                            V2=HO-S1_SO;
                            V1=V1/norm(V1);
                            V2=V2/norm(V2);
                            Res_PI=90-180/pi*acos(dot(V1,V2));
                        end
                        Obj.addConstant('S1APW',Res_APW);
                        Obj.addConstant('PI',Res_PI);




                    case 'lowerlimb'
                        SplineElement3 = varargin{2};
                        SphereElement1= varargin{3};

                        AN=SphereElement1.Centre;
                        HC=SplineElement3.AnchorPointSpheres(1).Centre;
                        HN=SplineElement3.AnchorPointSpheres(2).Centre;

                        % Add the SplineElement 
                        SPM=SplineElement('Pelvic Stick Model');
                        SPM.Type='Line';
                        SPM.InitialAnchorNumber=3; % Number of anchor points to start with
                        LBS={'1','2','3','4','5','6','7','8'};                        
                        SS1=SphereElement('Anchor 1');
                        SS1.Centre=HN; %Obj.OriginSphereElement.Centre;
                        SS1.Radius=1e-5;
                        SS1.Labels=LBS;
                        SS2=SphereElement('Anchor 2');
                        SS2.Centre=HC; %Obj.OriginSphereElement.Centre;
                        SS2.Radius=1e-5;
                        SS2.Labels=LBS;                    
                        SS3=SphereElement('Anchor 3');
                        SS3.Centre=AN; %Obj.OriginSphereElement.Centre;
                        SS3.Radius=1e-5;
                        SS3.Labels=LBS;                    
                        SPM.AnchorPointSpheres=[SS1;SS2;SS3];
                        Obj.ExtraElement=SPM;

                        RO=HC;
                        Vy=HC-AN;
                        Vy=Vy/norm(Vy);
                        Vt=HN-HC;
                        Vt=Vt/norm(Vt);
                        Vz=cross(Vt,Vy);
                        Vx=cross(Vz,Vy);

                        RefO=RO; % Place the origin at the posterior side of the Superior-posterior corner of the S1
                        RefX=RefO+50*Vx;
                        RefY=RefO+50*Vy;
                        RefZ=RefO+50*Vz;

                        % create the ExtraPlots data, here the anatomical
                        % features of pelvis on the sagittal view:
                        % Points 1 to 8 to form an polygon connecting the
                        % references of the spinopelvic anatomy:

                        P1=HN;
                        P2=HC;
                        P3=AN;

                        Obj.ExtraElement.AnchorPointSpheres(1).Centre=P1;
                        Obj.ExtraElement.AnchorPointSpheres(2).Centre=P2;
                        Obj.ExtraElement.AnchorPointSpheres(3).Centre=P3;

                        % Add LengLength to the constants of the
                        % element, so they can be displayed for the user:
                        LEG_LEN=norm(HC-AN);
                        Obj.addConstant('LegLen',LEG_LEN);

                    case 'femaxialrot'
                        SplineElement1 = varargin{2};
                        SplineElement2= varargin{3};

                        HC=SplineElement1.AnchorPointSpheres(1).Centre;
                        HCR=SplineElement1.AnchorPointSpheres(1).Radius;
                        HN=SplineElement1.AnchorPointSpheres(2).Centre;

                        K1=SplineElement2.AnchorPointSpheres(1).Centre;
                        K1R=SplineElement2.AnchorPointSpheres(1).Radius;
                        K2=SplineElement2.AnchorPointSpheres(2).Centre;
                        K2R=SplineElement2.AnchorPointSpheres(2).Radius;
                        VH=HC-HN;
                        VH=VH/norm(VH);
                        VK=K1-K2;
                        VK=VK/norm(VK);

                        if dot(VK,VH)>0 
                            Med=K1;
                            Lat=K2;
                            MR=K1R;
                            LR=K2R;
                        else
                            Med=K2;
                            Lat=K1;
                            MR=K2R;
                            LR=K1R;                        
                        end

                        VK_CR=Med-Lat; % corrected Knee axis
                        VK_CR=VK_CR/norm(VK_CR);

                        % Add the SplineElement 
                        SPM=SplineElement('Pelvic Stick Model');
                        SPM.Type='Line';
                        SPM.InitialAnchorNumber=3; % Number of anchor points to start with
                        LBS={'1','2','3','4','5','6','7','8'};                        
                        SS1=SphereElement('Anchor 1');
                        SS1.Centre=HC; %Obj.OriginSphereElement.Centre;
                        SS1.Radius=HCR;
                        SS1.Labels=LBS;
                        SS2=SphereElement('Anchor 2');
                        SS2.Centre=HN; %Obj.OriginSphereElement.Centre;
                        SS2.Radius=1e-5;
                        SS2.Labels=LBS;                    
                        SS3=SphereElement('Anchor 3');
                        SS3.Centre=(Med+Lat)/2; %Obj.OriginSphereElement.Centre;
                        SS3.Radius=1e-5;
                        SS3.Labels=LBS;  
                        SS4=SphereElement('Anchor 4');
                        SS4.Centre=Med; %Obj.OriginSphereElement.Centre;
                        SS4.Radius=MR;
                        SS4.Labels=LBS;                
                        SS5=SphereElement('Anchor 5');
                        SS5.Centre=Lat; %Obj.OriginSphereElement.Centre;
                        SS5.Radius=LR;
                        SS5.Labels=LBS;                     

                        SPM.AnchorPointSpheres=[SS1;SS2;SS3;SS4;SS5];
                        Obj.ExtraElement=SPM;

                        RO=HC;
                        Vy=HC-((Med+Lat)/2);
                        Vy=Vy/norm(Vy);
                        Vt=HN-HC;
                        Vt=Vt/norm(Vt);
                        Vz=cross(Vt,Vy);
                        Vx=cross(Vz,Vy);

                        RefO=RO; % Place the origin at the posterior side of the Superior-posterior corner of the S1
                        RefX=RefO+50*Vx;
                        RefY=RefO+50*Vy;
                        RefZ=RefO+50*Vz;

                        % create the ExtraPlots data, here the anatomical
                        % features of pelvis on the sagittal view:
                        % Points 1 to 8 to form an polygon connecting the
                        % references of the spinopelvic anatomy:

                        P1=HN;
                        P2=HC;
                        P3=(Med+Lat)/2;

                        Obj.ExtraElement.AnchorPointSpheres(1).Centre=P1;
                        Obj.ExtraElement.AnchorPointSpheres(2).Centre=P2;
                        Obj.ExtraElement.AnchorPointSpheres(3).Centre=P3;

                        % Add FEM)ROT to the constants of the
                        % element, so they can be displayed for the user:

                        % Knee Coordinates system
                        RO=((Med+Lat)/2);
                        Vx=Lat-RO;
                        Vx=Vx/norm(Vx);                    
                        Vt=HN-RO;
                        Vt=Vt/norm(Vt);
                        Vz=cross(Vx,Vt);
                        Vy=cross(Vz,Vx);                    
                        Vz=Vz/norm(Vz);
                        Vy=Vy/norm(Vy);

                        %projecte axis of the hip on the XZ plane:
                        VH_PR=dot(VH,Vx)*Vx+ dot(VH,Vz)*Vz;
                        VK_PR=dot(VK_CR,Vx)*Vx+ dot(VK_CR,Vz)*Vz;
                        VH_PR=VH_PR/norm(VH_PR);
                        VK_PR=VK_PR/norm(VK_PR);

                        FEM_ROT=dot(VH_PR,VK_PR)*180/pi;
                        Obj.addConstant('FEM_ROT',FEM_ROT);



                end            


                Obj.PlotType = 'GridElement';
                Obj.Colour = 'red'; % default color            
                Obj.PlaneType='XY';

                Obj.MajorTickSize=50;
                Obj.MinorTickSize=25;
                Obj.MajorTickColor='g';
                Obj.MinorTickColor='c';
                Obj.Guided=0;
                Obj.AutoPlaneSelect=0;

                %------------------------------------------------------------
                % Elements to be used for displaying the grid:
    %             OriginSphereElement % Origin of the Grid Element
    %             VxLineElement
    %             VyLineElement
    %             VzLineElement
    %             MeasureAnchorPoint
                %------------------------------------------------------------

                Obj.LayerName = [];
                Obj.SelLayer = {}; 
                Obj.InitialAnchorNumber = 0;

                %Obj.OriginSphereElement=SphereElement('Orig',OrgPoint,OrgRad); % Origin of the Grid Element
                Obj.MeasureAnchorPoint=Obj.OriginSphereElement;
                Obj.VxLineElement=LineElement('Vx');
                Obj.VyLineElement=LineElement('Vy');
                Obj.VzLineElement=LineElement('Vz');

                S1=SphereElement('Anchor 1');
                S1.Centre=RefO; %Obj.OriginSphereElement.Centre;
                S1.Radius=1e-5;

                S2=SphereElement('Anchor 2');
                S2.Centre=RefX; %Obj.XSphereElement.Centre;
                S2.Radius=1e-5;
                Obj.VxLineElement.AnchorPointSpheres=[S1,S2];
                V=S2.Centre-S1.Centre;
                V=V/norm(V);
                Obj.VxLineElement.Coords=S1.Centre;
                Obj.VxLineElement.Vector=V;

                S2=SphereElement('Anchor 2');
                S2.Centre=RefY; %Obj.YSphereElement.Centre;
                S2.Radius=1e-5;
                Obj.VyLineElement.AnchorPointSpheres=[S1,S2];
                V=S2.Centre-S1.Centre;
                V=V/norm(V);
                Obj.VyLineElement.Coords=S1.Centre;
                Obj.VyLineElement.Vector=V;

                %calculate the Vz orthogonal to Vx and Vy:

                S2=SphereElement('Anchor 3');
                S2.Centre=RefZ;
                S2.Radius=1e-5;
                Obj.VzLineElement.AnchorPointSpheres=[S1,S2];
                V=S2.Centre-S1.Centre;
                V=V/norm(V);
                Obj.VzLineElement.Coords=S1.Centre;
                Obj.VzLineElement.Vector=V;

                Obj.PlaneType='XY';            

                switch Obj.PlaneType
                    case 'XY'                   
                        Obj.PlaneNormal=Obj.VzLineElement.Vector;
                        Obj.PlaneOrigin=RefO;
                        Obj.Axis1Color='r';
                        Obj.Axis2Color='g';
                    case 'YZ'                   
                        Obj.PlaneNormal=Obj.VxLineElement.Vector;
                        Obj.PlaneOrigin=RefO;
                        Obj.Axis1Color='g';
                        Obj.Axis2Color='b';                    
                    case 'XZ'                   
                        Obj.PlaneNormal=Obj.VyLineElement.Vector;
                        Obj.PlaneOrigin=RefO;
                        Obj.Axis1Color='b';
                        Obj.Axis2Color='r';
                end
                
            end
            
        end
        
        
        
        function Obj = recreateElement(Obj,varargin)
            %Is expected to be run only for 'Simple' grid elements
%             Obj.Name = EName;
            Obj.ExtraPlots=[];
            Obj.Constants=[];
            Obj.StitchedViews=[];
            if ~isempty(varargin)
                Obj.TypeName=varargin{1};              
                TypeName=Obj.TypeName;           

                switch lower(TypeName)
                    case 'simple'
                        SphereElement1 = varargin{2};
                        SphereElement2= varargin{3};
                        SphereElement3 = varargin{4};
                        RefO=SphereElement1.Centre;
                        RefX=SphereElement2.Centre;
                        RefY=SphereElement3.Centre;
                        vx=RefX-RefO;
                        vy=RefY-RefO;
                        vz=cross(vx,vy);
                        vz=vz/norm(vz);
                        RefZ=RefO+50*vz;
                end            


                Obj.PlotType = 'GridElement';
                Obj.Colour = 'red'; % default color            
                Obj.PlaneType='XY';

                Obj.MajorTickSize=50;
                Obj.MinorTickSize=25;
                Obj.MajorTickColor='g';
                Obj.MinorTickColor='c';
                Obj.Guided=0;
                Obj.AutoPlaneSelect=0;

                Obj.LayerName = [];
                Obj.SelLayer = {}; 
                Obj.InitialAnchorNumber = 0;

                S1=SphereElement('Anchor 1');
                S1.Centre=RefO; %Obj.OriginSphereElement.Centre;
                S1.Radius=1e-5;

                S2=SphereElement('Anchor 2');
                S2.Centre=RefX; %Obj.XSphereElement.Centre;
                S2.Radius=1e-5;
                
                Obj.VxLineElement.AnchorPointSpheres(1).Centre=S1.Centre;
                Obj.VxLineElement.AnchorPointSpheres(2).Centre=S2.Centre;
                
                V=S2.Centre-S1.Centre;
                V=V/norm(V);
                Obj.VxLineElement.Coords=S1.Centre;
                Obj.VxLineElement.Vector=V;

                S2=SphereElement('Anchor 2');
                S2.Centre=RefY; %Obj.YSphereElement.Centre;
                S2.Radius=1e-5;
                % Obj.VyLineElement.AnchorPointSpheres=[S1,S2];
                Obj.VyLineElement.AnchorPointSpheres(1).Centre=S1.Centre;
                Obj.VyLineElement.AnchorPointSpheres(2).Centre=S2.Centre;
                
                
                V=S2.Centre-S1.Centre;
                V=V/norm(V);
                Obj.VyLineElement.Coords=S1.Centre;
                Obj.VyLineElement.Vector=V;

                %calculate the Vz orthogonal to Vx and Vy:
                S3=SphereElement('Anchor 3');
                S3.Centre=RefZ;
                S3.Radius=1e-5;
                % Obj.VzLineElement.AnchorPointSpheres=[S1,S2];
                Obj.VzLineElement.AnchorPointSpheres(1).Centre=S1.Centre;
                Obj.VzLineElement.AnchorPointSpheres(2).Centre=S3.Centre;
                
                V=S3.Centre-S1.Centre;
                V=V/norm(V);
                Obj.VzLineElement.Coords=S1.Centre;
                Obj.VzLineElement.Vector=V;

                Obj.PlaneType='XY';            

                switch Obj.PlaneType
                    case 'XY'                   
                        Obj.PlaneNormal=Obj.VzLineElement.Vector;
                        Obj.PlaneOrigin=RefO;
                        Obj.Axis1Color='r';
                        Obj.Axis2Color='g';
                    case 'YZ'                   
                        Obj.PlaneNormal=Obj.VxLineElement.Vector;
                        Obj.PlaneOrigin=RefO;
                        Obj.Axis1Color='g';
                        Obj.Axis2Color='b';                    
                    case 'XZ'                   
                        Obj.PlaneNormal=Obj.VyLineElement.Vector;
                        Obj.PlaneOrigin=RefO;
                        Obj.Axis1Color='b';
                        Obj.Axis2Color='r';
                end
                
            end
            
        end        
        
        
        
        function reconElement(Obj,IPs,featType)
            % Do nothing. This method is here just in case it is called
        end
        
        
        
        function updatePlaneType(Obj)
            switch Obj.PlaneType
                case 'XY'                   
                    Obj.PlaneNormal=Obj.VzLineElement.Vector;
%                     Obj.PlaneOrigin=RefO;
                    Obj.Axis1Color='r';
                    Obj.Axis2Color='g';
                case 'YZ'                   
                    Obj.PlaneNormal=Obj.VxLineElement.Vector;
%                     Obj.PlaneOrigin=RefO;
                    Obj.Axis1Color='g';
                    Obj.Axis2Color='b';                    
                case 'XZ'                   
                    Obj.PlaneNormal=Obj.VyLineElement.Vector;
%                     Obj.PlaneOrigin=RefO;
                    Obj.Axis1Color='b';
                    Obj.Axis2Color='r';
            end
        end
        
        function toggle(Obj,ToggleNum)            
            if ToggleNum % >0 then a specific toggle is selected
                switch ToggleNum
                    case 1
                        Obj.PlaneType='XY';
                    case 2
                        Obj.PlaneType='YZ';
                    case 3
                        Obj.PlaneType='XZ';
                end               
            else       % <0 cycle through possible PlaneTypes 
                switch Obj.PlaneType
                    case 'XY'
                        Obj.PlaneType='YZ';
                    case 'YZ'
                        Obj.PlaneType='XZ';
                    case 'XZ'
                        Obj.PlaneType='XY';
                end
            end   

           

            switch Obj.PlaneType
                case 'XY'                   
                    Obj.PlaneNormal=Obj.VzLineElement.Vector;
%                     Obj.PlaneOrigin=Obj.Origin; %SphereElement.Centre;
                    Obj.Axis1Color='r';
                    Obj.Axis2Color='g';
                case 'YZ'                   
                    Obj.PlaneNormal=Obj.VxLineElement.Vector;
%                     Obj.PlaneOrigin=Obj.Origin; %SphereElement.Centre;
                    Obj.Axis1Color='g';
                    Obj.Axis2Color='b';                    
                case 'XZ'                   
                    Obj.PlaneNormal=Obj.VyLineElement.Vector;
%                     Obj.PlaneOrigin=Obj.Origin; %SphereElement.Centre;
                    Obj.Axis1Color='r';
                    Obj.Axis2Color='b';                    
            end
            
            
        end
        
        function [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=calcGrid_old2(Obj,IP)
           Plot3D=1;
           Plot2D=1;
%            LX=[];
%            LY=[];
%            GX=[];
%            GY=[];           
%            CLR=[];
%            ALX=[];
%            ALY=[];
%            AGX=[];
%            AGY=[];
           
           GXX=[];
           GXY=[];
           GXA=[];
           GYX=[];
           GYY=[];
           GYA=[];
           
           LXX=[];
           LXY=[];
           LXA=[];
           LYX=[];
           LYY=[];
           LYA=[];


           if Plot3D
               h1=figure;
               IP.imshow3d(1,1);          
               hold on;
           end
           if Plot2D
               h2=figure;
               figure(h2);
               imshow(IP.Image);
               hold on;
           end
           
           MinX=-2000;
           MaxX=2000;

           
           %--------------------------------------------------------------
           % Calculate the intersection between the centre of image and the
           % XY plane:
           Pv=Obj.PlaneNormal;
           Pp=Obj.PlaneOrigin;
           Lv=IP.ImageNormal;
           P2=IP.SourcePosition;
           Pnt=linePlaneIntersect(Pv',Pp',Lv',P2');
           %--------------------------------------------------------------

%            XP=floor(Pnt(1)/Obj.MajorTickSize)*Obj.MajorTickSize;
%            YP=floor(Pnt(2)/Obj.MajorTickSize)*Obj.MajorTickSize;
%            ZP=floor(Pnt(3)/Obj.MajorTickSize)*Obj.MajorTickSize;
%            Org=[XP,YP,ZP]; %round point
           Org=Pp';
           
           
           Vx=Obj.VxLineElement.Vector';
           Vy=Obj.VyLineElement.Vector';
           Vz=Obj.VzLineElement.Vector';
           
%            switch Obj.PlaneType
%                case 'XY'                   
%                    Vx=Obj.VxLineElement.Vector';
%                    Vy=Obj.VyLineElement.Vector';
%                    Vz=Obj.VzLineElement.Vector';
%                case 'YZ'
%                    Vx=Obj.VxLineElement.Vector';
%                    Vy=Obj.VyLineElement.Vector';
%                    Vz=Obj.VzLineElement.Vector';                   
%                case 'ZX'
%                    Vx=Obj.VxLineElement.Vector';
%                    Vy=Obj.VyLineElement.Vector';
%                    Vz=Obj.VzLineElement.Vector';                   
%            end

           if Plot3D
               PP=[IP.SourcePosition';IP.CentrePoint'];
               plot3(PP(:,1),PP(:,2),PP(:,3),'m');
           end
           
           for x=MinX:Obj.MinorTickSize:MaxX
               P1=Org;
               P2=Org;
               switch Obj.PlaneType
                   case 'XY'
%                        P1=P1+[x,MinX,0];
%                        P2=P2+[x,MaxX,0];
                       P1=P1 + x*Vx + MinX*Vy + 0*Vz;
                       P2=P2 + x*Vx + MaxX*Vy + 0*Vz;                        
                   case 'YZ'
%                        P1=P1+[0,x,MinX];
%                        P2=P2+[0,x,MaxX];
                       P1=P1 + 0*Vx + x*Vy + MinX*Vz;
                       P2=P2 + 0*Vx + x*Vy + MaxX*Vz;                        
                   case 'XZ'
%                        P1=P1+[x,0,MinX];
%                        P2=P2+[x,0,MaxX];
                       P1=P1 + MinX*Vx + 0*Vy + x*Vz;
                       P2=P2 + MaxX*Vx + 0*Vy + x*Vz;                        
               end
                   
               PP=[P1;P2];
               if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);                   
               end
               Pv=IP.ImageNormal';
               Pp=IP.CentrePoint';
               Sp=IP.SourcePosition';
%                Lp=P1;
               Lv=(Sp-P1)/norm(Sp-P1);
               PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
               Lv=(Sp-P2)/norm(Sp-P2);
               PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
               PP=[PP1;PP2];
               if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               end
               R=size(IP.Image,1)*IP.PixelScale*1/2;
               PP=lineCircleIntersect(Pp,PP1,PP2,R);
               if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'r');%Obj.MinorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'r'); %Obj.MinorTickColor);
                    end
                    GXX=[GXX;PP2D(:,1)'];
                    GXY=[GXY;PP2D(:,2)'];
                    GXA=[GXA,x];
%                     CL=PP2D(:,1)*0+'r';
%                     CLR=[CLR;CL];
               end
           end
           for y=MinX:Obj.MinorTickSize:MaxX
               P1=Org;
               P2=Org;
               switch Obj.PlaneType
                   case 'XY'               
%                        P1=P1+[MinX,y,0];
%                        P2=P2+[MaxX,y,0];
                       P1=P1 + MinX*Vx + y*Vy + 0*Vz;
                       P2=P2 + MaxX*Vx + y*Vy + 0*Vz;                        
                   case 'YZ'
%                        P1=P1+[0,MinX,y];
%                        P2=P2+[0,MaxX,y];
                       P1=P1 + 0*Vx + MinX*Vy + y*Vz;
                       P2=P2 + 0*Vx + MaxX*Vy + y*Vz;                        
                   case 'XZ'
%                        P1=P1+[MinX,0,y];
%                        P2=P2+[MaxX,0,y];
                       P1=P1 + y*Vx + 0*Vy + MinX*Vz;
                       P2=P2 + y*Vx + 0*Vy + MaxX*Vz;                        
               end
               PP=[P1;P2];
               
               if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               end
               
               Pv=IP.ImageNormal';
               Pp=IP.CentrePoint';
               Sp=IP.SourcePosition';
%                Lp=P1;
               Lv=(Sp-P1)/norm(Sp-P1);
               PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
               Lv=(Sp-P2)/norm(Sp-P2);
               PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
               PP=[PP1;PP2];
               
               if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               end
               
               R=size(IP.Image,1)*IP.PixelScale*1/2;
               PP=lineCircleIntersect(Pp,PP1,PP2,R);
               if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'b'); %Obj.MinorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'b'); %Obj.MinorTickColor); 
                    end
                    GYX=[GYX;PP2D(:,1)'];
                    GYY=[GYY;PP2D(:,2)'];
                    GYA=[GYA,y];                    
%                     CL=PP2D(:,1)*0+'r';
%                     CLR=[CLR;CL];                    
               end
           end
           

           
           for x=MinX:Obj.MajorTickSize:MaxX
               P1=Org;
               P2=Org;
               switch Obj.PlaneType
                   case 'XY'
%                        P1=P1+[x,MinX,0];
%                        P2=P2+[x,MaxX,0];                       
                       P1=P1 + x*Vx + MinX*Vy + 0*Vz;
                       P2=P2 + x*Vx + MaxX*Vy + 0*Vz;                       
                       
                   case 'YZ'
%                        P1=P1+[0,x,MinX];
%                        P2=P2+[0,x,MaxX];
                       P1=P1 + 0*Vx + x*Vy + MinX*Vz;
                       P2=P2 + 0*Vx + x*Vy + MaxX*Vz;                       
                       
                   case 'XZ'
%                        P1=P1+[x,0,MinX];
%                        P2=P2+[x,0,MaxX];
                       P1=P1 + MinX*Vx + 0*Vy + x*Vz;
                       P2=P2 + MaxX*Vx + 0*Vy + x*Vz;                       
                       
               end              
               PP=[P1;P2];
%                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               Pv=IP.ImageNormal';
               Pp=IP.CentrePoint';
               Sp=IP.SourcePosition';
%                Lp=P1;
               Lv=(Sp-P1)/norm(Sp-P1);
               PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
               Lv=(Sp-P2)/norm(Sp-P2);
               PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
               PP=[PP1;PP2];
%                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               R=size(IP.Image,1)*IP.PixelScale*1/2;
               PP=lineCircleIntersect(Pp,PP1,PP2,R);
               if ~isempty(PP)
                   if Plot3D
                       figure(h1);
                       plot3(PP(:,1),PP(:,2),PP(:,3),'c'); %Obj.MajorTickColor);
                   end
                   PP2D=transTo2d(PP,IP);
                   if Plot2D
                       figure(h2);
                       plot(PP2D(:,1),PP2D(:,2),'c'); %Obj.MajorTickColor);
                   end
                   LXX=[LXX;PP2D(:,1)'];
                   LXY=[LXY;PP2D(:,2)'];
                   LXA=[LXA,x];
%                    CL=PP2D(:,1)*0+'r';
%                    CLR=[CLR;CL];
               end
           end
           for y=MinX:Obj.MajorTickSize:MaxX
               P1=Org;
               P2=Org;
               switch Obj.PlaneType
                   case 'XY'               
%                        P1=P1+[MinX,y,0];
%                        P2=P2+[MaxX,y,0];
                       P1=P1 + MinX*Vx + y*Vy + 0*Vz;
                       P2=P2 + MaxX*Vx + y*Vy + 0*Vz;                        
                   case 'YZ'
%                        P1=P1+[0,MinX,y];
%                        P2=P2+[0,MaxX,y];
                       P1=P1 + 0*Vx + MinX*Vy + y*Vz;
                       P2=P2 + 0*Vx + MaxX*Vy + y*Vz;                        
                       
                   case 'XZ'
%                        P1=P1+[MinX,0,y];
%                        P2=P2+[MaxX,0,y];
                       P1=P1 + y*Vx + 0*Vy + MinX*Vz;
                       P2=P2 + y*Vx + 0*Vy + MaxX*Vz;                                               
               end             
               PP=[P1;P2];
%                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               Pv=IP.ImageNormal';
               Pp=IP.CentrePoint';
               Sp=IP.SourcePosition';
%                Lp=P1;
               Lv=(Sp-P1)/norm(Sp-P1);
               PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
               Lv=(Sp-P2)/norm(Sp-P2);
               PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
               PP=[PP1;PP2];
%                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
               R=size(IP.Image,1)*IP.PixelScale*1/2;
               PP=lineCircleIntersect(Pp,PP1,PP2,R);
               if ~isempty(PP)
                   if Plot3D
                       figure(h1);
                       plot3(PP(:,1),PP(:,2),PP(:,3),'m'); %Obj.MajorTickColor);
                       axis equal;
                   end
                   PP2D=transTo2d(PP,IP);
                   if Plot2D                   
                       figure(h2);
                       plot(PP2D(:,1),PP2D(:,2),'m'); %Obj.MajorTickColor);
                       axis equal;
                   end
                   LYX=[LYX;PP2D(:,1)'];
                   LYY=[LYY;PP2D(:,2)'];
                   LYA=[LYA,y];
%                    CL=PP2D(:,1)*0+'r';
%                    CLR=[CLR;CL];                   
               end               
           end
           
           
           if Plot2D
               figure;
               title ('assem');
               hold on;
               plot(GXX',GXY',[Obj.MinorTickColor,':']);%,'Parent',hh);
               plot(GYX',GYY',[Obj.MinorTickColor,':']);%,'Parent',hh);               
               plot(LXX',LXY',[Obj.Axis1Color,'-']);%,'Parent',hh);
               plot(LYX',LYY',[Obj.Axis2Color,'-']);%,'Parent',hh);
           end

        end

        function [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=calcGrid_old(Obj,IP)
            Plot3D=0;
            Plot2D=0;
            %            LX=[];
            %            LY=[];
            %            GX=[];
            %            GY=[];
            %            CLR=[];
            %            ALX=[];
            %            ALY=[];
            %            AGX=[];
            %            AGY=[];
            
            GXX=[];
            GXY=[];
            GXA=[];
            GYX=[];
            GYY=[];
            GYA=[];
            
            LXX=[];
            LXY=[];
            LXA=[];
            LYX=[];
            LYY=[];
            LYA=[];
            
            
            if Plot3D
                h1=figure;
                IP.imshow3d(1,1);
                hold on;
            end
            if Plot2D
                h2=figure;
                figure(h2);
                imshow(IP.Image);
                hold on;
            end
            
            MinX=-2000;
            MaxX=2000;
            
            
            %--------------------------------------------------------------
            % Calculate the intersection between the centre of image and the
            % XY plane:
            Pv=Obj.PlaneNormal;
            Pp=Obj.PlaneOrigin;
            Lv=IP.ImageNormal;
            P2=IP.SourcePosition;
            Pnt=linePlaneIntersect(Pv',Pp',Lv',P2');
            %--------------------------------------------------------------
            
            %            XP=floor(Pnt(1)/Obj.MajorTickSize)*Obj.MajorTickSize;
            %            YP=floor(Pnt(2)/Obj.MajorTickSize)*Obj.MajorTickSize;
            %            ZP=floor(Pnt(3)/Obj.MajorTickSize)*Obj.MajorTickSize;
            %            Org=[XP,YP,ZP]; %round point
            Org=Pp';
            
            
            Vx=Obj.VxLineElement.Vector';
            Vy=Obj.VyLineElement.Vector';
            Vz=Obj.VzLineElement.Vector';
            
            %            switch Obj.PlaneType
            %                case 'XY'
            %                    Vx=Obj.VxLineElement.Vector';
            %                    Vy=Obj.VyLineElement.Vector';
            %                    Vz=Obj.VzLineElement.Vector';
            %                case 'YZ'
            %                    Vx=Obj.VxLineElement.Vector';
            %                    Vy=Obj.VyLineElement.Vector';
            %                    Vz=Obj.VzLineElement.Vector';
            %                case 'ZX'
            %                    Vx=Obj.VxLineElement.Vector';
            %                    Vy=Obj.VyLineElement.Vector';
            %                    Vz=Obj.VzLineElement.Vector';
            %            end
            
            if Plot3D
                PP=[IP.SourcePosition';IP.CentrePoint'];
                plot3(PP(:,1),PP(:,2),PP(:,3),'m');
            end
            
            for x=MinX:Obj.MinorTickSize:MaxX
                P1=Org;
                P2=Org;
                switch Obj.PlaneType
                    case 'XY'
                        %                        P1=P1+[x,MinX,0];
                        %                        P2=P2+[x,MaxX,0];
                        P1=P1 + x*Vx + MinX*Vy + 0*Vz;
                        P2=P2 + x*Vx + MaxX*Vy + 0*Vz;
                    case 'YZ'
                        %                        P1=P1+[0,x,MinX];
                        %                        P2=P2+[0,x,MaxX];
                        P1=P1 + 0*Vx + x*Vy + MinX*Vz;
                        P2=P2 + 0*Vx + x*Vy + MaxX*Vz;
                    case 'XZ'
                        %                        P1=P1+[x,0,MinX];
                        %                        P2=P2+[x,0,MaxX];
                        P1=P1 + MinX*Vx + 0*Vy + x*Vz;
                        P2=P2 + MaxX*Vx + 0*Vy + x*Vz;
                end
                
                PP=[P1;P2];
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'r');%Obj.MinorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'r'); %Obj.MinorTickColor);
                    end
                    GXX=[GXX;PP2D(:,1)'];
                    GXY=[GXY;PP2D(:,2)'];
                    GXA=[GXA,x];
                    %                     CL=PP2D(:,1)*0+'r';
                    %                     CLR=[CLR;CL];
                end
            end
            for y=MinX:Obj.MinorTickSize:MaxX
                P1=Org;
                P2=Org;
                switch Obj.PlaneType
                    case 'XY'
                        %                        P1=P1+[MinX,y,0];
                        %                        P2=P2+[MaxX,y,0];
                        P1=P1 + MinX*Vx + y*Vy + 0*Vz;
                        P2=P2 + MaxX*Vx + y*Vy + 0*Vz;
                    case 'YZ'
                        %                        P1=P1+[0,MinX,y];
                        %                        P2=P2+[0,MaxX,y];
                        P1=P1 + 0*Vx + MinX*Vy + y*Vz;
                        P2=P2 + 0*Vx + MaxX*Vy + y*Vz;
                    case 'XZ'
                        %                        P1=P1+[MinX,0,y];
                        %                        P2=P2+[MaxX,0,y];
                        P1=P1 + y*Vx + 0*Vy + MinX*Vz;
                        P2=P2 + y*Vx + 0*Vy + MaxX*Vz;
                end
                PP=[P1;P2];
                
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'b'); %Obj.MinorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'b'); %Obj.MinorTickColor);
                    end
                    GYX=[GYX;PP2D(:,1)'];
                    GYY=[GYY;PP2D(:,2)'];
                    GYA=[GYA,y];
                    %                     CL=PP2D(:,1)*0+'r';
                    %                     CLR=[CLR;CL];
                end
            end
                        
            
            for x=MinX:Obj.MajorTickSize:MaxX
                P1=Org;
                P2=Org;
                switch Obj.PlaneType
                    case 'XY'
                        %                        P1=P1+[x,MinX,0];
                        %                        P2=P2+[x,MaxX,0];
                        P1=P1 + x*Vx + MinX*Vy + 0*Vz;
                        P2=P2 + x*Vx + MaxX*Vy + 0*Vz;
                        
                    case 'YZ'
                        %                        P1=P1+[0,x,MinX];
                        %                        P2=P2+[0,x,MaxX];
                        P1=P1 + 0*Vx + x*Vy + MinX*Vz;
                        P2=P2 + 0*Vx + x*Vy + MaxX*Vz;
                        
                    case 'XZ'
                        %                        P1=P1+[x,0,MinX];
                        %                        P2=P2+[x,0,MaxX];
                        P1=P1 + MinX*Vx + 0*Vy + x*Vz;
                        P2=P2 + MaxX*Vx + 0*Vy + x*Vz;
                        
                end
                PP=[P1;P2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'c'); %Obj.MajorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'c'); %Obj.MajorTickColor);
                    end
                    LXX=[LXX;PP2D(:,1)'];
                    LXY=[LXY;PP2D(:,2)'];
                    LXA=[LXA,x];
                    %                    CL=PP2D(:,1)*0+'r';
                    %                    CLR=[CLR;CL];
                end
            end
            for y=MinX:Obj.MajorTickSize:MaxX
                P1=Org;
                P2=Org;
                switch Obj.PlaneType
                    case 'XY'
                        %                        P1=P1+[MinX,y,0];
                        %                        P2=P2+[MaxX,y,0];
                        P1=P1 + MinX*Vx + y*Vy + 0*Vz;
                        P2=P2 + MaxX*Vx + y*Vy + 0*Vz;
                    case 'YZ'
                        %                        P1=P1+[0,MinX,y];
                        %                        P2=P2+[0,MaxX,y];
                        P1=P1 + 0*Vx + MinX*Vy + y*Vz;
                        P2=P2 + 0*Vx + MaxX*Vy + y*Vz;
                        
                    case 'XZ'
                        %                        P1=P1+[MinX,0,y];
                        %                        P2=P2+[MaxX,0,y];
                        P1=P1 + y*Vx + 0*Vy + MinX*Vz;
                        P2=P2 + y*Vx + 0*Vy + MaxX*Vz;
                end
                PP=[P1;P2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'m'); %Obj.MajorTickColor);
                        axis equal;
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'m'); %Obj.MajorTickColor);
                        axis equal;
                    end
                    LYX=[LYX;PP2D(:,1)'];
                    LYY=[LYY;PP2D(:,2)'];
                    LYA=[LYA,y];
                    %                    CL=PP2D(:,1)*0+'r';
                    %                    CLR=[CLR;CL];
                end
            end
            
            
            if Plot2D
                figure;
                title ('assem');
                hold on;
                plot(GXX',GXY',[Obj.MinorTickColor,':']);%,'Parent',hh);
                plot(GYX',GYY',[Obj.MinorTickColor,':']);%,'Parent',hh);
                plot(LXX',LXY',[Obj.Axis1Color,'-']);%,'Parent',hh);
                plot(LYX',LYY',[Obj.Axis2Color,'-']);%,'Parent',hh);
            end
            
        end
        
        function [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=calcGrid(Obj,IP,PlaneType)
            Plot3D=0;
            Plot2D=0;
            %            LX=[];
            %            LY=[];
            %            GX=[];
            %            GY=[];
            %            CLR=[];
            %            ALX=[];
            %            ALY=[];
            %            AGX=[];
            %            AGY=[];
            
            GXX=[];
            GXY=[];
            GXA=[];
            GYX=[];
            GYY=[];
            GYA=[];
            
            LXX=[];
            LXY=[];
            LXA=[];
            LYX=[];
            LYY=[];
            LYA=[];
            
            
            if Plot3D
                h1=figure;
                IP.imshow3d(1,1);
                hold on;
            end
            if Plot2D
                h2=figure;
                figure(h2);
                imshow(IP.Image);
                hold on;
            end
            
            MinX=-2000;
            MaxX=2000;
            
            
            %--------------------------------------------------------------
            % Calculate the intersection between the centre of image and the
            % XY plane:
            Pv=Obj.PlaneNormal;
            Pp=Obj.PlaneOrigin;
            Lv=IP.ImageNormal;
            P2=IP.SourcePosition;
            Pnt=linePlaneIntersect(Pv',Pp',Lv',P2');
            %--------------------------------------------------------------
            
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            % Calculating the range of grid lines ": Feb 4, 2016
            % This avoids the slow rendering issues:
            
            % caculate the 4 conners of the Image (with a conservative 2
            % times of the size of the image: (Feb 4 2016)
            RV=cross(-IP.ImageNormal,IP.ImageUp);
            RV=RV/norm(RV);
            VectorUp=max(size(IP.Image))*IP.ImageUp*1.5; %*IP.PixelScale;
            VectorRight=max(size(IP.Image))*RV*1.5; %*IP.PixelScale;
            Pn1=IP.CentrePoint+VectorUp+VectorRight;
            Pn2=IP.CentrePoint+VectorUp-VectorRight;
            Pn3=IP.CentrePoint-VectorUp+VectorRight;
            Pn4=IP.CentrePoint-VectorUp-VectorRight;

            Lv1=IP.SourcePosition-Pn1;
            Lv1=Lv1/norm(Lv1);
            Lv2=IP.SourcePosition-Pn2;
            Lv2=Lv2/norm(Lv2);
            Lv3=IP.SourcePosition-Pn3;
            Lv3=Lv3/norm(Lv3);
            Lv4=IP.SourcePosition-Pn4;
            Lv4=Lv4/norm(Lv4);
            
            PPn1=linePlaneIntersect(Pv',Pp',-Lv1',P2');
            PPn2=linePlaneIntersect(Pv',Pp',-Lv2',P2');
            PPn3=linePlaneIntersect(Pv',Pp',-Lv3',P2');
            PPn4=linePlaneIntersect(Pv',Pp',-Lv4',P2');
            
            PPPn1=PPn1'-Obj.PlaneOrigin;
            PPPn2=PPn2'-Obj.PlaneOrigin;
            PPPn3=PPn3'-Obj.PlaneOrigin;
            PPPn4=PPn4'-Obj.PlaneOrigin;
            
            VVx1=dot(PPPn1,Obj.VxLineElement.Vector);
            VVx2=dot(PPPn2,Obj.VxLineElement.Vector);
            VVx3=dot(PPPn3,Obj.VxLineElement.Vector);
            VVx4=dot(PPPn4,Obj.VxLineElement.Vector);

            VVy1=dot(PPPn1,Obj.VyLineElement.Vector);
            VVy2=dot(PPPn2,Obj.VyLineElement.Vector);
            VVy3=dot(PPPn3,Obj.VyLineElement.Vector);
            VVy4=dot(PPPn4,Obj.VyLineElement.Vector);
            
            VVz1=dot(PPPn1,Obj.VzLineElement.Vector);
            VVz2=dot(PPPn2,Obj.VzLineElement.Vector);
            VVz3=dot(PPPn3,Obj.VzLineElement.Vector);
            VVz4=dot(PPPn4,Obj.VzLineElement.Vector);
                        
            XLims=[VVx1,VVx2,VVx3,VVx4];
            YLims=[VVy1,VVy2,VVy3,VVy4];
            ZLims=[VVz1,VVz2,VVz3,VVz4];
            
            XMinLim=min(XLims);
            XMaxLim=max(XLims);
            YMinLim=min(YLims);
            YMaxLim=max(YLims);
            ZMinLim=min(ZLims);
            ZMaxLim=max(ZLims);            
            
            %Round up to the closest MinorTickSize
%             Delt=Obj.MinorTickSize;
            Delt=Obj.MajorTickSize;
            XMinR=Delt*round(XMinLim/Delt);
            YMinR=Delt*round(YMinLim/Delt);
            ZMinR=Delt*round(ZMinLim/Delt);
            
            XMaxR=Delt*round(XMaxLim/Delt);
            YMaxR=Delt*round(YMaxLim/Delt);
            ZMaxR=Delt*round(ZMaxLim/Delt);
            
            switch PlaneType
                case 'XY'                   
                    MinX=XMinR;
                    MaxX=XMaxR;
                    MinY=YMinR;
                    MaxY=YMaxR;
                case 'YZ'                   
                    MinX=YMinR;
                    MaxX=YMaxR;
                    MinY=ZMinR;
                    MaxY=ZMaxR;                   
                case 'XZ'                   
                    MinX=ZMinR;
                    MaxX=ZMaxR;
                    MinY=XMinR;
                    MaxY=XMaxR;
            end

            % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            
            %            XP=floor(Pnt(1)/Obj.MajorTickSize)*Obj.MajorTickSize;
            %            YP=floor(Pnt(2)/Obj.MajorTickSize)*Obj.MajorTickSize;
            %            ZP=floor(Pnt(3)/Obj.MajorTickSize)*Obj.MajorTickSize;
            %            Org=[XP,YP,ZP]; %round point
            Org=Pp';
            
            
            Vx=Obj.VxLineElement.Vector';
            Vy=Obj.VyLineElement.Vector';
            Vz=Obj.VzLineElement.Vector';
            
            %            switch Obj.PlaneType
            %                case 'XY'
            %                    Vx=Obj.VxLineElement.Vector';
            %                    Vy=Obj.VyLineElement.Vector';
            %                    Vz=Obj.VzLineElement.Vector';
            %                case 'YZ'
            %                    Vx=Obj.VxLineElement.Vector';
            %                    Vy=Obj.VyLineElement.Vector';
            %                    Vz=Obj.VzLineElement.Vector';
            %                case 'ZX'
            %                    Vx=Obj.VxLineElement.Vector';
            %                    Vy=Obj.VyLineElement.Vector';
            %                    Vz=Obj.VzLineElement.Vector';
            %            end
            
            if Plot3D
                PP=[IP.SourcePosition';IP.CentrePoint'];
                plot3(PP(:,1),PP(:,2),PP(:,3),'m');
            end
            
            
            XMinors=setdiff(MinX:Obj.MinorTickSize:MaxX, MinX:Obj.MajorTickSize:MaxX);
            YMinors=setdiff(MinY:Obj.MinorTickSize:MaxY, MinY:Obj.MajorTickSize:MaxY);
            
            %for x=MinX:Obj.MinorTickSize:MaxX
            for x=XMinors %MinX:Obj.MinorTickSize:MaxX
                P1=Org;
                P2=Org;
                switch PlaneType
                    case 'XY'
                        %                        P1=P1+[x,MinX,0];
                        %                        P2=P2+[x,MaxX,0];
                        P1=P1 + x*Vx + MinX*Vy + 0*Vz;
                        P2=P2 + x*Vx + MaxX*Vy + 0*Vz;
                    case 'YZ'
                        %                        P1=P1+[0,x,MinX];
                        %                        P2=P2+[0,x,MaxX];
                        P1=P1 + 0*Vx + x*Vy + MinX*Vz;
                        P2=P2 + 0*Vx + x*Vy + MaxX*Vz;
                    case 'XZ'
                        %                        P1=P1+[x,0,MinX];
                        %                        P2=P2+[x,0,MaxX];
                        P1=P1 + MinX*Vx + 0*Vy + x*Vz;
                        P2=P2 + MaxX*Vx + 0*Vy + x*Vz;
                end
                
                PP=[P1;P2];
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'r');%Obj.MinorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'r'); %Obj.MinorTickColor);
                    end
                    GXX=[GXX;PP2D(:,1)'];
                    GXY=[GXY;PP2D(:,2)'];
                    GXA=[GXA,x];
                    %                     CL=PP2D(:,1)*0+'r';
                    %                     CLR=[CLR;CL];
                end
            end
            %for y=MinY:Obj.MinorTickSize:MaxY
            for y=YMinors %MinY:Obj.MinorTickSize:MaxY
                P1=Org;
                P2=Org;
                switch PlaneType
                    case 'XY'
                        %                        P1=P1+[MinX,y,0];
                        %                        P2=P2+[MaxX,y,0];
                        P1=P1 + MinX*Vx + y*Vy + 0*Vz;
                        P2=P2 + MaxX*Vx + y*Vy + 0*Vz;
                    case 'YZ'
                        %                        P1=P1+[0,MinX,y];
                        %                        P2=P2+[0,MaxX,y];
                        P1=P1 + 0*Vx + MinX*Vy + y*Vz;
                        P2=P2 + 0*Vx + MaxX*Vy + y*Vz;
                    case 'XZ'
                        %                        P1=P1+[MinX,0,y];
                        %                        P2=P2+[MaxX,0,y];
                        P1=P1 + y*Vx + 0*Vy + MinX*Vz;
                        P2=P2 + y*Vx + 0*Vy + MaxX*Vz;
                end
                PP=[P1;P2];
                
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                
                if Plot3D
                    plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                end
                
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'b'); %Obj.MinorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'b'); %Obj.MinorTickColor);
                    end
                    GYX=[GYX;PP2D(:,1)'];
                    GYY=[GYY;PP2D(:,2)'];
                    GYA=[GYA,y];
                    %                     CL=PP2D(:,1)*0+'r';
                    %                     CLR=[CLR;CL];
                end
            end
                        
            
            for x=MinX:Obj.MajorTickSize:MaxX
                P1=Org;
                P2=Org;
                switch PlaneType
                    case 'XY'
                        %                        P1=P1+[x,MinX,0];
                        %                        P2=P2+[x,MaxX,0];
                        P1=P1 + x*Vx + MinX*Vy + 0*Vz;
                        P2=P2 + x*Vx + MaxX*Vy + 0*Vz;
                        
                    case 'YZ'
                        %                        P1=P1+[0,x,MinX];
                        %                        P2=P2+[0,x,MaxX];
                        P1=P1 + 0*Vx + x*Vy + MinX*Vz;
                        P2=P2 + 0*Vx + x*Vy + MaxX*Vz;
                        
                    case 'XZ'
                        %                        P1=P1+[x,0,MinX];
                        %                        P2=P2+[x,0,MaxX];
                        P1=P1 + MinX*Vx + 0*Vy + x*Vz;
                        P2=P2 + MaxX*Vx + 0*Vy + x*Vz;
                        
                end
                PP=[P1;P2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'c'); %Obj.MajorTickColor);
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'c'); %Obj.MajorTickColor);
                    end
                    LXX=[LXX;PP2D(:,1)'];
                    LXY=[LXY;PP2D(:,2)'];
                    LXA=[LXA,x];
                    %                    CL=PP2D(:,1)*0+'r';
                    %                    CLR=[CLR;CL];
                end
            end
            for y=MinY:Obj.MajorTickSize:MaxY
                P1=Org;
                P2=Org;
                switch PlaneType
                    case 'XY'
                        %                        P1=P1+[MinX,y,0];
                        %                        P2=P2+[MaxX,y,0];
                        P1=P1 + MinX*Vx + y*Vy + 0*Vz;
                        P2=P2 + MaxX*Vx + y*Vy + 0*Vz;
                    case 'YZ'
                        %                        P1=P1+[0,MinX,y];
                        %                        P2=P2+[0,MaxX,y];
                        P1=P1 + 0*Vx + MinX*Vy + y*Vz;
                        P2=P2 + 0*Vx + MaxX*Vy + y*Vz;
                        
                    case 'XZ'
                        %                        P1=P1+[MinX,0,y];
                        %                        P2=P2+[MaxX,0,y];
                        P1=P1 + y*Vx + 0*Vy + MinX*Vz;
                        P2=P2 + y*Vx + 0*Vy + MaxX*Vz;
                end
                PP=[P1;P2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                Pv=IP.ImageNormal';
                Pp=IP.CentrePoint';
                Sp=IP.SourcePosition';
                %                Lp=P1;
                Lv=(Sp-P1)/norm(Sp-P1);
                PP1=linePlaneIntersect(Pv,Pp,Lv,P1);
                Lv=(Sp-P2)/norm(Sp-P2);
                PP2=linePlaneIntersect(Pv,Pp,Lv,P2);
                PP=[PP1;PP2];
                %                plot3(PP(:,1),PP(:,2),PP(:,3),Obj.MinorTickColor);
                R=size(IP.Image,1)*IP.PixelScale*1/2;
                PP=lineCircleIntersect(Pp,PP1,PP2,R);
                if ~isempty(PP)
                    if Plot3D
                        figure(h1);
                        plot3(PP(:,1),PP(:,2),PP(:,3),'m'); %Obj.MajorTickColor);
                        axis equal;
                    end
                    PP2D=transTo2d(PP,IP);
                    if Plot2D
                        figure(h2);
                        plot(PP2D(:,1),PP2D(:,2),'m'); %Obj.MajorTickColor);
                        axis equal;
                    end
                    LYX=[LYX;PP2D(:,1)'];
                    LYY=[LYY;PP2D(:,2)'];
                    LYA=[LYA,y];
                    %                    CL=PP2D(:,1)*0+'r';
                    %                    CLR=[CLR;CL];
                end
            end
            
            
            if Plot2D
                figure;
                title ('assem');
                hold on;
                plot(GXX',GXY',[Obj.MinorTickColor,':']);%,'Parent',hh);
                plot(GYX',GYY',[Obj.MinorTickColor,':']);%,'Parent',hh);
                plot(LXX',LXY',[Obj.Axis1Color,'-']);%,'Parent',hh);
                plot(LYX',LYY',[Obj.Axis2Color,'-']);%,'Parent',hh);
            end
            
        end
        
        
        
        function PO=createProj(Obj,IP,axesHandle)
            % Function: Creates two SphereElement projections, with one at
            %   the centre. Outputs the plot object, which contains a
            %   property that has all the handles for the plotted spheres
            %   and lines            
            PO = [];
            theHandles = Obj.createPlot(IP,axesHandle);
            extract = extractPlotProps(theHandles);
            PO = PlotObject(Obj,extract,theHandles,axesHandle);
        end
        
        function h=createPlot(Obj,IP,parent,varargin)
            
            [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=Obj.calcGrid(IP,Obj.PlaneType);
%             [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=Obj.calcGrid_old2(IP,Obj.PlaneType);
                        
            % extra segments to be added to the handles so any size can be
            % handles:

            n=100;
            repm=repmat([-999,-999],n,1);
            reps=repmat(-999,n,1);
            GXX(end+1:end+n,:)=repm;
            GXY(end+1:end+n,:)=repm;
            GXA(end+1:end+n)=reps;

            GYX(end+1:end+n,:)=repm;
            GYY(end+1:end+n,:)=repm;
            GYA(end+1:end+n)=reps;
            
            LXX(end+1:end+n,:)=repm;
            LXY(end+1:end+n,:)=repm;
            LXA(end+1:end+n)=reps;

            LYX(end+1:end+n,:)=repm;
            LYY(end+1:end+n,:)=repm;
            LYA(end+1:end+n)=reps;
                                   
            h={};
            h{1}=plot(GXX',GXY',[Obj.MinorTickColor,':'],'Parent',parent);
            h{2}=plot(GYX',GYY',[Obj.MinorTickColor,':'],'Parent',parent);

            h{3}=plot(LXX',LXY',[Obj.Axis1Color,'-'],'Parent',parent);
            h{4}=plot(LYX',LYY',[Obj.Axis2Color,'-'],'Parent',parent);
                        
            XX=mean(LXX');
            XY=mean(LXY');
            
            YX=mean(LYX');
            YY=mean(LYY');
            
            dx=-14;
            dy=-14;
            switch Obj.PlaneType
                case 'XY'
                    Axis1Label='X';
                    Axis2Label='Y';
                case 'YZ'
                    Axis1Label='Y';
                    Axis2Label='Z';                    
                case 'XZ'
                    Axis1Label='X';
                    Axis2Label='Z';                    
            end
            
            Ang=[];
            Txt={};
            for f=1:max(size(LXA))
                Txt{f}=sprintf('%s: %i',Axis1Label,LXA(f));
                A=[LXX(f,1),LXY(f,1)];
                B=[LXX(f,2),LXY(f,2)];
                V=B-A;
                Ang(f)=atan(V(2)/V(1))*180/pi;
            end
            hh=text(XX+dx,XY+dy,Txt,'color',Obj.Axis1Color,'rotation',mean(Ang),'FontSize',16,'Parent',parent);
            h{5}=hh;

            Ang=[];
            Txt={};
            for f=1:max(size(LYA))
                Txt{f}=sprintf('%s: %i',Axis2Label,LYA(f));
                A=[LYX(f,1),LYY(f,1)];
                B=[LYX(f,2),LYY(f,2)];
                V=B-A;
                Ang(f)=atan(V(2)/V(1))*180/pi;                
            end
            hh=text(YX+dx,YY+dy,Txt,'color',Obj.Axis2Color,'rotation',mean(Ang),'FontSize',16,'Parent',parent);
            h{6}=hh;

        end

        
        
        function h=createFullPlot(Obj,IP,parent,PlaneType,varargin)
            
            [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=Obj.calcGrid(IP,PlaneType);
            
            nLXA=numel(LXA);
            nLYA=numel(LYA);
            % extra segments to be added to the handles so any size can be
            % handles:

            n=100;
            repm=repmat([-999,-999],n,1);
            reps=repmat(-999,n,1);
            GXX(end+1:end+n,:)=repm;
            GXY(end+1:end+n,:)=repm;
            GXA(end+1:end+n)=reps;

            GYX(end+1:end+n,:)=repm;
            GYY(end+1:end+n,:)=repm;
            GYA(end+1:end+n)=reps;
            
            LXX(end+1:end+n,:)=repm;
            LXY(end+1:end+n,:)=repm;
            LXA(end+1:end+n)=reps;

            LYX(end+1:end+n,:)=repm;
            LYY(end+1:end+n,:)=repm;
            LYA(end+1:end+n)=reps;
            
            Axis1Color='w';
            Axis2Color='w';
            MajorTickColor='w';
            MinorTickColor='w';
            
                                   
            h={};
            h{1}=plot(GXX',GXY',[MinorTickColor,':'],'Parent',parent);
            h{2}=plot(GYX',GYY',[MinorTickColor,':'],'Parent',parent);

            h{3}=plot(LXX',LXY',[Axis1Color,'-'],'Parent',parent);
            h{4}=plot(LYX',LYY',[Axis2Color,'-'],'Parent',parent);
                        
            XX=mean(LXX');
            XY=mean(LXY');
            
            YX=mean(LYX');
            YY=mean(LYY');
            
            dx=-14;
            dy=-14;
            switch PlaneType
                case 'XY'
                    Axis1Label='X';
                    Axis2Label='Y';
                case 'YZ'
                    Axis1Label='Y';
                    Axis2Label='Z';                    
                case 'XZ'
                    Axis1Label='X';
                    Axis2Label='Z';                    
            end
            
            Ang=[];
            Txt={};
            % for f=1:max(size(LXA))
            for f=1:nLXA
                Txt{f}=sprintf('%s: %i',Axis1Label,LXA(f));
                A=[LXX(f,1),LXY(f,1)];
                B=[LXX(f,2),LXY(f,2)];
                V=B-A;
                Ang(f)=atan(V(2)/V(1))*180/pi;
            end
            hh=text(XX(1:nLXA)+dx,XY(1:nLXA)+dy,Txt,'color',Axis1Color,'rotation',mean(Ang),'FontSize',16,'Parent',parent);
            h{5}=hh;

            Ang=[];
            Txt={};
            % for f=1:max(size(LYA))
            for f=1:nLYA                
                Txt{f}=sprintf('%s: %i',Axis2Label,LYA(f));
                A=[LYX(f,1),LYY(f,1)];
                B=[LYX(f,2),LYY(f,2)];
                V=B-A;
                Ang(f)=atan(V(2)/V(1))*180/pi;                
            end
            hh=text(YX(1:nLYA)+dx,YY(1:nLYA)+dy,Txt,'color',Axis2Color,'rotation',mean(Ang),'FontSize',16,'Parent',parent);
            h{6}=hh;

        end
     
        
        
        function updateProj(Obj,plotObj,IP)                         
            Obj.updatePlot(IP,plotObj.theHandle);
            extract = extractPlotProps(plotObj.theHandle);
            plotObj.setProps(extract);            
        end

        function updatePlot(Obj,IP,h,varargin)
            
            [GXX,GXY,GXA,GYX,GYY,GYA,LXX,LXY,LXA,LYX,LYY,LYA]=Obj.calcGrid(IP,Obj.PlaneType);
            % set(h{1},'xdata',num2cell(GXX,2));
            % set(h{1},'ydata',num2cell(GXY,2));
                                               
            
            sethandlexydata(h{1},GXX,GXY,[Obj.MinorTickColor,':']);            
            sethandlexydata(h{2},GYX,GYY,[Obj.MinorTickColor,':']);            
            sethandlexydata(h{3},LXX,LXY,[Obj.Axis1Color,'-']);            
            sethandlexydata(h{4},LYX,LYY,[Obj.Axis2Color,'-']);            

            XX=mean(LXX');
            XY=mean(LXY');
            
            YX=mean(LYX');
            YY=mean(LYY');
            
            dx=-14;
            dy=-14;
            switch Obj.PlaneType
                case 'XY'
                    Axis1Label='X';
                    Axis2Label='Y';
                case 'YZ'
                    Axis1Label='Y';
                    Axis2Label='Z';                    
                case 'XZ'
                    Axis1Label='Z';
                    Axis2Label='X';                    
            end
            
            Ang=[];
            Txt={};
            for f=1:max(size(LXA))
                Txt{f}=sprintf('%s:%i',Axis1Label,LXA(f));
                A=[LXX(f,1),LXY(f,1)];
                B=[LXX(f,2),LXY(f,2)];
                V=B-A;
                Ang(f)=atan(V(2)/V(1))*180/pi;
            end
%             hh=text(XX+dx,XY+dy,Txt,'color',Obj.Axis1Color,'rotation',mean(Ang),'FontSize',12,'Parent',parent);
%             h{5}=hh;
            sethandletextdata(h{5},XX+dx,XY+dy,Txt,Obj.Axis1Color,mean(Ang));
%             set(h{5},'xdata',XX+dx);
%             set(h{5},'ydata',XY+dy);
%             set(h{5},'color',Obj.Axis1Color);
%             ser(h{5},'rotation',mean(Ang));
            
            Ang=[];
            Txt={};
            for f=1:max(size(LYA))
                Txt{f}=sprintf('%s:%i',Axis2Label,LYA(f));
                A=[LYX(f,1),LYY(f,1)];
                B=[LYX(f,2),LYY(f,2)];
                V=B-A;
                Ang(f)=atan(V(2)/V(1))*180/pi;                
            end
%             hh=text(YX+dx,YY+dy,Txt,'color',Obj.Axis2Color,'rotation',mean(Ang),'FontSize',12,'Parent',parent);
%             h{6}=hh;
            sethandletextdata(h{6},YX+dx,YY+dy,Txt,Obj.Axis2Color,mean(Ang));
%             set(h{6},'xdata',XX+dx);
%             set(h{6},'ydata',XY+dy);
%             set(h{6},'color',Obj.Axis2Color);
%             ser(h{6},'rotation',mean(Ang));
            
        end
        
        function projPoint = projPointOnPlane(Obj,Point,NormalVec,AnyPointOnPlane)
            % Function: calculates the coordinates of the projection of a
            %   3D point onto a plane
            % Output: projPoint - a 3D vector containing the projected
            %   point
            
            % Get x, y, and z of coordinate of interest
            x0 = Point(1);
            y0 = Point(2);
            z0 = Point(3);
            
            % Separate coordinates of any point on the plane
            c1 = AnyPointOnPlane(1);
            c2 = AnyPointOnPlane(2);
            c3 = AnyPointOnPlane(3);
            
            % get the x, y and z components of a normal vector to the plane
            n1 = NormalVec(1);
            n2 = NormalVec(2);
            n3 = NormalVec(3);
            
            % calculate the line parameter t needed to give a point on the
            % plane that is directly underneath the point of interest
            t = (n1*(c1 - x0) + n2*(c2 - y0) + n3*(c3 - z0))/(n1^2 + n2^2 + n3^2);
            
            projPoint = [x0 + t*n1, y0 + t*n2, z0 + t*n3];
            
        end
        
        function Colour = checkProximity(Obj,cent1,cent2)
            % Check if TargetElement circles are close enough. If they are,
            % output a parameter with the string 'green'. Otherwise, use
            % the default color stored in Obj.Colour
            diff = norm(cent1-cent2);
            
            if diff < Obj.ProximitySensitivity 
                Colour = 'green';
            else
                Colour = Obj.Colour;
            end     
        end
        
        function TargetCentre = calculateTargetPosition(Obj,IP,Centre)
            % Function: Calculates an appropriate location to plot the
            %   moving part of the TargetElement plot object. Uses the
            %   centre point of the ImagePlane, and the coordinates of the
            %   Selected elements to do this. Staring directly normal to
            %   the image plane, if the vector spanning the centers between
            %   the SelectedElements is parallel, the TargetElement spheres
            %   will be exactly lined ontop of each other. As the angle
            %   between the image plane normal and SelectedElements vector
            %   grows, up to a limit of 90 degrees, the TargetElement
            %   spheres get farther apart.
            % Output: the 3d coordinates that the moving part of the Target
            %   Element will be placed at
            
                % Get centre points of each of the elements that the target
                % element is pointing to and from
                Point1Centre = Obj.LineElement.Coords;
                Point2Centre = Obj.LineElement.Coords + Obj.LineElement.Vector;

                % Need the normal vector to the current image plane
                IPNormal = IP.ImageNormal;

                % vector that defines the direction of the target element
                TargetVector = Obj.LineElement.Vector;

                % the angle between the normal vector and target vector (rad)
                angle = atan2(norm(cross(IPNormal,TargetVector)),dot(IPNormal,TargetVector));
                
                if angle > pi/2
                    angle = pi - angle;
                end
                
                % calculate magnitude of the distance target will be away from
                % the centre of the IP
                ScalingFactor = 200;
                distMag = (angle/pi)*ScalingFactor; 
                
                % Get direction that target vector should be pointing. This
                % is done by projecting the centre points of the elements
                % situated at the head and tail of the target onto the
                % image plane. Subsequently, a vector is calculated from
                % these projected points and normalized. This vector is
                % then added to the centre point of the plane to determine
                % what the coordinates of the target should be placed at.
                ProjPoint1 = Obj.projPointOnPlane(Point1Centre,IPNormal,IP.CentrePoint);
                ProjPoint2 = Obj.projPointOnPlane(Point2Centre,IPNormal,IP.CentrePoint);
                ProjVector = ProjPoint2 - ProjPoint1;
                ProjUnitVector = ProjVector/norm(ProjVector);
                
                TargetCentre = distMag*ProjUnitVector' + Centre;
        end
        
        function hh = plotDirectionPointer(Obj,IP,parent)
            % Function: Plots an arrow near the edge of the screen that
            %   points from the centre of the image plane to the 1st point
            %   between the two SelectedElements
            % Output: hh - the handle to the plotted arrow
            
            % Get a 2x2 array of arrow coordinates projection onto the
            % image plane. x coordinates in column 1, y coordinates in col 2
            Proj = Obj.calculateArrowCoordinates(IP);
            
            Colour = 'green';
            
            % Plot arrow and assign handle of plot to output
            hh = plot(Proj(:,1),Proj(:,2),'-', 'Color', Colour,...
                'LineWidth', 2, 'Parent', parent,'MarkerSize',50);            
            
        end
        
        function Proj = calculateArrowCoordinates(Obj,IP)
            % Function: Looks at the centre of the ImagePlane and the first
            %   point of the two SelectedElements, and uses the point's
            %   coordinates to calculate an appropriate location to plot
            %   the base and tip of the arrow.
            % Output: A 2x2 array with the coordinates for the base and tip
            %   of the arrow. X coordinates are in column 1 and Y coord.
            %   are in col 2. 
            
            % Get the projection of the point 1 of the two landmarks that
            % were registered in previous steps
            point1 = Obj.LineElement.Coords;
            point1Proj = Obj.projPointOnPlane(point1,IP.ImageNormal,IP.CentrePoint);
            
            % Subtract center of ImagePlane from the point's projection.
            % Make this vector a unit vector.
            
            centre2point = point1Proj' - IP.CentrePoint; 
            unitVec = centre2point/norm(centre2point);
            
            % From the center of the image plane, calculate 2 distances in
            % the direction of the unit vector that will be used as the
            % coordinates for the base and tip of the arrow
            
            arrowBase = IP.CentrePoint + unitVec*70;
            arrowTip = IP.CentrePoint + unitVec*90;
            
            BaseProj = projSpherePoints(IP,arrowBase,0);
            TipProj = projSpherePoints(IP,arrowTip,0);
            
            Proj = [ BaseProj(1,:); TipProj(1,:) ];
        end
        
        function showArrowOn = LandmarksOnIP(Obj,IP)
            % Function: Calculate the distance between the centre of the IP and the
            %   projection onto the IP of the first point within the LineElement
            %   object. If the projected point is too close to the centre of the IP,
            %   then output variable to turn visibility of arrow off
            % Output: showArrowOn - 1 if the IP centre is sufficiently far
            %   from the point, otherwise 0
            point = Obj.projPointOnPlane(Obj.LineElement.Coords,IP.ImageNormal,IP.CentrePoint);
            diff = norm(IP.CentrePoint - point');
            
            if diff < 70
                % do not show arrow if IP centre is too close to
                % the first point of the landmarks
                showArrowOn = 0;   
            else
                showArrowOn = 1;
            end
        end
        
        function hh = plotImage(Obj,imageSrc,parent,arrowCoords,scale,varargin)
            % Function: Plots the image at the location specified in the
            %   'imageSrc' parameter, on the axes specified in 'parent', and at
            %   the location defined by 'center' (a 2x1 vector). Scaling is
            %   defined in 'scale'. The varargin is for optional rotation
            %   parameters. The first argument of varargin should be the
            %   current coordinates of the center target circle 
            %   the second argument should be 1 for rotation
            %   on and 0 for rotation off
            % Output: hh - the handle to the image object
            
            if ~isempty(varargin)
                IPCentre = varargin{1};
                rotationOn = varargin{2}; % either 0(off) or 1(on)
            else
                rotationOn = 0;
            end
            
            % get the image data, and resize as desired
            [im, ~, alpha] = imread(imageSrc);
            im = imresize(im,scale);
            alpha = imresize(alpha,scale);
            
            if rotationOn == 1
                
                if isempty(Obj.arrowImgData) % only do once to save time
                    Obj.arrowImgData = {im,alpha}; % store non-rotated arrow info here for use in updateProj()
                end
                
                [im,alpha] = Obj.rotateArrowObj(im,alpha,arrowCoords,IPCentre);    
                
            end
            
            f = imshow(im);
            
            % set the transparency of the image
            set(f,'AlphaData', alpha)
            % set the axes the image should be plotted on
            set(f,'Parent',parent);
            % Calculate the size of the image, then use that information to
            % plot the image in relation to its center
            imgObject = get(f);
            imgSize = size(imgObject.CData);
            set(f,'xdata',arrowCoords(1)-imgSize(1)/2,'ydata',arrowCoords(2)-imgSize(1)/2)
            
            hh = f;
        end
        
        function [im,alpha] = rotateArrowObj(Obj,im,alpha,arrowCoords,IPCentre)
            % Function: Takes in the image pixel data (im), alpha channel 
            %   data (alpha), coords for where the Arrow is supposed to be
            %   plotted (arrowCoords) and the coords for the center of the
            %   current IP. All coordinates are 2d, with x in the first
            %   column and y in the second column. This function is meant
            %   to rotate an arrow so that it's new x axis is aligned with
            %   a vector spanning the center of the IP to the tip of where
            %   the arrow should be plotted. In other words, the arrow is
            %   rotated to point from the center of the IP to where it is
            %   plotted.
            % Output: im - new, rotated, image pixel data
            %         alpha - new, rotated, alpha channel data
            
            % For rotating the arrow png
            IPCenter2Arrow = arrowCoords - IPCentre;
            angle = atan2(norm(cross([IPCenter2Arrow, 0],[1 0 0])),dot([IPCenter2Arrow, 0],[1 0 0]));

            if arrowCoords(2) >= IPCentre(2)
                % atan2(..) only calculates the magnitude of the angle
                % between the x axis and the IPCenter2Arrow vector.
                % Must figure out whether the arrow is above or beneath
                % the axis
                angle = -angle;
            end

            angle = angle*180/pi; % in degrees
            im = imrotate(im,angle);
            alpha = imrotate(alpha,angle);
        end
        
        function initializeLastIPInfo(Obj,axesHandle,cent1,cent2)
            % Function: Create a structure LastIPInfo assigned to the
            %   TargetElement object's properties that keeps track of the
            %   last colour for each axes. Important for updating the
            %   target png images properly. The structure is created
            %   anew if the number of IPs change and the reconElement
            %   button is pressed again.
            % Output: None - only assigns structure to object property
            if isempty(Obj.LastIPInfo) || (max(size(Obj.IPs)) ~= max(size(Obj.LastIPInfo.AllIPs))) || any(Obj.lastLineElementVector ~= Obj.LineElement.Vector)
                
                Obj.LastIPInfo = struct();
                Obj.LastIPInfo.AllIPs = {};

                for i = 1:max(size(Obj.IPs))
                    Obj.LastIPInfo.AllIPs = {Obj.LastIPInfo.AllIPs{:},Obj.IPs{i}.Name}; 
                end
                
                % get the current axes and IP pair
                Obj.LastIPInfo.Axes.Handle = axesHandle;
                Obj.LastIPInfo.Axes.LastColour = {Obj.checkProximity(cent1,cent2)};
                
                % save the last line element vector in the case that it
                % changed
                Obj.lastLineElementVector = Obj.LineElement.Vector;
                
            else if max(size(Obj.LastIPInfo.Axes.Handle)) < 2
                    % only want two axes added to the Axes struct, even if
                    % reconstruct is pressed again
                    Obj.LastIPInfo.Axes.Handle = [Obj.LastIPInfo.Axes.Handle,axesHandle];
                    Obj.LastIPInfo.Axes.LastColour = {Obj.LastIPInfo.Axes.LastColour{:},Obj.checkProximity(cent1,cent2)};
                end
            end
            
        end
        
        
        function [B, varargout]=getConstant(Obj,CName)
            B=[];
            for f=1:size(Obj.Constants,2)
                if strcmp(CName, Obj.Constants(f).Name)
                    B=Obj.Constants(f);
                    varargout{1} = f;
                end
            end
        end
        
        function addConstant(Obj,CName,Value)
            C=Constant(CName,Value);
            T=Obj.getConstant(CName);
            if isempty(T)
                Obj.Constants=[Obj.Constants,C];
            else
                for f=1:size(Obj.Constants,2)
                    if strcmp(CName, Obj.Constants(f).Name)
                        Obj.Constants(f).Value=Value;
                    end
                end
            end
        end
        
        
        function Obj = CopyfromGridElement(Obj,RGrid)
%             Obj.Name = RGrid.Name;
            %Obj.ExtraPlots=[];
            %Obj.Constants=[];

            %Obj.TypeName=RGrid.TypeName;            
                        
            for f=1:Obj.ExtraElement.InitialAnchorNumber
                Obj.ExtraElement.AnchorPointSpheres(f).Centre=RGrid.ExtraElement.AnchorPointSpheres(f).Centre;
                Obj.ExtraElement.AnchorPointSpheres(f).Radius=RGrid.ExtraElement.AnchorPointSpheres(f).Radius;
            end
            Obj.Constants=RGrid.Constants;
            
            
            %Obj.PlotType = RGrid.PlotType;
            Obj.Colour = RGrid.Colour;
            %Obj.PlaneType = RGrid.PlaneType;
            Obj.MajorTickSize=RGrid.MajorTickSize;
            Obj.MinorTickSize=RGrid.MinorTickSize;
            Obj.MajorTickColor=RGrid.MajorTickColor;
            Obj.MinorTickColor=RGrid.MinorTickColor;
            %Obj.Guided=RGrid.Guided;
            %Obj.AutoPlaneSelect=RGrid.AutoPlaneSelect;

            %Obj.LayerName = RGrid.LayerName;
            %Obj.SelLayer = RGrid.SelLayer;
            Obj.InitialAnchorNumber = RGrid.InitialAnchorNumber;
            
            %Obj.OriginSphereElement=SphereElement('Orig',OrgPoint,OrgRad); % Origin of the Grid Element
%             Obj.OriginSphereElement.Centre=RGrid.OriginSphereElement.Centre;
%            Obj.MeasureAnchorPoint.Centre=RGrid.MeasureAnchorPoint.Centre;
%             Obj.VxLineElement=RGrid.VxLineElement;
%             Obj.VyLineElement=RGrid.VyLineElement;
%             Obj.VzLineElement=RGrid.VzLineElement;
            
            for f=1:2
                Obj.VxLineElement.AnchorPointSpheres(f).Centre=RGrid.VxLineElement.AnchorPointSpheres(f).Centre;
                Obj.VyLineElement.AnchorPointSpheres(f).Centre=RGrid.VyLineElement.AnchorPointSpheres(f).Centre;
                Obj.VzLineElement.AnchorPointSpheres(f).Centre=RGrid.VzLineElement.AnchorPointSpheres(f).Centre;
                
                Obj.VxLineElement.AnchorPointSpheres(f).Radius=RGrid.VxLineElement.AnchorPointSpheres(f).Radius;
                Obj.VyLineElement.AnchorPointSpheres(f).Radius=RGrid.VyLineElement.AnchorPointSpheres(f).Radius;
                Obj.VzLineElement.AnchorPointSpheres(f).Radius=RGrid.VzLineElement.AnchorPointSpheres(f).Radius;                
                                                              
            end
            
            Obj.VxLineElement.Coords=RGrid.VxLineElement.Coords;
            Obj.VyLineElement.Coords=RGrid.VyLineElement.Coords;
            Obj.VzLineElement.Coords=RGrid.VzLineElement.Coords;
            
            Obj.VxLineElement.Vector=RGrid.VxLineElement.Vector;
            Obj.VyLineElement.Vector=RGrid.VyLineElement.Vector;
            Obj.VzLineElement.Vector=RGrid.VzLineElement.Vector;
            
            %Obj.PlaneType=RGrid.PlaneType;
            Obj.PlaneNormal=RGrid.PlaneNormal;
            Obj.PlaneOrigin=RGrid.PlaneOrigin;
            Obj.Axis1Color=RGrid.Axis1Color;
            Obj.Axis2Color=RGrid.Axis2Color;
            
            % Copy the ExtraElement property:
            Obj.ExtraElement.InitialAnchorNumber=RGrid.ExtraElement.InitialAnchorNumber;
            for f=1:RGrid.ExtraElement.InitialAnchorNumber
                Obj.ExtraElement.AnchorPointSpheres(f).Centre=RGrid.ExtraElement.AnchorPointSpheres(f).Centre;
                Obj.ExtraElement.AnchorPointSpheres(f).Radius=RGrid.ExtraElement.AnchorPointSpheres(f).Radius;
                Obj.ExtraElement.AnchorPointSpheres(f).Name=RGrid.ExtraElement.AnchorPointSpheres(f).Name;
            end
        end
        
        
        
        
    end
    
end


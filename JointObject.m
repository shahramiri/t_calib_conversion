classdef JointObject < handle
    %VISOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        Type
        Origin
        Dir    
        Trans
        ResetValue
        ScaleFactor
        SensorRead
        PoseFieldName
    end  
    
    methods
        
        function Obj=JointObject(BName,BType,FdName)
            Obj.Name=BName;
            Obj.Type=BType;
            Obj.PoseFieldName=FdName;
        end   
        
        function T=getTmat(Obj,varargin)
            switch Obj.Type
                
                case 'Translation' %obsolete
%                     V=(Obj.SensorRead-Obj.ResetValue)*Obj.ScaleFactor*Obj.Dir;
%                     T=eye(4);
%                     T(1,4)=V(1);
%                     T(2,4)=V(2);
%                     T(3,4)=V(3);
                    
                case 'Rotation' %obsolete
%                     V=(Obj.SensorRead-Obj.ResetValue)*Obj.ScaleFactor;
%                     T=AxelRot(V,Obj.Dir,Obj.Origin);
                    
                case 'SixDOF'                    
                    T_Ref_Imspace=Obj.Trans;
                    T_Marker_Camera=Obj.SensorRead;
                    T_Ref_Marker=varargin{1}; % Retreive from marker registration matrix from the input (supposed to be available for this case)
                    T=T_Marker_Camera*T_Ref_Marker/T_Ref_Imspace; % Only return the relevant part of the chain of matrices

            end
        end
                      
        function ImageBasedCalibSixDOF(Obj, TT_R, TT_O, TT_X, TT_Yaw_1, TT_Yaw_2, TT_Y1, TT_Y2) 
            % Finding the Vx and Vy Directions First:
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            % Vx:
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^            
            %Slab wrt Image Transformations:
            T_R=TT_R; %Pose_R.Tmat; %reference pose:
            T_O=TT_O; %Pose_O.Tmat;
            T_X=TT_X; %Pose_X.Tmat;
        
            T_S_I_O=T_O;
            T_S_I_X=T_X;
            T_X_O=T_S_I_X/T_S_I_O;
            
            Vx_Vector=T_X_O*[0,0,0,1]';
            
            Vx_Vector=Vx_Vector(1:3);
            Vx_Mag=norm(Vx_Vector);
            Vx_Vector=Vx_Vector/Vx_Mag;
            
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            % Vy:
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^              
            %Slab wrt Image Transformations:
            T_R=TT_R;  %Pose_R.Tmat; %reference pose:
            T_O=TT_Y1; %Pose_Y1.Tmat;
            T_X=TT_Y2; %Pose_Y2.Tmat;
            
            T_S_I_O=T_O;
            T_S_I_X=T_X;            
            T_X_O=T_S_I_X/T_S_I_O;
            
            Vy_Vector=T_X_O*[0,0,0,1]';
            
            Vy_Vector=Vy_Vector(1:3);
            Vy_Mag=norm(Vy_Vector);
            Vy_Vector=Vy_Vector/Vy_Mag;

            %=====================================================================
            Approx_Vector=Vy_Vector; % The Vertical vector to be input            

            %=====================================================================
            % Finding the WW rotation:
            T_O=TT_O;        %Pose_O.Tmat;
            T_X=TT_X;       %Pose_X.Tmat;
            T_Yaw_1=TT_Yaw_1; %Pose_Yaw_1.Tmat;
            T_Yaw_2=TT_Yaw_2; %Pose_Yaw_2.Tmat;
            T_Y1=TT_Y1;      %Pose_Y1.Tmat;                     

            T_S_I_O=T_O;
            T_S_I_X=T_X;
            T_S_I_Yaw_1=T_Yaw_1;
            T_S_I_Yaw_2=T_Yaw_2;            

            P1=T_S_I_X*[0,0,0,1]';
            P2=T_S_I_Yaw_1*[0,0,0,1]';
            P3=T_S_I_Yaw_2*[0,0,0,1]';
            
            P1=P1(1:3);
            P2=P2(1:3);
            P3=P3(1:3);
            
            %========================================================================
            % Projection of the points P1, P2, and P3 over the plane represented by
            % Vy_vector:
            % Using the equation:
            % U = V - (V dot N)N
            % where V is the position vector, and N is the normal vector
            Vy_Vector=Approx_Vector;
            N=Vy_Vector;
            V=P1;
            P1=V-N*(dot(V,N));
            V=P2;
            P2=V-N*(dot(V,N));
            V=P3;
            P3=V-N*(dot(V,N));
            %========================================================================
            
            %normal vector to where the rotation plane is:
            v1=P2-P1;
            v2=P3-P2;
            vn=cross(v1,v2);
            vn=vn/norm(vn);
            %mid-point between P1 and P2:
            P12m=(P1+P2)/2;
            P23m=(P3+P2)/2;
            P12v=cross(v1,vn);
            P23v=cross(v2,vn);
            pp1=P12m;
            pp2=pp1+P12v*100;
            pp3=P23m;
            pp4=pp3+P23v*100;
            % intersection the bisectors defines a point located on the pivot axis:
            XPnt=skewclosestpoint(pp1, pp2, pp3, pp4);
            % Angle between X and Yaw 2:
            ve1=P1-XPnt;
            ve2=P3-XPnt;
            ve1=ve1/norm(ve1);
            ve2=ve2/norm(ve2);
            Vyaw_Rot=acos(dot(ve1,ve2))*180/pi;
            Vyaw_Point=XPnt;
            Vyaw_Vector=vn;
            
            % Coordinate System defined with respect to the Image Space:
            Iorg=Vyaw_Point';
            % Convention is:
            %   x going upwards (along the length of the WW axis
            %   y going in-> in/out direction
            Ivx=Vx_Vector; % negative because Arm moves out in the first manuver so that direction needs to be reversed
            Ivy=Vy_Vector;
            vy=Ivy;  % following the save convention used for registering the axes in the TC-arm tracking module
            vz=cross(Ivx,vy);
            vx=cross(vy,vz);            
            
            Ix=vx'/norm(vx);
            Iy=vy'/norm(vy);
            Iz=vz'/norm(vz);
            P1=[0,0,0,1;10,0,0,1;0,20,0,1;0,0,30,1]';
            P2=[Iorg,1;Iorg+Ix*10,1;Iorg+Iy*20,1;Iorg+Iz*30,1]';
            T=P2/P1; % Carm Tracking Coordinates - to - Image Space Coordinates
            
            Obj.Trans=T;           
        end
              
        function ImageBasedCalibSixDOF_RefineXZ(Obj, CS, T_R, T_WW1, T_WW2, T_WW3, T_WW4)
            T={};
                      
            Pnt=[0,0,0,1;10,0,0,1;0,20,0,1;0,0,30,1]';
            CCPnts1=zeros(4,3);
            CCPnts2=zeros(4,3);
            CCPnts3=zeros(4,3);
            CCPnts4=zeros(4,3);
            CCPnts5=zeros(4,3);
            for f=1:4
                CPnts1=T_R/T_R*Pnt(:,f);
                CPnts2=T_WW1/T_R*Pnt(:,f);
                CPnts3=T_WW2/T_R*Pnt(:,f);
                CPnts4=T_WW3/T_R*Pnt(:,f);
                CPnts5=T_WW4/T_R*Pnt(:,f);            
                CCPnts1(f,:)=CPnts1(1:3)';
                CCPnts2(f,:)=CPnts2(1:3)';
                CCPnts3(f,:)=CPnts3(1:3)';
                CCPnts4(f,:)=CPnts4(1:3)';
                CCPnts5(f,:)=CPnts5(1:3)';
            end

            options = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',1000,'TolX',1e-6,'TolFun',1e-20);
            x0=[0,0];
            ref= lsqnonlin(@(x) refineXZObjFunc(CS,CCPnts1,CCPnts2,CCPnts3,CCPnts4,CCPnts5,x),x0,[],[],options);
            
            Pnts=[0,0,0,1;10,0,0,1;0,20,0,1;0,0,30,1]';
            P=CS*Pnts;
            PO=P(1:3,1)';
            PX=P(1:3,2)';
            PY=P(1:3,3)';
            PZ=P(1:3,4)';
            VX=PX-PO;
            VX=VX/norm(VX);
            VY=PY-PO;
            VY=VY/norm(VY);
            VZ=PZ-PO;
            VZ=VZ/norm(VZ);

            %centre in 3D
            Cent=PO+VX*ref(1)+VZ*ref(2);
            CS2=CS;
            CS2(1,4)=Cent(1);
            CS2(2,4)=Cent(2);
            CS2(3,4)=Cent(3);

            Obj.Trans=CS2;
        end
              
        function ImageBasedCalibSixDOF_RefineTilt(Obj, CS, T_R, T_WW1, T_WW2, T_WW3, T_WW4, T_Y1, T_Y2)
            % Acount for axis tilt:
                   
            % Two points on the WW1 and WW4
            T1=inv(T_R/T_R);
            T2=inv(T_WW4/T_R);
            T1_O=T1(1:3,4)';
            T2_O=T2(1:3,4)';
                       
            % Old Vx, Vy Vectors:
            Pnts=[0,0,0,1;10,0,0,1;0,20,0,1;0,0,30,1]';
            P=CS*Pnts;
            PO=P(1:3,1)';
            PX=P(1:3,2)';
            PY=P(1:3,3)';
            PZ=P(1:3,4)';
            VX=PX-PO;
            VX=VX/norm(VX);
            VY=PY-PO;
            VY=VY/norm(VY);
            VZ=PZ-PO;
            VZ=VZ/norm(VZ);

            % Calculate the Vy' orientation:
            Vec1=T1_O-PO;
            Vec2=T2_O-PO;
            Vy_Corr=cross (Vec1,Vec2);
            Vy_Corr=Vy_Corr/norm(Vy_Corr);
            
            % correct the Vy' sign:
            
            if dot(Vy_Corr,VY)<0
                Vy_Corr=-1*Vy_Corr;
            end
            
            % Recalculate the Vx, Vy;
            Vz_Corr=cross(VX,Vy_Corr);
            Vx_Corr=cross(Vy_Corr,Vz_Corr);
            

            P1=[0,0,0,1;10,0,0,1;0,20,0,1;0,0,30,1]';
            P2=[PO,1;PO+Vx_Corr*10,1;PO+Vy_Corr*20,1;PO+Vz_Corr*30,1]';
            T=P2*inv(P1); % Carm Tracking Coordinates - to - Image Space Coordinates
            
            Obj.Trans=T;            
                      
        end
        

        function V = getCalibratedValue(Obj)
            switch Obj.Type
                case 'Translation' %Obsolete
                    % V=(Obj.SensorRead-Obj.ResetValue)*Obj.ScaleFactor;
                case 'Rotation' %Obsolete
                    % V=(Obj.SensorRead-Obj.ResetValue)*Obj.ScaleFactor;
                case 'SixDOF'                    

                    V=Obj.getTmat;
            end
        end
                
        function T=getCorrectedTmat(Obj,varargin)
            if strcmpi(Obj.Type,'SixDOF')
                T=getTmat(Obj,varargin{1});
            else
                T=getTmat(Obj);
            end
        end
        
    end      
end

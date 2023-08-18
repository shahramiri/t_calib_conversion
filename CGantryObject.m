classdef CGantryObject < handle
    %VISOBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name
        LookupTable    
    end
    
    methods        
        function Obj=CGantryObject(BName)
            Obj.Name=BName;           
        end
               
        function save(Obj,FileName)            
            Name=Obj.Name;
            LookupTable=Obj.LookupTable;
            save(FileName,'Name','LookupTable','-mat');
        end
        
        function load(Obj,FileName)            
            load(FileName,'-mat');
            Obj.Name=Name;
            Obj.LookupTable=LookupTable;
        end
        

        %****************************************************************
        function Res=lookup(Obj,Tilti, Roti) %Interpolate and result calibration values

            AsmRot=Obj.LookupTable.Rot;
            AsmTilt=Obj.LookupTable.Tilt;
            AsmSource=Obj.LookupTable.Source;
            AsmPPoint=Obj.LookupTable.PPoint;
            AsmImageNormal=Obj.LookupTable.ImageNormal;
            AsmImageUp=Obj.LookupTable.ImageUp;
            AsmPrinDist=Obj.LookupTable.PrinDist;
            AsmXOffset=Obj.LookupTable.XOffset;
            AsmYOffset=Obj.LookupTable.YOffset;
            AsmPixelScale=Obj.LookupTable.PixelScale;
            
            %[Roti,Tilti,WWi,XXi,YYi]=readFullSensorDataFile(ObjectSensorDataFileName,PlotIt);                                                                                             % InterpM='linear';
            InterpM='pchip';
            % InterpM='cubic';
            % InterpM='spline';            
            MotionPath=[];
            % Interpolating the geometric parameters of the system based on Gantry
            % Angle tracked data:
            Ai=[];
            Ii=[];
            Xi=[];
            GOi=[];
            Ui=[];
            Vi=[];

            Sz=max(size(Roti));
            Xi=zeros(Sz,3);
            Ii=zeros(Sz,3);
            NVi=zeros(Sz,3);
            UVi=zeros(Sz,3);
            PDi=zeros(Sz,1);
            POXi=zeros(Sz,1);
            POYi=zeros(Sz,1);
            PSi=ones(Sz,1)*AsmPixelScale;

            for f=1:3
                xx=AsmSource(:,f);
                F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                XIin=[Roti;Tilti]';
                x=F(XIin);
                Xi(:,f)=x;
                %--------------------------------------------------------------------
                xx=AsmPPoint(:,f);
                F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                XIin=[Roti;Tilti]';
                x=F(XIin);
                Ii(:,f)=x;
                %--------------------------------------------------------------------
                xx=AsmImageNormal(:,f);
                F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                XIin=[Roti;Tilti]';
                x=F(XIin);
                NVi(:,f)=x;
                %--------------------------------------------------------------------
                xx=AsmImageUp(:,f);
                F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                XIin=[Roti;Tilti]';
                x=F(XIin);
                UVi(:,f)=x;
                %--------------------------------------------------------------------
                if f==1
                    xx=AsmPrinDist(:,f);       
                    F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                    XIin=[Roti;Tilti]';
                    x=F(XIin);
                    PDi(:,f)=x;
                    %--------------------------------------------------------------------
                    xx=AsmXOffset(:,f);       
                    F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                    XIin=[Roti;Tilti]';
                    x=F(XIin);
                    POXi(:,f)=x;
                    %--------------------------------------------------------------------
                    xx=AsmYOffset(:,f);       
                    F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
                    XIin=[Roti;Tilti]';
                    x=F(XIin);
                    POYi(:,f)=x;
                    %--------------------------------------------------------------------
                end
            end
                        
            PrincipalDist=PDi;
            Xoffset=POXi;
            Yoffset=POYi;
            PixelScale=PSi(1);
            SourcePosition=Xi;
            ImageNormal=NVi;
            ImageUp=UVi;
            CentrePoint=SourcePosition*0;
            
            % Calculate the CentrePoint:
            for g=1:max(size(PrincipalDist))
                j = ImageUp(g,:)/norm(ImageUp(g,:));
                k = ImageNormal(g,:)/norm(ImageNormal(g,:));
                i = cross(j,k);                
                % Convert the pixel points on the image into coordinates in 3D space.
                % We assume that the origin of the image coordinate system is at the
                % centre of the image                
                ctrpt = SourcePosition(g,:) + PrincipalDist(g,:)*-k - Xoffset(g,:)*i - Yoffset(g,:)*j;
                CentrePoint(g,:)=ctrpt;
            end            
            Res=[];
            Res.PrincipalDist=PrincipalDist;
            Res.Xoffset=Xoffset;
            Res.Yoffset=Yoffset;
            Res.PixelScale=PixelScale;
            Res.SourcePosition=SourcePosition;
            Res.ImageNormal=ImageNormal;
            Res.ImageUp=ImageUp;
            Res.CentrePoint=CentrePoint;                        
        end
        %****************************************************************
        
                %****************************************************************
        function Res=simpleLookup(Obj,Roti) %Interpolate and result calibration values

            AsmRot=Obj.LookupTable.Rot;
            AsmTilt=Obj.LookupTable.Tilt;
            AsmSource=Obj.LookupTable.Source;
            AsmPPoint=Obj.LookupTable.PPoint;
            AsmImageNormal=Obj.LookupTable.ImageNormal;
            AsmImageUp=Obj.LookupTable.ImageUp;
            AsmPrinDist=Obj.LookupTable.PrinDist;
            AsmXOffset=Obj.LookupTable.XOffset;
            AsmYOffset=Obj.LookupTable.YOffset;
            AsmPixelScale=Obj.LookupTable.PixelScale;
            
            %[Roti,Tilti,WWi,XXi,YYi]=readFullSensorDataFile(ObjectSensorDataFileName,PlotIt);                                                                                             % InterpM='linear';
            InterpM='pchip';
            % InterpM='cubic';
            % InterpM='spline';            
            MotionPath=[];
            % Interpolating the geometric parameters of the system based on Gantry
            % Angle tracked data:
            Ai=[];
            Ii=[];
            Xi=[];
            GOi=[];
            Ui=[];
            Vi=[];

            Sz=max(size(Roti));
            Xi=zeros(Sz,3);
            Ii=zeros(Sz,3);
            NVi=zeros(Sz,3);
            UVi=zeros(Sz,3);
            PDi=zeros(Sz,1);
            POXi=zeros(Sz,1);
            POYi=zeros(Sz,1);
            PSi=ones(Sz,1)*AsmPixelScale;

            for f=1:3
                xx=AsmSource(:,f);
%                 F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                 XIin=[Roti;Tilti]';
%                 x=F(XIin);
                x=interp1(AsmRot,xx,Roti,'pchip');
                Xi(:,f)=x;
                %--------------------------------------------------------------------
                xx=AsmPPoint(:,f);
%                 F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                 XIin=[Roti;Tilti]';
%                 x=F(XIin);
                x=interp1(AsmRot,xx,Roti,'pchip');
                Ii(:,f)=x;
                %--------------------------------------------------------------------
                xx=AsmImageNormal(:,f);
%                 F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                 XIin=[Roti;Tilti]';
%                 x=F(XIin);
                x=interp1(AsmRot,xx,Roti,'pchip');
                NVi(:,f)=x;
                %--------------------------------------------------------------------
                xx=AsmImageUp(:,f);
%                 F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                 XIin=[Roti;Tilti]';
%                 x=F(XIin);
                x=interp1(AsmRot,xx,Roti,'pchip');
                UVi(:,f)=x;
                %--------------------------------------------------------------------
                if f==1
                    xx=AsmPrinDist(:,f);       
%                     F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                     XIin=[Roti;Tilti]';
%                     x=F(XIin);
                    x=interp1(AsmRot,xx,Roti,'pchip');
                    PDi(:,f)=x;
                    %--------------------------------------------------------------------
                    xx=AsmXOffset(:,f);       
%                     F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                     XIin=[Roti;Tilti]';
%                     x=F(XIin);
                    x=interp1(AsmRot,xx,Roti,'pchip');
                    POXi(:,f)=x;
                    %--------------------------------------------------------------------
                    xx=AsmYOffset(:,f);       
%                     F=TriScatteredInterp(AsmRot,AsmTilt,xx,'natural');
%                     XIin=[Roti;Tilti]';
%                     x=F(XIin);
                    x=interp1(AsmRot,xx,Roti,'pchip');
                    POYi(:,f)=x;
                    %--------------------------------------------------------------------
                end
            end
                        
            PrincipalDist=PDi;
            Xoffset=POXi;
            Yoffset=POYi;
            PixelScale=PSi(1);
            SourcePosition=Xi;
            ImageNormal=NVi;
            ImageUp=UVi;
            CentrePoint=SourcePosition*0;
            
            % Calculate the CentrePoint:
            for g=1:max(size(PrincipalDist))
                j = ImageUp(g,:)/norm(ImageUp(g,:));
                k = ImageNormal(g,:)/norm(ImageNormal(g,:));
                i = cross(j,k);                
                % Convert the pixel points on the image into coordinates in 3D space.
                % We assume that the origin of the image coordinate system is at the
                % centre of the image                
                ctrpt = SourcePosition(g,:) + PrincipalDist(g,:)*-k - Xoffset(g,:)*i - Yoffset(g,:)*j;
                CentrePoint(g,:)=ctrpt;
            end            
            Res=[];
            Res.PrincipalDist=PrincipalDist;
            Res.Xoffset=Xoffset;
            Res.Yoffset=Yoffset;
            Res.PixelScale=PixelScale;
            Res.SourcePosition=SourcePosition;
            Res.ImageNormal=ImageNormal;
            Res.ImageUp=ImageUp;
            Res.CentrePoint=CentrePoint;                        
        end
        %
        
    end
    
end

% =========================================================================
% Set up the Data for loading calibration and creating an ImagePlane
% =========================================================================
Data=[];
Data.Cmd = 'TLoadCalibration';
Data.CmdID = '54';
Data.RegFile = 'M01_27-Jan-2022+15%3A02%3A47.mrkmat';
Data.CalibFolder = 'ATEC_OEC9900_Jan27_2022_01';

handles=[];
handles.tc=[];
handles.CalibsRootDataFolder='';
handles.SimExamDataFolder='';
handles.SaveOnDiskMode=0;
handles.W1=World('World!');
handles.W1.FrameGrabber=FrameGrabber(Data.CalibFolder,'FG1')
handles.W1.Carm=CarmObject;
handles.W1.Carm.World=handles.W1;
handles.W1.FrameGrabber.SimMode=1;
handles.W1.StitchMaster=StitchMaster('SM1');
handles.TestData=Exam('Exam_1','[Incomplete]')
handles.W1.Handles=handles;

GOrg=[0,0,0];
GVx=[1,0,0];
GVy=[0,1,0];
GVz=[0,0,1];
handles.W1.StitchMaster.Grid.PlaneOrigin=GOrg;
handles.W1.StitchMaster.Grid.VxLineElement.Vector=GVx;
handles.W1.StitchMaster.Grid.VyLineElement.Vector=GVy;
handles.W1.StitchMaster.Grid.VzLineElement.Vector=GVz;


FromWithinExam=0;
DoNOtWriteUsedCalibReg=1;
[handles,output_args] = TLoadCalibration(handles,Data, DoNOtWriteUsedCalibReg,FromWithinExam);

% =========================================================================
% Simulate Image Collection:
% =========================================================================
Data=[];
Data.Cmd = 'TCollectImageAndCoords';
Data.CmdID = '41';
Data.PlaneID = '1';
Data.TrackerData = '[0.300401,0.400098,0.500363,0.000000,0.000000,0.000000,1.000000]';

fprintf('Acquire Image..\n');
handles.W1.Carm.switchPlane(Data.PlaneID);
handles.W1.Carm.lookupCalib(Data.TrackerData);
CarmAngle=handles.W1.StitchMaster.returnCarmAngle(handles.W1.FrameGrabber.PlaneIndex); % display the physical C-arm angle:
[FG_thefile,FG_ErrorCode,FG_ErrorMessage]=handles.W1.FrameGrabber.framecap(0,1,CarmAngle, handles);

if FG_ErrorCode
    ImgPlane=[];
else
    [ImgPlane,~,~]=handles.W1.collectPosition(FG_thefile);
end

% figure;ImgPlane.imshow3d(1,1); axis equal;


% =========================================================================
% Camera Coordinate System:
% =========================================================================

Src=ImgPlane.SourcePosition;
Cnt=ImgPlane.CentrePoint;
Vz=ImgPlane.ImageNormal;
Vy=ImgPlane.ImageUp;
Vx=cross(Vy,Vz);
Vx=Vx/norm(Vx);
Vy=Vy/norm(Vy);
Vz=Vz/norm(Vz);

% Piercing Point:


% Org=Src;
Vorg=Cnt-Src;
Cx=dot(Vorg,Vx);
Cy=dot(Vorg,Vy);
Cz=dot(Vorg,Vz);

IP=ImagePlane('PlaneName');
IP.CentrePoint=[Cx, Cy, Cz]';
IP.SourcePosition =[0, 0, 0]';
IP.ImageNormal=[0,0,1]';
IP.ImageUp=[0,1,0]';
IP.PixelScale=ImgPlane.PixelScale;
IP.Image=ImgPlane.Image;

Pnt_UV=[0,0; 200,400; 512,512; 512,900];
cols={'r+','g+','b+','m+','y+'};
figure;imshow(IP.Image); hold on;

for f=1:size(Pnt_UV,1)
    plot(Pnt_UV(f,1),Pnt_UV(f,2),cols{f});
end
    
figure;IP.imshow3d(1,1);axis equal;

hold on;
CO=[0,0,0]';
Vx=[1,0,0]';
Vy=[0,-1,0]';
% Vy=[0,1,0]';
Vz=[0,0,-1]';
% Vz=[0,0,1]';

zdir=CO+(Vz*50);
xdir=CO+(Vx*50);
ydir=CO+(Vy*50);
plot3([CO(1),zdir(1)],[CO(2),zdir(2)],[CO(3),zdir(3)],'b');
plot3([CO(1),xdir(1)],[CO(2),xdir(2)],[CO(3),xdir(3)],'r');
plot3([CO(1),ydir(1)],[CO(2),ydir(2)],[CO(3),ydir(3)],'g');

PS=IP.PixelScale;
% Pnt_UV=[1024,512];
% Pnt_UV=[1024*0.75,512*0.5];

d_UV=[-512,-512];
% d_UV=[-512,512];

PrincipalDist= 978.6592;
Xoffset= -51.3001;
Yoffset= 63.9625;
PixelScale= 0.2109;

CntPoint=[-IP.Xoffset, -IP.Yoffset, -IP.PrincipalDist];


for f=1:size(Pnt_UV,1)
    % Pnt_XYZ=IP.CentrePoint + (Pnt_UV(f,1)+d_UV(1))*PS*Vx + (-Pnt_UV(f,2)-d_UV(2))*PS*Vy ;
    Pnt_X = -ImgPlane.Xoffset + (Pnt_UV(f,1)+d_UV(1))*PS;
    Pnt_Y = -ImgPlane.Yoffset + (Pnt_UV(f,2)+d_UV(2))*PS;
    Pnt_Z = -ImgPlane.PrincipalDist;
    Pnt_XYZ=[Pnt_X, Pnt_Y, Pnt_Z];
    plot3(Pnt_XYZ(1),Pnt_XYZ(2),Pnt_XYZ(3),cols{f});    
end


F=ImgPlane.PrincipalDist;
Ox=ImgPlane.Xoffset + 512*IP.PixelScale;
Oy=ImgPlane.Yoffset + 512*IP.PixelScale;
Mint=[F,0,Ox,0; 0,F,Oy,0; 0,0,1,0];
%figure;imshow(flipdim(IP.Image,1));

% Calculate the Extrinsic Matrix
Vz=[-ImgPlane.ImageNormal];
Vy=-ImgPlane.ImageUp;
Vx=cross(Vy,Vz);
Vy=cross(Vz,Vx);
Vx=Vx/norm(Vx);
Vy=Vy/norm(Vy);
Vz=Vz/norm(Vz);

Src=ImgPlane.SourcePosition;
Mext=[Vx',0; Vy',0; Vz',0; Src',1]';

PrjMat=Mint * inv(Mext);

% Construct sample points:
Pnts_3D=zeros(size(Pnt_UV,1),3);
Pnts_2D=zeros(size(Pnt_UV,1),2);

ty=ImgPlane.ImageUp;
tz=ImgPlane.ImageNormal;
tx=cross(ty,tz);
ty=cross(tz,tx);
tx=tx/norm(tx);
ty=ty/norm(ty);
tz=tz/norm(tz);

figure; ImgPlane.imshow3d(1,1); hold on; axis equal;
for f=1:size(Pnt_UV,1)
    Pnts_3D(f,:)= ImgPlane.CentrePoint + tx*(Pnt_UV(f,1)-512)*PS + ty*(Pnt_UV(f,2)-512)*PS;
    plot3(Pnts_3D(f,1),Pnts_3D(f,2),Pnts_3D(f,3),cols{f});
end

% Calculate the projection from the 3D points:
figure; imshow(ImgPlane.Image); hold on;
for f=1:size(Pnt_UV,1)
    Prj = PrjMat * [Pnts_3D(f,:),1]';
    Prj=Prj/Prj(3);
    Pnts_2D(f,:) = Prj(1:2)';
    plot(Pnts_2D(f,1),Pnts_2D(f,2),cols{f});
end

Pnts_2D
Pnt_UV
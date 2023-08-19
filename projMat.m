function PrjMat=projMat(ImgPlane)

% Take an ImagePlane and calcvulate the projection matrix:
PS=ImgPlane.PixelScale;
MidPoint=size(ImgPlane.Image,1)/2;
F=ImgPlane.PrincipalDist/PS;
Ox=ImgPlane.Xoffset/PS + MidPoint;
Oy=ImgPlane.Yoffset/PS + MidPoint;
Mint=[F, 0, Ox, 0; 0, -F, Oy, 0; 0, 0, 1, 0];

% Calculate the Extrinsic Matrix
Vz=-ImgPlane.ImageNormal;
Vy=-ImgPlane.ImageUp;
Vx=cross(Vy,Vz);
Vy=cross(Vz,Vx);
Vx=Vx/norm(Vx);
Vy=Vy/norm(Vy);
Vz=Vz/norm(Vz);

Src=ImgPlane.SourcePosition;
Mext=inv([Vx',0; Vy',0; Vz',0; Src',1]');

% Projection matrix:
PrjMat=Mint * Mext;
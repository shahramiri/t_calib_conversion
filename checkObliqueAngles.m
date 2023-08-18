function [ObliqueAngles]=checkObliqueAngles(CalibFolderName)
ObliqueAngles=[]; 
try
    FileName=dir(sprintf('%s\\*.CFG',CalibFolderName));
    % Read the cfg file and extract the carm name and description
    fid=fopen(sprintf('%s\\%s',CalibFolderName,FileName.name),'r');
    % skip the first two lines
    fgetl(fid);
    fgetl(fid);
    % skip three lines:
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    % read the tracking mountblockes and markers:
    if ~feof(fid)
        fgetl(fid); %skip the line of markers
    end
    
    % read the optional Oblique angles
    if ~feof(fid)
       OblStrCell=eval(fgetl(fid)); 
    end
    
    if (size(OblStrCell,1)==1) && (size(OblStrCell,2)==2)        
        if strcmpi(OblStrCell{1,1},'ObliqueAngles')
            ObliqueAngles=OblStrCell{1,2};
        end
    end
    
catch
   ObliqueAngles=[]; 
end
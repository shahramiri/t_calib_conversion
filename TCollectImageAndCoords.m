function [handles,output_args] = TCollectImageAndCoords(handles,Data)
% This command handles collecting image and adding it to the composite
% views

output_args=0;
fprintf('Acquire Image..\n');

% 1. DAQ> Update C-arm tracking based on raw TrackerData
handles.W1.Carm.switchPlane(Data.PlaneID);
handles.W1.Carm.lookupCalib(Data.TrackerData);

% 2. StitchMaster> Return Physical C-arm angle (PlaneID)
CarmAngle=handles.W1.StitchMaster.returnCarmAngle(handles.W1.FrameGrabber.PlaneIndex); % display the physical C-arm angle:

% 3. DAQ> Frame Capture
[FG_thefile,FG_ErrorCode,FG_ErrorMessage]=handles.W1.FrameGrabber.framecap(0,0,CarmAngle, handles);

% 4. StitchMaster> b. Add image to the corresponding view
if FG_ErrorCode
    ImgPlane=[];
else
    [ImgPlane,~,~]=handles.W1.collectPosition(FG_thefile);
end
sendCmdProgress(handles,0.25); %25 percent progress;
if (FG_ErrorCode==0) % if no error in getting an image from the gramegrabber then continue:
    
    % Add image to the stitcher:
    switch handles.W1.StitchMaster.CurrentPlane
        case 'AP'
            CP=1;
        case 'ML'
            CP=2;
    end
    handles.W1.StitchMaster.ShowCropped=1; % Always show the cropped view if possible

    % temprary file move
    movefile('Tcarmcalib\LastShot Data\Jointtrack Data\*.*', handles.TestData.ImageFolder,'f');

    ImagePlaneIndex=numel(handles.TestData.ImagePlanes)+1; 
    
    % addImage
    [AddedImage,AddedImageCoords,AddedImageSrcIndx,AddedImageDepthIndx,ErrCode]=handles.W1.StitchMaster.addImage(CP,handles.W1.StitchMaster.alignImagePlane(ImgPlane),0,[],0,[],ImagePlaneIndex);
    
    sendCmdProgress(handles,0.5); %50 percent progress;
    if ErrCode
        output_args=[];
        output_args.HasAttachment=0;         
        output_args.ErrorCode=2;
        output_args.ErrorMessage='Image Out of Bound! Reposition the Calibration Grid and Try again.';       
     
    else
        
        % Only add to the ImagePlanes if no error is returned:
        if numel(handles.TestData.ImagePlanes)==0
            handles.TestData.ImagePlanes=ImgPlane;
        else
            handles.TestData.ImagePlanes(end+1)=ImgPlane;
        end
        
        SrcPosition3D=handles.W1.StitchMaster.alignImagePlane(ImgPlane).SourcePosition';        
        CentrePoint3D=handles.W1.StitchMaster.alignImagePlane(ImgPlane).CentrePoint';
        ImageNormal3D=handles.W1.StitchMaster.alignImagePlane(ImgPlane).ImageNormal';

        if ~isempty(AddedImage)
            %Generate the filename:
            temp_name = tempname;
            [~, FName] = fileparts(temp_name);
            FileName=['SendingAttachments\',FName,'.png'];
            MemoryFileName=FName;

            % imwrite(AddedImage,FileName,'png');
            fprintf(['Added Image Portion is ready for retreive: ',FName,'\n']);
            handles.W1.StitchMaster.IsVoxelBasedCrop=0; % Set the cropping mode back to normal
            handles.W1.StitchMaster.cropViews(1);

            if handles.SaveOnDiskMode
                %imwrite(cat(3,AddedImage,AddedImageSrcIndx,AddedImageSrcIndx*0),FileName,'png');
                savepng(cat(3,AddedImage,AddedImageSrcIndx,AddedImageDepthIndx),FileName,1);
                ByteSize=checkFileBytes(FileName);
            else
                handles.Memory.MemoryFileName=FName;    
                %handles.Memory.Files.(FName)=img2PngByteArray(AddedImage); % ByteArray wo need for writting to the disk
                %imwrite(cat(3,AddedImage,AddedImageSrcIndx,AddedImageSrcIndx*0),FileName,'png');
                savepng(cat(3,AddedImage,AddedImageSrcIndx,AddedImageDepthIndx),FileName,1);
                fid=fopen(FileName);
                BytesRead=fread(fid);
                fclose(fid);
                delete(FileName); % delete the temporary;
                handles.Memory.Files.(FName)=BytesRead;
                ByteSize=numel(handles.Memory.Files.(FName));
            end

            sendCmdProgress(handles,0.7); %70 percent progress;
            
            %Corrected Coordinate according to Unity's convention:
            AddedImageCoords_Corrected=AddedImageCoords;
            switch CP
                case 1 %AP Plane:
                    AddedImageCoords_Corrected(1)=size(handles.W1.StitchMaster.APPlane.Image,1)-size(AddedImage,1)-AddedImageCoords(1);       
                case 2 %ML Plane
                    AddedImageCoords_Corrected(1)=size(handles.W1.StitchMaster.MLPlane.Image,1)-AddedImageCoords(1)-size(AddedImage,1); 
            end

            % 5. StitchMaster> c. crop views
            [handles,output_args_Crops] = TGetBothCroppedCoords(handles);
            
            sendCmdProgress(handles,0.8); %80 percent progress;

            output_args=[];
            output_args.HasAttachment=1;
            output_args.PlaneID=Data.PlaneID;            
            if handles.SaveOnDiskMode
                output_args.AttachmentFileName=FileName;
            else
                output_args.AttachmentFileName=MemoryFileName;
            end
            output_args.AttachmentType='Image';
            output_args.AttachmentSize=size(AddedImage);
            output_args.ImageCoordinates=AddedImageCoords_Corrected;
            output_args.ByteSize=ByteSize;

            % Replacement of outputs fed by the new faster function:
            output_args.CropAreaSize_APView=output_args_Crops.Imagesizes(1,:);
            output_args.CropAreaCoordinates_APView=output_args_Crops.ImageCoordinates(1,:);
            output_args.CropAreaSize_MLView=output_args_Crops.Imagesizes(2,:);
            output_args.CropAreaCoordinates_MLView=output_args_Crops.ImageCoordinates(2,:);

            % Added Field 2018 Oct 29; Source Positions 3D:
            output_args.SrcPosition3D=SrcPosition3D;
            
            % Added Field 2020 Jun 29; CentrePoint 3D:
            output_args.CentrePoint3D=CentrePoint3D;
            output_args.ImageNormal3D=ImageNormal3D;
            


        else % if out of bound or errors in collecting the images:

            temp_name = tempname;
            [~, FName] = fileparts(temp_name);
            FileName=['SendingAttachments\',FName,'.png'];  
            MemoryFileName=FName;
            AddedImage=zeros(2,2);

            if handles.SaveOnDiskMode        
                %imwrite(AddedImage,FileName,'png');
                savepng(AddedImage,FileName,1);
                % ByteSize=checkFileBytes(FileName);
            else
                handles.Memory.MemoryFileName=FName;
                %handles.Memory.Files.(FName)=img2PngByteArray(AddedImage); % ByteArray wo need for writting to the disk
                %imwrite(AddedImage,FileName,'png');
                savepng(AddedImage,FileName,1);
                fid=fopen(FileName);
                BytesRead=fread(fid);
                fclose(fid);
                delete(FileName); % delete the temporary;
                handles.Memory.Files.(FName)=BytesRead;
                % ByteSize=numel(handles.Memory.Files.(FName));
            end 
            
            sendCmdProgress(handles,0.9); %90 percent progress;
            
            %imwrite(AddedImage,FileName,'png');
            fprintf(['Added Image Portion is ready for retreive: ',FName,'\n']);        
            output_args=[];
            output_args.HasAttachment=1;
            output_args.PlaneID=Data.PlaneID;
            % output_args.AttachmentFileName=FileName;
            if handles.SaveOnDiskMode
                output_args.AttachmentFileName=FileName;
            else
                output_args.AttachmentFileName=MemoryFileName;
            end    
            output_args.AttachmentType='Empty';
            output_args.AttachmentSize=size(AddedImage);
            output_args.ImageCoordinates=[1,1];

            % Replacement of outputs fed by the new faster function:
            output_args.CropAreaSize_APView=[1,1];
            output_args.CropAreaCoordinates_APView=[1,1];
            output_args.CropAreaSize_MLView=[1,1];
            output_args.CropAreaCoordinates_MLView=[1,1];
            output_args.SrcPosition3D=0;
            output_args.CentrePoint3D=0;
            output_args.ImageNormal3D=0;
            fprintf('Undoing Last Collected Image - No UndoImage is found. Null Image Created\n');
        end

    end
else
        output_args=[];
        output_args.HasAttachment=0;                        
        output_args.ErrorCode=1;
        output_args.ErrorMessage='Error Receiving Image from the FrameGrabber Card.';        
end

% 7. DAQ> c. Reset AutoCollect Image

%% Supporting the AutoCollectImage mode:
if handles.W1.FrameGrabber.AutoCheck
    if handles.W1.FrameGrabber.SimMode
        handles.W1.FrameGrabber.resetAuto;    % If in simulation mode delete the sim file throough reset
    else
        handles.W1.FrameGrabber.checkAuto(1,handles); % Otherwise reset the Reference Frame of the frameGrabber
    end
end
    
end

function sendCmdProgress(handles, progressPerc)

if ~isempty(handles.tc)
    if strcmpi(handles.tc.status,'open')
        handles.CmdP.CmdProgress=progressPerc;
%         display(['>> Send update!! ', num2str(progressPerc)]);
        clientBytesReceived([],[],handles);
    end
end

end

classdef UserGroup < handle
    %BODY class to define bodies such as instruments, jig, patient, C-arm,
    %etc
    %   Detailed explanation goes here
    
    properties
        UserDatabase
        DataFolder='';
        RawDataFileName='Users\UserDatabase.txt';
        DataFileName='Users\UserDatabase.dat';
        EncryptKey = 'Measure2Perfect!';
        CurrentUser=[];
    end
    
    
    methods
        function Obj=UserGroup()
            Obj.UserDatabase=[];
            % Update the DataFileName paths:
        end
        
        function setFolder(Obj,FolderName)
            Obj.DataFolder=FolderName;
        end
        
        % https://github.com/dbossnirvana/AES_MATLAB
        function encryptFile(Obj)            
            try
            UserDB=[];
            UserDB.UserInfo=[];
            UserDB.PCName=[];
            preallocations;
            key = zerofill(Obj.EncryptKey);
            % Key Schedule
            round_keys = key_schedule(double(key));
            UserDB.PCName=aes_encryption(zerofill(getenv('COMPUTERNAME')),round_keys);
            fid=fopen(Obj.RawDataFileName);
            while (~feof(fid))
                %  [AA]=fscanf(fid,'%s;%s;%s;%s');
                AA=fgetl(fid);
                if ~(strcmpi(AA(1),'%')) %Skip comments
                    A=strsplit(AA,';','CollapseDelimiters',false);
                    fprintf('\n\r User: "%s" >',A{1});
                    Nametext = A{1}; %aes_encryption(zerofill(A{1}),round_keys);
                    Passtext = aes_encryption(zerofill(A{2}),round_keys);
                    Typetext = aes_encryption(zerofill(A{3}),round_keys);
                    ReportPasstext = aes_encryption(zerofill(A{4}),round_keys);
                    if isempty(UserDB.UserInfo)
                        UserDB.UserInfo.UserInfo.Name=Nametext;
                        UserDB.UserInfo.Pass=Passtext;
                        UserDB.UserInfo.Type=Typetext;
                        UserDB.UserInfo.ReportPass=ReportPasstext;
                    else
                        N=numel(UserDB.UserInfo);
                        UserDB.UserInfo(N+1).Name=Nametext;
                        UserDB.UserInfo(N+1).Pass=Passtext;
                        UserDB.UserInfo(N+1).Type=Typetext;
                        UserDB.UserInfo(N+1).ReportPass=ReportPasstext;
                    end
                end
            end
            fclose(fid);
            save([Obj.DataFolder, Obj.DataFileName],'UserDB','-mat');
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in Reading the User Database');
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end
        end
        
        function load(Obj) %update the QItems property based on the available Items            
            try
                load([Obj.DataFolder,Obj.DataFileName],'UserDB','-mat');
                Obj.UserDatabase=[];                                                
                Obj.UserDatabase=UserDB;
            catch e %e is an MException struct
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in Reading the User Database .resetAuto()');
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end            
        end
        
        function [Access,UserFound,UserIndex]=checkUser(Obj,User,Pass)
            try
                Access='0';
                UserFound='0';
                UserIndex=[];
                preallocations;
                key = zerofill(Obj.EncryptKey);
                round_keys = key_schedule(double(key));
                Dec_PCName = aes_decryption(Obj.UserDatabase.PCName, round_keys);
                if strcmpi(Dec_PCName,getenv('COMPUTERNAME'))
                    for f=1:numel(Obj.UserDatabase.UserInfo)
                        
                        if isempty(User) && isempty(Obj.UserDatabase.UserInfo(f).Name)
                            UserChecked=1;
                        else
                            UserChecked=strcmpi(User,Obj.UserDatabase.UserInfo(f).Name);
                        end
                        if (UserChecked)
                            UserFound='1';
                            UserIndex=f;
                            Dec_Pass = aes_decryption(Obj.UserDatabase.UserInfo(f).Pass, round_keys);
                            if strcmp(Pass,Dec_Pass)
                                Access='1';
                            end
                        end
                    end
                end
            catch e
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in Checking the user credentials');
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');
            end
        end
        
        function UserType=checkUserType(Obj,UserIndex)
            try
                UserType=[];
                key = zerofill(Obj.EncryptKey);
                round_keys = key_schedule(double(key));
                UserType = aes_decryption(Obj.UserDatabase.UserInfo(UserIndex).Type, round_keys);
            catch
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in Checking the usertype');
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');                
            end
        end
        
        function ReportPass=checkReportPass(Obj,UserIndex)
            try
                ReportPass=[];
                key = zerofill(Obj.EncryptKey);
                round_keys = key_schedule(double(key));
                ReportPass = aes_decryption(Obj.UserDatabase.UserInfo(UserIndex).ReportPass, round_keys);
            catch
                fprintf('\r\n---------------------------------------');
                fprintf('\r\nError in Checking the usertype');
                fprintf('\r\nThe identifier was:\n%s',e.identifier);
                fprintf('\r\nThere was an error! The message was:\n%s',e.message);
                fprintf('\r\n---------------------------------------');                
            end
        end      
        
        
        function setCurrentUser(Obj,Name,Type,ReportPass)
            Obj.CurrentUser=[];
            Obj.CurrentUser.Name=Name;
            Obj.CurrentUser.Type=Type;
            Obj.CurrentUser.ReportPass=ReportPass;
        end
        
        function resetCurrentUser(Obj)
            Obj.CurrentUser=[];
            Obj.CurrentUser.Name='';
            Obj.CurrentUser.Type='';
            Obj.CurrentUser.ReportPass='';
        end
    end
end
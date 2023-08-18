classdef List < handle
    %BODY class to define bodies such as instruments, jig, patient, C-arm,
    %etc
    %   Detailed explanation goes here
    
    properties
        Name
        Items
        QItems
    end
    
    
    properties (Dependent)
        Count
    end
        
    methods
        function Obj=List(LName)
            Obj.Name=LName;
            Obj.Items=[];
            Obj.QItems=[];
        end
        
        function updateQItems(Obj) %update the QItems property based on the available Items
            Obj.QItems=[];
            for f=1:Obj.Count
                Obj.QItems.(genvarname(Obj.Items(f).Name))=Obj.Items(f);
            end
        end
        
        function ListCount=get.Count(Obj)
            ListCount=size(Obj.Items,1);
        end
        
        function addItem(Obj,Item)
            Obj.Items=[Obj.Items;Item];
            Obj.updateQItems;
        end
        
        function emptyList(Obj,Item)
            Obj.Items=[];
            Obj.QItems=[];
        end
        
        function Item=getByName(Obj,Name) %new version is for faster processing
            Item=[];
            if isfield(Obj.QItems,genvarname(Name))
                Item=Obj.QItems.(genvarname(Name));
            end
        end
        
        function Item=getByName_Old(Obj,Name) %the old version uses a loop and so it is slow
            Item=[];
            if Obj.Count>0
                for f=1:Obj.Count
                    if strcmp(Obj.Items(f).Name,Name)
                        Item=[Item;Obj.Items(f)];
                    end
                end
            end
        end

        function removeByName(Obj,Name)
            cc=Obj.Count;
            cnt=0;
            i=0;
            if cc>0
                while (i<cc)
                    i=i+1;
                    if strcmp(Obj.Items(i-cnt).Name,Name)
                        Obj.Items(i-cnt)=[];
                        cnt=cnt+1;
                    end
                end
            end
            Obj.updateQItems;
        end
                  
    end
    
end
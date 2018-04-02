% Written by Zhen Liu. 



clear
close all
clc
%--------------------------------------------------------------------------
% Extract bacteria information under the same folder
%--------------------------------------------------------------------------
dataPaths=uigetdir;
cd(dataPaths);
files_bacteria=dir([dataPaths '\*.mat']);
cell_temp=cell(1,1);
for i=1:length(files_bacteria)
    data=importdata([dataPaths '\' files_bacteria(i).name]);   
    cell_small=data.cellList{1,1};
    cell_temp=[cell_temp,cell_small];
end
emptycells=find(~cell2mat(cellfun(@length, cell_temp, 'UniformOutput',false)));
if ~isempty(emptycells)
    cell_temp(:,emptycells)=[];
end
%--------------------------------------------------------------------------
% Extracted the fluorescence signals and plot the signal profile 
%--------------------------------------------------------------------------
x_all=[];
y_all=[];
figure
interest_all=[];
temp=cell(length(cell_temp),3);  
calculated=[];
exception=[];
cell_length=[];
for  i=1:length(cell_temp)
   try
        x_temp=cell_temp{1,i};
        x1=x_temp.steplength;
        x2=cell_temp{1,i}.length;  % x2 as the output
        cell_length=[cell_length,x2];
        %x3=x1/x2; % comment it 
        %x=cumsum(x3);
        x=cumsum(x1); % x1
        y=x_temp.signal1;
        %interest(1,1)=i; % comment it 
        k=1;
        for j=2:(length(y)-1)
            if y(j)>y(j-1)&&y(j)>y(j+1) %identify the peak of signal profile
               temp{i,k}=x(j);  % without normalization     
               k=k+1; 
            end
        end
             
        %y=y/max(y);
        %x_all=[x_all,x];
        %y_all=[y_all,y];
        hold on
        %plot(x,y,'-')
        clear x_temp x1 x2 x3 x y 
        calculated=[calculated,i];     
    catch
        exception=[exception,i];
   end
    
end
hold off

% Revised by yongliang%%

%find and eliminate the cells in which only one smMapZ ring was identified%%
%
emptycells2=find(~cell2mat(cellfun(@length, temp(:,2), 'UniformOutput',false)))
temp(emptycells2,:)=[]
cell_length(:,emptycells2)=[]
%%%---------------------------------------------------



%%------------------------------------------------------------%%
%calcaulate the peak_distance between two outer smMapZ rings
%-----------------------------------------------------------
peak_distance=zeros(1,length(temp));

for i=1:length(temp)
    if ~isempty(temp{i,4})
        peak_distance(1,i)=temp{i,4}-temp{i,1};
    elseif ~isempty(temp{i,3})
        peak_distance(1,i)=temp{i,3}-temp{i,1};
    elseif ~isempty(temp{i,2})
        peak_distance(1,i)=temp{i,2}-temp{i,1};
    else
        peak_distance(1,i)=0;
    end
end
hold on
figure
plot(cell_length,peak_distance,'*')

% -------------------------------
%calcaulate the peak_to_pole
%-------------------------------------

peak_to_pole=zeros(1,length(temp));

for i=1:length(temp)
    if ~isempty(temp{i,4})
        peak_to_pole(1,i)=min((cell_length(1,i)-temp{i,4}),temp{i,1});
    elseif ~isempty(temp{i,3})
        peak_to_pole(1,i)=min((cell_length(1,i)-temp{i,3}),temp{i,1});
    elseif ~isempty(temp{i,2})
        peak_to_pole(1,i)=min((cell_length(1,i)-temp{i,2}),temp{i,1});
    else
        peak_to_pole(1,i)=min((cell_length(1,i)-temp{i,1}),temp{i,1});
    end
end

figure
plot(cell_length,peak_to_pole,'*')
%---------------------------------------------





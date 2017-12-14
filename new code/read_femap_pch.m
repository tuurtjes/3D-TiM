%<manager>
%<version>1.0.0.0</version>
%<date>17/06/2014</date>
%<type>Matlab</type>
%</manager>

%function [NASTRAN_DATA]=read_femap_pch()
clear all
% NMa, 17/06/2014
% READS STIFFNESS AND MASS MATRIX FROM NASTRAN PUNCH FILE
% USE IN FEMAP FOLLOWING COMMAND (BEFORE BULK)
% EXTSEOUT(EXTID=10,STIFFNESS,MASS,DMIGPCH)
% AND PRESCRIBE REQUIRED NODES AND DOFS VIA ASET
%
% Returns structure NASTRAN_DATA with 
% Kred=reduced stiffness matrix
% Mred=reduced mass matrix
% dof= overview with nodes and dofs
% grid= coordinates of used nodes [node nr, coord sys, x,y,z]
%
%EXAMPLE NASTRAN FILE
% INIT MASTER(S)
% ASSIGN,OUTPUT4='test.op4',UNIT=36
% NASTRAN SYSTEM(442)=-1,SYSTEM(319)=1
% ID husky,Femap
% SOL SEMODES
% TIME 10000
% CEND
%   TITLE = modal
%   ECHO = NONE
%   DISPLACEMENT(PLOT) = ALL
%   SPCFORCE(PLOT) = ALL
%   ESE(PLOT) = ALL
%   EXTSEOUT(EXTID=10,STIFFNESS,MASS,DMIGPCH)
%   METHOD = 1
%   SPC = 2
% BEGIN BULK
% $ ***************************************************************************
% $   Written by : Femap with NX Nastran
% $   Version    : 11.0.0
% $   Translator : NX Nastran
% $   From Model : C:\Users\nma\Desktop\New folder\husky_01.modfem
% $   Date       : Tue Jun 17 17:15:40 2014
% $ ***************************************************************************
% $
% PARAM,POST,-1
% PARAM,OGEOM,NO
% PARAM,AUTOSPC,YES
% PARAM,GRDPNT,0
% EIGRL          1                      10       0                    MASS
% CORD2C         1       0      0.      0.      0.      0.      0.      1.+FEMAPC1
% +FEMAPC1      1.      0.      1.
% CORD2S         2       0      0.      0.      0.      0.      0.      1.+FEMAPC2
% +FEMAPC2      1.      0.      1.
% $ Femap with NX Nastran Constraint Set 2 : pinned
% ASET        5933     123
% ASET        5937     123
% ASET       12120     123
% ASET       12121     123
%


[file1, path1] = uigetfile('*.pch','Select *.pch file');
[~, fname,~] = fileparts(file1);
[file2, path2] = uiputfile('*.mat', 'Save NASTRAN_DATA as', sprintf('%s%s.mat', path1, fname));

fid=fopen(strcat(path1,file1));
tline=fgetl(fid);

% super element id
tline=fgetl(fid);
superid = str2num(tline(63:end));
Kname = sprintf('K%d',superid);
Mname = sprintf('M%d',superid);

while isempty(findstr(tline,'GRID*'))
    tline=fgetl(fid);
end
grid=[];
i=0;
%read node coordinates
while tline(1)~='$'
    i=i+1;
    nr=str2num(tline(6:25));
    cs=str2num(tline(26:41));
    x=str2num(tline(41:56));
    y=str2num(tline(57:end));
    tline=fgetl(fid);
    z=str2num(tline(2:25));
    %grid(i,:)=[nr cs x y z]
    grid(i,:)=[nr x y z];
    tline=fgetl(fid);
end

% find all dofs
while isempty(findstr(tline,'KAAX')) && isempty(findstr(tline,Kname))
    tline=fgetl(fid);
end
%start of KAAX found
dof=[];
tline=fgetl(fid);
i = 0;
while isempty(findstr(tline,'MAAX')) && isempty(findstr(tline,Mname))
    i=i+1; %row
    dof(i,:)=str2num(tline(15:end));
    tline=fgetl(fid);
    while tline(1)=='*'
        tline=fgetl(fid);
    end
end

% re-open
fclose(fid);
fid=fopen(strcat(path1,file1));

while isempty(findstr(tline,'KAAX')) && isempty(findstr(tline,Kname))
    tline=fgetl(fid);
end
%start of KAAX found
i=0;
Kred=[];
tline=fgetl(fid);
while isempty(findstr(tline,'MAAX')) && isempty(findstr(tline,Mname))
    i=i+1; %row
    j=0;
    dof1 = str2num(tline(15:end));
    ind1 = find(dof(:,1)==dof1(1) & dof(:,2)==dof1(2));
    tline=fgetl(fid);
    while tline(1)=='*'
        j=j+1;
        dof2 = str2num(tline(17:40));
        ind2 = find(dof(:,1)==dof2(1) & dof(:,2)==dof2(2));
        Kred(ind1,ind2)=str2num(tline(41:end));
        Kred(ind2,ind1)=str2num(tline(41:end));
        tline=fgetl(fid);
    end
end

fprintf(1,'Stiffness matrix done\n')

%Read mass matrix
i=0;
Mred=[];
tline=fgetl(fid);
while isempty(findstr(tline,'VAX'))
    i=i+1; %row
    j=0;
    dof1 = str2num(tline(15:end));
    ind1 = find(dof(:,1)==dof1(1) & dof(:,2)==dof1(2));
    tline=fgetl(fid);
    while tline(1)=='*'
        j=j+1;
        dof2 = str2num(tline(17:40));
        ind2 = find(dof(:,1)==dof2(1) & dof(:,2)==dof2(2));
        Mred(ind1,ind2)=str2num(tline(41:end));
        Mred(ind2,ind1)=str2num(tline(41:end));
        tline=fgetl(fid);
    end
end

fprintf(1,'Mass matrix done\n')

%Store data in struc
NASTRAN_DATA.Mred=Mred;
NASTRAN_DATA.Kred=Kred;
NASTRAN_DATA.dof=dof;
NASTRAN_DATA.grid=grid;
NASTRAN_DATA.file=strcat(path1,file1);

clear Mred Kred dof grid

fclose(fid)

if ~isnumeric(file2)
    save([path2, file2], 'NASTRAN_DATA');
end

disp('Done!')

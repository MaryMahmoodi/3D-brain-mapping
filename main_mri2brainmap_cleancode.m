close all
clear 
clc
tic,
%% initial parameters
method= 'cubic';           % cubic % RBF % linear % avg 
fprintf('your electrode labels should be aligned with uppercase or lower case no difference....  ');
sprintf(' Cz,Cpz,Pz,Poz,Oz,Iz,Fpz,AFz,Fz,Fcz,O2,Po4,P2,Cp2,P4,Po8,P10,T6,P6,Tp8,Cp6,Cp4,C2,C4,C6,T4,Ft8,Fc6,Fc4,F4,F6,F8,AF8,Fc2,F2,AF4,Fp2,P1,P3,Po7,O1,Po3,Fc1,F1,AF3,Fp1,AF7,F3,F5,F7,Ft7,Fc5,Fc3,C1,Cp1,C3,C5,T3,Tp7,Cp5,Cp3,P5,T5,P9,nasion,lpa,rpa')

Channels= {'Cz' 'C4' 'T4' 'C3' 'T3' 'Fz' 'F4' 'F8' 'F3' 'F7' 'Fpz' 'Fp2' 'Fp1' 'Pz' 'P3' 'T5' 'P4' 'T6' 'Oz' 'O1' 'O2'}';
Values=[0.42  0.75  0.85  0.45  0.43  0.4  0.43  0.55  0.35  0.37  0.25  0.39  0.35  0.78  0.6  0.6  1  0.95  0.8  0.6  0.87]';


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Channels={'Cz'; 'C4'; 'T4'; 'C3'; 'T3'; 'Fz'; 'F4'; 'F8'; 'F3'; 'F7'; 'Fpz'; 'Fp2'; 'Fp1'; 'Pz'; 'P3'; 'T5'; 'P4'; 'T6'; 'Oz'; 'O1'; 'O2'};
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Values=[0.7,0.9,0.95,0.55,0.14,0.15,0.26,0.84,0.25,0.8,0.24,0.93,0.35,0.2,0.25,0.6,0.47,0.35,0.8,0.59,0.55];

%  Channels={'Cz';'Cpz';'Pz';'Poz';'Oz';'Iz';'Fpz';'AFz';'Fz';'Fcz';'O2';'Po4';'P2';'Cp2';'P4';'Po8';'P10';'T6';'P6';'Tp8';'Cp6';'Cp4';'C2';'C4';'C6';'T4';'Ft8';'Fc6';'Fc4';'F4';'F6';'F8';'AF8';'Fc2';'F2';'AF4';'Fp2';'P1';'P3';'Po7';'O1';'Po3';'Fc1';'F1';'AF3';'Fp1';'AF7';'F3';'F5';'F7';'Ft7';'Fc5';'Fc3';'C1';'Cp1';'C3';'C5';'T3';'Tp7';'Cp5';'Cp3';'P5';'T5';'P9';'nasion';'lpa';'rpa'};
%  Values=[0.5 0.3 1.0 0.5 0.7 0.1 1.0 0.8 0.5 0.3 1 1 0.7 0.1 1.0 0.0 0.1 0.3 0.0 0.2 1 0.5 0.3 1.0 0.5 0.7 0.1 1.0 0.8 0.5 0.3 1 1 0.7 0.1 1.0 0.0 0.1 0.3 0.0 0.2 1 0.5 0.3 1.0 0.5 0.7 0.1 1.0 0.8 0.5 0.3 1 1 0.7 0.1 1.0 0.0 0.1 0.3 0.0 0.2 1 1 1 1 1];

%  Values=randn(length(Channels),1);
% % % % % % % % % %    Values=[0.7,0.9,0.95,0.55,0.14,0.15,0.26,0.84,0.25,0.8,0.24,0.93,0.35,0.2,0.25,0.6,0.47,0.35,0.8,0.59,0.55 0.7,0.9,0.95,0.55,0.14,0.15,0.26,0.84,0.25,0.8,0.24,0.93,0.35,0.2,0.25,0.6,0.47,0.35,0.8,0.59,0.55 0.7,0.9,0.95,0.55,0.14,0.15,0.26,0.84,0.25,0.8,0.24,0.93,0.35,0.2,0.25,0.6,0.47,0.35,0.8,0.59,0.55 ,0.55,0.59,0.55 ,0.55];

%% main brain map

if exist('scalpmesh.mat','file') && exist('elec_realigned_67.mat','file') %&& exist( 'elec_realigned_25.mat','file')% elec_realigned or elec_realigned_25 
load ('scalpmesh');
load elec_realigned_67

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % load bnd
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % scalpmesh.pos= bnd(3).pos;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % scalpmesh.tri= bnd(3).tri;
else
    
% inputs for electrode registation and generation of scalp mesh  and other boundaries
    %1. build 67 electrode locations
    elec_67= [];
    %2. read in the anatomical data
    mri  = ft_read_mri('Subject01.mri');
    
% make scalpmesh and elec_realigned
[ scalpmesh, elec_realigned, bnd,ER1, rms_ER1 ] = make_scalpmesh_elec_realigned( elec_67, mri );
save scalpmesh scalpmesh
% brainmesh=bnd(3);
% save brainmesh  brainmesh
save bnd bnd %3 boundaries
save elec_realigned_67  elec_realigned

end

% show_3Dbrainmap

%%%comment one of the below lines

%  show_3Dbrainmap( method, Channels, Values,scalpmesh,elec_realigned)
show_3Dbrainmap2( method, Channels, Values,scalpmesh,elec_realigned)

toc,

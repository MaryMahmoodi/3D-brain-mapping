clear
%% from MRI acquisition to brain mapping %%
%%     initial parameters 
ICPmethod=1;%using ICP for electrode regisration
electrodealign_fieldtrip=1;% using fieldtrip (concluding ICP and warping methods)
showmri=1;% to see MRI (3 slices in 2 diemensional)
show_3Orthogonalslices=1; % plot 3 orthoginal MRI slices (3D)(to localize fiducial points)
showsegmri=1;% to see segmented MRI
%% inputs

method= 'cubic' ;          % cubic % RBF % linear % avg 
fprintf('your electrode labels should be aligned with following labels uppercase or lower case are not important....  ');
sprintf(' Cz,Cpz,Pz,Poz,Oz,Iz,Fpz,AFz,Fz,Fcz,O2,Po4,P2,Cp2,P4,Po8,P10,T6,P6,Tp8,Cp6,Cp4,C2,C4,C6,T4,Ft8,Fc6,Fc4,F4,F6,F8,AF8,Fc2,F2,AF4,Fp2,P1,P3,Po7,O1,Po3,Fc1,F1,AF3,Fp1,AF7,F3,F5,F7,Ft7,Fc5,Fc3,C1,Cp1,C3,C5,T3,Tp7,Cp5,Cp3,P5,T5,P9,nasion,lpa,rpa')

Channels= {'Cz' 'C4' 'T4' 'C3' 'T3' 'Fz' 'F4' 'F8' 'F3' 'F7' 'Fpz' 'Fp2' 'Fp1' 'Pz' 'P3' 'T5' 'P4' 'T6' 'Oz' 'O1' 'O2'}';
Values=[0.42  0.75  0.85  0.45  0.43  0.4  0.43  0.55  0.35  0.37  0.25  0.39  0.35  0.78  0.6  0.6  1  0.95  0.8  0.6  0.87]';


%1. build 67 electrode locations
elec_67= makelec_67(  );

%2. read in the anatomical data
mri  = ft_read_mri('Subject01.mri');


%% main from Mri (segmentation, mesh, electrode registartion) 2 Brainmap    %%%%%%%%%%%

% ft_volumereslice-> interpolates and reslices a volume along the
%   principal axes of the coordinate system according to a specified resolution.
cfg              = [];
cfg.dim          = [256 256 256];                 % original dimension
mri              = ft_volumereslice(cfg,mri);

% plot 3  slices of MRI (2D)
cfg = [];
cfg.method        =  'ortho';%ortho %slice
cfg.anaparameter='anatomy';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2];
cfg.opacitymap    = 'rampup';

if showmri
    ft_sourceplot(cfg, mri);
end

% plot 3 orthoginal MRI slices (3D)
if show_3Orthogonalslices
    fprintf (' write %s and determine coordsys for +x+y+z directions respectively...%s', 'y');
    sprintf('directions: (a)nterior/ (p)osterior, (s)uperior/ (i)nferior, (r)ight/(l)eft')
   
    mri=ft_determine_coordsys(mri);
end
elec_67=[];
elec_67= makelec_67( mri.coordsys);

% 3. segmentation -> obtain the gray, skull and scalp tissue probability maps(tpm)
cfg              = [];
cfg.output       = {'gray','skull','scalp'};
tpm              = ft_volumesegment(cfg,mri);
i=3;%scalp
cfg.funparameter = cfg.output{1,i};
cfg.location     = 'center';
if showsegmri
    ft_sourceplot(cfg, tpm);title (cfg.output{1,i});
end

% 4. mesh generation (3D surface triangulars) from the segmented MRI
cfg             = [];
cfg.tissue      = {'scalp', 'skull', 'gray'};
cfg.numvertices = [10000, 1600, 2400];
bnd             = ft_prepare_mesh(cfg, tpm);
for i=1:length(cfg.tissue);bnd(i).label=cfg.tissue{i} ;end

for i=1:size(bnd,2)
    
    figure; ft_plot_mesh(bnd(i), 'facecolor', 'skin','edgecolor', 'none', 'facealpha',0.4,'edgealpha',0.4);camlight;title(cfg.tissue{i});
     lighting gouraud
      material dull
%       lightangle(0, 90);
      alpha 0.9
      if i==1
          figure; ft_plot_mesh(bnd(i), 'facecolor', 'skin','edgecolor', 'k', 'facealpha',0.4,'edgealpha',0.4);camlight;title(cfg.tissue{i});
     lighting gouraud
      material dull
%       lightangle(0, 90);
      alpha 0.9
      end
end


% 5. Register electrodes on scalp mesh (ICP)
if ICPmethod==1
    [elec_new, scalpmesh,index,ER,rms_ER]=ICP_electrode_realign(cfg,mri,bnd,elec_67);
    p=elec_new.pos;
    q=scalpmesh.pos(index,:);
    ER  = distance_error( p,q );%"mm"
    rms_ER = rms_error(p,q);%"mm"
    
    [skincolor,electrodecolor, fontcolor]=requestcolor();
    fprintf('the rms and distance error with ICP in "mm" is:');ER(end)
    rms_ER(end)
    pause(5);
    % figure(8);clf ;ft_plot_mesh(scalpmesh, 'facecolor', skincolor,'edgecolor', 'none', 'facealpha',0.4,'edgealpha',0.4);camlight;title(cfg.tissue{i});
    % hold on; ft_plot_sens(elec_new,'label','on','fontcolor', fontcolor,'facecolor',electrodecolor);title ('electrode registration with ICP');
end

% 6. Register electrodes on scalp mesh  by fieldtrip

if electrodealign_fieldtrip
    [elec_realigned, scalpmesh,index1,ER1,rms_ER1]=fieldtrip_electrode_realign(mri,bnd,elec_67,elec_new);
    p1=elec_realigned.elecpos;
    q1=scalpmesh.pos(index1,:);
    ER1  = distance_error( p1,q1 );%"mm"
    rms_ER1 = rms_error(p1,q1);%"mm"
    
    [skincolor,electrodecolor, fontcolor]=requestcolor();
    fprintf('the rms and distance error with fieldtrip(with ICP and warping method) in "mm" is:');ER1(end)
    rms_ER1(end)
    pause(5);
    figure(9);clf ;ft_plot_mesh(scalpmesh, 'facecolor', skincolor,'edgecolor', 'k', 'facealpha',0.4,'edgealpha',0.4);camlight;
    title(cfg.tissue{i});
    hold on; ft_plot_sens(elec_realigned,'label','on','fontcolor', fontcolor,'facecolor',electrodecolor);title ('final electrode registration');
end

save scalpmesh scalpmesh
save elec_realigned_67  elec_realigned
save bnd bnd
%   if 0


R=rotz(-2); % after a final registration it seems this electrode file needs a 5 degree rotation around +y axis
for i=1:size (elec_realigned.elecpos,1)
    
    elec_realigned.elecpos(i,:)=(R*elec_realigned.elecpos(i,:)')';
    elec_realigned.chanpos(i,:)=(R*elec_realigned.elecpos(i,:)')';
    elec_realigned.pos(i,:)=(R*elec_realigned.elecpos(i,:)')';
    elec_realigned.pnt(i,:)=(R*elec_realigned.elecpos(i,:)')';
end

[elec_realigned, scalpmesh,index1,ER1,rms_ER1]=fieldtrip_electrode_realign(mri,bnd,elec_67,elec_realigned);
p1=elec_realigned.elecpos;
q1=scalpmesh.pos(index1,:);
ER1  = distance_error( p1,q1 );%"mm"
rms_ER1 = rms_error(p1,q1);%"mm"


save elec_realigned_67 elec_realigned
[skincolor,electrodecolor, fontcolor]=requestcolor();
figure(9);clf ;ft_plot_mesh(scalpmesh, 'facecolor', skincolor,'edgecolor', 'k', 'facealpha',0.4,'edgealpha',0.4);camlight;
hold on; ft_plot_sens(elec_realigned,'label','on','fontcolor', fontcolor,'facecolor',electrodecolor);title ('final electrode registration');

%   end % if 0


if size(elec_realigned.label,1)>1;elec_realigned.label=elec_realigned.label';end
[elec_realigned,scalpmesh, elec_realigned4, index_scalpmesh]=apply_preliminaries_brainmap ( elec_realigned,scalpmesh);

%%% adjust electrodes to desired channels defined by following labels and brain map

elec_realigned.label  ={'Cz';'Cpz';'Pz';'Poz';'Oz';'Iz';'Fpz';'AFz';'Fz';'Fcz';'O2';'Po4';'P2';'Cp2';'P4';'Po8';'P10';'T6';'P6';'Tp8';'Cp6';'Cp4';'C2';'C4';'C6';'T4';'Ft8';'Fc6';'Fc4';'F4';'F6';'F8';'AF8';'Fc2';'F2';'AF4';'Fp2';'P1';'P3';'Po7';'O1';'Po3';'Fc1';'F1';'AF3';'Fp1';'AF7';'F3';'F5';'F7';'Ft7';'Fc5';'Fc3';'C1';'Cp1';'C3';'C5';'T3';'Tp7';'Cp5';'Cp3';'P5';'T5';'P9';'nasion';'lpa';'rpa'};

if length(Channels)> 67
    error('more than %s electrodes is not supported for brain mapping','67')
end

elec_realigned3=[];
elec_realigned3.cfg=elec_realigned.cfg;
counter=1;
for i=1:length(elec_realigned.label  )
    j=[];
    j = find(strcmp(lower(elec_realigned.label{i}),lower (Channels) ));
    if ~isempty(j)
        elec_realigned3.chanpos(counter,:)=elec_realigned.chanpos(i,:);
        elec_realigned3.chantype{counter}=elec_realigned.label{i};
        elec_realigned3.chanunit{counter}=elec_realigned.chanunit {i};
        elec_realigned3.coordsys=elec_realigned.coordsys;
        elec_realigned3.elecpos(counter,:)= elec_realigned.elecpos(i,:);
        elec_realigned3.label{counter}=elec_realigned.label{i};
        elec_realigned3.type=elec_realigned.type;
        elec_realigned3.unit=elec_realigned.unit;
        EEG(counter)=Values(j);
        counter=counter+1;
        
    else
        counter=counter;
    end
end
elec_realigned=elec_realigned3;
elec_realigned.tra=eye(size(elec_realigned.elecpos,1));

[elec_realigned,scalpmesh, elec_realigned2, index_scalpmesh]=apply_preliminaries_brainmap ( elec_realigned,scalpmesh);



elec=elec_realigned2;

[scalpmesh,elec]=apply_brainmap(EEG,elec,scalpmesh, index_scalpmesh,method);%or elec_realigned %scalpmesh.power



scalpmesh2=[];counter=1;counter1=1;scalpmesh3=[];
for i=1:length(index_scalpmesh)
    if ~isempty(index_scalpmesh{i})
        scalpmesh2.pos(counter,:)=scalpmesh.pos(i,:);
        scalpmesh2.power(counter,1)=scalpmesh.power(i);
        counter=counter+1;
    else
        scalpmesh3.pos  (counter1,:)=scalpmesh.pos(i,:);
        counter1=counter1+1;
        
    end
end
scalpmesh2.tri=delaunay(scalpmesh2.pos(:,1), scalpmesh2.pos(:,2),scalpmesh2.pos(:,3));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   [scalpmesh2,neighbours]= make_tri_neighbors(scalpmesh2 );scalpmesh2.tri=scalpmesh2.tri2;

scalpmesh3.tri=delaunay(scalpmesh3.pos(:,1), scalpmesh3.pos(:,2),scalpmesh3.pos(:,3));


elec_realigned2=elec_realigned4;

[electrodecolor, active_electrodecolor]=requestcolor2();


if strcmp(method,'cubic') ;  a1=14;figure(a1);clf;end
if strcmp(method,'linear') ; b1=13; figure(b1);clf;end
if strcmp(method,'RBF') ;  figure(12);clf;end
if strcmp(method,'avg') ;  figure(11);clf;end

ft_plot_mesh(scalpmesh,'facecolor', 'k','edgecolor', 'none', 'facealpha',0.1,'edgealpha',0.1,'maskstyle','opacity');%'maskstyle','opacity'

%  if using scalpmesh.power and 'facealpha',0.4,'edgealpha',0.4
% lighting gouraud
%       material dull
%       lightangle(0, 90);
%       alpha 0.9

hold on; ft_plot_mesh(scalpmesh2,'vertexcolor',scalpmesh2.power ,'edgecolor','none','facealpha',0.4,'edgealpha',0.4,'maskstyle','opacity');title('BrainMap- random signal');%'edgecolor','none or [0 0 0]'
view(2);

hold on; ft_plot_sens(elec_realigned2,'label','on','facecolor',electrodecolor,'fontsize',14);title (method)%'facecolor',electrodecolor

% if size(elec_realigned2.elecpos,1)==size(elec.elecpos,1)
hold on; ft_plot_sens(elec,'label','off','facecolor',active_electrodecolor,'fontsize',14);title (method,'fontsize',14)% plot active electrodes
% else
%     hold on; ft_plot_sens(elec,'label','on','facecolor',active_electrodecolor,'fontsize',14);title (method,'fontsize',14)% plot active electrodes
% end

legend('scalp','electrode cap','deactive electrodes','active electrodes','location','NorthEast')

c=1;
for i=1:length(index_scalpmesh)
    if isempty(index_scalpmesh{i})
        c=c+1;
    end
end

fprintf(' BrainMap - %s electrodes- resolution %s  points...',num2str(length(Channels)),num2str(c));

cM=colorMap([0:0.01:1]');
colormap(cM);
cc1=colorbar('Location','eastoutside');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % figure(13);clf;ft_plot_topo3d(scalpmesh.pos, scalpmesh.power);title('topoplot');colorbar;title(method);  hold on; ft_plot_sens(elec,'label','on','facecolor','g','fontsize',14);title (method)


% if using  [scalpmesh2,neighbours]= make_tri_neighbors(scalpmesh2 );scalpmesh2.tri=scalpmesh2.tri2;
% lighting gouraud
%       material dull
%       lightangle(0, 90);
%       alpha 0.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% plot with mesh triangles
if strcmp(method,'cubic') ;  a1=14;figure(a1+4);clf;end
if strcmp(method,'linear') ; b1=13; figure(b1+4);clf;end
if strcmp(method,'RBF') ;  figure(12+4);clf;end
if strcmp(method,'avg') ;  figure(11+4);clf;end

ft_plot_mesh(scalpmesh,'facecolor', 'k','edgecolor', 'k', 'facealpha',0.09,'edgealpha',1);%camlight
view(3);


hold on; ft_plot_mesh(scalpmesh2,'vertexcolor',scalpmesh2.power ,'edgecolor','none','facealpha',0.4,'edgealpha',0.4);title('BrainMap- random signal');%'edgecolor','none or [0 0 0]'
view(2);

hold on; ft_plot_sens(elec_realigned2,'label','on','facecolor',electrodecolor,'fontsize',14);title (method)%'facecolor',electrodecolor


hold on; ft_plot_sens(elec,'label','off','facecolor',active_electrodecolor,'fontsize',14);title (method,'fontsize',14)% plot active electrodes
legend('scalp','electrode cap','deactive electrodes','active electrodes','location','NorthEast')

c=1;
for i=1:length(index_scalpmesh)
    if isempty(index_scalpmesh{i})
        c=c+1;
    end
end

fprintf(' BrainMap - %s electrodes- electrode cap resolution %s  points...',num2str(length(Channels)),num2str(c));


cM=colorMap([0:0.01:1]');
colormap(cM);

cc1=colorbar('Location','eastoutside');


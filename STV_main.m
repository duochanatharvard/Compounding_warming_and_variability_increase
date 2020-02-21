% *************************************************************************
% Prepare for data
% *************************************************************************
clear;

clearvars -except exp_num dir

dir = [pwd,'/'];

load([dir,'Subset_summertime_tas_all_models_20N_poleward.mat'],'tas');
tas_detrend = CDC_detrend(tas,4,3);

lon   = repmat([1.25:2.5:360]',1,28);
lat   = repmat([21.25:2.5:90],144,1);

clear('model_reg','model_reg_member','model_var','model_var_std')

for md = 1:23
    for prd = 1:2
        
        disp(['Model: ',num2str(md),' Period: ',num2str(prd)])
        
        N = 100;
        N_block = 6;
        dim = 3;
        input_abs     = squeeze(tas(:,:,md,:,prd));
        input_detrend = squeeze(tas_detrend(:,:,md,:,prd));
        
        disp('compute variance begins')
        [output_mean, out_mean_member, ~] = CDC_mean_bt(input_abs,N,dim,N_block);
        model_mean(:,:,md,prd) = output_mean;
        model_mean_std(:,:,md,prd) = CDC_std(out_mean_member,3);

        [output_var,~,out_var_member,out_std] = CDC_var_bt(input_detrend,N,dim,[],N_block);
        model_var(:,:,md,prd) = output_var;
        model_var_std(:,:,md,prd) = CDC_std(out_var_member,3);
        disp('compute variance completes')
        
        % compute regional mean
        for ct_reg = 1:3
            
            switch ct_reg
                case 1   % Europe mean
                    MASK = zeros(size(lon,1),size(lat,2));
                    MASK(lon>0 & lon<60 & lat>45 & lat<60) = 1;
                case 2   % Canada mean
                    MASK = zeros(size(lon,1),size(lat,2));
                    MASK(lon>240 & lon<270 & lat>50 & lat<60) = 1;
                case 3   % High latitude
                    MASK = zeros(size(lon,1),size(lat,2));
                    MASK(lat>60 & lat<80) = 1;
            end
            
            lat_in   = 21.25:2.5:90;
            
            model_var_reg(md,prd,ct_reg)          = CDC_mask_mean(output_var,lat_in,MASK);
            model_var_reg_member(md,prd,:,ct_reg) = CDC_mask_mean(out_var_member,lat_in,MASK);

            model_mean_reg(md,prd,ct_reg)          = CDC_mask_mean(output_mean,lat_in,MASK);
            model_mean_reg_member(md,prd,:,ct_reg) = CDC_mask_mean(out_mean_member,lat_in,MASK);
   
        end  
    end
end

file_save = [dir,'ACDC_paper_revisit_2020.mat'];
save(file_save,'model_var','model_var_std','model_var_reg','model_var_reg_member',...
    'model_mean','model_mean_std','model_mean_reg','model_mean_reg_member','lon','lat','-v7.3')
 

%% *************************************************************************
% Generating figures
% *************************************************************************
clear;

dir  = [pwd,'/'];
load([dir,'ACDC_paper_revisit_2020.mat'])

model_use  = {'ACCESS1-0'  ,  'ACCESS1-3' ,   'CCSM4'    ,'CESM1-BGC' ,   'CESM1-CAM5'   , 'CNRM-CM5'  ,  'CSIRO-Mk3-6-0',...
    'GFDL-CM3','GFDL-ESM2G',    'GFDL-ESM2M'  ,  'GISS-E2-H-CC' ,   'GISS-E2-H' ,   'GISS-E2-R-CC' ,   'GISS-E2-R',...
    'HadGEM2-CC','HadGEM2-ES' ,   'inmcm4'  ,  'IPSL-CM5A-MR' ,   'MIROC5'  ,  'MRI-CGCM3'  ,  'MRI-ESM1',...
    'NorESM1-ME','NorESM1-M'};

l_use = true(23,1);

% Set the style of all maps in the paper
load('mask.mat')
land_mask_0 = [ones(144,44) mask ];
st_full_subcoast = {'region',[-180 180 20 80],'mask',land_mask_0,...
    'fontsize',24,'barloc','southoutside',...
    'xtick',[-90 0 90],'ytick',[10:20:70],'bckgrd',[1 1 1]*.9,...
    'subcoast',1,'docoast',1,'coastwi',1,'daspect',[1 0.6 1]};


%% *************************************************************************
% Figure 1: Changes in mean and variability and cross-model uncertainty of changes
% *************************************************************************
figure(1);  clf;
pic = [nan(144,44) nanmean(model_mean(:,:,l_use,2) - model_mean(:,:,l_use,1),3)];
pic_std = [nan(144,44) CDC_std(model_mean(:,:,l_use,2) - model_mean(:,:,l_use,1),3)];
clear('sig')
% sig(:,:,1) = abs(pic./pic_std) > 1.96;
sig(:,:,1) = [nan(144,44) sum((model_mean(:,:,l_use,2) - model_mean(:,:,l_use,1))>0,3) > (23*.75)]; 
sig(:,:,2) = 0;
CDF_plot_map('pcolor',pic,...
    [st_full_subcoast,'crange',6,'cmap',b2rCD(12),...
    'bartit','CMIP5 mean change in temperature [^oC]','plabel',' ','sig',sig,'sigtype','marker']);

figure(2);  clf;
pic = [nan(144,44) nanmean(model_var(:,:,l_use,2) - model_var(:,:,l_use,1),3)];
pic_std = [nan(144,44) CDC_std(model_var(:,:,l_use,2) - model_var(:,:,l_use,1),3)];
clear('sig')
% sig(:,:,1) = abs(pic./pic_std) > 1.96;
sig(:,:,1) = [nan(144,44) sum((model_var(:,:,l_use,2) - model_var(:,:,l_use,1))>0,3) > (23*.75)]; 
sig(:,:,2) = 0;
CDF_plot_map('pcolor',pic,...
    [st_full_subcoast,'crange',2,'cmap',b2rCD(12),...
    'bartit','CMIP5 mean change in variance [^oC^2]','plabel',' ','sig',sig,'sigtype','marker']);

figure(3);  clf;
aa  = model_mean(:,:,l_use,2) - model_mean(:,:,l_use,1);
q   = quantile(aa,[0.25 0.75],3);
qr  = q(:,:,2) - q(:,:,1);
pic = [nan(144,44) qr];
CDF_plot_map('pcolor',pic,...
    [st_full_subcoast,'crange',6,'cmap',hotCD(12),...
    'bartit','Inter-quantile range of change in temperature [^oC]','plabel',' ']);

figure(4);  clf;
aa  = model_var(:,:,l_use,2) - model_var(:,:,l_use,1);
q   = quantile(aa,[0.25 0.75],3);
qr  = q(:,:,2) - q(:,:,1);
pic = [nan(144,44) qr];
CDF_plot_map('pcolor',pic,...
    [st_full_subcoast,'crange',2,'cmap',hotCD(12),...
    'bartit','Inter-quantile range of change in variance [^oC^2]','plabel',' ']);


%% *************************************************************************
% Figure 2 - map: variability change and warming
% *************************************************************************
[slope, ~, slope_member, ~] = ...
    CDC_yorkfit_bt(model_var(:,:,l_use,2)-model_var(:,:,l_use,1),model_mean(:,:,l_use,2)-model_mean(:,:,l_use,1),...
                   sqrt(model_var_std(:,:,l_use,2).^2 + model_var_std(:,:,l_use,1).^2),...
                   sqrt(model_mean_std(:,:,l_use,2).^2 + model_mean_std(:,:,l_use,1).^2),0,3,1000);
                
figure(5); clf;

ref = 0;
sig_data = (quantile(slope_member,0.025,3) > ref & slope > ref) | ...
    (quantile(slope_member,0.975,3) < ref & slope < ref);

pic = [nan(144,44) slope];
clear('sig')
sig(:,:,1) = [zeros(144,44) sig_data];
sig(:,:,2) = 0;
CDF_plot_map('pcolor',pic,...
    [st_full_subcoast,'crange',0.8,'cmap',b2rCD(12),...
    'bartit','Scaling between changes in mean temperature and variance [^oC^2 / ^oC]','plabel',' ','sig',sig,'sigtype','marker']);
caxis([-1 1]*.6)

reg_list = [0 60 45 60; -120 -90 50 60];
for ct = 1
m_plot([reg_list(ct,1) reg_list(ct,2) reg_list(ct,2) reg_list(ct,1) reg_list(ct,1)],...
    [reg_list(ct,3) reg_list(ct,3) reg_list(ct,4) reg_list(ct,4) reg_list(ct,3)],...
    '-','color',[0 0 1],'linewi',3);
end

set(gcf,'position',[.1 10 15 9],'unit','inches')
set(gcf,'position',[.1 14 15 9],'unit','inches')


%% *************************************************************************
% Figure 2 - scatter: variability change and warming
% *************************************************************************
figure(6); clf; hold on

ct_reg = 1;

x       = model_mean_reg(l_use,2,ct_reg) - model_mean_reg(l_use,1,ct_reg);
y       = model_var_reg(l_use,2,ct_reg)  - model_var_reg(l_use,1,ct_reg);
x_std   = CDC_std(model_mean_reg_member(l_use,2,:,ct_reg) - model_mean_reg_member(l_use,1,:,ct_reg),3);
y_std   = CDC_std(model_var_reg_member(l_use,2,:,ct_reg) - model_var_reg_member(l_use,1,:,ct_reg),3);
r       = 0;

N       = 10000;
P.mute_output = 1;

l = 1:23;

[slope, inter, slope_member, inter_member] = CDC_yorkfit_bt(y(l),x(l),y_std(l),x_std(l),r,1,N,P);

[slope_member,I] = sort(slope_member);
inter_member = inter_member(I);

% Plot the shading for 2 s.d.
clear('yy')
x1_pic = 2;
x2_pic = 9;
x_temp = x1_pic:0.01:x2_pic;

for ct = 1:numel(x_temp)
    yy(ct,:) = quantile(x_temp(ct) * slope_member + inter_member,[0.025 0.975]);
end
patch([x_temp,fliplr(x_temp)],[yy(:,1)' fliplr(yy(:,2)')],[1 1 1]*.6,'facealpha',0.3,'linest','none')

plot([x1_pic x2_pic],[x1_pic x2_pic]*slope+inter,'w-','linewi',4)
plot([x1_pic x2_pic],[x1_pic x2_pic]*slope+inter,'r-','linewi',2)
disp([num2str(slope,'%6.2f'),' ',num2str(quantile(slope_member,[0.025 0.5 0.975]),'%6.2f')])
[out_st,out_col] = CDF_scatter(x,y,5,'x_std',x_std,'y_std',y_std,'mksize',13);
daspect([7 4 1])
CDF_panel([x1_pic x2_pic -1 3],'',{},'changes in mean temperature [^oC]','changes in variance [^oC^2]','fontsize',22)
set(gcf,'position',[.1 10 7 7]*0.9,'unit','inches')
set(gcf,'position',[.1 10 7 7]*1,'unit','inches')

% plot legend
figure(13); clf; hold on;
CDF_scatter_legend(model_use,out_st,out_col,[1,4,1],'fontsize',15,'mksize',13)
axis([1 7.3 -24 0])
set(gcf,'color','w')
set(gcf,'position',[.1 1 3 6],'unit','inches')
set(gcf,'position',[.1 1 18 10],'unit','inches')


%% *************************************************************************
% Figure 3 - risk analysis
% *************************************************************************
clear;

target    = 32;
hist_clim = 19.38;
hist_var  = 2.54;

warming   = 3.3;  % warming rate such that the variance does not change
rcp_var   = hist_var + 0.40 * warming - 1.32;

t = 1:0.01:40;

figure(1); clf; hold on
h(1) = plot(t,normpdf(t,hist_clim,sqrt(hist_var)),'k','linewi',2);
warming = 7.5;
h(2) = plot(t,normpdf(t,hist_clim + warming,sqrt(hist_var + 0.40 .* warming - 1.32)),'r','linewi',2);
h(3) = plot(t,normpdf(t,hist_clim + warming,sqrt(rcp_var)),'r -- ','linewi',2);
h(4) = plot([1 1]*target,[0 0.30],'k--','linewi',2); 

a = normpdf(t,hist_clim + warming,sqrt(hist_var + 0.40 .* warming - 1.32));
a(t < target) = 0;
nansum(a)
h(5) = bar(t,a,'linest','none','facecolor',[1 .4 .4]);

a    = normpdf(t,hist_clim + warming,sqrt(rcp_var));
a(t < target) = 0;
nansum(a)
h(6) = bar(t,a,'linest','none','facecolor',[1 .7 .7]);

CDF_panel([10 37 0 0.03],'','','Monthly mean temperature [^oC]','pdf','fontsize',20)

set(gcf,'position',[.1 10 15 9],'unit','inches')
set(gcf,'position',[.1 14 11 5],'unit','inches')


clear('Tab')
ct_1 = 0;
intv = 0.05;
l1   = 28: intv :35;
l2   =  0: intv :8;
list = [1/10000000 1/200 1/100 1/50 1/30 1/20 1/15 1/10 1/5 1/3 1/2 1 2 3 5 10 15 20 30 50 100 200 1000000];
for target = l1
    ct_1 = ct_1 + 1;
    ct_2 = 0;
    for warming = l2
        ct_2 = ct_2 + 1;
        Tab(ct_1,ct_2,1) = 1 - normcdf(target,hist_clim + warming,sqrt(rcp_var));
        Tab(ct_1,ct_2,2) = 1 - normcdf(target,hist_clim + warming,sqrt(hist_var + 0.4 .* warming - 1.32));
    end
end


figure(2); clf; hold on;
a = (Tab(:,:,2)) ./ Tab(:,:,1);
a = discretize(a,list);
CDF_pcolor(l2,l1,a'+0.1);

for i = 3:1:8    plot([i i],[28 40],'color',[1 1 1]*.75); end
for i = 29:1:34  plot([2 8],[i i],'color',[1 1 1]*.75,'linewi',1); end

a = 1 ./ Tab(:,:,1);
[C,h] =  contour(l2,l1,a,[1 10 100 1000 10000 100000],'k--','linewi',2);
clabel(C,h,'LabelSpacing',800,'Color','k','FontWeight','bold','FontSize',13)

a = 1 ./ Tab(:,:,2);
[C,h] =  contour(l2,l1,a,[1 10 100 1000 10000 100000],'b','linewi',2);
clabel(C,h,'LabelSpacing',2000,'Color','b','FontWeight','bold','FontSize',13)


b2rCD(11,0);
CDF_panel([3 8 28 35],'','','Warming [^oC]','Threshold [^oC]','fontsize',20)
h = colorbar;
ylabel(h,'Probability ratio')
a = {'','1/200','1/100','1/50','1/30','1/20','1/15','1/10','1/5','1/3','1/2',...
    '1','2','3','5','10','15','20','30','50','100','200',''};
set(h,'ytick',[2:2:numel(list)],'yticklabel',a(2:2:end));
set(gca,'xtick',[3 4 5 6 7 8])
set(gca,'ytick',28:1:35)
caxis([1 numel(list)])
daspect([1 3 1])

set(gcf,'position',[.1 10 15 9],'unit','inches')
set(gcf,'position',[.1 10 11 9],'unit','inches')
%Compare results from several conditions output from the script:
%"Track_Results_MSD_V2.n" or "SMT_Motion_Classifier.m"
%Pablo Aurelio Gomez Garcia, 2018.

clear all
close all
clc


%% Inputs
%------------------------------------------------------------------
%Data---------------------------------------------------
%------------------------------------------------------------------
%Set the data type description (Maximum 8 types)
type.t1 = 'H2B NPC';
type.t2 = 'H2B ESC';
type.t3 = 'H2B ESC H1tKO';
type.t4 = 'H2B Fixed';
% type.t5 = 'FOXA1 Directed';
% type.t6 = 'H3 Directed';
% type.t7 = 'FOXA1 Butterfly';
% type.t8 = 'H3 Butterfly';

% type.t1 = 'ESC butt direct';
% type.t2 = 'ESC brownian';
% type.t3 = 'Oct4';
% type.t4 = '';

% type.t1 = 'ESC H2B PairedTracks';
% type.t2 = 'ESC H2B UnpairedTracks';
% type.t3 = 'ESC H2B Euchromatin';

% type.t1 = 'ESC H2B';
% type.t2 = 'ESC H1tKO H2B';
% type.t3 = 'ESC Centromeres';
% type.t4 = 'ESC Telomeres';
% type.t5 = 'Fixed Cells';



%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)
% startlocn = {'C:\Users\pgomez.CRG\Desktop\Temporal\FoxA1\JON SMT data\FOXA1 WT\MSD_Results_6min_frames_subROIs0_Inside0'};
startlocn = {'C:\Users\pgomez.CRG\Documents\PhD\My papers\SMT\Data\H2B 15ms'};

Save_Results = 0;
File_name = 'MSD_Results';


%------------------------------------------------------------------
%Numerical Inputs---------------------------------------------------
%------------------------------------------------------------------
n_dim = 2; %Dimensionality of the movement (2 for 2D and 3 for 3D).
Frame_interval = 0.015; %Exposure time in seconds

TMSD_fitting_points = 3;   %Minimum 3 point for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 3;  %Minimum 3 point for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; 

R2LIMIT = 0.7; %R-squared minimum value from the D fits on each T-MSD Curve.

%-----------------------------------------------------------------
%For the Confined Circle Difussion Model Fitting---------------------------------------------------
%---------------------------------------------------------------------
num_points = 6; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
level = 0.01; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in ?m]. 
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (um^2/s)
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (um)

Radius_interval = [30 60]; %Radius interval (in nm) in order to divide tracks into several populations. (log10(60) = 1.78)

%-----------------------------------------------------------------
%For the Histogram of Difussion Coefficients---------------------------------------------------
%---------------------------------------------------------------------
binWidth = 0.1; %Bin width
bin_start = -3; %Bin start
bin_end = 1.5; %Bin end

%-----------------------------------------------------------------
%For the angles analysis---------------------------------------------------
%---------------------------------------------------------------------
AC_arc = 30; %Arc in degrees in order to calculate the assimetry coefficient. 30 will correspond to 150º-210º versus 30º-33º.
AC_bin_size = 50; %Bin size in nm for grouping the pair of jumps and calculating the assimetry coefficient
AC_num_groups = 10; %Number of groups to make. Notice that maximum jump would be AC_num_groups*AC_bin_size

max_lag = 50; %Maximum number of lag times to study on the assimetry coefficient of angles distribution.

%-----------------------------------------------------------------
%For the Survival distribution fitting---------------------------------------------------
%---------------------------------------------------------------------
k1_init = 2; %Starting value for k1. (1/s). 2 means 0.5seconds.
k2_init = 0.02; %Starting value for k2. (1/s). 0.02 means 50 seconds.
fraction_init = 0.5; %Initial value for the fraction belonging to each population.
range = 2:300; %Point range in order to perform the fitting.

k1_init_bleach = 4; %Starting value for k1. (1/s). 2 means 0.5seconds.
k2_init_bleach = 0.01; %Starting value for k2. (1/s). 0.02 means 50 seconds.
fraction_init_bleach = 0.5; %Initial value for the fraction belonging to each population.
range_bleach = 1:300; %Point range in order to perform the fitting.


%-----------------------------------------------------------------
%For the Radius of gyration analysis---------------------------------------------------
%---------------------------------------------------------------------
k_init_gyration = 1; %This is the initial half life value for the exponential fitting of the radius of gyration time evolution.
range_gyration = 1:100; %Range for the fitting of the radius of gyration evolution over time.
edges_gyration_density_plot = {0:5:100 0:5:100}; %This is for the density plot of radius of gyration vs radius of confinement.

%------------------------------------------------------------------------
%Initialize variables-------------------------------------------------
%----------------------------------------------------------------------
Dcoeffs = {};
datatypes = fieldnames(type);
ntypes = size(datatypes,1);         
str_rad = {strcat('>',num2str(Radius_interval(2)),'nm') strcat('[',num2str(Radius_interval(1)),'nm-',num2str(Radius_interval(2)),'nm]') strcat('<',num2str(Radius_interval(1)),'nm')};


%% Load the Data
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%Ask the user to select the files
files = cell(ntypes,1);
for t = 1:ntypes;
    files{t}= Select1DataGroup(strcat(type.(datatypes{t}),'_MSD_Results'),'mat',startlocn{1});
    files_Reference{t}= Select1DataGroup(strcat(type.(datatypes{t}),'_CellReference'),'mat',startlocn{1});
    files_CenterOfMass{t}= Select1DataGroup(strcat(type.(datatypes{t}),'_CellCenterOfMass'),'mat',startlocn{1});
end

%Load the .mat files with the Results of the Tracking Analysis
%It will ONLY Load one file per category!!!
% for t=1:ntypes;
%     MSD_Results{t}=load(strcat(files{t}.data{1,2},files{t}.data{1,1}));
% end

%Load the .mat files with the Results of the Tracking Analysis
%You can load multiple files per category.
for t=1:ntypes;

MSD_Results{t}.ma = msdanalyzer(2, 'um', 's');  %Initialize MSD analizer
  
    for c=1:size(files{t}.data,1);
        
    MSD_Results_prov{t}{c} = load(strcat(files{t}.data{c,2},files{t}.data{c,1}));
    Cell_Reference_prov{t}{c} = struct2cell(load(strcat(files_Reference{t}.data{c,2},files_Reference{t}.data{c,1})));

    
    MSD_Results{t}.ma = MSD_Results{t}.ma.addAll(MSD_Results_prov{t}{c}.ma.tracks);

    end
    
    MSD_Results{t}.ma = MSD_Results{t}.ma.computeMSD;
%     MSD_Results{t}.ma = MSD_Results{t}.ma.LogTMSD(TLOGLOG_fitting_points);
%     MSD_Results{t}.ma = MSD_Results{t}.ma.TMSD(TMSD_fitting_points); 

  Cell_Reference{t} = cat(1,Cell_Reference_prov{t}{:});
  Cell_Reference{t} = cat(1,Cell_Reference{t}{:});
  
  
  Cell_CenterOfMass{t} = struct2cell(load(strcat(files_CenterOfMass{t}.data{1,2},files_CenterOfMass{t}.data{1,1})));

end




%% Calculate some numbers
for t=1:ntypes;
    
Num_Tracks(t) = size(MSD_Results{t}.ma.tracks,1);

for j=1:Num_Tracks(t);
Track_Length(j,t) = size(MSD_Results{t}.ma.tracks{j},1);
Residence_Times(j,t)=Track_Length(j,t)*Frame_interval;
end

end

%Compute the Residence times and the empirical Cumulative Density
%Distribution to obtain the Survival curve of the Tracks
% [f1,x1] = ecdf(Residence_Times(:,1));




%% Calculating and plotting

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure()
for t=1:ntypes;
      fprintf('\n');
      fprintf('------------------------------------');
      fprintf('\n');
    fprintf('TE-MSD , Global Behaviour');
    fprintf('\n');
    fprintf(files{t}.dataname(1:end-9));
    fprintf('\n');
    
    Color=lines(ntypes);
    msmsd = MSD_Results{t}.ma.getMeanMSD;
    %MSD_Results{t}.ma.plotMeanMSD(gca, true);

     h(t) = plot(msmsd(:,1), msmsd(:,2),'Color',Color(t,:),'LineWidth', 2);

    
    hold on
%     xlim([0 tmax])
%     ylim([0 ymax]);

    %Fitting curve
    [fo, gof,D_lsqFit,lsq_conf_int] = MSD_Results{t}.ma.TEMSD(TEMSD_fitting_points);
    ci = confint(fo);
%     D_ensemble(t) = fo.p1/2/n_dim; 
%     D_emsemble_upper_ci(t) = ci(1)/2/n_dim;
%     D_emsemble_lower_ci(t) = ci(2)/2/n_dim;
    D_ensemble(t) = D_lsqFit(1)/2/n_dim; 
    D_emsemble_upper_ci(t) = lsq_conf_int(1,1)/2/n_dim;
    D_emsemble_lower_ci(t) = lsq_conf_int(1,2)/2/n_dim;
    
%     plot(fo);
    
    if ntypes == 1;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))));
    end
    
    hold on
    if ntypes == 2;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))));
    end
    if ntypes == 3;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))));
    end
    if ntypes ==4;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))));
    end
    
   if ntypes ==5;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))));
   end
    
     if ntypes ==6;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))));
    end
    
 if ntypes ==7;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))));
 end
    
     if ntypes ==8;
    l=legend([h],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))));
    end
    
    set(l, 'Interpreter', 'none');
    
    xlabel('Time lag (seconds)','FontSize',20,'FontWeight','bold');
    ylabel('MSD (µm^2)','FontSize',20,'FontWeight','bold');
   
    hold on
    title('Time Ensembled - Mean Square Displacement','FontSize',20,'FontWeight','bold');
    
set(l, 'FontSize',20) 
set(gca,'FontSize',20,'FontWeight','bold');

    
fprintf('\n');
fprintf('----------------------------------------------');
fprintf('\n');
 
   
end



%% --------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure()
%Plot the diffusion coefficient for the different conditions fitted from
%the TE-MSD. Fitted value +- confidence intervals
for t=1:ntypes;

ha{t} = errorbar(t, D_ensemble(t),D_ensemble(t)-D_emsemble_lower_ci(t), D_emsemble_upper_ci(t)-D_ensemble(t),'o','MarkerSize',5,'MarkerFaceColor',Color(t,:),'Linewidth',1.4);
% xlim([0 0.1]);
% ylim([-0.02 0.05]);

%   xlabel('');
%   ylabel('');
%   title(''); 
  
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
     end
      
title('Errorbar plot of Diffusion Coefficients obtained from linear fitting of the TE-MSD','FontSize',20,'FontWeight','bold');   
set(l, 'Interpreter', 'none')  
set(l, 'FontSize',20) 
set(gca, 'XTick', 0:ntypes)
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Diffusion Coefficient (um2/s)','FontSize',20,'FontWeight','bold');
hold on
drawnow
errorbar_tick(ha{t},50);    

    
end



%Plot the Empirical Survival curve of the Tracks Normalized to ymax
% scatter(x1,(1-f1)*ymax);



%% Histogram of diffusion Coefficients
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure()
%This other approach fits the every single MSD curve and then the histogram
%of D is very informative
for t=1:ntypes;
  fprintf('\n');
fprintf('----------------------------------------------');
fprintf('\n');
fprintf('T-MSD , Fitting each MSD curve');
fprintf('\n');
    fprintf(files{t}.dataname(1:end-9));
    fprintf('\n');
    
    
Color=lines(ntypes);
ma{t} = MSD_Results{t}.ma.TMSD(TMSD_fitting_points);
good_enough_fit = ma{t}.lfit.r2fit > R2LIMIT;
Dmean = mean( ma{t}.lfit.a(good_enough_fit) ) / 2 / ma{t}.n_dim;
Dstd  =  std( ma{t}.lfit.a(good_enough_fit) ) / 2 / ma{t}.n_dim;



Ds = ma{t}.lfit.a(good_enough_fit)/2/ma{t}.n_dim;
%Take out negative values from the Diffusion Coefficients list
idx2 = find(Ds>0);
Ds = Ds(idx2);

Num_D(t) = length(Ds);
D_individual{t} = Ds;


%Plot Histogram specifying bin sizes
%histogram(log(Ds),bins,'FaceColor',Color(t,:));

%Plot histogram of the Log(Ds) specifying bin centers and normalizing
%The Diffusion Coefficients distribution is a Log-normal, so that the
%Log(Ds) follows a Normal distribution.
   binCtrs = bin_start:binWidth:bin_end; %Bin centers, depends on your data
   n=length(Ds);
   counts = hist(log10(Ds),binCtrs);
   prob = counts / (n * binWidth);
   H = bar(binCtrs,prob,'hist');
   set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
   hold on;
   
   
 %Fit a Normal distribution to each condition
 %------------------------------------------------------------------
 %Test Bimodality of the distribution:
probability_unimodal_D(t) = test_bimodality(log10(Ds));
 
% %Fit a Normal distribution
% 
% if probability_unimodal(t) > 0.9;
%    pd(t) = fitdist(log10(Ds),'normal');
%    y = pdf(pd(t),binCtrs);
%    plot(binCtrs,y,'LineWidth',2,'Color',Color(t,:));
%    xlim([bin_start bin_end]);
%    mu_normal(t) = pd(t).mu;
%    sigma_normal(t) = pd(t).sigma;
%    D_est(t) = exp(mu_normal(t));
%    D_sigma_est(t) = exp(mu_normal(t)+sigma_normal(t))-exp(mu_normal(t));
%    
% else

%Fit a Bimodal Normal distribution   
%Define the function. (5 parameters)
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2)p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
                     
pStart = 0.5;
muStart = quantile(log10(Ds),[0.25 0.75]);
sigmaStart = sqrt(var(log10(Ds)) - 0.25*diff(muStart).^2);
start = [pStart muStart sigmaStart sigmaStart];

lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];
paramEsts{t} = mle(log10(Ds), 'pdf',pdf_normmixture, 'start',start,'lower',lb, 'upper',ub);

xgrid = linspace(1.1*min(log10(Ds)),1.1*max(log10(Ds)),200);
pdfgrid{t} = pdf_normmixture(xgrid,paramEsts{t}(1),paramEsts{t}(2),paramEsts{t}(3),paramEsts{t}(4),paramEsts{t}(5));
hold on
plot(xgrid,pdfgrid{t},'-','LineWidth',2,'Color',Color(t,:));
xlabel('x')
ylabel('Probability Density')
xlim([bin_start bin_end]);


% end  

Dcoeffs{t} = Ds;

fprintf('\n');
fprintf(strcat('Total number of Tracks = ',num2str(Num_Tracks(t))));
fprintf('\n');
fprintf('**Estimation of the diffusion coefficient from linear fit of the MSD curves (Fitting every MSD curve)**:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));
fprintf('\n');

end


if ntypes == 1;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting');    
end


if ntypes == 2;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting');
end

if ntypes == 3;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting',strcat(files{3}.dataname(1:end-9),' Num of D = ',num2str(Num_D(3))),'Gaussian Fitting');
    
end

if ntypes == 4;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting',strcat(files{3}.dataname(1:end-9),' Num of D = ',num2str(Num_D(3))),'Gaussian Fitting',strcat(files{4}.dataname(1:end-9),' Num of D = ',num2str(Num_D(4))),'Gaussian Fitting');
    
end

if ntypes ==5;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting',strcat(files{3}.dataname(1:end-9),' Num of D = ',num2str(Num_D(3))),'Gaussian Fitting',strcat(files{4}.dataname(1:end-9),' Num of D = ',num2str(Num_D(4))),'Gaussian Fitting',strcat(files{5}.dataname(1:end-9),' Num of D = ',num2str(Num_D(5))),'Gaussian Fitting');
end

 if ntypes ==6;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting',strcat(files{3}.dataname(1:end-9),' Num of D = ',num2str(Num_D(3))),'Gaussian Fitting',strcat(files{4}.dataname(1:end-9),' Num of D = ',num2str(Num_D(4))),'Gaussian Fitting',strcat(files{5}.dataname(1:end-9),' Num of D = ',num2str(Num_D(5))),'Gaussian Fitting',strcat(files{6}.dataname(1:end-9),' Num of D = ',num2str(Num_D(6))),'Gaussian Fitting');   
 end
 
 if ntypes ==7;
 l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting',strcat(files{3}.dataname(1:end-9),' Num of D = ',num2str(Num_D(3))),'Gaussian Fitting',strcat(files{4}.dataname(1:end-9),' Num of D = ',num2str(Num_D(4))),'Gaussian Fitting',strcat(files{5}.dataname(1:end-9),' Num of D = ',num2str(Num_D(5))),'Gaussian Fitting',strcat(files{6}.dataname(1:end-9),' Num of D = ',num2str(Num_D(6))),'Gaussian Fitting',strcat(files{7}.dataname(1:end-9),' Num of D = ',num2str(Num_D(7))),'Gaussian Fitting');
 end
 
if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9),' Num of D = ',num2str(Num_D(1))),'Gaussian Fitting',strcat(files{2}.dataname(1:end-9),' Num of D =',num2str(Num_D(2))),'Gaussian Fitting',strcat(files{3}.dataname(1:end-9),' Num of D = ',num2str(Num_D(3))),'Gaussian Fitting',strcat(files{4}.dataname(1:end-9),' Num of D = ',num2str(Num_D(4))),'Gaussian Fitting',strcat(files{5}.dataname(1:end-9),' Num of D = ',num2str(Num_D(5))),'Gaussian Fitting',strcat(files{6}.dataname(1:end-9),' Num of D = ',num2str(Num_D(6))),'Gaussian Fitting',strcat(files{7}.dataname(1:end-9),' Num of D = ',num2str(Num_D(7))),'Gaussian Fitting',strcat(files{8}.dataname(1:end-9),' Num of D = ',num2str(Num_D(8))),'Gaussian Fitting');    
 end
set(l, 'Interpreter', 'none');

xlabel('Log10(D)/µm^2/s','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontSize',20,'FontWeight','bold');
title('Log10(D) Histogram and Gaussian Fitting','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');


% fprintf('\n');
% 
% fprintf('**Estimation of the diffusion coefficient from Gaussian fitting of the Logaritmic distribution of Ds:\n')
% fprintf('D = %.3g ± %.3g (N = %d)\n', ...
%     D_est(t), D_sigma_est(t)/2, Num_Tracks(t));
% fprintf('\n');
% fprintf('----------------------------------------------');
% fprintf('\n');

%% -------------------------------------------------------------------------------
%---------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure()
%Boxplot comparing the diffusion coefficients from the T-MSD from different
%conditions
if ntypes == 2;

C = [D_individual{1}' D_individual{2}'];
grp = [zeros(1,length(D_individual{1})),ones(1,length(D_individual{2}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9))});
title('Difusion coeficients')
% ylim([-0.02 0.6]);

 end
 
 if ntypes == 3;    

C = [D_individual{1}' D_individual{2}' D_individual{3}'];
grp = [zeros(1,length(D_individual{1})),ones(1,length(D_individual{2})),2*ones(1,length(D_individual{3}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9))});
title('Difusion coeficients')
% ylim([-0.02 0.6]);

 end
 
 if ntypes ==4;    
     
C = [D_individual{1}' D_individual{2}' D_individual{3}' D_individual{4}'];
grp = [zeros(1,length(D_individual{1})),ones(1,length(D_individual{2})),2*ones(1,length(D_individual{3})),3*ones(1,length(D_individual{4}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9))});
title('Difusion coeficients')
% ylim([-0.02 0.6]);

 end

  
 if ntypes ==5;
     
C = [D_individual{1}' D_individual{2}' D_individual{3}' D_individual{4}' D_individual{5}'];
grp = [zeros(1,length(D_individual{1})),ones(1,length(D_individual{2})),2*ones(1,length(D_individual{3})),3*ones(1,length(D_individual{4})),4*ones(1,length(D_individual{5}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9))});
title('Difusion coeficients')
% ylim([-0.02 0.6]);

 end

 
 if ntypes ==6;

C = [D_individual{1}' D_individual{2}' D_individual{3}' D_individual{4}' D_individual{5}' D_individual{6}'];
grp = [zeros(1,length(D_individual{1})),ones(1,length(D_individual{2})),2*ones(1,length(D_individual{3})),3*ones(1,length(D_individual{4})),4*ones(1,length(D_individual{5})),5*ones(1,length(D_individual{6}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9))});
title('Difusion coeficients')
% ylim([-0.02 0.6]);

end


if ntypes ==7;
         
C = [D_individual{1}' D_individual{2}' D_individual{3}' D_individual{4}' D_individual{5}' D_individual{6}' D_individual{7}'];
grp = [zeros(1,length(D_individual{1})),ones(1,length(D_individual{2})),2*ones(1,length(D_individual{3})),3*ones(1,length(D_individual{4})),4*ones(1,length(D_individual{5})),5*ones(1,length(D_individual{6})),6*ones(1,length(D_individual{7}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9))});
title('Difusion coeficients')
% ylim([-0.02 0.6]);
end
 

%% --------------------------------------------------------------
%----------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Compare the different conditions throught a Kolmogorov-Smirnof Statistical Test

% if ntypes == 2;
% [h_ks,p_ks] = kstest2(Dcoeffs{1},Dcoeffs{2});
% [h_ttest,p_ttest] = ttest2(Dcoeffs{1},Dcoeffs{2});
% 
% 
% fprintf('\n');
% fprintf('----------------------------------------------');
% fprintf('\n');
% fprintf('**Kolmogorov-Smirnof Test comparing both populations of Diffusion Coefficients:\n')
% fprintf('h = %.3g ,  p = %d\n', ...
%     h_ks, p_ks);
% 
% fprintf('\n');
% fprintf('h=1 rejects the null hypothesis that the two samples are drawn from the same population\n')
% fprintf('----------------------------------------------');
% fprintf('\n');
% 
% fprintf('\n');
% fprintf('----------------------------------------------');
% fprintf('\n');
% fprintf('**t-Test comparing both populations of Diffusion Coefficients:\n')
% fprintf('h = %.3g ,  p = %d\n', ...
%     h_ttest, p_ttest);
% 
% fprintf('\n');
% fprintf('h=1 rejects the null hypothesis that the two samples are drawn from the same population\n')
% fprintf('----------------------------------------------');
% fprintf('\n');
% end


%% --------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%Analysis of motion type by fitting an anomalous diffusion model (gammas and alphas)
% 
% % figure()
% for t=1:ntypes;
%     
% fprintf('\n');
% fprintf('----------------------------------------------');
% fprintf('\n');
% fprintf('Anomalous Diffusion Analysis');
% fprintf('\n');
% fprintf(files{t}.dataname(1:end-9));
% fprintf('\n');   
%     
% 
% ma = MSD_Results{t}.ma.LogTMSD(TLOGLOG_fitting_points);
% r2fits = MSD_Results{t}.ma.loglogfit.r2fit;
% alphas = MSD_Results{t}.ma.loglogfit.alpha;
% 
% % Remove bad fits
% bad_fits = r2fits < R2LIMIT;
% fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
% alphas(bad_fits) = [];
% %Remove NaNs
% alphas(find(isnan(alphas))) = [];
% 
% gammas = MSD_Results{t}.ma.loglogfit.gamma;
% gammas(bad_fits) = []; % discard bad fits
% gammas(find(isnan(gammas))) = []; %Remove NaNs
% 
% Dmean = mean( gammas ) / 2 / MSD_Results{t}.ma.n_dim;
% Dstd  =  std( gammas ) / 2 / MSD_Results{t}.ma.n_dim;
% 
% fprintf('Estimation of the diffusion coefficient from log-log fit of the MSD curves (Anomalous diffusion):\n')
% fprintf('D = %.2e ± %.2e (mean ± std, N = %d)\n', ...
%     Dmean, Dstd, numel(gammas));
% fprintf('\n');
% fprintf('----------------------------------------------');
% fprintf('\n');


% msmsd = MSD_Results{t}.ma.getMeanMSD;
% gr(t) = plot(msmsd(:,1), msmsd(:,2),'Color',Color(t,:),'LineWidth', 2);
% hold on
% time=0:0.01:tmax;
% f3 = (median(gammas))*(time.^(median(alphas)));
% plot(time,(median(gammas))*(time.^(median(alphas))));
% if ntypes == 1;
%  l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))));   
% end
% 
%  if ntypes == 2;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))));
%     end
%     if ntypes == 3;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))));
%     end
%     if ntypes ==4;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))));
%     end
%    
%    if ntypes ==5;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))));
%    end
%     
%      if ntypes ==6;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))));
%     end
%     
%  if ntypes ==7;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))));
%  end
%     
%      if ntypes ==8;
%     l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))));
%     end    
%   
% 
%     set(l, 'Interpreter', 'none');
%     
%     xlabel('Time lag (seconds)');
%     ylabel('MSD (µm^2)');
    
%     xlim([0 tmax])
%     ylim([0 ymax]);
% end
%Plot the Empirical Survival curve of the Tracks Normalized to ymax
% scatter(x1,(1-f1)*ymax);



%% --------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%Plot MSD of all Tracks
% for t=1:ntypes;
% 
% figure()
% MSD_Results{t}.ma.plotTracks;
% MSD_Results{t}.ma.labelPlotTracks;
% jo=legend(strcat(files{t}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(t))));
% set(jo, 'Interpreter', 'none');
% 
% figure()
% MSD_Results{t}.ma.plotMSD
% joo=legend(strcat(files{t}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(t))));
% set(joo, 'Interpreter', 'none');
% 
% end


%Plot and Analyse Residence Times
% figure()
% for t=1:ntypes;
%     Res_times = Residence_Times(:,t);
%     Res_times(Res_times==0) = [];
% 
%    binCtrs = 0:0.2:1000; %Bin centers, depends on your data
%    n=length(Res_times);
%    counts = hist(Res_times,binCtrs);
%    prob = counts / (n * 0.2);
%    H = bar(binCtrs,prob,'hist');
%    set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
%    hold on;
%    xlim([0 100]);
% 
% mean_Res_Times(t) = mean(Res_times);
% median_Res_Times(t) = median(Res_times);
% 
% hold on
% [f1,x1] = ecdf(Res_times);
% scatter(x1,(1-f1));
% 
%    
% clear Res_times
% 
% 
% 
% end


%% -------------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Plot and analyse the Cumulative Density distribution of Tracks
figure()
for t=1:ntypes;
    Res_times = Residence_Times(:,t);
    Res_times(Res_times==0) = [];

[f1,x1] = ecdf(Res_times);
f2 = fit(x1,(1-f1),'exp2');
x1=x1(2:end);
f1=f1(2:end);
gr(t) = scatter(x1,(1-f1)*Num_Tracks(t),'MarkerFaceColor',Color(t,:));
ylim([0 150]);
xlim([0 20]);
hold on
% plot(f2);
% hold on
   

if ntypes == 1;
 l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))));   
end

 if ntypes == 2;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))));
    end
    if ntypes == 3;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))));
    end
    if ntypes ==4;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))));
    end
    
     
   if ntypes ==5;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))));
   end
    
     if ntypes ==6;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))));
    end
    
 if ntypes ==7;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))));
 end
    
     if ntypes ==8;
    l=legend([gr],strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))));
    end

    set(l, 'Interpreter', 'none');
    
    xlabel('Time (seconds)','FontSize',20,'FontWeight','bold');
    ylabel('Absolute Survival Tracks','FontSize',20,'FontWeight','bold');
    title('Absolute Survival Tracks');
    

    clear Res_times

end


%% ------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Compare the TE-MSD with errorbars for looking at the radius of confinement
figure()
for t=1:ntypes;

mmsd = MSD_Results{t}.ma.getMeanMSD;

temps = mmsd(:,1);
xs = mmsd(:,2);
%plot(temps,xs,'Color',Color(t,:),'LineWidth', 1);
%hold on
dx_plot = mmsd(:,3) ./ sqrt(mmsd(:,4));

ha{t} = errorbar(temps(2:1:end), xs(2:1:end), dx_plot(2:1:end),'o','MarkerSize',5,'Color',Color(t,:),'MarkerFaceColor',Color(t,:),'Linewidth',1.4);

  xlabel('Time (seconds)','FontSize',20,'FontWeight','bold');
  ylabel('MSD (µm^2)','FontSize',20,'FontWeight','bold');
  title('TE-MSD Curves','FontSize',24,'FontWeight','bold'); 
  set(gca,'FontSize',20,'FontWeight','bold');
  
  if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model','Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model','Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Circle Confined Difussion Model','Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Circle Confined Difussion Model',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Circle Confined Difussion Model','Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Circle Confined Difussion Model',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Circle Confined Difussion Model',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Circle Confined Difussion Model','Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Circle Confined Difussion Model',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Circle Confined Difussion Model',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Circle Confined Difussion Model',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Circle Confined Difussion Model','Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Circle Confined Difussion Model',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Circle Confined Difussion Model',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Circle Confined Difussion Model',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Circle Confined Difussion Model',strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Circle Confined Difussion Model','Interpreter','none');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Circle Confined Difussion Model',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Circle Confined Difussion Model',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Circle Confined Difussion Model',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Circle Confined Difussion Model',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Circle Confined Difussion Model',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Circle Confined Difussion Model',strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Circle Confined Difussion Model',strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Circle Confined Difussion Model','Interpreter','none');
     end
     
set(l, 'Interpreter', 'none')  
set(l, 'FontSize',16) 

  

%% Fitting the confined difusion model to the TE-MSD-----------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
[Val_fitted{t}, residual{t},jacobian{t},resnorm{t}] = confined_diffusion_fit_V2(mmsd(:,:),level,num_points,D0,R0);
conf{t} = nlparci(Val_fitted{t},residual{t},'jacobian',jacobian{t});
f_confined_circ_diff_D{t} = Val_fitted{t}(2);
f_confined_circ_diff_R{t} = abs(Val_fitted{t}(1));
f_confined_circ_diff_offset{t} = Val_fitted{t}(3);
Localization_precision_estimation{t} = sqrt(f_confined_circ_diff_offset{t}/(n_dim^2));

fprintf('\n');
fprintf('----------------------------------------------');
fprintf('\n');
fprintf('Circle Confined Model Diffusion Analysis');
fprintf('\n');
fprintf(files{t}.dataname(1:end-9));
fprintf('\n');   

fprintf('\n')
fprintf('**Estimation of the diffusion coefficient from circle confined diffusion model fitted to the TE-MSD)**:\n')
fprintf('D = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n',Val_fitted{t}(2), conf{t}(2,1), conf{t}(2,2));
fprintf('\n')

fprintf('\n')
fprintf('**Estimation of the Radius of confinement from circle confined diffusion model fitted to the TE-MSD)**:\n')
fprintf('R = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n',abs(Val_fitted{t}(1)), abs(conf{t}(1,1)), abs(conf{t}(1,2)));
fprintf('\n')

fprintf('\n')
fprintf('**Localization precision estimation:**:\n')
fprintf('Loc_precision = %.3f nm.\n',Localization_precision_estimation{t}*1000);
fprintf('\n')


xlim([0 num_points*Frame_interval+Frame_interval/2]);
ylim([0 0.03]);
drawnow
errorbar_tick(ha{t},20);
hold on
msd_confined_circle_diffusion = @(X,Xdata)(X(1)^2)*(1-exp(-4*X(2)*Xdata/(X(1)^2))) + X(3);
times = linspace(0,20,100000); %From 0 to 20 seconds. Just for plotting
plot(times,msd_confined_circle_diffusion(Val_fitted{t},times),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
% hold on
% plot(times,msd_confined_circle_diffusion(conf{t}(:,1)',times),'color',Color(t,:),'LineWidth',1,'LineStyle','-');
% hold on
% plot(times,msd_confined_circle_diffusion(conf{t}(:,2)',times),'color',Color(t,:),'LineWidth',1,'LineStyle','-');

end

clear ha

%% Radius of confinement comparison---------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Plot Radius of confinement for the different conditions with 95% CI
%From Confined Circle Difusion Model
figure()
for t=1:ntypes;

    

ha{t} = errorbar(t, f_confined_circ_diff_R{t}*1000,(f_confined_circ_diff_R{t}-conf{t}(1,1))*1000, (conf{t}(1,2)-f_confined_circ_diff_R{t})*1000,'o','MarkerSize',5,'MarkerFaceColor',Color(t,:),'Linewidth',1.4);
% xlim([0 0.1]);
%ylim([-0.02 0.05]);

%   xlabel('');
%   ylabel('');
%   title(''); 
  
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
     end
 
title('Radius of confinement (nm) with 95% Confidence Intervals','FontSize',24,'FontWeight','bold'); 
set(l, 'Interpreter', 'none')  
set(l, 'FontSize',16) 
set(gca,'FontSize',20,'FontWeight','bold');
set(gca, 'XTick', 0:ntypes)
ylabel('Radius of confinement (nm)','FontSize',20,'FontWeight','bold');
hold on
drawnow
errorbar_tick(ha{t},50);    

    
end

clear ha


%% Diffusion Coefficient Comparison-----------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Plot Diffusion Coefficient for the different conditions with 95% CI from
%Confined Circle Diffusion Model
figure()
for t=1:ntypes;


ha{t} = errorbar(t, f_confined_circ_diff_D{t},f_confined_circ_diff_D{t}-conf{t}(2,1), conf{t}(2,2)-f_confined_circ_diff_D{t},'o','MarkerSize',5,'MarkerFaceColor',Color(t,:),'Linewidth',1.4);
% xlim([0 0.1]);
%ylim([-0.02 0.05]);

%   xlabel('');
%   ylabel('');
%   title(''); 
  
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
     end
 
title('Diffusion Coefficients (um2/s) with 95% Confidence Intervals','FontSize',24,'FontWeight','bold'); 
set(l, 'Interpreter', 'none')  
set(l, 'FontSize',16) 
set(gca,'FontSize',20,'FontWeight','bold');
set(gca, 'XTick', 0:ntypes)
ylabel('Diffusion Coefficient (um2/s)','FontSize',20,'FontWeight','bold');
hold on
drawnow
errorbar_tick(ha{t},50);    

    
end

clear ha



%% - Calculate the distance traveled by Directed Tracks
%This would be something equivalent to the Radius of confimenent for
%Confined Tracks
figure()
for t=1:ntypes;

for ggg = 1:size(MSD_Results{t}.ma.tracks,1);
points_coord = MSD_Results{t}.ma.tracks{ggg}(:,2:3);
[max_dist{t}(ggg), min_dist{t}(ggg), avg_dist{t}(ggg)] = distance_scatter(points_coord);
end


binCtrs = 50:10:700; %Bin centers, depends on your data
   n=length(max_dist{t});
   counts = hist(max_dist{t}*1000,binCtrs);
   prob = counts / (n * 10);
   H = bar(binCtrs,prob,'hist');
   set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
   hold on;
   
xlabel('Maximum Distance travelled by the Track (nm)','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontSize',20,'FontWeight','bold');
title('Histogram of Distance Travelled by each Track','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlim([0 500]);

end

figure()
%Boxplot comparing the maximum distances travelled by each track for each
%condition.

 if ntypes == 2;    
C = [max_dist{1} max_dist{2}];
grp = [zeros(1,length(max_dist{1})),ones(1,length(max_dist{2}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9))});
title('Maximum distance traveled by the molecule(um)','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Distance (um)','FontSize',20,'FontWeight','bold');
% ylim([-0.02 0.6]);
 end
 
 if ntypes == 3;    

C = [max_dist{1} max_dist{2} max_dist{3}];
grp = [zeros(1,length(max_dist{1})),ones(1,length(max_dist{2})),2*ones(1,length(max_dist{3}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9))});
title('Maximum distance traveled by the molecule(um)','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Distance (um)','FontSize',20,'FontWeight','bold');
% ylim([-0.02 0.6]);

 end
 
 if ntypes ==4;    
     
C = [max_dist{1} max_dist{2} max_dist{3} max_dist{4}];
grp = [zeros(1,length(max_dist{1})),ones(1,length(max_dist{2})),2*ones(1,length(max_dist{3})),3*ones(1,length(max_dist{4}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9))});
title('Maximum distance traveled by the molecule(um)','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Distance (um)','FontSize',20,'FontWeight','bold');
% ylim([-0.02 0.6]);

 end

  
 if ntypes ==5;
     
C = [max_dist{1} max_dist{2} max_dist{3} max_dist{4} max_dist{5}];
grp = [zeros(1,length(max_dist{1})),ones(1,length(max_dist{2})),2*ones(1,length(max_dist{3})),3*ones(1,length(max_dist{4})),4*ones(1,length(max_dist{5}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9))});
title('Maximum distance traveled by the molecule(um)','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Distance (um)','FontSize',20,'FontWeight','bold');
% ylim([-0.02 0.6]);

 end

 
 if ntypes ==6;

C = [max_dist{1} max_dist{2} max_dist{3} max_dist{4} max_dist{5} max_dist{6}];
grp = [zeros(1,length(max_dist{1})),ones(1,length(max_dist{2})),2*ones(1,length(max_dist{3})),3*ones(1,length(max_dist{4})),4*ones(1,length(max_dist{5})),5*ones(1,length(max_dist{6}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9))});
title('Maximum distance traveled by the molecule(um)','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Distance (um)','FontSize',20,'FontWeight','bold');
% ylim([-0.02 0.6]);

end


if ntypes ==7;
         
C = [max_dist{1} max_dist{2} max_dist{3} max_dist{4} max_dist{5} max_dist{6} max_dist{7}];
grp = [zeros(1,length(max_dist{1})),ones(1,length(max_dist{2})),2*ones(1,length(max_dist{3})),3*ones(1,length(max_dist{4})),4*ones(1,length(max_dist{5})),5*ones(1,length(max_dist{6})),6*ones(1,length(max_dist{7}))];
boxplot(C,grp,'Labels',{strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9))});
title('Maximum distance traveled by the molecule(um)','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Distance (um)','FontSize',20,'FontWeight','bold');
% ylim([-0.02 0.6]);
end









%% Fitting the confined difusion model to each trajectory (T-MSD)-----------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure()
strB = '';
for t=1:ntypes;
fprintf('\n')
fprintf('Circle Confined Diffusion Model')
fprintf('\n')
fprintf('-----------------------------------------------------')
fprintf('\n')
fprintf('fitting each individual T-MSD Curve')
fprintf('\n')
fprintf('\n')
fprintf(strcat('Step_',num2str(t),'_of_',num2str(ntypes)))
fprintf('\n')   


for i=1:size(MSD_Results{t}.ma.msd,1);
    
local_msd = MSD_Results{t}.ma.msd{i};
nonnan = ~isnan(local_msd(:,2));  
local_msd = local_msd(nonnan,1:2);

[Val_fitted_all{t}{i}, residual_all{t}{i},jacobian_all{t}{i},resnorm_all{t}{i}] = confined_diffusion_fit_V2(local_msd,level,num_points,D0,R0);
conf_all{t}{i} = nlparci(Val_fitted_all{t}{i},residual_all{t}{i},'jacobian',jacobian_all{t}{i});
f_confined_circ_diff_D_all{t}(i) = Val_fitted_all{t}{i}(2);
f_confined_circ_diff_R_all{t}(i) = abs(Val_fitted_all{t}{i}(1));
f_confined_circ_diff_offset_all{t}(i) = Val_fitted_all{t}{i}(3);
Localization_precision_estimation_all{t}(i) = sqrt(f_confined_circ_diff_offset_all{t}(i)/(n_dim^2));

Track_length_radius{t}(i) = size(local_msd,1);


% fprintf('\n');
% fprintf('----------------------------------------------');
% fprintf('\n');
% fprintf('Circle Confined Model Diffusion Analysis');
% fprintf('\n');
% fprintf(files{t}.dataname(1:end-9));
% fprintf('\n');   
% 
% fprintf('\n')
% fprintf('**Estimation of the diffusion coefficient from circle confined diffusion model fitted to the TE-MSD)**:\n')
% fprintf('D = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n',Val_fitted{t}(2), conf{t}(2,1), conf{t}(2,2));
% fprintf('\n')
% 
% fprintf('\n')
% fprintf('**Estimation of the Radius of confinement from circle confined diffusion model fitted to the TE-MSD)**:\n')
% fprintf('R = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n',abs(Val_fitted{t}(1)), abs(conf{t}(1,1)), abs(conf{t}(1,2)));
% fprintf('\n')

% xlim([0 num_points*Frame_interval+Frame_interval/2]);
% % ylim([0 0.01]);
% drawnow
% errorbar_tick(ha{t},20);
% hold on
% msd_confined_circle_diffusion = @(X,Xdata)(X(1)^2)*(1-exp(-4*X(2)*Xdata/(X(1)^2))) + X(3);
% times = linspace(0,20,100000); %From 0 to 20 seconds. Just for plotting
% plot(times,msd_confined_circle_diffusion(Val_fitted{t},times),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');

  strOut = sprintf('Completed: % 4.1f',100*i/(size(MSD_Results{t}.ma.msd,1)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);
end

fprintf('\n')

%Filter. Just take the good fits. 
iidx{t} = find(cell2mat(resnorm_all{t}) < 0.00001);
bin_size = 0.05;
binCtrs = 0:bin_size:4; %Bin centers. (depends on your data)
n=length(f_confined_circ_diff_R_all{t}(iidx{t}));
counts = hist(log10(f_confined_circ_diff_R_all{t}(iidx{t})*1000),binCtrs);

   prob = counts / (n*bin_size);
   H = bar(binCtrs,prob,'hist');
   set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
   hold on;
   
xlabel('Log10(Radius of confinement) (nm)','FontSize',20,'FontWeight','bold');
% xlabel('Radius of confinement (nm)','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontSize',20,'FontWeight','bold');
title('Histogram of Radius of confinement','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlim([0 3]);


%Fit a Normal distribution to each condition
 %------------------------------------------------------------------
 %Test Bimodality of the distribution:
% probability_unimodal_Radius(t) = test_bimodality(log10(f_confined_circ_diff_R_all{t}(iidx)*1000));
%  
% %Fit a Normal distribution
% 
% % if probability_unimodal(t) > 0.9;
% %    pd(t) = fitdist(log10(f_confined_circ_diff_R_all{t}(iidx)*1000),'normal');
% %    y = pdf(pd(t),binCtrs);
% %    plot(binCtrs,y,'LineWidth',2,'Color',Color(t,:));
% %    xlim([bin_start bin_end]);
% %    mu_normal(t) = pd(t).mu;
% %    sigma_normal(t) = pd(t).sigma;
% %    D_est(t) = exp(mu_normal(t));
% %    D_sigma_est(t) = exp(mu_normal(t)+sigma_normal(t))-exp(mu_normal(t));
% %    
% % else
% 
% %Fit a Bimodal Normal distribution   
% %Define the function. (5 parameters)
%                   
% pStart = 0.5;
% muStart = quantile(log10(f_confined_circ_diff_R_all{t}(iidx)*1000),[0.25 0.75]);
% sigmaStart = sqrt(var(log10(f_confined_circ_diff_R_all{t}(iidx)*1000)) - 0.25*diff(muStart).^2);
% start = [pStart muStart sigmaStart sigmaStart];
% 
% lb = [0 0 0 0 0];
% ub = [1 10 10 Inf Inf];
% %%%pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2)p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2)
% paramEsts_radius{t} = mle(log10(f_confined_circ_diff_R_all{t}(iidx)*1000), 'pdf',pdf_normmixture, 'start',start,'lower',lb, 'upper',ub);
% 
% xgrid = linspace(0.9*min(log10(f_confined_circ_diff_R_all{t}(iidx)*1000)),1.1*max(log10(f_confined_circ_diff_R_all{t}(iidx)*1000)),200);
% pdfgrid_radius{t} = pdf_normmixture(xgrid,paramEsts_radius{t}(1),paramEsts_radius{t}(2),paramEsts_radius{t}(3),paramEsts_radius{t}(4),paramEsts_radius{t}(5));
% hold on
% plot(xgrid,pdfgrid_radius{t},'-','LineWidth',2,'Color',Color(t,:));
% xlim([0 3.5]);


fprintf('\n')
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end
    
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end

end

hold on
plot([log10(mean(Localization_precision_estimation_all{1}(iidx{1})*1000)),log10(mean(Localization_precision_estimation_all{1}(iidx{1})*1000))],[0,3],'LineWidth',4,'Color',[0 0 0]);

fprintf('\n')
fprintf('**Localization precision estimation calculated from all trajectories:**:\n')
fprintf('Loc_precision = %.3f nm ± %.3f nm.\n',mean(Localization_precision_estimation_all{t}(iidx{t})*1000),std(Localization_precision_estimation_all{t}(iidx{t})*1000));
fprintf('\n')

%% Plot the radius of confinement histogram

figure()
for t=1:ntypes;
     
bin_size = 5;
binCtrs = 0:bin_size:1000; %Bin centers. (depends on your data)
n=length(f_confined_circ_diff_R_all{t}(iidx{t}));
counts = hist((f_confined_circ_diff_R_all{t}(iidx{t})*1000),binCtrs);

   prob = counts / (n*bin_size);
   H = bar(binCtrs,prob,'hist');
   set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
   hold on;
   
xlabel('Radius of confinement (nm)','FontSize',20,'FontWeight','bold');
% xlabel('Radius of confinement (nm)','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontSize',20,'FontWeight','bold');
title('Histogram of Radius of confinement','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlim([0 300]);


%Fit a Normal distribution to each condition
 %------------------------------------------------------------------
 %Test Bimodality of the distribution:
% probability_unimodal_Radius(t) = test_bimodality(log10(f_confined_circ_diff_R_all{t}(iidx)*1000));
%  
% %Fit a Normal distribution
% 
% % if probability_unimodal(t) > 0.9;
% %    pd(t) = fitdist(log10(f_confined_circ_diff_R_all{t}(iidx)*1000),'normal');
% %    y = pdf(pd(t),binCtrs);
% %    plot(binCtrs,y,'LineWidth',2,'Color',Color(t,:));
% %    xlim([bin_start bin_end]);
% %    mu_normal(t) = pd(t).mu;
% %    sigma_normal(t) = pd(t).sigma;
% %    D_est(t) = exp(mu_normal(t));
% %    D_sigma_est(t) = exp(mu_normal(t)+sigma_normal(t))-exp(mu_normal(t));
% %    
% % else
% 
% %Fit a Bimodal Normal distribution   
% %Define the function. (5 parameters)
%                   
% pStart = 0.5;
% muStart = quantile(log10(f_confined_circ_diff_R_all{t}(iidx)*1000),[0.25 0.75]);
% sigmaStart = sqrt(var(log10(f_confined_circ_diff_R_all{t}(iidx)*1000)) - 0.25*diff(muStart).^2);
% start = [pStart muStart sigmaStart sigmaStart];
% 
% lb = [0 0 0 0 0];
% ub = [1 10 10 Inf Inf];
% %%%pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2)p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2)
% paramEsts_radius{t} = mle(log10(f_confined_circ_diff_R_all{t}(iidx)*1000), 'pdf',pdf_normmixture, 'start',start,'lower',lb, 'upper',ub);
% 
% xgrid = linspace(0.9*min(log10(f_confined_circ_diff_R_all{t}(iidx)*1000)),1.1*max(log10(f_confined_circ_diff_R_all{t}(iidx)*1000)),200);
% pdfgrid_radius{t} = pdf_normmixture(xgrid,paramEsts_radius{t}(1),paramEsts_radius{t}(2),paramEsts_radius{t}(3),paramEsts_radius{t}(4),paramEsts_radius{t}(5));
% hold on
% plot(xgrid,pdfgrid_radius{t},'-','LineWidth',2,'Color',Color(t,:));
% xlim([0 3.5]);


fprintf('\n')
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end
    
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end

 
end
 
hold on
plot([mean(Localization_precision_estimation_all{1}(iidx{1})*1000),mean(Localization_precision_estimation_all{1}(iidx{1})*1000)],[0,0.03],'LineWidth',4,'Color',[0 0 0]);

%% Plot Radius of confinement vs Tracklength of each Track to check if there
%is correlation:
figure()
for t=1:ntypes;
scatter(f_confined_circ_diff_R_all{t}(iidx{t})*1000,Track_length_radius{t}(iidx{t}));
hold on
end

xlabel('Radius of confinement (nm)','FontSize',20,'FontWeight','bold');
ylabel('Track length (frames)','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');  
xlim([0 200]);
ylim([0 400]);

if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end
    
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end
 
 
 %% Plot a histogram showing the distribution of tracks lengths
 figure()
 hist(Track_length_radius{t}(iidx{t}),50);
xlabel('Track length (frames)','FontSize',20,'FontWeight','bold');
ylabel('Counts','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');  


 
%% Divide tracks into THREE populations based on their radius of confinement
%Also study if the Tracklength distribution changes between these different
%radius of confinement groups of tracks.
%This makes sense for studying the confined trajectories

%1. High R
for t=1:ntypes;
idxr_high{t} = find(f_confined_circ_diff_R_all{t}*1000 > Radius_interval(2));

figure()
plotTracks(MSD_Results{t}.ma,gca,idxr_high{t});
MSD_Results{t}.ma.labelPlotTracks;

title(strcat('Trayectories with Radius of confinement ',str_rad{1}),'FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
set(gca,'Ydir','reverse');
l = legend(strcat(files{t}.dataname(1:end-9))); 
set(l, 'Interpreter', 'none')

for i=1:length(idxr_high{t});
Track_length_radius_interval{t}{1}(i) = size(MSD_Results{t}.ma.tracks{idxr_high{t}(i)},1);
end

Trajectories_High_R{t} = MSD_Results{t}.ma.tracks(idxr_high{t});

end


%2. Mid R
for t=1:ntypes;
idxr_between{t} = find(f_confined_circ_diff_R_all{t}*1000 >= Radius_interval(1) & f_confined_circ_diff_R_all{t}*1000 <= Radius_interval(2));

figure()
plotTracks(MSD_Results{t}.ma,gca,idxr_between{t});
MSD_Results{t}.ma.labelPlotTracks;

title(strcat('Trayectories with Radius of confinement ',str_rad{2}),'FontSize',24,'FontWeight','bold'); 
l = legend(strcat(files{t}.dataname(1:end-9))); 
set(l, 'Interpreter', 'none')

for i=1:length(idxr_between{t});
Track_length_radius_interval{t}{2}(i) = size(MSD_Results{t}.ma.tracks{idxr_between{t}(i)},1);
end


Trajectories_Mid_R{t} = MSD_Results{t}.ma.tracks(idxr_between{t});
end


%3. Low R
for t=1:ntypes;
idxr_low{t} = find(f_confined_circ_diff_R_all{t}*1000 < Radius_interval(1));

figure()
plotTracks(MSD_Results{t}.ma,gca,idxr_low{t});
MSD_Results{t}.ma.labelPlotTracks;

title(strcat('Trayectories with Radius of confinement ',str_rad{3}),'FontSize',24,'FontWeight','bold'); 
l = legend(strcat(files{t}.dataname(1:end-9))); 
set(l, 'Interpreter', 'none')

for i=1:length(idxr_low{t});
Track_length_radius_interval{t}{3}(i) = size(MSD_Results{t}.ma.tracks{idxr_low{t}(i)},1);
end

Trajectories_Low_R{t} = MSD_Results{t}.ma.tracks(idxr_low{t});

end

%% Calculate the distance of the tracks to the center of mass of each cell

%1. High R
for t=1:ntypes;
   
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      Cell_Reference_High_R{t} =  Cell_Reference{t}(idxr_high{t});
      idx_cell = find(cell2mat(Cell_Reference_High_R{t}) == i);
       
      Tracks_Cell_High_R{t}{i} = Trajectories_High_R{t}(idx_cell);
      
          for j = 1:size(Tracks_Cell_High_R{t}{i},1);
          Tracks_CenterOfMass_High_R{t}{i}(j,1:2) = mean(Tracks_Cell_High_R{t}{i}{j}(:,2:3));
          
          Track_Distance_to_CellCenter_High_R{t}{i}(j) = sqrt( (Tracks_CenterOfMass_High_R{t}{i}(j,1) - Cell_CenterOfMass{t}{1}{i}(1))^2 + (Tracks_CenterOfMass_High_R{t}{i}(j,2) - Cell_CenterOfMass{t}{1}{i}(2))^2);

          end
      
        

  end
  
  Track_Distance_to_CellCenter_High_R_all{t} = [Track_Distance_to_CellCenter_High_R{t}{:}];
  
end


%1. Mid R
for t=1:ntypes;
   
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      Cell_Reference_Mid_R{t} =  Cell_Reference{t}(idxr_between{t});
      idx_cell = find(cell2mat(Cell_Reference_Mid_R{t}) == i);
       
      Tracks_Cell_Mid_R{t}{i} = Trajectories_Mid_R{t}(idx_cell);
      
          for j = 1:size(Tracks_Cell_Mid_R{t}{i},1);
          Tracks_CenterOfMass_Mid_R{t}{i}(j,1:2) = mean(Tracks_Cell_Mid_R{t}{i}{j}(:,2:3));
          
          Track_Distance_to_CellCenter_Mid_R{t}{i}(j) = sqrt( (Tracks_CenterOfMass_Mid_R{t}{i}(j,1) - Cell_CenterOfMass{t}{1}{i}(1))^2 + (Tracks_CenterOfMass_Mid_R{t}{i}(j,2) - Cell_CenterOfMass{t}{1}{i}(2))^2);

          end
      
        

  end
  
  
    Track_Distance_to_CellCenter_Mid_R_all{t} = [Track_Distance_to_CellCenter_Mid_R{t}{:}];

end


%1. Low R
for t=1:ntypes;
   
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      Cell_Reference_Low_R{t} =  Cell_Reference{t}(idxr_low{t});
      idx_cell = find(cell2mat(Cell_Reference_Low_R{t}) == i);
       
      Tracks_Cell_Low_R{t}{i} = Trajectories_Low_R{t}(idx_cell);
      
          for j = 1:size(Tracks_Cell_Low_R{t}{i},1);
          Tracks_CenterOfMass_Low_R{t}{i}(j,1:2) = mean(Tracks_Cell_Low_R{t}{i}{j}(:,2:3));
          
          Track_Distance_to_CellCenter_Low_R{t}{i}(j) = sqrt( (Tracks_CenterOfMass_Low_R{t}{i}(j,1) - Cell_CenterOfMass{t}{1}{i}(1))^2 + (Tracks_CenterOfMass_Low_R{t}{i}(j,2) - Cell_CenterOfMass{t}{1}{i}(2))^2);

          end
      
        

  end
  
    Track_Distance_to_CellCenter_Low_R_all{t} = [Track_Distance_to_CellCenter_Low_R{t}{:}];
  
end


%Plot the Results
for t=1:ntypes;
   
   figure()
        Color_R = lines(3);
         
        bin_size = 0.5;
        binCtrs = 0:bin_size:16; %Bin centers. (depends on your data)
        n=length(Track_Distance_to_CellCenter_High_R_all{t});
        counts = hist(Track_Distance_to_CellCenter_High_R_all{t},binCtrs); %In um

        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(1,:),'FaceAlpha',0.6);
        hold on;
       
        n=length(Track_Distance_to_CellCenter_Mid_R_all{t});
        counts = hist(Track_Distance_to_CellCenter_Mid_R_all{t},binCtrs); %In um

        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(2,:),'FaceAlpha',0.6);
        hold on;
        
         n=length(Track_Distance_to_CellCenter_Low_R_all{t});
        counts = hist(Track_Distance_to_CellCenter_Low_R_all{t},binCtrs); %In um

        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(3,:),'FaceAlpha',0.6);

        xlabel('Distance (um)','FontSize',20,'FontWeight','bold');
        ylabel('Frequency','FontSize',20,'FontWeight','bold');
        title('Histogram of Track distance to the center of the cell','FontSize',24,'FontWeight','bold');
        legend('High Radius Conf','Medium Radius Conf','Low Radius Conf');
        set(gca,'FontSize',20,'FontWeight','bold');  
       

    
     
end

%% Calculate the NND of the center of mass of each track splitted into the Radius of confinement intervals

%1. High R
for t=1:ntypes;
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      for j = 1:size(Tracks_Cell_High_R{t}{i},1);
          if size(Tracks_Cell_High_R{t}{i},1)>1;
      NND_High_R{t}{i}(j) = real(min(nonzeros(pdist2(Tracks_CenterOfMass_High_R{t}{i}(j,1:2),Tracks_CenterOfMass_High_R{t}{i},'euclidean'))));
          else
       NND_High_R{t}{i}(j) = real(min(pdist2(Tracks_CenterOfMass_High_R{t}{i}(j,1:2),Tracks_CenterOfMass_High_R{t}{i},'euclidean')));
          end
      end
      
  end
  
    NND_High_R_all{t} = [NND_High_R{t}{:}];
    NND_High_R_all{t} = nonzeros(NND_High_R_all{t});

end

%2. Mid R
for t=1:ntypes;
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      for j = 1:size(Tracks_Cell_Mid_R{t}{i},1);
          if size(Tracks_Cell_Mid_R{t}{i},1)>1;
      NND_Mid_R{t}{i}(j) = real(min(nonzeros(pdist2(Tracks_CenterOfMass_Mid_R{t}{i}(j,1:2),Tracks_CenterOfMass_Mid_R{t}{i},'euclidean'))));
          else
       NND_Mid_R{t}{i}(j) = real(min(pdist2(Tracks_CenterOfMass_Mid_R{t}{i}(j,1:2),Tracks_CenterOfMass_Mid_R{t}{i},'euclidean'))); 
          end
      end
      
  end
    
    NND_Mid_R_all{t} = [NND_Mid_R{t}{:}];
    NND_Mid_R_all{t} = nonzeros(NND_Mid_R_all{t});
end

%3. Low R
for t=1:ntypes;
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      for j = 1:size(Tracks_Cell_Low_R{t}{i},1);
          if size(Tracks_Cell_Low_R{t}{i},1)>1;
      NND_Low_R{t}{i}(j) = real(min(nonzeros(pdist2(Tracks_CenterOfMass_Low_R{t}{i}(j,1:2),Tracks_CenterOfMass_Low_R{t}{i},'euclidean'))));
          else
       NND_Low_R{t}{i}(j) = real(min(pdist2(Tracks_CenterOfMass_Low_R{t}{i}(j,1:2),Tracks_CenterOfMass_Low_R{t}{i},'euclidean')));
          end
      end
      
  end
    
    NND_Low_R_all{t} = [NND_Low_R{t}{:}];
    NND_Low_R_all{t} = nonzeros(NND_Low_R_all{t});
end


%Plot the Results
for t=1:ntypes;
   
   figure()         
        bin_size = 500;
        binCtrs = 0:bin_size:17000; %Bin centers. (depends on your data)
        n=length(NND_High_R_all{t});
        counts = hist(NND_High_R_all{t}*1000,binCtrs); %In nm

        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(1,:),'FaceAlpha',0.6);
        hold on;
       
        n=length(NND_Mid_R_all{t});
        counts = hist(NND_Mid_R_all{t}*1000,binCtrs); %In nm

        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(2,:),'FaceAlpha',0.6);
        hold on;
        
         n=length(NND_Low_R_all{t});
        counts = hist(NND_Low_R_all{t}*1000,binCtrs); %In nm

        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(3,:),'FaceAlpha',0.6);

        xlabel('NND (nm)','FontSize',20,'FontWeight','bold');
        ylabel('Frequency','FontSize',20,'FontWeight','bold');
        title('Histogram of NND between Tracks','FontSize',24,'FontWeight','bold');
        legend('High Radius Conf','Medium Radius Conf','Low Radius Conf');
        set(gca,'FontSize',20,'FontWeight','bold');  
        xlim([-300 18000]);

    
     
end

%% Calculate the NND of the center of mass of each track between Tracks of different Radius of confinement intervals

%1. Between High-R and Mid-R
for t=1:ntypes;
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      for j = 1:size(Tracks_Cell_High_R{t}{i},1);
          if size(Tracks_Cell_Mid_R{t}{i},1)>1;
      NND_High_Mid_R{t}{i}(j) = real(min(nonzeros(pdist2(Tracks_CenterOfMass_High_R{t}{i}(j,1:2),Tracks_CenterOfMass_Mid_R{t}{i},'euclidean'))));
          else if size(Tracks_Cell_Mid_R{t}{i},1)>0;
      NND_High_Mid_R{t}{i}(j) = real(min(pdist2(Tracks_CenterOfMass_High_R{t}{i}(j,1:2),Tracks_CenterOfMass_Mid_R{t}{i},'euclidean')));
              else
      NND_High_Mid_R{t}{i}(j) = 0;           
              end
          end
      end
      
  end
  
    NND_High_Mid_R_all{t} = [NND_High_Mid_R{t}{:}];
    NND_High_Mid_R_all{t} = nonzeros(NND_High_Mid_R_all{t});

end

%2. Between High-R and Low-R
for t=1:ntypes;
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      for j = 1:size(Tracks_Cell_High_R{t}{i},1);
          if size(Tracks_Cell_Low_R{t}{i},1)>1;
      NND_High_Low_R{t}{i}(j) = real(min(nonzeros(pdist2(Tracks_CenterOfMass_High_R{t}{i}(j,1:2),Tracks_CenterOfMass_Low_R{t}{i},'euclidean'))));
          else if size(Tracks_Cell_Low_R{t}{i},1)>0;
      NND_High_Low_R{t}{i}(j) = real(min(pdist2(Tracks_CenterOfMass_High_R{t}{i}(j,1:2),Tracks_CenterOfMass_Low_R{t}{i},'euclidean'))); 
              else
      NND_High_Low_R{t}{i}(j) = 0; 
              end
          end
      end
      
  end
    
    NND_High_Low_R_all{t} = [NND_High_Low_R{t}{:}];
    NND_High_Low_R_all{t} = nonzeros(NND_High_Low_R_all{t});
end

%3. Between Mid-R and Low-R
for t=1:ntypes;
    %Iterate on each Cell (1 cell per image)
  for i=1:size(Cell_CenterOfMass{t}{1},2);
      
      for j = 1:size(Tracks_Cell_Mid_R{t}{i},1);
          if size(Tracks_Cell_Low_R{t}{i},1)>1;
      NND_Mid_Low_R{t}{i}(j) = real(min(nonzeros(pdist2(Tracks_CenterOfMass_Mid_R{t}{i}(j,1:2),Tracks_CenterOfMass_Low_R{t}{i},'euclidean'))));
          else if size(Tracks_Cell_Low_R{t}{i},1)>0;
      NND_Mid_Low_R{t}{i}(j) = real(min(pdist2(Tracks_CenterOfMass_Mid_R{t}{i}(j,1:2),Tracks_CenterOfMass_Low_R{t}{i},'euclidean')));
              else
      NND_Mid_Low_R{t}{i}(j) = 0;          
              end
          end
      end
      
  end
    
    NND_Mid_Low_R_all{t} = [NND_Mid_Low_R{t}{:}];
    NND_Mid_Low_R_all{t} = nonzeros(NND_Mid_Low_R_all{t});
end


%Plot the Results
for t=1:ntypes;
   
   figure()         
        bin_size = 100;
        binCtrs = 0:bin_size:17000; %Bin centers. (depends on your data)
        
        n=length(NND_High_Mid_R_all{t});
        counts = hist(NND_High_Mid_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(1,:),'FaceAlpha',0.6);
        hold on;
       
        n=length(NND_High_Low_R_all{t});
        counts = hist(NND_High_Low_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(2,:),'FaceAlpha',0.6);
        hold on;
        
        n=length(NND_Mid_Low_R_all{t});
        counts = hist(NND_Mid_Low_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(3,:),'FaceAlpha',0.6);

        xlabel('NND (nm)','FontSize',20,'FontWeight','bold');
        ylabel('Frequency','FontSize',20,'FontWeight','bold');
        title('Histogram of NND between Tracks','FontSize',24,'FontWeight','bold');
        legend('High-Mid Radius Conf','High-Low Radius Conf','Mid-Low Radius Conf');
        set(gca,'FontSize',20,'FontWeight','bold');  
        xlim([-300 6000]);
    
     
end



%Plot the Results
for t=1:ntypes;
   
   figure()         
        bin_size = 100;
        binCtrs = 0:bin_size:17000; %Bin centers. (depends on your data)
        
        n=length(NND_High_R_all{t});
        counts = hist(NND_High_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(1,:),'FaceAlpha',0.6);
        hold on;
       
        n=length(NND_High_Mid_R_all{t});
        counts = hist(NND_High_Mid_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(2,:),'FaceAlpha',0.6);
        hold on;
        
        n=length(NND_High_Low_R_all{t});
        counts = hist(NND_High_Low_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(3,:),'FaceAlpha',0.6);

        xlabel('NND (nm)','FontSize',20,'FontWeight','bold');
        ylabel('Frequency','FontSize',20,'FontWeight','bold');
        title('Histogram of NND between Tracks','FontSize',24,'FontWeight','bold');
        legend('High-High Radius Conf','High-Mid Radius Conf','High-Low Radius Conf');
        set(gca,'FontSize',20,'FontWeight','bold');  
        xlim([-300 6000]);
    
     
end


%Plot the Results
for t=1:ntypes;
   
   figure()         
        bin_size = 500;
        binCtrs = 0:bin_size:17000; %Bin centers. (depends on your data)
        
        n=length(NND_Mid_R_all{t});
        counts = hist(NND_Mid_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(1,:),'FaceAlpha',0.6);
        hold on;
       
        n=length(NND_High_Mid_R_all{t});
        counts = hist(NND_High_Mid_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(2,:),'FaceAlpha',0.6);
        hold on;
        
        n=length(NND_Mid_Low_R_all{t});
        counts = hist(NND_Mid_Low_R_all{t}*1000,binCtrs); %In nm
        prob = counts / (n*bin_size);
        H = bar(binCtrs,prob,'hist');
        set(H,'Facecolor',Color_R(3,:),'FaceAlpha',0.6);

        xlabel('NND (nm)','FontSize',20,'FontWeight','bold');
        ylabel('Frequency','FontSize',20,'FontWeight','bold');
        title('Histogram of NND between Tracks','FontSize',24,'FontWeight','bold');
        legend('Mid-Mid Radius Conf','Mid-High Radius Conf','Mid-Low Radius Conf');
        set(gca,'FontSize',20,'FontWeight','bold');  
        xlim([-300 6000]);
    
     
end


%% Voronoi Analysis of Tracks center of mass splitted into different radius of confinement intervals

%3. Between Mid-R and Low-R
% for t=1:ntypes;
%     %Iterate on each Cell (1 cell per image)
%   for i=1:size(Cell_CenterOfMass{t}{1},2);
%       
%       Tracks_CenterOfMass_High_R{t}{i};
%       Tracks_CenterOfMass_Mid_R{t}{i};
%       Tracks_CenterOfMass_Low_R{t}{i};
% 
%      
%       
%   end
%     
%   
% end





%% Plot an histogram of Residence times for each radius interval
figure() 
for t=1:ntypes;
subplot(ntypes,1,t);   
    for gy=1:3;
Color_rad = lines(3);
binCtrs = 0:0.04:2.5; %Bin centers, depends on your data    
n=length(Track_length_radius_interval{t}{gy});
counts = hist(Track_length_radius_interval{t}{gy}*Frame_interval,binCtrs);
prob = counts / (n*0.05);
H = bar(binCtrs,prob,'hist');
set(H,'Facecolor',Color_rad(gy,:),'FaceAlpha',0.5);
hold on;
    end

xlabel('Residence times (s)');
ylabel('frequency');

l=legend(strcat(files{t}.dataname(1:end-9),'_Radius interval',str_rad{1}),strcat(files{t}.dataname(1:end-9),'_Radius interval',str_rad{2}),strcat(files{t}.dataname(1:end-9),'_Radius interval',str_rad{3}));    
set(l, 'Interpreter', 'none');
xlim([0 2]);
end
suptitle('Histogram of Residence times of the tracks from several Radius of confinement intervals'); 




%Analyze and compare the residence time distribution of the Tracks
%belonging to different radius of confinement range
Residence_Times_radius = {};

for t=1:ntypes;
    for gy=1:length(Radius_interval)+1;
        
    Num_particles_rad{t}{gy} = size(Track_length_radius_interval{t}{gy},2);
    
        for i=1:Num_particles_rad{t}{gy};   
        Residence_Times_radius{t}{gy}(i) = Track_length_radius_interval{t}{gy}(i)*Frame_interval;
        end
    
    end
end

%Fit 
figure()
Esp_coeff_rad = {};
Esp_Sigma_rad = {};
Esp_fit_rad = {};

for t=1:ntypes;
    for gy=1:3;
        
survival_rad{t}{gy} = [];       

max_time{t}{gy} = max(Residence_Times_radius{t}{gy});
min_time{t}{gy} = min(Residence_Times_radius{t}{gy});

        for r_time = min_time{t}{gy}:Frame_interval:max_time{t}{gy};
            
            survival_rad{t}{gy} = [survival_rad{t}{gy}; r_time length(find(Residence_Times_radius{t}{gy} >= r_time))/Num_particles_rad{t}{gy}];
            
        end

        
scatter(survival_rad{t}{gy}(:,1),survival_rad{t}{gy}(:,2),'MarkerFaceColor',Color(t,:));
hold on
ylabel('Normalized survival fraction (Counts)','FontSize',20,'FontWeight','bold');
xlabel('Time (s)','FontSize',20,'FontWeight','bold');
title('Normalized survival fraction for several Rc intervals','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  

%Fit a 2-component exponential decay to the survival probability function
[Esp_coeff_rad{t}{gy},Esp_Sigma_rad{t}{gy}, Esp_fit_rad{t}{gy}] = ExpDecay_2Cmp_fit(survival_rad{t}{gy}, [k1_init k2_init], fraction_init,range);
plot(Esp_fit_rad{t}{gy}(:,1),Esp_fit_rad{t}{gy}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
xlim([0 5]);

%Display results
fprintf('\n')
fprintf('Residence Times')
fprintf('\n')
fprintf('-----------------------------------------------------')
fprintf('\n');
fprintf(files{t}.dataname(1:end-9));
fprintf('\n');
fprintf(strcat('Radius interval',str_rad{gy}));
fprintf('\n')
fprintf(strcat('k1 = ',num2str(Esp_coeff_rad{t}{gy}(1)),' ± ',num2str(Esp_Sigma_rad{t}{gy}(1))))
fprintf('\n')   
fprintf(strcat('t1 = ',num2str(1/Esp_coeff_rad{t}{gy}(1)),' ± ',num2str( -1/(Esp_coeff_rad{t}{gy}(1) + Esp_Sigma_rad{t}{gy}(1)) + 1/(Esp_coeff_rad{t}{gy}(1)) )));
fprintf('\n')  
fprintf(strcat('k2 = ',num2str(Esp_coeff_rad{t}{gy}(2)),' ± ',num2str(Esp_Sigma_rad{t}{gy}(2))))
fprintf('\n')   
fprintf(strcat('t2 = ',num2str(1/Esp_coeff_rad{t}{gy}(2)),' ± ',num2str( -1/(Esp_coeff_rad{t}{gy}(2) + Esp_Sigma_rad{t}{gy}(2)) + 1/(Esp_coeff_rad{t}{gy}(2)) )));
fprintf('\n')  
fprintf(strcat('f = ',num2str(Esp_coeff_rad{t}{gy}(3)),' ± ',num2str(Esp_Sigma_rad{t}{gy}(3))))
fprintf('\n')   

    end
end



%% Analyze and compare the angle distributions and the Jumps distributions.
ThetaSign = {};
figure()
for t=1:ntypes;
    
    for i=1:size(MSD_Results{t}.ma.tracks,1);
        
        dvector = [];
        SignTan = [];
        TanTheta_1 = [];
        TanTheta_2 = [];
        translocation = [];
        jump_per_angle = [];
        tracktemp = MSD_Results{t}.ma.tracks{i}; %In um.
        
 
        dvector = [ ((tracktemp(2:1:end,2) - tracktemp(1:1:end-1,2)))*1000 ((tracktemp(2:1:end,3) - tracktemp(1:1:end-1,3)))*1000 ]; %In nm.
        dvector(:,3) = 0; %Make the vector 3D for the cross product. 
        
        
            for hh=1:1:length(dvector)-1;        
                  TanTheta_1 = [TanTheta_1; norm(cross(dvector(hh,:),dvector(hh+1,:)))];
                  TanTheta_2 = [TanTheta_2; dot(dvector(hh,:),dvector(hh+1,:))];
                  signo = sign(cross(dvector(hh,:),dvector(hh+1,:)));
                  SignTan = [SignTan; signo(3)];  

                  translocation(hh,1) = mean([ norm(dvector(hh,:)), norm(dvector(hh+1,:)) ]); 
                  jump_per_angle(hh,1) = norm(dvector(hh,:) + dvector(hh+1,:));
            end  
            
            for hh=1:length(dvector);
            jumps_vector{t}{i}(hh) = norm(dvector(hh,:));     
            end
            total_track_displacement{t}(i) = norm( [ (tracktemp(end,2) - tracktemp(1,2))*1000 (tracktemp(end,3) - tracktemp(1,3))*1000 0] );
      
        angle{t}{i} = SignTan.*atan2(TanTheta_1,TanTheta_2); %Angle (sign included)
        mean_translocation{t}{i} = translocation;
        jumps_per_angle{t}{i} = jump_per_angle;
        
        mean_jump_vector{t}(i) = mean(jumps_vector{t}{i});
        
    
    end
        
Angle_all{t} = angle{t}(:);
Angle_all{t} = cat(1,Angle_all{t}{:});

f_asim{t} = (length(find(abs(Angle_all{t}*180/pi) >= 180-AC_arc)))/(length(find(abs(Angle_all{t}*180/pi) <= AC_arc)));%Calculate asimetry coefficient for each track.

Mean_translocation_all{t} = mean_translocation{t}(:);
Mean_translocation_all{t} = cat(1,Mean_translocation_all{t}{:});

Jumps_per_angle_all{t} = jumps_per_angle{t}(:);
Jumps_per_angle_all{t} = cat(1,Jumps_per_angle_all{t}{:});

[t_an, r_an] = rose(Angle_all{t},50);
r_an = r_an./numel(Angle_all{t}); % normalize
polar(t_an, r_an) % polar plot
hold on    
end

title('Histogram of angles between consecutive jumps','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==6;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end  
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end

%  
%  figure()
%  for t=1:ntypes;
%  scatter(Mean_translocation_all{t},Angle_all{t});
%  hold on
%  title('Angles vs Jump Displacement','FontSize',24,'FontWeight','bold');
%  xlabel('Mean jump (nm)');
% ylabel('Angle (radians)');
% set(gca,'FontSize',20,'FontWeight','bold');  
%  end
% 
%  
%  figure()
%  for t=1:ntypes;
%      
%    for i=1:10;
%  bin_size = 30; %Bin size in nanometers.    
%  range_local = [bin_size*i - bin_size; bin_size*i];
%  
%  jumps_local_idx = find(mean_jump_vector{t} >=  range_local(1) & mean_jump_vector{t} <  range_local(2)); 
%  angles_local{t} = angle{t}(jumps_local_idx);
%  angles_local{t} = cat(1,angles_local{t}{:});
%  
%  f_asim_local{t}(i) = (length(find(abs(angles_local{t}*180/pi) >= 180-AC_arc)))/(length(find(abs(angles_local{t}*180/pi) <= AC_arc)));%Calculate asimetry coefficient for each track.
% 
%  
%    end  
%      
% plot((bin_size/2):bin_size:bin_size*i,f_asim_local{t});
% hold on
% title('f180/0 vs Average Jump Displacement for each track','FontSize',24,'FontWeight','bold');
% xlabel('Average jump displacement (nm)');
% ylabel('Assimetry coefficient (f180/0)');
% set(gca,'FontSize',20,'FontWeight','bold');
%  end
%  
%  if ntypes == 1;
%  l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
%  end
%  if ntypes == 2;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes == 3;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes ==4;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes ==5;
%     l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
%  end   
%  if ntypes ==6;
% l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
%  end   
%  if ntypes ==7;
%     l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
%  end  
%  if ntypes ==8;
% l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
%  end
%  
% clear angles_local f_asim_local 
% 
%  figure()
%  for t=1:ntypes;
%      
%    for i=1:10;
%  bin_size = 30; %Bin size in nanometers.    
%  range_local = [bin_size*i - bin_size; bin_size*i];
%  
%  jumps_local_idx = find(total_track_displacement{t} >=  range_local(1) & total_track_displacement{t} <  range_local(2)); 
%  angles_local{t} = angle{t}(jumps_local_idx);
%  angles_local{t} = cat(1,angles_local{t}{:});
%  
%  f_asim_local{t}(i) = (length(find(abs(angles_local{t}*180/pi) >= 180-AC_arc)))/(length(find(abs(angles_local{t}*180/pi) <= AC_arc)));%Calculate asimetry coefficient for each track.
% 
%  
%    end  
%      
% plot((bin_size/2):bin_size:bin_size*i,f_asim_local{t});
% hold on
% title('f180/0 vs Total Track Displacement','FontSize',24,'FontWeight','bold'); 
% xlabel('Total track displacement (nm)');
% ylabel('Assimetry coefficient (f180/0)');
% set(gca,'FontSize',20,'FontWeight','bold');
%  end
%  
%  if ntypes == 1;
%  l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
%  end
%  if ntypes == 2;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes == 3;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes ==4;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes ==5;
%     l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
%  end   
%  if ntypes ==6;
% l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
%  end   
%  if ntypes ==7;
%     l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
%  end  
%  if ntypes ==8;
% l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
%  end
%   
%  
%  clear angles_local f_asim_local 

 figure()
 for t=1:ntypes;
     
   for i=1:AC_num_groups;
 range_local = [AC_bin_size*i - AC_bin_size; AC_bin_size*i];
 
 jumps_local_idx = find(Mean_translocation_all{t} >=  range_local(1) & Mean_translocation_all{t} <  range_local(2)); 
 angles_local{t} = Angle_all{t}(jumps_local_idx);
 
 f_asim_local{t}(i) = (length(find(abs(angles_local{t}*180/pi) >= 180-AC_arc)))/(length(find(abs(angles_local{t}*180/pi) <= AC_arc)));%Calculate asimetry coefficient for each track.

 AC_num_angles{t}(1,i) = length(angles_local{t});
 AC_num_angles{t}(2,i) = length(find(abs(angles_local{t}*180/pi) >= 180-AC_arc));
 AC_num_angles{t}(3,i) = length(find(abs(angles_local{t}*180/pi) <= AC_arc));
 AC_num_angles{t}(4,i) = length(find( abs(angles_local{t}*180/pi) > AC_arc & abs(angles_local{t}*180/pi) < 180-AC_arc));

   end  
     
plot((AC_bin_size/2):AC_bin_size:AC_bin_size*i,f_asim_local{t});
hold on
title('f180/0 vs mean displacement of each two corresponding jumps','FontSize',24,'FontWeight','bold'); 
xlabel('Mean displacement (nm)');
ylabel('Assimetry coefficient (f180/0)');
set(gca,'FontSize',20,'FontWeight','bold');

 end
 
 if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
 end
 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==6;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end  
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end
 
%  
%   clear angles_local f_asim_local 
%  
% 
%  figure()
%  for t=1:ntypes;
%      
%    for i=1:15;
%  range_local = [AC_bin_size*i - AC_bin_size; AC_bin_size*i];
%  
%  jumps_local_idx = find(Jumps_per_angle_all{t} >=  range_local(1) & Jumps_per_angle_all{t} <  range_local(2)); 
%  angles_local{t} = Angle_all{t}(jumps_local_idx);
%  
%  f_asim_local{t}(i) = (length(find(abs(angles_local{t}*180/pi) >= 180-AC_arc)))/(length(find(abs(angles_local{t}*180/pi) <= AC_arc)));%Calculate asimetry coefficient for each track.
% 
%  
%    end  
%      
% plot((bin_size/2):bin_size:bin_size*i,f_asim_local{t});
% hold on
% title('f180/0 vs absolute displacement after each two corresponding jumps','FontSize',24,'FontWeight','bold'); 
% xlabel('Displacement after two jumps (nm)');
% ylabel('Assimetry coefficient (f180/0)');
% set(gca,'FontSize',20,'FontWeight','bold');
% 
%  end
%  
%  if ntypes == 1;
%  l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
%  end
%  if ntypes == 2;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes == 3;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes ==4;
%  l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
%  end
%  if ntypes ==5;
%     l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
%  end   
%  if ntypes ==6;
% l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
%  end   
%  if ntypes ==7;
%     l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
%  end  
%  if ntypes ==8;
% l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
%  end  
  



%% Analyze and compare the angle distributions and the time lag distributions.
for t=1:ntypes;
    
    for i=1:size(MSD_Results{t}.ma.tracks,1);
        
        tracktemp = MSD_Results{t}.ma.tracks{i}; %In um.
        
        [Angles_lagtimes{t}{i}, Lag_times{t}{i}] = SMT_Angles_vs_LagTime_V2(tracktemp,Frame_interval,Localization_precision_estimation{t});

    end
 
end


for t=1:ntypes;
       for lg = 1:max_lag;
            Angle_lagtimes_idx{t}{lg} = cellfun(@(x) find(size(x,2)>=lg), Angles_lagtimes{t},'UniformOutput',0);
            Angle_lagtimes_idx_logi{t}{lg} = ~cellfun('isempty', Angle_lagtimes_idx{t}{lg});
            Angle_lagtimes_temp{t}{lg} = cellfun(@(x) x{lg}, Angles_lagtimes{t}(find(Angle_lagtimes_idx_logi{t}{lg}==1)),'UniformOutput',0);
            Angle_lagtimes_all{t}{lg} = cat(1,Angle_lagtimes_temp{t}{lg}{:});
            f_asim_lagtimes{t}(lg) = (length(find(abs(Angle_lagtimes_all{t}{lg}*180/pi) >= 180-AC_arc)))/(length(find(abs(Angle_lagtimes_all{t}{lg}*180/pi) <= AC_arc)));%Calculate asimetry coefficient for each track.
       end
end



 figure()
 for t=1:ntypes;
          
plot(Frame_interval:Frame_interval:max_lag*Frame_interval,f_asim_lagtimes{t});
hold on
title('f180/0 vs Time lag of the jumps','FontSize',24,'FontWeight','bold'); 
ylabel('AC');
xlabel('Lag time (s)');
set(gca,'FontSize',20,'FontWeight','bold');
 end
 
 if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
 end
 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==6;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end  
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end




 
 
%% Analyze and compare the radius of gyration evolution over time
figure()
strB = '';
 for t=1:ntypes;
    
     
fprintf('\n')
fprintf('Radius of gyration exponential model')
fprintf('\n')
fprintf('-----------------------------------------------------')
fprintf('\n')
fprintf('fitting each individual Radius of gyration time evolution curve')
fprintf('\n')
fprintf('\n')
fprintf(strcat('Step_',num2str(t),'_of_',num2str(ntypes)))
fprintf('\n')   

     
    for i=1:size(MSD_Results{t}.ma.tracks,1);
        
        %Loop through each trajectory
        for k = 2:length(MSD_Results{t}.ma.tracks{i});
              
            x_seg = MSD_Results{t}.ma.tracks{i}(1:k,2); 
            y_seg = MSD_Results{t}.ma.tracks{i}(1:k,3);
            R_Tensor{t}{i}(1,1) = mean(x_seg.^2) - mean(x_seg).^2;
            R_Tensor{t}{i}(1,2) = mean(x_seg.*y_seg) - mean(x_seg).*mean(y_seg);
            R_Tensor{t}{i}(2,1) = R_Tensor{t}{i}(1,2);
            R_Tensor{t}{i}(2,2) = mean(y_seg.^2) - mean(y_seg).^2;

            %Determine the eigenvalues of R_Tensor and calculate Rg in microns
            R_values = eig(R_Tensor{t}{i});
            R1sq = R_values(1);
            R2sq = R_values(2);
            
            Radius_of_gyration{t}{i}(k-1) = (R1sq + R2sq).^(0.5); 
        end
        
     %Fit an exponential function to each individual Rg curve to extract the
    %HalfTime (how fast it reaches the plateau) and the asimptotic value:
    
    [Esp_coeff_gyration{t}{i},Esp_Sigma_gyration{t}{i}, Esp_fit_gyration{t}{i},resnorm_gyration_all{t}{i}] = ExpGrow_1Cmp_fit([(MSD_Results{t}.ma.tracks{i}(2:end,1)-min(MSD_Results{t}.ma.tracks{i}(:,1))) Radius_of_gyration{t}{i}'], k_init_gyration, range_gyration);
    Half_life_gyration_all{t}(i) = Esp_coeff_gyration{t}{i}(1);
    Plateau_value_gyration_all{t}(i) = Esp_coeff_gyration{t}{i}(2);
    
    Max_gyration_all{t}(i) = max(Radius_of_gyration{t}{i});
    
    strOut = sprintf('Completed: % 4.1f',100*i/(size(MSD_Results{t}.ma.msd,1)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);           
            
    end
    fprintf('\n') 
    
            %Calculate ensemble average of Rg from all particles for N time steps.
            for n=1:20;
            idx = find( (cellfun(@(x) length(x), Radius_of_gyration{t})) >= n);    
            Radius_of_gyration_mean{t}(n) = mean(cellfun(@(x) x(n), Radius_of_gyration{t}(idx)));
            end
    
    plot(Frame_interval:Frame_interval:20*Frame_interval,Radius_of_gyration_mean{t},'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
    hold on

  iidx_gyration{t} = find(cell2mat(resnorm_gyration_all{t}) < 0.005);

    
 end
fprintf('\n') 

title('Time evolution of emsemble average of radius of gyration','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Time (s)');
ylabel('Radius of gyration (um)');
 



%Show an example of an individual trayectory with the exponential fitting
figure()
suptitle('Some examples of the fitting of radius of gyration time evolution'); 

subplot(4,4,1)
plot(Esp_fit_gyration{t}{i}(:,1),Esp_fit_gyration{t}{i}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i}(2:end,1)-min(MSD_Results{t}.ma.tracks{i}(:,1))),Radius_of_gyration{t}{i}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,2)
plot(Esp_fit_gyration{t}{i-1}(:,1),Esp_fit_gyration{t}{i-1}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-1}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-1}(:,1))),Radius_of_gyration{t}{i-1}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,3)
plot(Esp_fit_gyration{t}{i-2}(:,1),Esp_fit_gyration{t}{i-2}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-2}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-2}(:,1))),Radius_of_gyration{t}{i-2}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,4)
plot(Esp_fit_gyration{t}{i-3}(:,1),Esp_fit_gyration{t}{i-3}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-3}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-3}(:,1))),Radius_of_gyration{t}{i-3}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,5)
plot(Esp_fit_gyration{t}{i-4}(:,1),Esp_fit_gyration{t}{i-4}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-4}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-4}(:,1))),Radius_of_gyration{t}{i-4}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,6)
plot(Esp_fit_gyration{t}{i-5}(:,1),Esp_fit_gyration{t}{i-5}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-5}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-5}(:,1))),Radius_of_gyration{t}{i-5}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,7)
plot(Esp_fit_gyration{t}{i-6}(:,1),Esp_fit_gyration{t}{i-6}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-6}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-6}(:,1))),Radius_of_gyration{t}{i-6}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,8)
plot(Esp_fit_gyration{t}{i-7}(:,1),Esp_fit_gyration{t}{i-7}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-7}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-7}(:,1))),Radius_of_gyration{t}{i-7}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,9)
plot(Esp_fit_gyration{t}{i-8}(:,1),Esp_fit_gyration{t}{i-8}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-8}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-8}(:,1))),Radius_of_gyration{t}{i-8}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,10)
plot(Esp_fit_gyration{t}{i-9}(:,1),Esp_fit_gyration{t}{i-9}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-9}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-9}(:,1))),Radius_of_gyration{t}{i-9}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,11)
plot(Esp_fit_gyration{t}{i-10}(:,1),Esp_fit_gyration{t}{i-10}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-10}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-10}(:,1))),Radius_of_gyration{t}{i-10}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,12)
plot(Esp_fit_gyration{t}{i-11}(:,1),Esp_fit_gyration{t}{i-11}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-11}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-11}(:,1))),Radius_of_gyration{t}{i-11}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,13)
plot(Esp_fit_gyration{t}{i-12}(:,1),Esp_fit_gyration{t}{i-12}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-12}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-12}(:,1))),Radius_of_gyration{t}{i-12}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,14)
plot(Esp_fit_gyration{t}{i-13}(:,1),Esp_fit_gyration{t}{i-13}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-13}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-13}(:,1))),Radius_of_gyration{t}{i-13}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,15)
plot(Esp_fit_gyration{t}{i-14}(:,1),Esp_fit_gyration{t}{i-14}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-14}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-14}(:,1))),Radius_of_gyration{t}{i-14}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

subplot(4,4,16)
plot(Esp_fit_gyration{t}{i-15}(:,1),Esp_fit_gyration{t}{i-15}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
hold on
plot((MSD_Results{t}.ma.tracks{i-15}(2:end,1)-min(MSD_Results{t}.ma.tracks{i-15}(:,1))),Radius_of_gyration{t}{i-15}');
xlabel('Time (s)');
ylabel('Radius of gyration (um)');

%% Analyze and compare the residence time distribution of the Tracks
Residence_Times = {};

figure()
for t=1:ntypes;
Num_particles{t} = size(MSD_Results{t}.ma.tracks,1);
    
    for i=1:Num_particles{t};
        
    Residence_Times{t}(i) = MSD_Results{t}.ma.tracks{i}(end,1) - MSD_Results{t}.ma.tracks{i}(1,1);
    
    end
    
    
    
   binCtrs = 0:5:350; %Bin centers. (depends on your data)
   n = length(Residence_Times{t});
   counts = hist(Residence_Times{t},binCtrs);
   prob = counts / (n * 5);
   H = bar(binCtrs,prob,'hist');
   set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
   hold on;
   
xlabel('Residence time (s)','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontSize',20,'FontWeight','bold');
title('Histogram of Residence Times','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  

   
end


%Fit 
figure()
Esp_coeff = {};
Esp_Sigma = {};
Esp_fit = {};
for t=1:ntypes;
survival{t} = [];       

max_time{t} = max(Residence_Times{t});
min_time{t} = min(Residence_Times{t});

        for r_time = min_time{t}:Frame_interval:max_time{t};
            
            survival{t} = [survival{t}; r_time length(find(Residence_Times{t} >= r_time))/Num_particles{t}];
            
        end

        
scatter(survival{t}(:,1),survival{t}(:,2),'MarkerFaceColor',Color(t,:));
hold on
ylabel('Normalized survival fraction (Counts)','FontSize',20,'FontWeight','bold');
xlabel('Time (s)','FontSize',20,'FontWeight','bold');
title('Normalized survival fraction','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  


%Fit a 2-component exponential decay to the survival probability function
[Esp_coeff{t},Esp_Sigma{t}, Esp_fit{t}] = ExpDecay_2Cmp_fit(survival{t}, [k1_init k2_init], fraction_init,range);
plot(Esp_fit{t}(:,1),Esp_fit{t}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');

xlim([0 60]);

% if ntypes == 1;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
% end
% 
%  if ntypes == 2;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
%  end
%  if ntypes == 3;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
%  end
%  if ntypes ==4;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
%  end
%  if ntypes ==5;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
%    end
%     
%      if ntypes ==6;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
%     end
%     
%  if ntypes ==7;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
%  end
%     
%      if ntypes ==8;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
%      end


%Fit a 2-component exponential decay to the number of molecules per frame
%evolution over time. This is to estimate the effect of the photobleaching.
counter = 1;
for uu = min(MSD_Results{t}.ma.tracks{1}(:,1)):Frame_interval:max(MSD_Results{t}.ma.tracks{end}(:,1)); 
Molecules{t} = cellfun(@(x) find(x(:,1) == uu),MSD_Results{t}.ma.tracks,'UniformOutput',false);
Molecules_evolution{t}(counter,1) = uu;
Molecules_evolution{t}(counter,2) = length(find(~cellfun('isempty',Molecules{t})));
counter = counter + 1;
end


[Esp_coeff_bleach{t},Esp_Sigma_bleach{t}, Esp_fit_bleach{t}] = ExpDecay_2Cmp_fit(Molecules_evolution{t}, [k1_init_bleach k2_init_bleach], fraction_init_bleach, range_bleach);



%Display results
fprintf('\n')
fprintf('Residence Times')
fprintf('\n')
fprintf('-----------------------------------------------------')
fprintf('\n');
fprintf(files{t}.dataname(1:end-9));
fprintf('\n'); 
fprintf('\n')
fprintf(strcat('k1 = ',num2str(Esp_coeff{t}(1)),' ± ',num2str(Esp_Sigma{t}(1))))
fprintf('\n')   
fprintf(strcat('t1 = ',num2str(1/Esp_coeff{t}(1)),' ± ',num2str( -1/(Esp_coeff{t}(1) + Esp_Sigma{t}(1)) + 1/(Esp_coeff{t}(1)) )));
fprintf('\n')  
fprintf(strcat('k2 = ',num2str(Esp_coeff{t}(2)),' ± ',num2str(Esp_Sigma{t}(2))))
fprintf('\n')   
fprintf(strcat('t2 = ',num2str(1/Esp_coeff{t}(2)),' ± ',num2str( -1/(Esp_coeff{t}(2) + Esp_Sigma{t}(2)) + 1/(Esp_coeff{t}(2)) )));
fprintf('\n')  
fprintf(strcat('f = ',num2str(Esp_coeff{t}(3)),' ± ',num2str(Esp_Sigma{t}(3))))
fprintf('\n')   

fprintf('\n')
fprintf('Residence Times with Photobleaching correction')
fprintf('\n')
fprintf('-----------------------------------------------------')
fprintf('\n');
fprintf(files{t}.dataname(1:end-9));
fprintf('\n'); 
fprintf('\n')
fprintf(strcat('k1 = ',num2str(Esp_coeff{t}(1) - Esp_coeff_bleach{t}(1)),' ± ',num2str(Esp_Sigma{t}(1))))
fprintf('\n')   
fprintf(strcat('t1 = ',num2str(1/(Esp_coeff{t}(1) - Esp_coeff_bleach{t}(1))),' ± ',num2str( -1/( (Esp_coeff{t}(1) - Esp_coeff_bleach{t}(1)) + Esp_Sigma{t}(1)) + 1/( (Esp_coeff{t}(1)) - Esp_coeff_bleach{t}(1)) )));
fprintf('\n')  
fprintf(strcat('k2 = ',num2str(Esp_coeff{t}(2) - Esp_coeff_bleach{t}(2)),' ± ',num2str(Esp_Sigma{t}(2))))
fprintf('\n')   
fprintf(strcat('t2 = ',num2str(1/(Esp_coeff{t}(2) - Esp_coeff_bleach{t}(2))),' ± ',num2str( -1/( (Esp_coeff{t}(2) - Esp_coeff_bleach{t}(2)) + Esp_Sigma{t}(2)) + 1/( (Esp_coeff{t}(2)) - Esp_coeff_bleach{t}(2)) )));
fprintf('\n')  
fprintf(strcat('f = ',num2str(Esp_coeff{t}(3)),' ± ',num2str(Esp_Sigma{t}(3))))
fprintf('\n') 

end

figure()
for t=1:ntypes;
scatter(Molecules_evolution{t}(:,1),Molecules_evolution{t}(:,2));
ylabel('Number of molecules','FontSize',20,'FontWeight','bold');
xlabel('Time (s)','FontSize',20,'FontWeight','bold');
title('Number of molecules evolution (Photobleaching kinetics','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlim([0 60]);
hold on
plot(Esp_fit_bleach{t}(:,1),Esp_fit_bleach{t}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
end

% if ntypes == 1;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
% end
% 
%  if ntypes == 2;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
%  end
%  if ntypes == 3;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
%  end
%  if ntypes ==4;
%  l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
%  end
%  if ntypes ==5;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
%    end
%     
%      if ntypes ==6;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
%     end
%     
%  if ntypes ==7;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
%  end
%     
%      if ntypes ==8;
%     l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
%      end

% figure()
% for t=1:ntypes;
% scatter(survival{t}(:,1),survival{t}(:,2),'MarkerFaceColor',Color(t,:));
% hold on
% ylabel('Normalized survival fraction (Counts) with Photobleaching correction','FontSize',20,'FontWeight','bold');
% xlabel('Time (s)','FontSize',20,'FontWeight','bold');
% title('Normalized survival fraction','FontSize',24,'FontWeight','bold'); 
% set(gca,'FontSize',20,'FontWeight','bold');  
% 
% yg{t} = ExpDecay_2Cmp_fun([Esp_coeff{t}(1)-Esp_coeff_bleach{t}(1);  Esp_coeff{t}(2)-Esp_coeff_bleach{t}(2); Esp_coeff{t}(3); Esp_coeff{t}(4)], 0.5:0.5:100);
% plot(0.5:0.5:100,yg{t},'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
% end

%% Plot an error bar with the residence times results

figure()
for t=1:ntypes;

   t2 = 1/Esp_coeff{t}(2);
   t2_sigma = -1/(Esp_coeff{t}(2) + Esp_Sigma{t}(2)) + 1/(Esp_coeff{t}(2));
    
    
ha{t} = errorbar(t, t2, t2_sigma, t2_sigma,'o','MarkerSize',5,'MarkerFaceColor',Color(t,:),'Linewidth',1.4);
% xlim([0 0.1]);
% ylim([-0.02 0.05]);
  
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
     end
      
title('Errorbar plot of Residence times','FontSize',20,'FontWeight','bold');   
set(l, 'Interpreter', 'none')  
set(l, 'FontSize',20) 
set(gca, 'XTick', 0:ntypes)
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Residence times (s)','FontSize',20,'FontWeight','bold');
hold on
drawnow
errorbar_tick(ha{t},50);    

    
end



figure()
for t=1:ntypes;

   t2 = 1/(Esp_coeff{t}(2) - Esp_coeff_bleach{t}(2));
   t2_sigma = -1/( (Esp_coeff{t}(2) - Esp_coeff_bleach{t}(2)) + Esp_Sigma{t}(2)) + 1/( (Esp_coeff{t}(2)) - Esp_coeff_bleach{t}(2));
    
    
ha{t} = errorbar(t, t2, t2_sigma, t2_sigma,'o','MarkerSize',5,'MarkerFaceColor',Color(t,:),'Linewidth',1.4);
% xlim([0 0.1]);
% ylim([-0.02 0.05]);
  
if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Interpreter','none');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Interpreter','none');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Interpreter','none');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Interpreter','none');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Interpreter','none');
     end
      
title('Errorbar plot of Residence times after photobleaching correction','FontSize',20,'FontWeight','bold');   
set(l, 'Interpreter', 'none')  
set(l, 'FontSize',20) 
set(gca, 'XTick', 0:ntypes)
set(gca,'FontSize',20,'FontWeight','bold');
ylabel('Residence times (s)','FontSize',20,'FontWeight','bold');
hold on
drawnow
errorbar_tick(ha{t},50);    

    
end



%% Analise and compare the jump distributions of the tracks

figure()



for t=1:ntypes;

jumps_vector_cat{t} = horzcat(jumps_vector{t}{:});    

binCtrs = 0:5:500; %Bin centers, depends on your data    
n=length(jumps_vector_cat{t});
counts = hist(jumps_vector_cat{t},binCtrs);
prob = counts / (n * 10);
H = bar(binCtrs,prob,'hist');
set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
hold on;
end

xlim([-10 200]);

title('Histogram of jumps length','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Jump length (nm)');
ylabel('frequency');

 if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9)),'Interpreter','none');   
 end
 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),'Interpreter','none');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==6;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),'Interpreter','none');
 end   
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),'Interpreter','none');
 end  
 if ntypes ==8;
l=legend(strcat(files{1}.dataname(1:end-9)),strcat(files{2}.dataname(1:end-9)),strcat(files{3}.dataname(1:end-9)),strcat(files{4}.dataname(1:end-9)),strcat(files{5}.dataname(1:end-9)),strcat(files{6}.dataname(1:end-9)),strcat(files{7}.dataname(1:end-9)),strcat(files{8}.dataname(1:end-9)),'Interpreter','none');
 end


 %% Plot the Radius of confinement againt the radius of gyration fit for each track
 
 for t=1:ntypes;     
densityplot(f_confined_circ_diff_R_all{t}(iidx{t})*1000,Half_life_gyration_all{t}(iidx{t}),edges_gyration_density_plot);    
title(strcat(files{t}.dataname(1:end-9),'_Filtered'),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Radius of confinement (nm)');
ylabel('K radius of gyration (1/s)');
 end

 for t=1:ntypes;    
densityplot(f_confined_circ_diff_R_all{t}(iidx{t})*1000,Plateau_value_gyration_all{t}(iidx{t})*1000,edges_gyration_density_plot); 
title(strcat(files{t}.dataname(1:end-9),'_Filtered'),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Radius of confinement (nm)');
ylabel('Radius of gyration (nm)');
 end
 
 
for t=1:ntypes;     
densityplot(f_confined_circ_diff_R_all{t}*1000,Half_life_gyration_all{t},edges_gyration_density_plot);    
title(strcat(files{t}.dataname(1:end-9)),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Radius of confinement (nm)');
ylabel('K radius of gyration (1/s)');
axis('square');
end

for t=1:ntypes;    
densityplot(f_confined_circ_diff_R_all{t}*1000,Plateau_value_gyration_all{t}*1000,edges_gyration_density_plot); 
title(strcat(files{t}.dataname(1:end-9)),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Radius of confinement (nm)');
ylabel('Radius of gyration (nm)');
axis('square');
end
 
for t=1:ntypes;  
figure()
out = scatplot(f_confined_circ_diff_R_all{t}*1000,Plateau_value_gyration_all{t}*1000,'voronoi',5,10,5,1,4);
title(strcat(files{t}.dataname(1:end-9)),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
set(gca,'FontSize',20,'FontWeight','bold');  
xlabel('Radius of confinement (nm)');
ylabel('Radius of gyration (nm)');
axis('square');
xlim([0 100]);
ylim([0 100]);
end
 
%% Save the results from the MSD Analysis

Results = {};
Results{1,2} = 'Radius of confinement (Confined circle model) ensemble (um)';
Results{1,3} = 'Difussion Coefficient (Confined circle model) ensemble (um2/s)';
Results{1,4} = 'Confidence Intervals (Confined circle model) ensemble (um and s)';
Results{1,5} = 'Radius of confinement (Confined circle model) all tracks (um)';
Results{1,6} = 'Difussion Coefficient (Confined circle model) all tracks (um2/s)';
Results{1,7} = 'Confidence Intervals (Confined circle model) all tracks (um and s)';
Results{1,8} = 'Track length of the trajectories with good enough confined difussion model fit (frames)';
Results{1,9} = 'Radius of gyration (Ensemble Mean) (um)';
Results{1,10} = 'Radius of gyration (All Tracks) (um)';
Results{1,11} = 'Difussion Coefficient (TE-MSD emsemble fit) (um2/s)';
Results{1,12} = 'Difussion Coefficient CI (TE-MSD emsemble fit) (um2/s)';
Results{1,13} = 'Difussion Coefficient (T-MSD fit) (um2/s)';
Results{1,14} = 'Probability Unimodal distribution of Diffusion Coefficients';
Results{1,15} = 'Residence times (s)';
Results{1,16} = 'Maximum distance travelled (um)';
Results{1,17} = 'Angle between segments (radians)';
Results{1,18} = 'Jumps length distribution (um)';
Results{1,19} = 'Half life radius of gyration exponential fit all tracks (1/s)';
Results{1,20} = 'Asymptote radius of gyration exponential fit all tracks (um)';
Results{1,21} = 'Half life radius of gyration exponential fit [Only tracks with good circle model fit] (1/s)';
Results{1,22} = 'Asymptote radius of gyration exponential fit [Only tracks with good circle model fit] (um)';
Results{1,23} = 'Maximum value of the radius of gyration time evolution (um)';
% Results{1,24} = '';


for t=1:ntypes;
Results{t+1,1} = strcat(files{t}.dataname(1:end-9));

Results{t+1,2} = f_confined_circ_diff_R{t};
Results{t+1,3} = f_confined_circ_diff_D{t};
Results{t+1,4} = conf{t};
Results{t+1,5} = f_confined_circ_diff_R_all{t}(iidx{t});
Results{t+1,6} = f_confined_circ_diff_D_all{t}(iidx{t});
Results{t+1,7} = conf_all{t}(iidx{t});
Results{t+1,8} = Track_length_radius{t}(iidx{t});
Results{t+1,9} = Radius_of_gyration_mean{t};
Results{t+1,10} = Radius_of_gyration{t};
Results{t+1,11} = D_ensemble(t);
Results{t+1,12} = [D_emsemble_lower_ci(t) D_emsemble_lower_ci(t)];
Results{t+1,13} = D_individual{t};
Results{t+1,14} = probability_unimodal_D(t);
Results{t+1,15} = Residence_Times{t};
Results{t+1,16} = max_dist{t};
Results{t+1,17} = Angle_all{t};
Results{t+1,18} = jumps_vector_cat{t}/1000;
Results{t+1,19} = Half_life_gyration_all{t};
Results{t+1,20} = Plateau_value_gyration_all{t};
Results{t+1,21} = Half_life_gyration_all{t}(iidx{t});
Results{t+1,22} = Plateau_value_gyration_all{t}(iidx{t});
Results{t+1,23} = Max_gyration_all{t}(iidx{t});
% Results{t+1,24} = ;


end

if Save_Results == 1;
fprintf('\n');     
fprintf('Saving...\n');

dir = uigetdir(startlocn{1});
% mkdir(strcat(files{1}.data{1,2},'\Compare_MSD_Results_'));
% save(fullfile(strcat(files{1}.data{1,2},'\MSD_Results_')),'Results');
save(fullfile(dir,strcat(File_name,'.mat')),'Results');

end



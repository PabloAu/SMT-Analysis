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
type.t1 = 'H2B ESC';
type.t2 = 'H2B NPC';
% type.t5 = '';
% type.t6 = '';


%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)
startlocn = {'C:\Users\pgomez.CRG\Documents\PhD\My papers\SMT\Data'};
Save_Results = 1;
File_name = 'MSD_Results';


%------------------------------------------------------------------
%Numerical Inputs---------------------------------------------------
%------------------------------------------------------------------
n_dim = 2; %Dimensionality of the movement (2 for 2D and 3 for 3D).
Frame_interval = 0.5; %Exposure time in seconds

TMSD_fitting_points = 3;   %Minimum 3 point for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 3;  %Minimum 3 point for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; 

R2LIMIT = 0.8; %R-squared minimum value from the D fits on each T-MSD Curve.

%-----------------------------------------------------------------
%For the Confined Circle Difussion Model Fitting---------------------------------------------------
%---------------------------------------------------------------------
num_points = 20; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
level = 0.01; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in um]. 
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (um^2/s)
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (um)

Radius_limit = 60; %Radius limit in order to divide tracks into two populations. (log10(60) = 1.78)

%-----------------------------------------------------------------
%For the Histogram of Difussion Coefficients---------------------------------------------------
%---------------------------------------------------------------------
binWidth = 0.1; %Bin width
bin_start = -5; %Bin start
bin_end = -0.5; %Bin end

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

%------------------------------------------------------------------------
%Initialize variables-------------------------------------------------
%----------------------------------------------------------------------
Dcoeffs = {};
datatypes = fieldnames(type);
ntypes = size(datatypes,1);         


%% Load the Data
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%Ask the user to select the files
files = cell(ntypes,1);
for t = 1:ntypes;
    files{t}= Select1DataGroup(type.(datatypes{t}),'mat',startlocn{1});
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
    
    MSD_Results{t}.ma = MSD_Results{t}.ma.addAll(MSD_Results_prov{t}{c}.ma.tracks);

    end
    
    MSD_Results{t}.ma = MSD_Results{t}.ma.computeMSD;
%     MSD_Results{t}.ma = MSD_Results{t}.ma.LogTMSD(TLOGLOG_fitting_points);
%     MSD_Results{t}.ma = MSD_Results{t}.ma.TMSD(TMSD_fitting_points); 

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
% ylim([0 1000]);
% xlim([0 200]);
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

xlim([0 num_points*Frame_interval+Frame_interval/2]);
ylim([0 0.12]);
drawnow
errorbar_tick(ha{t},20);
hold on
msd_confined_circle_diffusion = @(X,Xdata)(X(1)^2)*(1-exp(-4*X(2)*Xdata/(X(1)^2))) + X(3);
times = linspace(0,20,100000); %From 0 to 20 seconds. Just for plotting
plot(times,msd_confined_circle_diffusion(Val_fitted{t},times),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');


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


binCtrs = 0:10:1000; %Bin centers, depends on your data
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
xlim([0 1000]);

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

        
scatter(survival{t}(2:end,1),survival{t}(2:end,2),'MarkerFaceColor',Color(t,:));
hold on
ylabel('Normalized survival fraction (Counts)','FontSize',20,'FontWeight','bold');
xlabel('Time (s)','FontSize',20,'FontWeight','bold');
title('Normalized survival fraction','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  


%Fit a 2-component exponential decay to the survival probability function
[Esp_coeff{t},Esp_Sigma{t}, Esp_fit{t}] = ExpDecay_2Cmp_fit(survival{t}, [k1_init k2_init], fraction_init,range);
plot(Esp_fit{t}(:,1),Esp_fit{t}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.')

xlim([0 60]);

if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Double exponential fitting');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Double exponential fitting',strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Double exponential fitting');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Double exponential fitting',strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Double exponential fitting',strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Double exponential fitting');
     end
set(l, 'Interpreter', 'none');  



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

if ntypes == 1;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting');   
end

 if ntypes == 2;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting');
 end
 if ntypes == 3;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting');
 end
 if ntypes ==4;
 l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting');
 end
 if ntypes ==5;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting');
   end
    
     if ntypes ==6;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Double exponential fitting');
    end
    
 if ntypes ==7;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Double exponential fitting',strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Double exponential fitting');
 end
    
     if ntypes ==8;
    l=legend(strcat(files{1}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(1))),'Double exponential fitting',strcat(files{2}.dataname(1:end-9),' Tracks =',num2str(Num_Tracks(2))),'Double exponential fitting',strcat(files{3}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(3))),'Double exponential fitting',strcat(files{4}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(4))),'Double exponential fitting',strcat(files{5}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(5))),'Double exponential fitting',strcat(files{6}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(6))),'Double exponential fitting',strcat(files{7}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(7))),'Double exponential fitting',strcat(files{8}.dataname(1:end-9),' Tracks = ',num2str(Num_Tracks(8))),'Double exponential fitting');
     end
set(l, 'Interpreter', 'none');  


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
ylabel('Residence Times (s)','FontSize',20,'FontWeight','bold');
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
ylabel('Residence Times (s)','FontSize',20,'FontWeight','bold');
hold on
drawnow
errorbar_tick(ha{t},50);    

    
end


%% Save the results from the MSD Analysis

Results = {};
Results{1,2} = 'Radius of confinement (Confined circle model) ensemble (nm)';
Results{1,3} = 'Difussion Coefficient (Confined circle model) ensemble (um2/s)';
Results{1,4} = 'Difussion Coefficient CI (Confined circle model) ensemble (um2/s)';
Results{1,5} = 'Difussion Coefficient (TE-MSD emsemble fit) (um2/s)';
Results{1,6} = 'Difussion Coefficient CI (TE-MSD emsemble fit) (um2/s)';
Results{1,7} = 'Difussion Coefficient (T-MSD fit) (um2/s)';
Results{1,8} = 'Probability Unimodal distribution of Diffusion Coefficients';
Results{1,9} = 'Residence times (s)';
Results{1,10} = 'Maximum distance travelled (nm)';
% Results{1,15} = '';
% Results{1,16} = '';
% Results{1,17} = '';

for t=1:ntypes;
Results{t+1,1} = strcat(files{t}.dataname(1:end-9));

Results{t+1,2} = f_confined_circ_diff_R{t};
Results{t+1,3} = f_confined_circ_diff_D{t};
Results{t+1,4} = conf{t};
Results{t+1,5} = D_ensemble(t);
Results{t+1,6} = [D_emsemble_lower_ci(t) D_emsemble_lower_ci(t)];
Results{t+1,7} = D_individual{t};
Results{t+1,8} = probability_unimodal_D(t);
Results{t+1,9} = Residence_Times{t};
Results{t+1,10} = max_dist{t};
% Results{t+1,15} = ;
% Results{t+1,16} = ;
% Results{t+1,17} = ;

end

if Save_Results == 1;
fprintf('\n');     
fprintf('Saving...\n');

dir = uigetdir(startlocn{1});
% mkdir(strcat(files{1}.data{1,2},'\Compare_MSD_Results_'));
% save(fullfile(strcat(files{1}.data{1,2},'\MSD_Results_')),'Results');
save(fullfile(dir,strcat(File_name,'.mat')),'Results');

end



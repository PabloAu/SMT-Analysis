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
% type.t1 = 'FOXA1 all';
% type.t2 = 'H3 all';
% type.t3 = 'FOXA1 Confined';
% type.t4 = 'H3 Confined';
% type.t5 = 'FOXA1 Directed';
% type.t6 = 'H3 Directed';
% type.t7 = 'FOXA1 Butterfly';
% type.t8 = 'H3 Butterfly';

% type.t1 = 'ESC butt direct';
% type.t2 = 'ESC brownian';
% type.t3 = 'Oct4';
% type.t4 = '';

type.t1 = 'ESC H2B PairedTracks';
type.t2 = 'ESC H2B UnpairedTracks';
% type.t3 = 'ESC H2B Euchromatin';

% type.t1 = 'ESC H2B';
% type.t2 = 'ESC H2B H1tKO';
% type.t3 = 'ESC Telomeres';
% type.t4 = 'Fixed Cells';

%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)
% startlocn = {'E:\CRG\Data\190201_2ColorSMT_TelomereCentromere and H2B\TRF1-GFP H2B-JF646\JF646 20pM\TrackMate Analysis\Pairing 700nm 50%'};
startlocn = {'E:\CRG\Data\190201_2ColorSMT_TelomereCentromere and H2B\TRF1-GFP H2B-JF646\JF646 20pM\TrackMate Analysis\Pairing 300nm 40%'};

Save_Results = 1;
File_name = 'MSD_Results';


%------------------------------------------------------------------
%Numerical Inputs---------------------------------------------------
%------------------------------------------------------------------
n_dim = 2; %Dimensionality of the movement (2 for 2D and 3 for 3D).
Frame_interval = 0.015; %Exposure time in seconds

TMSD_fitting_points = 3;   %Minimum 3 point for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 3;  %Minimum 3 point for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; 

R2LIMIT = 0.8; %R-squared minimum value from the D fits on each T-MSD Curve.

%-----------------------------------------------------------------
%For the Confined Circle Difussion Model Fitting---------------------------------------------------
%---------------------------------------------------------------------
num_points = 10; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
level = 0.001; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in ?m]. 
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (?m^2/s)
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (?m)

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
AC_bin_size = 20; %Bin size in nm for grouping the pair of jumps and calculating the assimetry coefficient
AC_num_groups = 20; %Number of groups to make. Notice that maximum jump would be AC_num_groups*AC_bin_size

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
scatter(Track_length_radius{t}(iidx{t}),f_confined_circ_diff_R_all{t}(iidx{t})*1000);
hold on
end

xlabel('Track length (frames)','FontSize',20,'FontWeight','bold');
ylabel('Radius of confinement (nm)','FontSize',20,'FontWeight','bold');
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
 
 
 
%% Divide tracks into two populations based on their radius of confinement
for t=1:ntypes;
idxr_high = find(f_confined_circ_diff_R_all{t}*1000 > Radius_interval(2));

figure()
plotTracks(MSD_Results{t}.ma,gca,idxr_high);
MSD_Results{t}.ma.labelPlotTracks;

title('Trayectories with Radius of confinement higher than Rlimit','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
set(gca,'Ydir','reverse');
l = legend(strcat(files{t}.dataname(1:end-9))); 
set(l, 'Interpreter', 'none')
end


for t=1:ntypes;
idxr_between = find(f_confined_circ_diff_R_all{t}*1000 >= Radius_interval(1) & f_confined_circ_diff_R_all{t}*1000 <= Radius_interval(2));

figure()
plotTracks(MSD_Results{t}.ma,gca,idxr_between);
MSD_Results{t}.ma.labelPlotTracks;

title('Trayectories with Radius of confinement in the specified interval','FontSize',24,'FontWeight','bold'); 
l = legend(strcat(files{t}.dataname(1:end-9))); 
set(l, 'Interpreter', 'none')

end



for t=1:ntypes;
idxr_low = find(f_confined_circ_diff_R_all{t}*1000 < Radius_interval(1));

figure()
plotTracks(MSD_Results{t}.ma,gca,idxr_low);
MSD_Results{t}.ma.labelPlotTracks;

title('Trayectories with Radius of confinement lower than Rlimit','FontSize',24,'FontWeight','bold'); 
l = legend(strcat(files{t}.dataname(1:end-9))); 
set(l, 'Interpreter', 'none')
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
% ThetaSign = {};
% figure()
% 
% for t=1:ntypes;
%     
%     for i=1:size(MSD_Results{t}.ma.tracks,1);
%         
%         tracktemp = MSD_Results{t}.ma.tracks{i}; %In um.
%         
%         [Angles{t}{i}, Lag_times{t}{i}] = SMT_Angles_vs_LagTime(tracktemp,Frame_interval);
%     
%     end
%         
% f_asim{t} = (length(find(abs(Angle_all{t}*180/pi) >= 150)))/(length(find(abs(Angle_all{t}*180/pi) <= 30)));%Calculate asimetry coefficient for each track.
%  
% end
% 
% 
%  figure()
%  for t=1:ntypes;
%      
%    for i=1:12;
%  bin_size = 15; %Bin size in miliseconds.    
%  range_local = [bin_size*i - bin_size; bin_size*i];
%  
%  jumps_local_idx = find(mean_jump_vector{t} >=  range_local(1) & mean_jump_vector{t} <  range_local(2)); 
%  angles_local{t} = angle{t}(jumps_local_idx);
%  angles_local{t} = cat(1,angles_local{t}{:});
%  
%  f_asim_local{t}(i) = (length(find(abs(angles_local{t}*180/pi) >= 150)))/(length(find(abs(angles_local{t}*180/pi) <= 30)));%Calculate asimetry coefficient for each track.
% 
%  
%    end  
%      
% plot((bin_size/2):bin_size:bin_size*i,f_asim_local{t});
% hold on
% title('f180/0 vs Time lag of the jumps','FontSize',24,'FontWeight','bold'); 
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




 
 
%% Analyze and compare the radius of gyration evolution over time
figure()
 for t=1:ntypes;
    
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
    end
    
            %Calculate ensemble average of Rg from all particles for N time steps.
            for n=1:20;
            idx = find( (cellfun(@(x) length(x), Radius_of_gyration{t})) >= n);    
            Radius_of_gyration_mean{t}(n) = mean(cellfun(@(x) x(n), Radius_of_gyration{t}(idx)));
            end
    
    plot(Frame_interval:Frame_interval:20*Frame_interval,Radius_of_gyration_mean{t},'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.');
    hold on

 end

title('Time evolution of emsemble average of radius of gyration','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  
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
    
    
    
   binCtrs = 0:0.05:4; %Bin centers. (depends on your data)
   n = length(Residence_Times{t});
   counts = hist(Residence_Times{t},binCtrs);
   prob = counts / (n * 0.05);
   H = bar(binCtrs,prob,'hist');
   set(H,'Facecolor',Color(t,:),'FaceAlpha',0.6);
   hold on;
   
xlabel('Residence time (s)','FontSize',20,'FontWeight','bold');
ylabel('Frequency','FontSize',20,'FontWeight','bold');
title('Histogram of Residence Times','FontSize',24,'FontWeight','bold'); 
set(gca,'FontSize',20,'FontWeight','bold');  

   
end



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
[Esp_coeff{t},Esp_Sigma{t}, Esp_fit{t}] = ExpDecay_2Cmp_fit(survival{t}, [1 0.1]);
plot(Esp_fit{t}(:,1),Esp_fit{t}(:,2),'color',Color(t,:),'LineWidth',2.2,'LineStyle','-.')

xlim([0 60]);

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


%% Save the results from the MSD Analysis

Results = {};
Results{1,2} = 'Radius of confinement (Confined circle model) ensemble (?m)';
Results{1,3} = 'Difussion Coefficient (Confined circle model) ensemble (?m2/s)';
Results{1,4} = 'Confidence Intervals (Confined circle model) ensemble (?m and s)';
Results{1,5} = 'Radius of confinement (Confined circle model) all tracks (?m)';
Results{1,6} = 'Difussion Coefficient (Confined circle model) all tracks (?m2/s)';
Results{1,7} = 'Confidence Intervals (Confined circle model) all tracks (?m and s)';
Results{1,8} = 'Radius of gyration (Ensemble Mean) (?m)';
Results{1,9} = 'Radius of gyration (All Tracks) (?m)';
Results{1,10} = 'Difussion Coefficient (TE-MSD emsemble fit) (?m2/s)';
Results{1,11} = 'Difussion Coefficient CI (TE-MSD emsemble fit) (?m2/s)';
Results{1,12} = 'Difussion Coefficient (T-MSD fit) (?m2/s)';
Results{1,13} = 'Probability Unimodal distribution of Diffusion Coefficients';
Results{1,14} = 'Residence times (s)';
Results{1,15} = 'Maximum distance travelled (?m)';
Results{1,16} = 'Angle between segments (radians)';
Results{1,17} = 'Jumps length distribution (?m)';
% Results{1,18} = '';
% Results{1,19} = '';
% Results{1,20} = '';

for t=1:ntypes;
Results{t+1,1} = strcat(files{t}.dataname(1:end-9));

Results{t+1,2} = f_confined_circ_diff_R{t};
Results{t+1,3} = f_confined_circ_diff_D{t};
Results{t+1,4} = conf{t};
Results{t+1,5} = f_confined_circ_diff_R_all{t}(iidx{t});
Results{t+1,6} = f_confined_circ_diff_D_all{t}(iidx{t});
Results{t+1,7} = conf_all{t}(iidx{t});
Results{t+1,8} = Radius_of_gyration_mean{t};
Results{t+1,9} = Radius_of_gyration{t};
Results{t+1,10} = D_ensemble(t);
Results{t+1,11} = [D_emsemble_lower_ci(t) D_emsemble_lower_ci(t)];
Results{t+1,12} = D_individual{t};
Results{t+1,13} = probability_unimodal_D(t);
Results{t+1,14} = Residence_Times{t};
Results{t+1,15} = max_dist{t};
Results{t+1,16} = Angle_all{t};
Results{t+1,17} = jumps_vector_cat{t}/1000;
% Results{t+1,18} = ;
% Results{t+1,19} = ;
% Results{t+1,20} = ;

end

if Save_Results == 1;
fprintf('\n');     
fprintf('Saving...\n');

dir = uigetdir(startlocn{1});
% mkdir(strcat(files{1}.data{1,2},'\Compare_MSD_Results_'));
% save(fullfile(strcat(files{1}.data{1,2},'\MSD_Results_')),'Results');
save(fullfile(dir,strcat(File_name,'.mat')),'Results');

end



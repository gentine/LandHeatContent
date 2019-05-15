
% plot final selected FLUXNET files in higher quality for publication

clear all;
close all;
DIR = '../../../Datasets/FLUXNET/';
Nmin_sites = 40; % number of sites considered
Gmean_threshold = 0; % for p values
N = 100000; % bootstrap number

name_output_file    =   [DIR 'Stats_G_trend'];
% load([name_output_file '.mat'],'S')
% load([name_output_file '_yearly.mat'],'S')
header = 'SiteID	lon	lat	Gtable	Gtable_std	Gtable_closure	Gtable_closure_std	Gtable_CV	years_ 1	years_ 2	years_ 3	years_ 4	years_ 5	years_ 6	years_ 7	years_ 8	years_ 9	years_10	years_11	years_12	years_13	years_14	years_15	years_16	years_17	years_18	years_19	years_20	years_21	years_22	years_23	years_24	years_25	years_26	Gseas_amp'; 
fileID = fopen([name_output_file '.csv'],'rt'); % 35 columns
C = textscan(fileID, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter',',', 'HeaderLines',1, ...
    'MultipleDelimsAsOne',true, 'CollectOutput',false);
fclose(fileID);

FLUXNET2015_years = 1990:2015;
Nyears  = length(FLUXNET2015_years);
Nsites  = length(C{2});
G_years = NaN*ones(Nsites,Nyears);
for i=1:Nsites% loop on sites
    for j=1:Nyears
        G_years(i,j) = C{8+j}(i); 
    end
end
G_years(G_years==0) = NaN;

%Bootstrap data

random = int64(fix(Nsites*rand(fix(1*Nsites),N)))+1;
for i=1:Nyears
    G_vect      = G_years(:,i);
    Grand_mean  = nanmean(G_vect(random),1);
    Grand_std(i)= nanstd(Grand_mean);
end
  

random = int64(fix(Nsites*rand(fix(1*Nsites),N)))+1;
Grand = zeros(size(random,1),size(random,2),Nyears);
for i=1:Nyears
    G_vect         = G_years(:,i);
    Grand(:,:,i)   = G_vect(random);
end
Grand_mean  = squeeze(nanmean(Grand,1));
Grand_std   = nanstd(Grand_mean,0,1); 


% now copmute long-term estaimte when there were enough sites
number_sites = sum(~isnan(G_years),1);
years = FLUXNET2015_years;
ind = find(number_sites>=Nmin_sites);%find(years>=2004);
G_years_mean = nanmean(nanmean(G_years(:,ind)))
Grand = permute(Grand(:,:,ind),[2 1 3]);
Grand = reshape(Grand,[size(Grand,1) size(Grand,2)*size(Grand,3)]);
Grand_meantmp  = nanmean(Grand,2);  
Grand_mean_std = nanstd(Grand_meantmp,0,1)  
Grand_mean_mean= nanmean(Grand_meantmp,1);  

%p_value_Gpositive = sum(Grand_mean(:,ind)<Gmean_threshold,1)/N;

    

Gmean = nanmean(G_years,1);


figure 
subplot(2,1,1)
shadedplot(years, Gmean-Grand_std, Gmean+Grand_std,'b');
alpha(.15);
hold on
plot(years,Gmean,'Color',[0.4 0.4 0.7],'linewidth',2);
plot(years,zeros(size(years)));
xlabel('Years','Interpreter','latex')
ylabel('Across-site $\overline{G} \ {\rm [W/m^2]}$','Interpreter','latex')
set(gca,'Fontsize',16);
grid on
xlim([1997 2014])
subplot(2,1,2)
plot(years,number_sites,'Color',[0.4 0.4 0.7],'linewidth',2);
set(gca,'Fontsize',16);
xlabel('Years','Interpreter','latex')
ylabel('Number of measuring sites ','Interpreter','latex')
grid on 
xlim([1997 2014])





'end'
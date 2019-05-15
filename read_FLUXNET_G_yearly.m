
% plot final selected FLUXNET files in higher quality for publication

clear all;
close all;
DIR = '../../../Datasets/FLUXNET/';
system(['rm ' DIR 'tempdir/*'])
name_output_file    =   [DIR 'Stats_G_trend'];
use_hourly_daily_monthly = 0; % 0-> hourly, 1->daily, 2-> monthly
 
% Local Fluxnet measurements:
% Fluxnet stations information
Fluxnet_info = readtable([DIR 'Fluxnet_names.xlsx']);  
max_CV = 300/100;
exclude_single_year_sites = 0; % true or fals exclude single year sites
% min_CV = 50/100; 
FLUXNET2015_years = 1990:2015;

% statistics table
T = table;
index_table = 0;           
S = {};
year_table = [];

list_files = dir([DIR 'FULLSET/*.zip']);
for i=1:length(list_files)
    tempdir = [DIR 'tempdir/'];
    mkdir(tempdir);

    zipFilename = list_files(i).name;
    PREFIX   = zipFilename(1:3);
    
    if(strcmp(PREFIX,'FLX')) 
        
        %% Unzip files into "tempdir" directory (will create the directory if needed)
        unzip([DIR 'FULLSET/' zipFilename], tempdir)

        % read csv and xls files
        filename_csv = dir([DIR 'tempdir/*.csv']);
        if(~isempty(filename_csv))
            % Fluxnet 2015 dataset
            title_name = filename_csv(1).name(1:end-4);
            for ii=1:length(title_name) 
              if(strcmp(title_name(ii),'_'))
                  title_name(ii) = ' ';
              end
            end 
            sitecode     = title_name(5:10)
           if(strcmp(sitecode,'AU-Tum'))
              'test'
           end

           if(use_hourly_daily_monthly==0) % use hourly data
                filename_csv = dir([DIR 'tempdir/*FULLSET_HH*.csv']);
           end
           if(use_hourly_daily_monthly==1) % daily data
                filename_csv = dir([DIR 'tempdir/*FULLSET_DD*.csv']);
           end
           if(use_hourly_daily_monthly==2) % monthly data
                filename_csv = dir([DIR 'tempdir/*FULLSET_MM*.csv']);
           end
           if(~isempty(filename_csv))
               csv_file = readtable([DIR 'tempdir/' filename_csv(1).name]); 


               if(~any(strcmp('G_F_MDS',fieldnames(csv_file))) || ~any(strcmp('NETRAD',fieldnames(csv_file)))) % check if G exists in csv file
                    % do nothing
               else
                   G            = csv_file.G_F_MDS;
                   G_quality    = csv_file.G_F_MDS_QC;
                   G(abs(G_quality-1)<1e-6 | abs(G_quality-2)<1e-6 | abs(G_quality-3)<1e-6 | G<-900 | G>900) = nan; %  keep measured and god quality gap filled,	0 = measured; 1 = good quality gapfill; 2 = medium; 3 = poor
                   H            = csv_file.H_CORR;
                   H_quality    = csv_file.H_F_MDS_QC;
                   LE           = csv_file.LE_CORR; 
                   LE_quality   = csv_file.LE_F_MDS_QC; 
                   LE_Rn        = nansum(LE)./nansum(H+LE);
                   Rn           = csv_file.NETRAD;
                   G_closure    = Rn-H-LE;
                   G_closure_quality = max(LE_quality,H_quality);
                   G_closure(abs(G_closure_quality-1)<1e-6 | G_closure<-900 | G_closure>900) = nan;
                   if(use_hourly_daily_monthly==0) % use hourly data
                        time         = num2str(csv_file.TIMESTAMP_START);
                   else % use daily or montlhy output
                        time         = num2str(csv_file.TIMESTAMP);
                   end
                   GEP         = max(csv_file.GPP_NT_VUT_REF,0)/2.1; % warning some sites have 2.1 scaling becasue of weight

                   year     = str2num(time(:,1:4));
                   month    = str2num(time(:,5:6));
                   if(use_hourly_daily_monthly==2)
                        day      = zeros(size(month));    
                   else
                        day      = str2num(time(:,7:8));
                   end
                   if(use_hourly_daily_monthly==0)
                       hours    = str2num(time(:,9:10));
                       minutes  = str2num(time(:,11:12));
                   else
                      hours     = 0; 
                      minutes   = 0;
                   end

                   DateNumber   = datenum(year,month,day,hours,minutes,0);
                   DateRef      = datenum(0,0,0,0,0,0);
                   OneYear      = (datenum(12,0,0,0,0,0)-datenum(0,0,0,0,0,0))/12;
                   DateNumber   = (DateNumber - DateRef)/OneYear;
                   if(round(DateNumber(end))>=1996) % we need to haev FLuxnet data during the retrieval
                       % find Site code anme and latitude/longitude information
                       test     = 0;
                       n_sites  = 0;
                       while n_sites<length(Fluxnet_info.SiteCode) && test==0
                            n_sites = n_sites + 1;
                            test = strcmp(Fluxnet_info.SiteCode(n_sites),sitecode);
                       end
                       lat = Fluxnet_info.Latitude(n_sites);
                       lon = Fluxnet_info.Longitude(n_sites);

                       nb_years = round(DateNumber(end))-round(DateNumber(1))+1;

                       % monthly average values 
                       G_monthly    = nan(12,nb_years);
                       G_closure_monthly    = nan(12,nb_years);
                       GEP_monthly  = nan(12,nb_years);  

                       for yy=1:nb_years
                           for mm=1:12 
                               ind = find(year==year(1)+yy-1 & month==month(1)+mm-1); 
                               % firt copmute mean diurnal cycle
                               if(~isempty(ind))
                                   if(use_hourly_daily_monthly==0)% hourly time scale
                                       G_daily = nan(1,24);
                                       G_closure_daily = nan(1,24);
                                       GEP_daily = nan(1,24);
                                       for hh=0:23
                                           ind_hour = find(hours(ind) == hh);
                                           G_daily(hh+1)   = nanmean(G(ind(ind_hour)));
                                           G_closure_daily(hh+1)   = nanmean(G_closure(ind(ind_hour)));
                                           GEP_daily(hh+1) = nanmean(GEP(ind(ind_hour)));
                                       end
                                       G_monthly(mm,yy)         = nanmean(G_daily); % only then take the average diurnal cycle
                                       G_closure_monthly(mm,yy) = nanmean(G_closure_daily); % only then take the average diurnal cycle
                                       GEP_monthly(mm,yy)       = nanmean(GEP_daily);
                                   else
                                       G_monthly(mm,yy)         = nanmean(G(ind)); % only then take the average diurnal cycle
                                       G_closure_monthly(mm,yy) = nanmean(G_closure(ind)); % only then take the average diurnal cycle
                                       GEP_monthly(mm,yy)       = nanmean(GEP(ind));
                                   end
                               end
                           end
                       end

                       GEP_monthly(GEP_monthly<1e-6) = NaN;
                       mean_G  = nanmean(G_monthly,2)'; % mean seasonal cycle
                       sigma_G = nanstd(G_monthly,1,2)'; 
                       mean_G_closure  = nanmean(G_closure_monthly,2)';
                       sigma_G_closure = nanstd(G_closure_monthly,1,2)'; 
    %                    mean_G  = [mean_G(end) mean_G];
                       sigma_G = [sigma_G(end) sigma_G]; 
                       sigma_G_closure = [sigma_G_closure(end) sigma_G_closure]; 
                       mean_GEP  = nanmean(GEP_monthly,2)';
                       sigma_GEP = nanstd(GEP_monthly,1,2)'; 
    %                    mean_GEP  = [mean_GEP(end) mean_GEP];
                       sigma_GEP = [sigma_GEP(end) sigma_GEP]; 
                       max_GEP  = max(mean_GEP);
                       months_summer = find(mean_G>0.6*max_GEP);
                       months = round(DateNumber(1))+1/12:1/12:round(DateNumber(end))+1;
                       G_yearly = mean(G_monthly,1); % yearly averages do not accept nan values!

                       % copmute mean and std deviation over current months
                       years = round(DateNumber(1)):round(DateNumber(end));

                       % mean sesaonal values across years
                       G_seasonal_mean = mean(G_monthly(months_summer,:),1);
                       G_seasonal_mean_std = std(G_monthly(months_summer,:),0,1);

                       title_name = title_name(5:11); 


                       if(~isnan(nanmean(G_monthly(:))))

                           h = figure(1)
                           subplot(2,1,1)
                           plot(DateNumber,G)
                           xlabel('Time');
                           ylabel('G (W/m^2)');
                           title(['Hourly ' title_name])
                           set(gca,'FontSize',18)
                           xlim([year(1) year(end)+1])
                           subplot(2,1,2)
                           plot(year(1)+(month(1)-1)/12:1/12:year(1)+(length(G_monthly(:))-1)/12,G_monthly(:))
                           xlabel('Months');
                           ylabel('G (W/m^2)');
                           title(['Monthly ' title_name])
                           set(gca,'FontSize',18)
                           xlim([year(1) year(end)+1])
                           saveas(h,[DIR '/DeltaG/' title_name],'jpg')



                %            T.lon_table          = lon_table'; 
                %            T.lat_table          = lat_table'; 

                           % dirst copmute coefficient of vairaition to find
                           % outlier years
                           mean_Gs             = mean(G_monthly,1);
                           Gtablemean          = nanmean(mean_Gs,2); % use mean or median, nanmedian
                           nonans              = ~isnan(mean(G_monthly,1)); % need the entire year to assess
                           mean_Gs             = mean(G_monthly,1); 
                           std_mean_Gs         = nanstd(mean_Gs);
                           mean_Gs_closure     = mean(G_closure_monthly,1);
                           std_mean_Gs_closure = nanstd(mean_Gs_closure);
                           CV                  = std_mean_Gs./Gtablemean;
                           variations          = abs(mean_Gs/Gtablemean);  
                           % exclude yearly outliers 
                           test                = ~(variations>max_CV | isnan(variations) | isnan(mean_Gs));

                           % now we have a good interannual one
                           if(~isnan(nanstd(mean_Gs(test))) && nanmean(mean_Gs(test))~=0)
                                if(mean(CV)~=0 || ~exclude_single_year_sites) % exclude sites with only one year
%                                     if(nanstd(mean_Gs(test))==0) % if only one data point
%                                        index_table = 1;
%                                        Gtable(index_table) = nan;
%                                        Gtable_CV(index_table) = 0;
%                                        Gtable_closure(index_table)    = nan;
%                                     else 
                                       % TODO years are not correct only
                                       % beginning and end of site but not
                                       % actual measurements
                                       index_table = index_table+1;
                                       SiteID{index_table}   = sitecode;
%                                        year_beg(index_table) = year(1);
%                                        year_end(index_table) = year(end);
%                                        month_beg(index_table)= month(1);
%                                        month_end(index_table)= month(end);
%                                        day_beg(index_table) = day(1);
%                                        day_end(index_table) = day(end);
                                       year_beg             = year(1);
                                       year_end             = year(end);
%                                        month_beg(index_table)= month(1);
%                                        month_end(index_table)= month(end);
%                                        day_beg(index_table) = day(1);
%                                        day_end(index_table) = day(end);
                                       lat_fluxnet(index_table) = lat;
                                       lon_fluxnet(index_table) = lon;
                                       year_vector           = zeros(size(FLUXNET2015_years));
                                       years_site            = year_beg:year_end;
                                       good_years            = years_site(nonans(1:end-1));
                                       year_vector(intersect(FLUXNET2015_years,good_years)-FLUXNET2015_years(1)+1)=mean_Gs(~isnan(mean_Gs)); % save value
                                       % save statistics into file 
                                       clear T;
                                       T            = table;
                                       T.SiteID     = SiteID'; 
                                       T.lon        = lon_fluxnet';
                                       T.lat        = lat_fluxnet';
                                       
%                                        T.year_beg   = year_beg'; 
%                                        T.year_end   = year_end';
%                                        T.month_beg  = month_beg';
%                                        T.month_end  = month_end';
%                                        T.day_beg  = day_beg';
%                                        T.day_end  = day_end';
                                       Gtable(index_table)      = nanmean(mean_Gs(test));
                                       Gtable_std(index_table)  = nanstd(mean_Gs(test));
                                       Gtable_CV(index_table)   = Gtable_std(index_table)/Gtable(index_table);
                                       Gtable_closure(index_table)      = nanmean(mean_Gs_closure(test));
                                       Gtable_closure_std(index_table)  = nanstd(mean_Gs_closure(test));
                                       Nyears(index_table)      = sum(~isnan(mean_Gs));
                                       LE_Rn_table(index_table) = LE_Rn;
                                       T.Gtable                 = Gtable';
                                       T.Gtable_std             = Gtable_std'; 
                                       T.Gtable_closure         = Gtable_closure';
                                       T.Gtable_closure_std     = Gtable_closure_std'; 
                                       T.Gtable_CV              = Gtable_CV';
    %                                    T.Gtable_CV              = LE_Rn_table';
%                                        T.Nyears                 = sum(year_vector)';
                                       year_table   = [year_table;year_vector];
                                       T.years      = year_table;


                                       % seasonal cycle analysis
                                       vect = G_monthly(:);
                                       time = 1:1:length(vect);
                                       [Plomb,flomb] = plomb(vect(~isnan(vect)),time(~isnan(vect)));

        %                                 figure
        %                                 plot(flomb,Plomb)
        %                                 grid
        %                                 xlabel('Frequency (cycles/months)')
        %                                 ylabel('Power (dBW)')
                                       [M,ind_best] = min(abs(flomb-1/12)); % seasonal cycle one cycle per 12 months
                                       ind_best     = ind_best(1);
                                       Plomb_best   = Plomb(ind_best);
                                       Plomb_best   = 2*log10(Plomb_best);
                                       Gseas_amp(index_table) = Plomb_best;
                                       T.Gseas_amp  = Gseas_amp';
        %                                figure
        %                                plot(time,vect);
        %                                hold all;
        %                                plot(time,2*log10(Plomb_best)*sin(2*pi*1/12*time-pi/2))


                                       % save yearly values
                                       s.SiteID     = sitecode;
                                       s.year_beg   = year(1);
                                       s.year_end   = year(end);
                                       s.lat        = lat;
                                       s.lon        = lon;
                                       s.G_yearly   = G_yearly;
                                       nn           = length(S);
                                       S{nn+1,1}    = s;
                                       clear s;
%                                     end

                                    writetable(T,[name_output_file '.csv']);
                                    writetable(T,[name_output_file '.xls']);
                                    
                                    save([name_output_file '.mat'],'T')
                                    save([name_output_file '_yearly.mat'],'S')
                                end
                           end
                        end


                   end 
               end
            %     filename_xls = dir('tempdir/*.xls');
            %     if(~isempty(filename_xls)) 
            %        xls_file = readtable(['tempdir/' filename_xls.name]); 
            %     end 

            end

            close all;
            %% Delete the "tempdir" directory and all of its subdirectories
        %     cd ..;
            rmdir(tempdir, 's')
        end
    end
end
 
 





'end'
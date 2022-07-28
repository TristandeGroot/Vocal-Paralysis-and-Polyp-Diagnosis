%ELEN3014 Assignment
%Group Members: Keri-Lee Carstens, Zunaid Valodia, Tristan de Groot
%Signal Analysis of Polyps and Vocal Paralysis

close all;
clear all;

Start_Stop = 1;

disp('This program can be used to analyse the vocal behaviour of patients with Paralysis and Polyps Individually.');
disp('Follow the prompts to access the analysis for each patient.');

namesParalysis = {'P1-A.wav', 'P1-B.wav', 'P1-C.wav','P2-A.wav', 'P2-B.wav', 'P2-C.wav','P3-A.wav', 'P3-B.wav', 'P3-C.wav','P4-A.wav', 'P4-B.wav', 'P4-C.wav','P5-A.wav', 'P5-B.wav', 'P5-C.wav','P6-A.wav', 'P6-B.wav', 'P6-C.wav','P7-A.wav', 'P7-B.wav', 'P7-C.wav','P8-A.wav', 'P8-B.wav', 'P8-C.wav','P9-A.wav', 'P9-B.wav', 'P9-C.wav','P10-A.wav', 'P10-B.wav', 'P10-C.wav','P11-A.wav', 'P11-B.wav', 'P11-C.wav'};
namesPolyps = {'1A.wav', '1B.wav', '1C.wav','2A.wav', '2B.wav', '2C.wav','3A.wav', '3B.wav', '3C.wav','4A.wav', '4B.wav', '4C.wav','5A.wav', '5B.wav', '5C.wav','6A.wav', '6B.wav', '6C.wav','7A.wav', '7B.wav', '7C.wav','8A.wav', '8B.wav', '8C.wav','9A.wav', '9B.wav', '9C.wav','10A.wav', '10B.wav', '10C.wav','11A.wav', '11B.wav', '11C.wav'};
subplot_cols = 2;
subplot_rows = 3;
k2=1;

while Start_Stop == 1;

%promt = 'Would you like to analyse the patients with Polyps (input = 1) or Paralysis (input = 2)?';
PolorPar = input('Would you like to analyse the patients with Polyps or Paralysis ?','s')

if strcmp(PolorPar,'Paralysis')
    names = namesParalysis;
    wave = cell(size(namesParalysis));
    fs = cell(size(namesParalysis));
    PeakThreshold = 0.02;
    
    Patient = input('Which patient (1-11) would you like to analyse?','s')
    
    if strcmp(Patient,'1')
    position=1;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'2')
    position=4;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'3')
    position=7;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'4')
    position=10;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'5')
    position=13;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'6')
    position=16;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'7')
    position=19;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'8')
    position=22;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'9')
    position=25;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'10')
    position=28;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'11')
    position=31;
    PatientSex='Female';
    end 
end   

if strcmp(PolorPar,'Polyps')
    names = namesPolyps;
    wave = cell(size(namesPolyps));
    fs = cell(size(namesPolyps));
    PeakThreshold = 0.05;
    Patient = input('Which patient (1-11) would you like to analyse?','s')
    
    if strcmp(Patient,'1')
    position=1;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'2')
    position=4;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'3')
    position=7;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'4')
    position=10;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'5')
    position=13;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'6')
    position=16;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'7')
    position=19;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'8')
    position=22;
    PatientSex=' Male ';
    end 
    
    if strcmp(Patient,'9')
    position=25;
    PatientSex=' Male '
    end 
    
    if strcmp(Patient,'10')
    position=28;
    PatientSex='Female';
    end 
    
    if strcmp(Patient,'11')
    position=31;
    PatientSex=' Male ';
    end 
end    
   
j=2; % position in column
k=1; % position in row
i2=1; %control of storage position

for i= position:position+2
    
 
    [wave{i},fs{i}]=audioread(names{i});
    yOrigin = wave{i}(:,1);
    
    lentest=round(length(yOrigin)/2);
    lower=lentest-2500;
    higher=lentest+2500;
    
    q4=1;
    for q3=lower:higher
    y(q4)=yOrigin(q3,1);
    q4=q4+1;
    end
    
    filenamecut=[names{i} ' Cut Section.wav'];
    audiowrite(filenamecut,y,fs{i});
    
    dt = 1/fs{i};
    tOrigin = 0:dt:(length(yOrigin)*dt)-dt;
    t = 0:dt:(length(y)*dt)-dt;
    figure(1)
    subplot(subplot_rows,subplot_cols,k);
    plot(t,y); xlabel('Seconds'); ylabel('Amplitude');
    title(['File: ' names{i}]);
    grid on;
    
    y2 = fftshift(fft(y));
    n = length(y2);          % number of samples
    f = (-n/2:n/2-1)*(fs{i}/n);     % frequency range
    power = abs(y2).^2/n;    % power of the DFT
    dbs=10.*log10(power/(10)^(-12));
    figure(1)
    subplot(subplot_rows,subplot_cols,j);
     plot(f,power);
        ylabel('Power'); xlabel('Frequency (Hz)');
        title(['File: ' names{i}]);
        grid on;
    
        
       Sample{i2}=y;
       SampleFrequency{i2}=dbs;
       
       Stddevtime{i}=std(y);
       stddevfreq{i}=std(dbs);
       
       Normalizedy=normalize(yOrigin,'range',[-1 1]);
       
       
       %finding Shimmer
           peaks=findpeaks(Normalizedy,fs{i},'MinPeakHeight', 0.6);
            tem2=0;
            tem3=0;
            nlen=length(peaks);
            for j1=1:nlen-1
                tem1=abs(abs(peaks(j1))-abs(peaks(j1+1)));
                tem2=tem2+tem1;
            end
            tem3=sum(peaks(:,1));
       
            shim1=(1/(nlen-1))*tem2;
            shim2=(1/nlen)*tem3;
            shimmer(i)=abs(shim1/shim2)*100;
            
            %Finding Intensity
            peaksIntensity=findpeaks(yOrigin,'MinPeakHeight', PeakThreshold);
            Intensity(i)=mean(peaksIntensity(:,1));
            
            %Finding Jitter
            [pks,locs] = findpeaks(Normalizedy,fs{i},'MinPeakDistance',0.00286);
            %TimePerPos=tOrigin(locs(:,1));
            period = diff(locs);
            meanPeriod = mean(period);
            DiffPeriod = abs(diff(period));
            SumPeriodiff = sum(DiffPeriod);
            NoDiff = numel(DiffPeriod);
            meanJitter=(SumPeriodiff)/NoDiff;
            Jitter(i)= (meanJitter/meanPeriod)*100;
        
       
      fundementalFarray{i} = pitch(y',fs{i});
      
      MeanAmplitude(i2)=sum(y)/length(y);
       
      
      bw(i)=powerbw(y,fs{i});
      

            
            j=j+2;
            k=k+2;
            i2=i2+1;
             

end
    


   q3=1;
   fundementalF{q3}=mean(fundementalFarray{:,position});                               %Averaging Pitch function Frequency Values
   fundementalF{q3+1}=mean(fundementalFarray{:,position+1});
   fundementalF{q3+2}=mean(fundementalFarray{:,position+2});
   fundementalFMean=(fundementalF{q3}+fundementalF{q3+1}+fundementalF{q3+2})/3; %mean for the patient Fundemental Frequency
   bwave=(bw(position)+bw(position+1)+bw(position+2))/3;                                          %mean for the patient Bandwidth
   Stddevtimeave=(Stddevtime{position}+Stddevtime{position+1}+Stddevtime{position+2})/3;          %mean for the patient Standard Deviation
   ShimmerMean=(shimmer(position)+shimmer(position+1)+shimmer(position+2))/3;                       %mean of Shimmer
   FinalJitterMean=(Jitter(position)+Jitter(position+1)+Jitter(position+2))/3;                      %mean of Jitter
   IntensityMean=(Intensity(position)+Intensity(position+1)+Intensity(position+2))/3;               %mean of Intensity
   
   if PatientSex==' Male ';
       FrePolThres='<100';
       FreParaThres='>150';
   
   else PatientSex=='Female';
       FrePolThres='<180';
       FreParaThres='>250';
   end
 
   
   
    KeyFeatures={'Sex of Patient';'Fundamental Frequency (Hz)';'Standard Deviation Time Domain';'Bandwidth (Hz)';'Shimmer(%)';'Jitter(%)';'Intensity'};
    CalculatedValue={PatientSex;fundementalFMean;Stddevtimeave;bwave;ShimmerMean;FinalJitterMean;IntensityMean};%;jittC;shimC};
    PolypsValues={'_';FrePolThres;'>0.2';'<4.8';'>3';'>2';'>0.11'};      %Values that lie with range for Polyps
    ParalysisValues={'_';FreParaThres;'<0.04';'>7.5';'>3';'>2';'<0.08'};        %Values that lie with range for Paralysis
    Results = table(KeyFeatures,CalculatedValue,PolypsValues,ParalysisValues)
    
    %indicesConfirmationOfDiagnosis
   SCOREPOLYPS =0;
   SCOREPARALYSIS =0;
   if PatientSex==' Male '
       if fundementalFMean<100
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
       if Stddevtimeave>0.2
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
       if bwave<4.8
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
       if IntensityMean>0.11
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
   end
   
   if PatientSex=='Female' 
       if fundementalFMean<180 
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
       if Stddevtimeave>0.2
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
       if bwave<4.8
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
       if IntensityMean>0.11
           SCOREPOLYPS=SCOREPOLYPS+1;
       end
   end
    
    if PatientSex==' Male ' 
       if fundementalFMean>150
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
       if Stddevtimeave<0.04
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
       if bwave>7.5
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
       if IntensityMean<0.08
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
   end
    
    if PatientSex=='Female' 
       if fundementalFMean>250
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
       if Stddevtimeave<0.04
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
       if bwave>7.5
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
       if IntensityMean<0.08
           SCOREPARALYSIS=SCOREPARALYSIS+1;
       end
   end
    
   if SCOREPOLYPS >= 3;
       disp('At least half of indices confirm a positive diagnosis for Polyps');
   end
   if SCOREPARALYSIS >= 3;
       disp('At least half of indices confirm a positive diagnosis for Paralysis');
   end
   SCOREPOLYPS
   SCOREPARALYSIS


Continue = input('Would you like to continue (Yes or No)? ','s')

if strcmp(Continue,'Yes')
    Start_Stop = 1;
end
if strcmp(Continue,'No')
    Start_Stop = 0;
    
end
end
    
disp('Signal Analysis Program has ended.');    
    
    
   
   
   k1=1;
    i=1;
    
    
    subplot_cols = 2;
    subplot_rows = 3;
    
for k3=1:numel(namesPolyps)/3
    k=1;
    j=2;
    
for k2=1:numel(namesPolyps)/11
    [wave{i},fs{i}]=audioread(namesPolyps{i});
    yOriginPolyp = wave{i}(:,1);
    
    lentest=round(length(yOriginPolyp)/2);
    lower=lentest-2500;
    higher=lentest+2500;
    
    q4=1;
    for q3=lower:higher
    yPolyp(q4)=yOriginPolyp(q3,1);
    q4=q4+1;
    end
    
    
    dt = 1/fs{i};
    t = 0:dt:(length(yPolyp)*dt)-dt;
    tOriginPolyp = 0:dt:(length(yOriginPolyp)*dt)-dt;
    
    
    y2Polyp = fftshift(fft(yPolyp));
    n = length(y2Polyp);          % number of samples
    f = (-n/2:n/2-1)*(fs{i}/n);     % frequency range
    powerPolyp = abs(y2Polyp).^2/n;    % power of the DFT
    dbsPolyp=10.*log10(powerPolyp/(10)^(-12));
    
    
        k=k+2;
        j=j+2;
        
       SamplePolyp{i}=yPolyp;
       SampleFrequencyPolyp{i}=dbsPolyp;
       
       StddevtimePolyp{i}=std(yPolyp);
       stddevfreqPolyp{i}=std(dbsPolyp);
       
       
       
       %finding Shimmer
            
        NormalizedyPolyp=normalize(yOriginPolyp,'range',[-1 1]);
       
           peaksPolyp=findpeaks(NormalizedyPolyp,fs{i},'MinPeakHeight', 0.6);
            tem2=0;
            tem3=0;
            nlenPolyp=length(peaksPolyp);
            for j1=1:nlenPolyp-1
                tem1=abs(abs(peaksPolyp(j1))-abs(peaksPolyp(j1+1)));
                tem2=tem2+tem1;
            end
            tem3=sum(peaksPolyp(:,1));
       
            shim1Polyp=(1/(nlenPolyp-1))*tem2;
            shim2Polyp=(1/nlenPolyp)*tem3;
            ShimmerPolyp(i)=abs(shim1Polyp/shim2Polyp)*100;
       
                  %Finding Jitter
            [pksPolyp,locsPolyp] = findpeaks(NormalizedyPolyp,fs{i},'MinPeakDistance',0.00286);
            %TimePerPos=tOrigin(locs(:,1));
            periodPolyp = diff(locsPolyp);
            meanPeriodPolyp = mean(periodPolyp);
            DiffPeriodPolyp = abs(diff(periodPolyp));
            SumPeriodiffPolyp = sum(DiffPeriodPolyp);
            NoDiffPolyp = numel(DiffPeriodPolyp);
            meanJitterPolyp=(SumPeriodiffPolyp)/NoDiffPolyp;
            JitterPolyp(i)= (meanJitterPolyp/meanPeriodPolyp)*100; 
            
      %Finding Intensity
      peaksPolypI=findpeaks(yOriginPolyp,'MINPEAKHEIGHT', 0.05);
      IntensityPolyp(i)=mean(peaksPolypI(:,1));
       
       
      fundementalfPolyp{i} = pitch(yPolyp',fs{i});
      
      MeanAmplitudePolyp(i)=sum(yPolyp)/length(yPolyp);
            
        i=i+1;
        
    
end
k1=k1+1;
end    
    
    k1=1;
    i=1;
    

for k3=1:numel(namesParalysis)/3
    k=1;
    j=2;
    
for k2=1:numel(namesParalysis)/11
    [wave{i},fs{i}]=audioread(namesParalysis{i});
    yOriginParalysis = wave{i}(:,1);
    
    lentest=round(length(yOriginParalysis)/2);
    lower=lentest-2500;
    higher=lentest+2500;
    
    q4=1;
    for q3=lower:higher
    yParalysis(q4)=yOriginParalysis(q3,1);
    q4=q4+1;
    end
    
    
    dt = 1/fs{i};
    t = 0:dt:(length(yParalysis)*dt)-dt;
    tOrigin = 0:dt:(length(yOriginParalysis)*dt)-dt;
    
    y2Paralysis = fftshift(fft(yParalysis));
    n = length(y2Paralysis);          % number of samples
    f = (-n/2:n/2-1)*(fs{i}/n);     % frequency range
    powerParalysis = abs(y2Paralysis).^2/n;    % power of the DFT
    dbsParalysis=10.*log10(powerParalysis/(10)^(-12));

    
        k=k+2;
        j=j+2;
        
       SampleParalysis{i}=yParalysis;
       SampleFrequencyParalysis{i}=dbsParalysis;
       
       StddevtimeParalysis{i}=std(yParalysis);
       stddevfreqParalysis{i}=std(dbsParalysis);
       
       %finding Shimmer
            
        NormalizedyParalysis=normalize(yOriginParalysis,'range',[-1 1]);
       
           peaksParalysis=findpeaks(NormalizedyParalysis,fs{i},'MinPeakHeight', 0.6);
            tem2=0;
            tem3=0;
            nlenParalysis=length(peaksParalysis);
            for j1=1:nlenParalysis-1
                tem1=abs(abs(peaksParalysis(j1))-abs(peaksParalysis(j1+1)));
                tem2=tem2+tem1;
            end
            tem3=sum(peaksParalysis(:,1));
       
            shim1Paralysis=(1/(nlenParalysis-1))*tem2;
            shim2Paralysis=(1/nlenParalysis)*tem3;
            ShimmerParalysis(i)=abs(shim1Paralysis/shim2Paralysis)*100;
       
            %Finding Jitter
            [pksParalysis,locsParalysis] = findpeaks(NormalizedyParalysis,fs{i},'MinPeakDistance',0.00286);
            %TimePerPos=tOrigin(locs(:,1));
            periodParalysis = diff(locsParalysis);
            meanPeriodParalysis = mean(periodParalysis);
            DiffPeriodParalysis = abs(diff(periodParalysis));
            SumPeriodiffParalysis = sum(DiffPeriodParalysis);
            NoDiffParalysis = numel(DiffPeriodParalysis);
            meanJitterParalysis=(SumPeriodiffParalysis)/NoDiffParalysis;
            JitterParalysis(i)= (meanJitterParalysis/meanPeriodParalysis)*100;
       
      %Finding Intensity
      peaksParalysisI=findpeaks(yOriginParalysis,'MINPEAKHEIGHT', 0.02);
      IntensityParalysis(i)=mean(peaksParalysisI(:,1));
       
       
      fundementalfParalysis{i} = pitch(yParalysis',fs{i});
      
      MeanAmplitudeParalysis(i)=sum(yParalysis)/length(yParalysis);
             
        i=i+1;
        
    
end
k1=k1+1;
end

for q=1:33
    bwParalysis{q}=powerbw(SampleParalysis{q},fs{q});
    bwPolyp{q}=powerbw(SamplePolyp{q},fs{q});

end

q5=1;
for q4=1:33-2
   fundementalFPolyp{q4}=mean(fundementalfPolyp{:,q4});                               
   fundementalFPolyp{q4+1}=mean(fundementalfPolyp{:,q4+1});
   fundementalFPolyp{q4+2}=mean(fundementalfPolyp{:,q4+2});
   
   fundementalFParalysis{q4}=mean(fundementalfParalysis{:,q4});                               
   fundementalFParalysis{q4+1}=mean(fundementalfParalysis{:,q4+1});
   fundementalFParalysis{q4+2}=mean(fundementalfParalysis{:,q4+2});
   
   q4=q4+2;
end




q2=1;
for q1=1:33/3
   bwaveParalysis(q1)=(bwParalysis{q2}+bwParalysis{q2+1}+bwParalysis{q2+2})/3;
   bwavePolyp(q1)=(bwPolyp{q2}+bwPolyp{q2+1}+bwPolyp{q2+2})/3;
   
   StddevtimeaveParalysis(q1)=(StddevtimeParalysis{q2}+StddevtimeParalysis{q2+1}+StddevtimeParalysis{q2+2})/3;
   StddevtimeavePolyp(q1)= (StddevtimePolyp{q2}+StddevtimePolyp{q2+1}+StddevtimePolyp{q2+2})/3;
   
   ShimmerMeanPolyp(q1)=(ShimmerPolyp(q2)+ShimmerPolyp(q2+1)+ShimmerPolyp(q2+2))/3;
   ShimmerMeanParalysis(q1)=(ShimmerParalysis(q2)+ShimmerParalysis(q2+1)+ShimmerParalysis(q2+2))/3;
   
   IntensityMeanPolyp(q1)=(IntensityPolyp(q2)+IntensityPolyp(q2+1)+IntensityPolyp(q2+2))/3;
   IntensityMeanParalysis(q1)=(IntensityParalysis(q2)+IntensityParalysis(q2+1)+IntensityParalysis(q2+2))/3;
   
   fundementalFMeanParalysis(q1)=(fundementalFParalysis{q2}+fundementalFParalysis{q2+1}+fundementalFParalysis{q2+2})/3;
   fundementalFMeanPolyp(q1)=(fundementalFPolyp{q2}+fundementalFPolyp{q2+1}+fundementalFPolyp{q2+2})/3;
   
   JitterMeanPolyp(q1)=(JitterPolyp(q2)+JitterPolyp(q2+1)+JitterPolyp(q2+2))/3;
   JitterMeanParalysis(q1)=(JitterParalysis(q2)+JitterParalysis(q2+1)+JitterParalysis(q2+2))/3;
   
   q2=q2+2;
end


%figure(13)
%plot(IntensityMean);

%figure(14)
%plot(shimmer);

    figure(12)
    plot(StddevtimeavePolyp);
    hold on;
    plot(StddevtimeaveParalysis);
    title(['Standard Deviation Plot of all Patients-Paralysis:Red Polyps:Blue']);

    figure(13)
    plot(IntensityMeanPolyp);
    hold on;
    plot(IntensityMeanParalysis);
    title(['Intensity Plot of all Patients-Paralysis:Red Polyps:Blue']);
    
    figure(14)
    plot(bwavePolyp);
    hold on;
    plot(bwaveParalysis);
    title(['Bandwidth Plot of all Patients-Paralysis:Red Polyps:Blue']);
    
    figure(15)
    plot(fundementalFMeanPolyp);
    hold on;
    plot(fundementalFMeanParalysis);
    title(['Fundemental Frequency Plot of all Patients-Paralysis:Red Polyps:Blue']);
    
    figure(16)
    plot(ShimmerMeanPolyp);
    hold on;
    plot(ShimmerMeanParalysis);
    title(['Shimmer Plot of all Patients-Paralysis:Red Polyps:Blue']);
    
    figure(17)
    plot(JitterMeanPolyp);
    hold on;
    plot(JitterMeanParalysis);
    title(['Jitter Plot of all Patients-Paralysis:Red Polyps:Blue']);
  
  
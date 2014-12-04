
%% ECG Program, edited from ECG-LAb, Maria L Stavrinou, MArch 2008 %%
% Detects the RR interval and calculates the RR frequency 
% Input, the ECG Signal,in one-dimensional array. It asks the Sampling frequency %%
%% doulevei 21 sEPTEMBRIOU 2008/.  added the total power at 14 June 2009
%% 20 Iouliou 2009
%% 22 Augoustou 2009 
%% 1st September 2009
%% 8 November 2009  %% it removes the ectopic beat !!!!!! PSD based on
%% welch periodogram
%% used it 22 December 2009 for the patients group
%%% the ex ECGs_new16TXTB_suc.m
%%% Now its automated!!!!!!!!

clear all
global samplerate_ecg
%repeats=input('For how many times? /datasets? ');
samplerate_ecg=5000; %#ok<NASGU> %%input( ' Write Sampling Frequency:   ');

files=uipickfiles('num', [], 'out', 'str');
a7=size(files, 1);
for q=1:a7
    disp(q)
    current_name=files(q).name(end-13:end-10);
    pathname1=files(q).name(1:end-18);
    dataecg=load(files(q).name);
    name=current_name;
    ecg_sinal=dataecg;
   
    
     cd(pathname1)
        mkdir(name)
        cd(name)
        
        figure; plot(ecg_sinal); title(['primary ECG ' name]); 
        saveas(gcf, 'ECGsinalprimary', 'fig')
        ecg_sinal1=ecg_sinal;
        save ecg_sinal1 ecg_sinal1
%     %filtering the  data
%     h=fdesign.lowpass('Fp,Fst,Ap,Ast',0.0180,0.02,1,60);
%     d=design(h,'equiripple'); %Lowpass FIR filter
%     y=filtfilt(d.Numerator,1,ecg_sinal); %zero-phase filtering
% %     %y1=filter(d.Numerator,1,ecg_sinal); %conventional filtering
%     subplot(211);
%     plot(y);
%     title('Filtered Waveforms');
%     legend('Zero-phase Filtering','Conventional Filtering');
%     subplot(212);
%     plot(ecg_sinal);
%     title('Original Waveform'); 
%     saveas(gcf, 'Filtered_ECG', 'fig')
%     % end of fitering 
%     ecg_sinal=[]; 
%     ecg_sinal=y;
        
    clear current_name
        %samplerate_ecg=input( ' Write Sampling Frequency:   ');

        % repeats=input('For how many times? /datasets? ')
        % for q=1:repeats
        %     

        %% ECG Ondar 2
        ecg_filtrado=ecg_sinal;
       

        %regiao de busca do pico da onda R
        %(valor maior da mais imunidade a ruido)
        regiao_de_busca=round(0.070*samplerate_ecg);  % MLS Search region 

        %nivel em que o ecg filtrado é considerado onda R
        ganho_limiar=0.15; % MLS profit threshold

        %periodo de comparacao: 2 segundos
        limiar_comparacao=round(2*samplerate_ecg); % threshold comparison

        %intervalo minimo entre marcacoes *markings*
        %antes estava trabalhando com 350ms **MLS before it was working
        %with 350ms
        salto_entre_pulsos=round(0.350*samplerate_ecg);  %jump between pulses

        %busca ondas R de 10 em 10 msecs **msecs searchs waves R of 10 in
        %10
        salto_10msec=round(0.01*samplerate_ecg); 

        %depois que acha o limiar, volta esse tanto *later that it finds
        %the threshold, comes back this in such a way
        volta_amostras=round(0.030*samplerate_ecg); % return samples

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %inicializacoes de variaveis**inicializacoes of 0 variable
        tamanho=length(ecg_sinal); % tamanho=size
        n=1;

        %inicia ondasr com 'empty'**it initiates ondasr with ' empty'
        ondasr=[]; % waves
        ecg_ondar=[];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Localiza as regioes onde estao as ondas R
        %It locates the regions where estao waves R

        %repete para cada onda R **it repeats for each wave R
        while n<tamanho,

           %calcula o limiar para a onda R atual, baseado nos proximos 2 segundos
           %com ganho de 0.2 nao sao marcadas algumas ondas R no
           %sinal 124, posicao 300-336 ****it calculates the threshold for
           %current wave R, based in next the 2 seconds with 0.2 profit are not marked some waves R in signal 124, position 300-336
           if (n+limiar_comparacao)<=tamanho,
                limiar=ganho_limiar*max(abs(ecg_filtrado(n:n+limiar_comparacao)));
           else
               limiar=ganho_limiar*max(abs(ecg_filtrado(tamanho-limiar_comparacao:tamanho))); 
           end

           %procura a regiao da onda R, se encontrar...**it looks the
           %region of wave R, if to find…
           if ( (ecg_filtrado(n)>limiar) & (n<tamanho) ),

              %guarda o indice inicial da regiao da onda R **it keeps the
              %initial index of the region of wave R
              ondasr=[ondasr;n];

              %salta uns 200 msecs para evitar novos	 msecs jumps ones 200 to prevent new detonations in this exactly pulse.     
              %disparos neste mesmo pulso.
              n=n+salto_entre_pulsos;    

           %se nao foi encontrado, continua a busca de 10 em 10 msecs **if
           %it was not found, continues the search of 10 in 10 msecs 
           else
              n=n+salto_10msec;
              %n=n+1;
           end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % Localiza, na regiao, o pico da onda R
        %It locates, in the region, the peak of wave R

        %inicia constantes e variaveis **it initiates constants and 0
        %variable
        n=1;
        num_marcas=length(ondasr);
        marcas=[];

        %filtra o sinal com um passa altas para retirar a oscilação da linha de base
        %it filters the signal with one passes high to remove oscilação of the base line
        [B,A]= butter(4,1/(samplerate_ecg/2),'high'); %cria um butterworth de w0 = 1 Hz;
        ecg_filtrado=filtfilt(B,A,ecg_sinal); %filtfilt é uma filtragem de fase zero **a filtering of zero stage

        %volta um trecho do sinal para buscar a onda R**a stretch comes
        %back of the signal to search wave R
        ondasr=ondasr-volta_amostras;
        if ondasr(1)<1,ondasr(1)=1;end;

        %localiza os picos das ondas R**it locates the peaks of waves R
        for i=1:num_marcas,

           %acha o pico da regiao**it finds the peak of the region
           if tamanho>=ondasr(i)+regiao_de_busca,
              [valor,marca]=max(abs(ecg_filtrado(ondasr(i):ondasr(i)+regiao_de_busca)));   
           else
              [valor,marca]=max(abs(ecg_filtrado(ondasr(i):tamanho)));
           end

           %calcula a posição da marca e guarda o valor**it calculates
           %posição of the mark and keeps the value
           marca=marca+ondasr(i)-1;
           ecg_ondar=[ecg_ondar;marca];

        end

        if isempty(ondasr),
           ecg_ondar=-1;
        end

%         cd(pathname1)
%         mkdir(name)
%         cd(name)
%         save ecg_sinal ecg_sinal
%         figure; plot(ecg_sinal); title(['ECG ' name])
%         saveas(gcf, 'ECG', 'fig')

        %% Calcula IRR  **segundos=seconds
        segundos_sinal=(length(ecg_sinal)+1)/samplerate_ecg; %% Maria add in
        global samplerate_ecg

        tempoRR=0;
        intervaloRR=0;

        if length(ecg_ondar)>1,

            %cria o eixo do tempo**it creates the axle of the time
           tempo=0:1/samplerate_ecg:segundos_sinal-1/samplerate_ecg;

           %transforma as marcacoes em msec**it transforms the markings
           %into msec
           tempoRR=round(1000*tempo(ecg_ondar));

           %calcula os intervalos R-R**it calculates intervals R-R
           RR=tempoRR(1:length(tempoRR)-1);
           RR1=tempoRR(2:length(tempoRR));
           intervaloRR=RR1-RR;

           %transforma em segundos (a primeira marcacao nao é intervalo)
           %it transforms into seconds (the first marking not é interval)
            tempoRR=tempoRR(2:length(tempoRR))/1000;

           %arredonda os intervalos R-R
           intervaloRR=round(intervaloRR);

        else   
           intervaloRR=-1;
           tempoRR=-1;
        end

        %% End Calcula RR
        
        %% Ectopic Beat removal by Maria Stavrinou
        x=intervaloRR;
        stdX=std(intervaloRR);
        for i=2
            :(length(intervaloRR)-1) % xekiname apo to 2 gia g;itwsoume to x(-1)
            if x(i)>(mean(x)+3*std(x)) || x(i)<(mean(x)-3*std(x))
            x(i)=(x(i-1)+x(i+1))*0.5;
            end
        end
        figure; plot(intervaloRR); hold on; plot(x,'r')
        % trexoyme akoma mia fora ton algorithmo gia na beltiwsoume to meso oro kai
        % to std
        stdX=std(x);
        for i=2:(length(x)-1)
            if x(i)>(mean(x)+3*std(x)) || x(i)<(mean(x)-3*std(x))
                x(i)=(x(i-1)+x(i+1))*0.5;
            end
        end
        figure; plot(intervaloRR); ylim([290 1600]); hold on; plot(x,'r');ylim([400 1400]); title(['RR interval (blue original)' name])
        saveas(gcf, 'Corrected_RR', 'fig')
        intervaloRR=x;%% replace the intervaloRR with the new cleaned one
        
        %% End of Ectopic Beat removal 
        figure; plot(intervaloRR); title(['Interval RR' name]);ylabel('Duration of RR interval in msec'); xlabel('Number of RR Interval'); 
        save intervaloRR intervaloRR
        save tempoRR tempoRR %%% saved in the same directory as before 

	 saveas(gcf, 'RR_interpolated', 'fig')        
%% Some temporal aspects  
        min1=min(intervaloRR(1,:));
        save min1 min1;
        max1=max(intervaloRR(1,:));
        save max1 max1;
        mean1=mean(intervaloRR(1,:));
        save mean1 mean1
        std1=std(intervaloRR(1,:));
        save std1 std1
        numberRR=length(intervaloRR);
        save numberRR numberRR


        %% Frequency Domain Analysis
        % Interpolation
        % Interpolation
        Fsnew=4; 
        xi=1:(1/Fsnew):length(intervaloRR);
        tbefore=1:length(intervaloRR); 
        interpol_intervaloRR=interp1(tbefore, intervaloRR,xi, 'cubic');
        meanINTRR=mean(interpol_intervaloRR); %% It differs from mean1 LOOK! meanINTRR =  743.5458, mean1 =  743.6776
        %interpol2=detrend(interpol_intervaloRR);  % we detrend the dataset
        %so the mean appears to zero. 
        NN=length(interpol_intervaloRR);
        Freqresnew=Fsnew/NN %Fs= F(2)*length(interpol_intervaloRR) % H nea sampling frequency basei tis interpolation  2000 Hz
        save interpol_intervaloRR interpol_intervaloRR
        %save interpol2 interpol2 
        save Fsnew Fsnew
      
        
        % Power spectrum estimation 
        buffer=(interpol_intervaloRR-meanINTRR)/std(interpol_intervaloRR);
        Fs_buffer=Fsnew;
        % calculate length of window
        % calculate length of window
               
        p2=nextpow2(NN); % calculate length of window
        
        NFFT=2^p2;
        xnew=[buffer zeros(1,NFFT-length(buffer))];
        buffer=xnew; 
        clear xnew;
        F = Fs_buffer/2*linspace(0,1,NFFT/2);%rv(2*j),i); dati(1:131073,i);
        NN=[];
        NN =length(buffer);
        X=fft(buffer, NFFT)/NN;
        figure;plot(F,2*abs(X(1:NFFT/2))) 
        saveas(gcf, 'FFT','fig')
        
        P=(2*abs(X(1:NFFT/2)));
        hpsd=dspdata.psd(P, 'Fs', Fs_buffer);
        meanINTRR=mean(interpol_intervaloRR);
        Fw = hpsd.Frequencies;                 
        power=avgpower(hpsd)/(Fsnew); figure(23); plot(hpsd); title([' Total Power ' name]) 
        saveas(gcf, 'PSD','fig'); 
        Parsevalaki=(sum(buffer.^2))/(length(buffer)*Fsnew); % From Parseval 
  
        %% (Very) Low frequency LF 0.003 - 0.04 Hz
        freqVLF1=max(find(Fw<0.003));  %%%% indexfind1=max(find(F<freqrange(1)));
        freqVLF2=max(find(Fw<0.04));
        %freqVLF2=findLF2(end);
        NNLF=length(freqVLF1:freqVLF2);
        freqrangeVLF=[Fw(freqVLF1) Fw(freqVLF2)];
        %powerbandLF=sum(P(1:freqLF))/NNLF;
        powerbandVLF=avgpower(hpsd, freqrangeVLF)/(Fs_buffer);
        
        save powerbandVLF powerbandVLF

         %%  Low frequency LF 0.04 - 0.15 Hz
        findLF=find(Fw<0.15);
        freqLF=findLF(end);
        NNLF=length(1:freqLF);
        freqrangeLF=[Fw(freqVLF2) Fw(freqLF)];
        %powerbandLF=sum(P(1:freqLF))/NNLF;
        powerbandLF=avgpower(hpsd, freqrangeLF)/(Fs_buffer);
        save powerbandLF powerbandLF
       
        
        %% hIGH frequency hF 0.15 - 0.4 Hz
        freqHF=max(find(Fw<0.4));
        %freqHF=findHF(end);
        NNHF=length(freqLF:freqHF);
        freqrangeHF=[Fw(freqLF) Fw(freqHF)];
        powerbandHF=avgpower(hpsd, freqrangeHF)/(Fs_buffer);
        %powerbandHF=sum(P(freqLF:freqHF))/NNHF;
        %powerbandHF_psd=sum(psd(freqLF:freqHF));
%         figure; plot(F(freqLF:freqHF),P(freqLF:freqHF,1));%title([' EEG_channel ' s{i}]) %%% (1:469) (1:469)
%         set(gca, 'YLim', ylim);grid on; zoom on;title(['RR Frequency Spectrum' name ]); xlabel('Frequency (Hz)');
        save powerbandHF powerbandHF

        %% total power HF LF 0 - 0.15 Hz
        %% total power 0.003 - 0.4 Hz
        NNTF=length(freqVLF1:freqHF);
        %powerbandTotal=sum(P(1:freqHF))/NNTF;
        freqrangeT=[Fw(freqVLF1) Fw(freqHF)];
        powerbandTotal=avgpower(hpsd,freqrangeT)/(Fs_buffer);
        save powerbandTotal powerbandTotal

        
        %%% RATIO Chroni
        ratioch=(max1-min1)/mean1;
        
        %% save results on a text file
        fid=fopen('RR_Temporal_analysis_results.txt', 'wt');
        fprintf(fid, '      ÁíÜëõóç ôùí äéáóôçìÜôùí  RR\n\n');
        fprintf(fid, 'Number of RR intervals   %6.0f \n', numberRR);
		fprintf(fid, 'Time in minutes of ECG signal   %6.0f \n', length(ecg_sinal)/samplerate_ecg);
        fprintf(fid, 'min RR interval (msec) %6.2f\n', min1);
        fprintf(fid, 'max RR interval (msec)  %6.2f\n', max1);
        fprintf(fid, 'mean RR interval (msec)  %6.2f\n', mean1);
        fprintf(fid, 'STD RR interval (msec)  %6.2f\n', std1);
        fprintf(fid, 'ratio (max-min)/mean  %6.2f\n', ratioch);
        fprintf(fid, ' Very Low frequency VLF 0.003 - 0.04 Hz  %6.6f\n', powerbandVLF);
        fprintf(fid, ' Low frequency LF  0.04 - 0.15 Hz  %6.6f\n', powerbandLF);
        fprintf(fid, 'High frequency 0.15<f<0.4  Hz power  %6.6f\n', powerbandHF);
        fprintf(fid, 'Total band power 0 <f< 0.4 Hz % 6.6f\n', powerbandTotal);
        
        ALL_local(q,1)=str2num(name);
        ALL_local(q,2)=powerbandVLF;
        ALL_local(q,3)=powerbandLF;
        ALL_local(q,4)=powerbandHF;
        ALL_local(q,5)=powerbandTotal;
        ALL_local(q,6)=mean1;
        ALL_local(q,7)=std1;
        ALL_local(q,8)=ratioch;
        ALL_local(q,9)=min1;
        ALL_local(q,10)=max1;
        ALL_local(q,11)=numberRR;
        ALL_local(q,12)=power;
        ALL_local(q,13)=Parsevalaki;
        
        save ALL_local ALL_local
        cd ..
        ALL_controls(q,1)=str2num(name);
        ALL_controls(q,2)=powerbandVLF;
        ALL_controls(q,3)=powerbandLF;
        ALL_controls(q,4)=powerbandHF;
        ALL_controls(q,5)=powerbandTotal;
        ALL_controls(q,6)=mean1;
        ALL_controls(q,7)=std1;
        ALL_controls(q,8)=ratioch;
        ALL_controls(q,9)=min1;
        ALL_controls(q,10)=max1;
        ALL_controls(q,11)=numberRR;
        ALL_controls(q,12)=power;
        ALL_controls(q,13)=Parsevalaki;
        save ALL_controls ALL_controls
        close all
%         ALL_local(q,1)=str2num(name);
%         ALL_local(q,2)=powerbandVLF;
%         ALL_local(q,3)=powerbandLF;
%         ALL_local(q,4)=powerbandHF;
%         ALL_local(q,5)=powerbandTotal;
%         ALL_local(q,6)=mean1;
%         ALL_local(q,7)=std1;
%         ALL_local(q,8)=ratioch;
%         
%         save ALL_local ALL_local
%         cd ..
%         ALL_controls2(q,1)=str2num(name);
%         ALL_controls2(q,2)=powerbandVLF;
%         ALL_controls2(q,3)=powerbandLF;
%         ALL_controls2(q,4)=powerbandHF;
%         ALL_controls2(q,5)=powerbandTotal;
%         ALL_controls2(q,6)=mean1;
%         ALL_controls2(q,7)=std1;
%         ALL_controls2(q,8)=ratioch;
%         ALL_controls2(q,9)=parseval;
%         ALL_controls2(q,10)=power;
%         
%         save ALL_controls2 ALL_controls2
        
        close all
end

% function valuta_trained_NN(model)

% Fs=500;
Fs=H_Fs;

n_inputs=model.Layers(1).InputSize;
if(min(n_inputs(1:2)==224)<1),fprintf(' **ERR** NN _inputs: %6.0f %6.0f\n',n_inputs(1:2));end
ECG_mean=[ I_mean II_mean III_mean V1_mean V2_mean V3_mean V4_mean V5_mean V6_mean ];
ECG_TMP=[ I;II;III; V1; V2; V3; V4; V5; V6];
if(size(ECG_TMP,2)<5000), ECG_TMP(1,5000)=0;end

 ECG_PRO=ECG_TMP(1:9,1:5000)'; 
if(K_OPT_NN==2), ECG_PRO=ECG_TMP(1:9,end-5000+1:end)';end
%ECG_PRO=ECG_TMP(1:9,end-5000+1:end)';
%ECG_CINC.dati=ECG_PRO(:);

ECG_CINC.dati=[ECG_PRO(:)' ECG_mean];

    data=ECG_CINC.dati(:);
%     fprintf('ECG: size:%7.0f   sum:%10.0f   abs:%12.0f\n',numel(data),sum(data),sum(abs(data)));
    
    [~,signalLength] = size(data);
    signalLength = numel(data);
%     fprintf('size data:%6.0f%6.0f,  signal length ->  %6.0f \n',size(data),signalLength);
    fb = cwtfilterbank('SignalLength',signalLength, 'SamplingFrequency',Fs,'VoicesPerOctave',12);

    cfs = abs(fb.wt(data));
    im = ind2rgb(im2uint8(rescale(cfs)),jet(128));
    imwrite(imresize(im,[224 224]),'NEW_IMAGE.jpg');
    ECG_image=imread('NEW_IMAGE.jpg');
    
    if(K_OPT_NN==3), ECG_image=im2uint8(imresize(im,[224 224])); end
    
    fprintf('Image- size:');fprintf('%6.0f',size(ECG_image));
    fprintf(' sum:%10.0f%10.0f%10.0f \n',sum(sum(ECG_image)));
    
    [YPred,probs] = classify(model,ECG_image);
    iii=find(probs>0.4);
 %   if(numel(iii)<1),[iitmp,iii]=max(probs);end            % max of NN
    if(numel(iii)<1), iii=find(out_labels>0);end            % out_label of "do"
    out_labels=[];  out_labels(1:9)=0;
    out_labels(iii)=1;
    scores=probs;

    
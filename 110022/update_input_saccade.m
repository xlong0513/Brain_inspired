% update_input_saccade updates the spatial visual input to the network
% according to the saccade.
%
% created: Jakob Heinzle 01/07

for k=1:inputs.ninputs
    if inputs.external{k}.retinotopic
        middle=fov+(nfac+1)/2;
        if middle>21
            helparray=zeros(nfac,1);
        elseif middle>11
            helparray=[inputs.external{k}.inarray(middle-10:nfac) zeros(1,middle-11)];
        elseif middle==11
            helparray=inputs.external{k}.inarray;
        elseif middle<0
            helparray=zeros(nfac,1);
        elseif middle<11
            helparray=[zeros(1,11-middle) inputs.external{k}.inarray(1:middle+10)];

        end
        % auxiliary variables.
        nperpop=pops.population{n_to}.poolsize;
        InpH=zeros(nperpop*nfac,1);
        for oo=1:nfac
            InpH((oo-1)*nperpop+1:oo*nperpop)=helparray(oo);
        end
        inputs.external{k}.ExtInp=inputs.external{k}.MeanInp*InpH;
        inputs.external{k}.NoiseExtInp=sqrt(gmaxE_ext*InpH/2);
        if strcmp(inputs.external{k}.name,'VisualInput') % change input timing for visual input only.
            inputs.external{k}.t_on=t+50;
            inputs.external{k}.t_trans_off=t+90;
            if ((middle<21.5)&(middle>0.5))
                stimatfovea=helparray(middle);
            else
                stimatfovea=0;
            end
        end
    elseif strcmp(inputs.external{k}.name,'FeatureSpace')  % check if the input is the input to white spaces detectors.
        nperpop=pops.population{n_to}.poolsize;
        InpH=zeros(nperpop,1)+(stimatfovea==0);
        inputs.external{k}.ExtInp=inputs.external{k}.MeanInp*InpH;
        inputs.external{k}.NoiseExtInp=sqrt(gmaxE_ext*InpH/2);
    end
end

disp('Visual input changed')
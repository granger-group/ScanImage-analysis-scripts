function analysis = phAnalyze_IntrinsicProperties(wName, varargin)
%phAnalyzeIntrinsicProperties 

    if ~iswave(wName) || ~ischar(wName)
        error('need to provide a wave name in a string for analysis');
    end
    
	% some globals for posthoc analysis, if you want
	analysis=[];

	
%% set up the variables to process data

    hString=getWaveUserDataField(wName, 'headerString');
    if isempty(hString)
        error([wName ' UserData does not contain a header string']);
    end
    
    aiChan=getWaveUserDataField(wName, 'ai');
    if isempty(aiChan)
        disp([wName '.UserData.ai does not exist.  Assuming headstage 1 (i.e. ai0)']);
    end
    aiChanStr=num2str(aiChan);
    % get the output pulse pattern
    pulseString=valueFromHeaderString(['state.phys.internal.pulseString_ao' aiChanStr], hString);
    % get the RC check pulse pattern
    pulseStringRC=valueFromHeaderString('state.phys.internal.pulseString_RCCheck', hString);

    
	% Flag to determine if only the first time occurence of a pulse amplitude
	% should be used or if all repeated occurences should be processed 
	firstOnly=1; 

	% what are the start and end times of the variable current pulse
	% injections
	pulseStart=phUtil_parsePulsePatternString(pulseString, 'delay');
    numPulses=phUtil_parsePulsePatternString(pulseString, 'numPulses');
    pulseWidth=phUtil_parsePulsePatternString(pulseString, 'pulseWidth');
	pulseEnd=pulseStart+pulseWidth;
	pulseAmplitude=phUtil_parsePulsePatternString(pulseString, 'amplitude');
    if (numPulses>1)
        disp('Analyzing first pulse only');
    end
    
	% where the RC check occurs
	% set checkPulseStart=[] if there is no standard RC check pulse in
	% every trace
	checkPulseSize=phUtil_parsePulsePatternString(pulseStringRC, 'amplitude'); % pA we are in current clamp
	checkPulseStart=phUtil_parsePulsePatternString(pulseStringRC, 'delay');
	checkPulseEnd=checkPulseStart+phUtil_parsePulsePatternString(pulseStringRC, 'pulseWidth');

	% values for QC inclusion of individual sweeps
	maxRestSD=5; % max SD of the resting voltage to pass QC
 	minRm=50; % min Rm to pass QC
 	maxRm=1000; % max Rm to pass QC
	
	maxVm=-50; % max Vm for inclusion and to pass QC
	minVm=-90; % min Vm for inclusion and to pass QC
	
	initialDeltaVmSweeps=5; % How many of the first sweeps should we use to 
							% determine the starting Vm
	maxDeltaVm=5; % How far can the Vm move from the above value and still be included
	maxFractionalDeltaRm=0.2; % What fractional change in the Rm will be accept
	
	medianFilterSize=1; % filter the data?
	
%% Set up some parameters to format the output .CSV file and Overwrite the 
% default parameters set above using any parameters values specified 
% by the user in the varargin

	for c=1:2:length(varargin)
		vv=varargin{c+1};
		if isempty(vv)
			vStr='[]';
		elseif isnumeric(vv)
			vStr=num2str(v);
		elseif ischar(vv)
			vStr=['''' vv ''''];
		end
		
		disp(['Override: ' varargin{c} '=' vStr]);
		eval([varargin{c} '=' vStr ';']);
	end
%%	
%% nested function to return a subrange of the data 
    function rang=SR(startS, endS)
        rang=acqData(floor(startS*acqRate):floor(endS*acqRate));
	end
%%
%% nested function to extract a number from string dealing with possible multiple formats
    function ns=extractNum(s)
        if isnumeric(s) 
            ns=s;
        elseif ischar(s)
			si=strfind(s, '_');
            if isempty(si)
                ns=str2double(s);
            else
                ns=str2double(s(si+1:end));
            end
        end
	end
%%
%% nested function to return a value from the headerstring
    function hv=headerValue(sString, conv)
        if nargin<2
            conv=0;
        end
        hv=phUtil_HeaderValue(a.(['AD0_' num2str(acqNum)]).UserData.headerString, ...
            sString, conv);
	end
%%	
%% nested function that returns if a value falls within limits
	function ww=within(x, lo, hi)
		ww=(x>=lo) & (x<=hi);
		return
	end
%%
%% nested function to return only non-nan entries
	function ap=nonNan(a)
		ap=a(~isnan(a));
		return
    end
%%

	%% run through the acquisitions and calculate passive parameters
    % use to do a first pass QC 
	% examine resting potential, RC

    acqData=getWave(wName, 'data');
    acqRate=getWave(wName, 'xscale');
    acqRate=1/acqRate(2); % samples per ms

    acqEndPt=length(acqData)-1;
    acqLen=length(acqData)/acqRate;
    
    % define periods that are "baseline" and anylyze them 
    if isempty(checkPulseStart)
        notPulse=[SR(1, pulseStart-10) SR(pulseEnd+150, acqLen-1)]; 
        rPeak=NaN;
        rEnd=NaN;
        tau=NaN;
    else
        if checkPulseStart>pulseStart % the RC check comes late
            notPulse=[SR(1, pulseStart-10) SR(pulseEnd+150, checkPulseStart-10)]; 
        else
            notPulse=[SR(1, checkPulseStart-1) SR(checkPulseEnd+50, pulseStart-10) SR(pulseEnd+150, acqLen-1)]; 
        end
        [rPeak, rEnd, tau]=phUtil_CurrentClampPulseAnalysis(SR(checkPulseStart, checkPulseEnd)- mean(notPulse), acqRate, checkPulseSize);
    end

    analysis.restMean=mean(notPulse);
    analysis.restMedian=median(notPulse);
    analysis.restSD=std(notPulse);
    analysis.restMin=min(notPulse);
    analysis.restMax=max(notPulse);
    analysis.checkPulseRpeak=rPeak;
    analysis.checkPulseRend=rEnd;
    analysis.checkPulseTau=tau;
    analysis.deltaI=pulseAmplitude;
	analysis.pulseV=mode(round(SR(pulseStart, pulseEnd)));
    if isempty(analysis.pulseV)
        analysis.pulseV=mean(round(SR(pulseStart, pulseEnd)));
    end
    
	deltaI=pulseAmplitude;
    analysis.deltaI=deltaI;

    if deltaI~=0
        analysis.pulseRm=1000*...       % Rm in MOhm
            (analysis.pulseV-analysis.restMean)/deltaI;
        if deltaI<0
            minHyp=min(SR(pulseStart, pulseEnd));
            endHyp=mean(SR(pulseEnd-20,pulseEnd-1));
            analysis.sagV=endHyp-minHyp;
            analysis.reboundV=mean(SR(pulseEnd+20,pulseEnd+70)) ...
                -analysis.restMean;
        end
    end

    if pulseWidth>0
        analysis.pulseAP=phAnalyze_AP(SR(pulseStart, pulseEnd), acqRate);
        if isempty(analysis.pulseAP)
            analysis.nAP=0;
        else
            analysis.nAP=analysis.pulseAP.nAP;
        end

        analysis.postAP=phAnalyze_AP(SR(pulseEnd+1, floor(length(acqData)/acqRate)), acqRate);
        if isempty(analysis.postAP)
            analysis.reboundAP=0;
        else
            analysis.reboundAP=analysis.postAP.nAP;
        end
        analysis.pulseAHP=min(SR(pulseEnd+1, floor(length(acqData)/acqRate)))-...
            analysis.restMean;

    else
        analysis.pulseAP=[];
        analysis.postAP=[];
        analysis.pulseAHP=0;
    end        

    
    setWaveUserDataField(wName, 'analysis', analysis)

end



function [ results ] = phAnalyzeAP(dData, acqRate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here 
%Edited by Adam Granger on 6/6/2023 to detect APs using the 1st derivative
%instead of the zero crossing at 0 mV
    
    g2=gradient(dData);
    g3=gradient(g2);
    
    [gUp, ~]=phUtil_FindXings(g2, 1, 1); % find the up-crossings where the dV/dt goes over 1 - equilavent to 10 mV/ms
    
    if isempty(gUp)
        results=[];
%        disp([   'No action potentials found']);
        return
    end
    
    for i = 1:length(gUp)
        
        [~, gDowns] = phUtil_FindXings(g2(floor(gUp(i)):ceil(gUp(i)+30)),-0.3,1); %find the next crossing past -0.3
        
        
        if isempty(gDowns)
            
            gDown(i) = nan;
        
        else
            temp_gDown = gUp(i)+gDowns(1);
            if i<length(gUp) && gUp(i+1)>temp_gDown %don't record gDowns if the next gUp comes before it
                gDown(i) = temp_gDown; %get the first down crossing past -0.1 for each gUp    
            elseif i == length(gUp)
                gDown(i) = temp_gDown;
            else
                gDown(i) = nan;
            end
            
        end
    end
    
    gUp = gUp(~isnan(gDown)); %only keep gUps corresponding to crossings below -0.1
    gDown = gDown(~isnan(gDown));
    
    %[gUp, gDown]=phUtil_FindXings(dData, 0, 1); % find the zero crossing - this is the old implementation
    gUp=floor(gUp);
    gDown=ceil(gDown);
    
	if nargin<2
		acqRate=10;
	end
	

    
    xNum=min(length(gUp), length(gDown));
%    disp([   num2str(xNum) ' action potentials found']);

    results.nAP=xNum;
    
    results.AP_peak_V=zeros(1, xNum);
    results.AP_peak_time=zeros(1, xNum);
    results.AP_AHP_V=zeros(1,xNum);
    results.AP_thresh_V=zeros(1, xNum);
    results.AP_thresh_time=zeros(1, xNum);   
    results.AP_HW_V=zeros(1, xNum);
    results.AP_0W=zeros(1, xNum);
	results.AP_HW=zeros(1, xNum);
    results.AP_max_dVdT=zeros(1, xNum);
	
%     g2=gradient(dData);
%     g3=gradient(g2);
    
	if length(gDown)>length(gUp)
		if gDown(1)<gUp(1)
			gDown=gDown(2:(length(gUp)+1));
		else
			disp('problem with gDown');
		end
	end
    gDown(end+1)=length(dData);
    gUp(end+1)=length(dData);	
    

    lastMin=1;
    for counter=1:xNum
        [results.AP_peak_V(counter), Imax]=max(dData(gUp(counter):gDown(counter)));
        Imax=Imax+gUp(counter)-1;
        results.AP_peak_time(counter)=Imax;
        [results.AP_AHP_V(counter), Imin]=...
			min(...
			dData(gDown(counter):...
			min(gDown(counter)+10*acqRate,gUp(counter+1)))); % find the min between this AP and the next or 10 ms later, whichever comes first
        Imin=Imin+gDown(counter)-1;

        results.AP_max_dVdT(counter)=max(g2(lastMin:Imax));

		[~, I]=max(g3(lastMin:Imax));
		I=lastMin+I-1-1; % the extra -1 is because of a shift in points taking the derivative
        I=max(I,1);
        results.AP_thresh_V(counter)=dData(I); 
        results.AP_thresh_time(counter)=I;

		
		if (results.AP_peak_V(counter)-results.AP_thresh_V(counter))>20 % if there is not at least 20 mV between threshold and peak, then skip
			HW_V=(results.AP_peak_V(counter)-results.AP_thresh_V(counter))/2+results.AP_thresh_V(counter);
			[ggUp, ggDown]=phUtil_FindXings(dData(lastMin:Imin), HW_V, 1);
			if length(ggDown)>length(ggUp)
				if ggDown(1)<ggUp(1)
					ggDown=ggDown(2:(length(ggUp)+1));
				else
					disp('problem with ggDown');
				end
			end

			results.AP_HW(counter)=ggDown(1)-ggUp(1);
			results.AP_HW_V(counter)=HW_V;

			[ggUp, ggDown]=phUtil_FindXings(dData(lastMin:Imin), 0, 1);
			if ~isempty(ggUp) %Adam added this condition to handle detected APs that do not cross 0 mV
                if length(ggDown)>length(ggUp)
				    if ggDown(1)<ggUp(1)
					    ggDown=ggDown(2:(length(ggUp)+1));
				    else
					    disp('problem with ggDown');
				    end
			    end
			    results.AP_0W(counter)=ggDown(1)-ggUp(1);
            end
		end

		lastMin=Imin;
	end
    results.AP_thresh_time=results.AP_thresh_time/acqRate;
    results.AP_HW=results.AP_HW/acqRate;
    results.AP_0W=results.AP_0W/acqRate;
	results.AP_peak_time=results.AP_peak_time/acqRate;
    results.AP_max_dVdT=results.AP_max_dVdT*acqRate;
end


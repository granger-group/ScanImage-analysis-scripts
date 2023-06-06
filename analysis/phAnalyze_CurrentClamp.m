function analysis=phAnalyze_CurrentClamp(wName)
    global state phAnalysis
    analysis=phAnalyze_IntrinsicProperties(wName);
    
    fNames=fieldnames(analysis);
    
    eString=['e' num2str(state.epoch)];
    if ~isfield(phAnalysis, eString)
        phAnalysis.(eString)=[];
    end
    if isempty(phAnalysis.(eString))
        phAnalysis.(eString).acqs={wName};
        newCounter=1;
    else
        phAnalysis.(eString).acqs{end+1}=wName;
    end
    newCounter=length(phAnalysis.(eString).acqs);
    
    for fCounter=1:length(fNames)
        fName=fNames{fCounter};
        if isnumeric(analysis.(fName)) && ~isempty(analysis.(fName))
            phAnalysis.(eString).(fName)(newCounter)=analysis.(fName);
        else
            phAnalysis.(eString).(fName)(newCounter)={analysis.(fName)};
        end            
    end
    
    if state.files.autoSave
       save([state.files.baseName '_phAnalysis.mat'], 'phAnalysis');
    end
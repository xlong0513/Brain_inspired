%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of CONFIGR,
% as described in:
% Carpenter, G. A., Gaddam, C. S., & Mingolla, E. (2007). 
% CONFIGR: A vision-based system for long-range figure completion. Neural Networks, xx(x) xxx-xxx. 
%  Technical Report CAS/CNS TR-2007-016, Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(iSPPura{1})
    DeepSearch=0;
    
    searchlight=0;%checks if any corners have been found
    step_search=0; %increments area of search
    span_max=1; %checks if bounds of image have been reached
    
    
    
    
    N=FillInLimit+1; 
    iMIN=max(iSPPura{1}(i)-N,1);
    iMAX=max(iSPPura{1}(i),1);
    jMIN=min(jSPPura{1}(i),ImSize(2));
    jMAX=min(jSPPura{1}(i)+N,ImSize(2));
    
    
    if sum(sum(sign02(IsparsePura{3}(iMIN:iMAX-1,jMIN+1:jMAX))))~=0 %Offset in iMAX and jMIN is needed to not count 
        %                                                           the position of the initial spanning corner
        [iPart,jPart,sPart]=find(IsparsePura{3}(iMIN:iMAX-1,jMIN+1:jMAX));
        KernSizeTempI=iMAX-iMIN+1;
        KernSizeTempJ=jMAX-jMIN+1;
        
        [YshortestDist,indexK]=sort(abs(KernSizeTempI-iPart)+jPart); 
        
        TestPatchI=iPart(indexK(1))+iMIN-1;
        TestPatchJ=jPart(indexK(1))+jMIN;
        checkval=0;
        stepval=2;
        while checkval==0;
            if abs(TestPatchI-iMAX)<2 | abs(TestPatchJ-jMIN)<2%;
                try
                    TestPatchI=iPart(indexK(stepval))+iMIN-1;
                    TestPatchJ=jPart(indexK(stepval))+jMIN;
                    stepval=stepval+1;
                catch
                    checkval=1;
                end
                
            else
                checkval=1;
            end
        end
        
        
        SearchDepthDistance=YshortestDist(stepval-1);
        ChkNxtCorner=stepval-1;
        
        
        Bounded_tst_=CheckBound(TestPatchI,iMAX,jMIN,TestPatchJ,SubPixRes,IcrawlT);
        
        while (sum(sum(IpixelsTemp(TestPatchI+1:iMAX-1,jMIN+1:TestPatchJ-1)))~=0 | (Bounded_tst_==0)) && (ChkNxtCorner<length(indexK)) ...
                && YshortestDist(ChkNxtCorner)<FillInLimit+1
            ChkNxtCorner=ChkNxtCorner+1;
            TestPatchI=iPart(indexK(ChkNxtCorner))+iMIN-1;
            TestPatchJ=jPart(indexK(ChkNxtCorner))+jMIN;
            Bounded_tst_=CheckBound(TestPatchI,iMAX,jMIN,TestPatchJ,SubPixRes,IcrawlT);
        end
        
        
        if sum(sum(IpixelsTemp(TestPatchI+1:iMAX-1,jMIN+1:TestPatchJ-1)))==0 && (Bounded_tst_==1); 
            if sum(sum(IsparsePura{1}(TestPatchI+1:iMAX,jMIN:TestPatchJ-1)))>1
                ItempCorners(:,:)=0;
                ItempCorners(TestPatchI+1:iMAX,jMIN:TestPatchJ-1)=1;
                
                [iSpanners,jSpanners,sSpanners]=find(IsparsePura{1}.*ItempCorners);
                
                [YshortestDist,indexSpanners]=sort(abs(iSpanners-TestPatchI)+abs(-jSpanners+TestPatchJ));
                
                iMAX=iSpanners(indexSpanners(1));
                jMIN=jSpanners(indexSpanners(1));
            end
            
            IallRects(rectsAll,:)=[TestPatchI iMAX jMIN TestPatchJ 1];
            rectsAll=rectsAll+1;
            
        end
        
        
    end
    
        
end

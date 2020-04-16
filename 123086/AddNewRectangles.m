for i=1:length(iSPNewWalls{1})
    DeepSearch=0;
    
    searchlight=0;%checks if any corners have been found
    step_search=0; %increments area of search
    span_max=1; %checks if bounds of image have been reached
    
    N=FillInLimit+1; 
    iMIN=max(iSPNewWalls{1}(i)-N,1);
    iMAX=max(iSPNewWalls{1}(i),1);
    jMIN=min(jSPNewWalls{1}(i),1);
    jMAX=min(jSPNewWalls{1}(i)+N,ImSize(2));
    
    
    if sum(sum(sign02(IsparsePura{3}(iMIN:iMAX-1,jMIN+1:jMAX)+IsparseNewWalls{3}(iMIN:iMAX-1,jMIN+1:jMAX))))~=0 
        [iPart,jPart,sPart]=find(IsparsePura{3}(iMIN:iMAX-1,jMIN+1:jMAX)+IsparseNewWalls{3}(iMIN:iMAX-1,jMIN+1:jMAX));
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
            if sum(sum(IsparsePura{1}(TestPatchI+1:iMAX,jMIN:TestPatchJ-1)+IsparseNewWalls{1}(TestPatchI+1:iMAX,jMIN:TestPatchJ-1)))<=1
                
                
                
                IallRectsNew(rectsAllNew,:)=[TestPatchI iMAX jMIN TestPatchJ 1];
                WallPosts(rectsAllNew,:)=sWALL;
                rectsAllNew=rectsAllNew+1;
                
            end
        end
        
        
    end
    
    
end


for i=1:length(iSPNewWalls{2})
    
    
    searchlight=0;%checks if any corners have been found
    step_search=0; %increments area of search
    span_max=1; %checks if bounds of image have been reached
    
    
    
    N=FillInLimit+1; 
    iMIN=max(iSPNewWalls{2}(i)-N,1);
    iMAX=max(iSPNewWalls{2}(i),1);
    jMAX=max(jSPNewWalls{2}(i),1);
    jMIN=max(jSPNewWalls{2}(i)-N,1);
    sWALL=[0 sSPAll{2}(i) 0 0];
    
    
    
    
    if sum(sum(sign02(IsparsePura{4}(iMIN:iMAX-1,jMIN:jMAX-1)+IsparseNewWalls{4}(iMIN:iMAX-1,jMIN:jMAX-1))))~=0 
        [iPart,jPart,sPart]=find(IsparsePura{4}(iMIN:iMAX-1,jMIN:jMAX-1)+IsparseNewWalls{4}(iMIN:iMAX-1,jMIN:jMAX-1));
        KernSizeTempI=iMAX-iMIN+1;
        KernSizeTempJ=jMAX-jMIN+1;
        [YshortestDist,indexK]=sort(abs(KernSizeTempI-iPart)+abs(KernSizeTempJ-jPart));
        
        TestPatchI=iPart(indexK(1))+iMIN-1;
        TestPatchJ=jPart(indexK(1))+jMIN-1;
        
        checkval=0;
        stepval=2;
        while checkval==0;
            if abs(TestPatchI-iMAX)<2 | abs(TestPatchJ-jMAX)<2 ;
                try
                    TestPatchI=iPart(indexK(stepval))+iMIN-1;
                    TestPatchJ=jPart(indexK(stepval))+jMIN-1;
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
        
        
        Bounded_tst_=CheckBound(TestPatchI,iMAX,TestPatchJ,jMAX,SubPixRes,IcrawlT);
        
        
        while (sum(sum(IpixelsTemp(TestPatchI+1:iMAX-1,TestPatchJ+1:jMAX-1)))~=0 | (Bounded_tst_==0)) && (ChkNxtCorner<length(indexK)) ...
                && YshortestDist(ChkNxtCorner)<FillInLimit+1
            ChkNxtCorner=ChkNxtCorner+1;
            TestPatchI=iPart(indexK(ChkNxtCorner))+iMIN-1;
            TestPatchJ=jPart(indexK(ChkNxtCorner))+jMIN-1;
            Bounded_tst_=CheckBound(TestPatchI,iMAX,TestPatchJ,jMAX,SubPixRes,IcrawlT);
        end 
        
        
        
        if sum(sum(IpixelsTemp(TestPatchI+1:iMAX-1,TestPatchJ+1:jMAX-1)))==0 && (Bounded_tst_==1); %&& (Inonsides(TestPatchI,TestPatchJ)+Inonsides(iMAX+1,jMAX+1)>0);%...
            
            if sum(sum(IsparsePura{2}(TestPatchI+1:iMAX,TestPatchJ+1:jMAX)+IsparseNewWalls{2}(TestPatchI+1:iMAX,TestPatchJ+1:jMAX)))<=1
                
                
                
                IallRectsNew(rectsAllNew,:)=[TestPatchI iMAX TestPatchJ jMAX 2];
                
                rectsAllNew=rectsAllNew+1;
                
            end
            
        end
        
    end
    
end



%=============================================================



for i=1:length(iSPNewWalls{3})
    DeepSearch=0;
    
    searchlight=0;%checks if any corners have been found
    step_search=0; %increments area of search
    span_max=1; %checks if bounds of image have been reached
    
    N=FillInLimit+1; 
    iMIN=max(iSPNewWalls{3}(i),1);
    iMAX=min(iSPNewWalls{3}(i)+N,ImSize(1));
    jMIN=max(jSPNewWalls{3}(i)-N,1);
    jMAX=min(jSPNewWalls{3}(i),ImSize(2));
    
    
    if sum(sum(sign02(IsparsePura{1}(iMIN+1:iMAX,jMIN:jMAX-1)+IsparseNewWalls{1}(iMIN+1:iMAX,jMIN:jMAX-1))))~=0 
        [iPart,jPart,sPart]=find(IsparsePura{1}(iMIN+1:iMAX,jMIN:jMAX-1)+IsparseNewWalls{1}(iMIN+1:iMAX,jMIN:jMAX-1));
        KernSizeTempI=iMAX-iMIN+1;
        KernSizeTempJ=jMAX-jMIN+1;
        
        [YshortestDist,indexK]=sort(abs(iPart)+abs(KernSizeTempJ-jPart)); 
        
        TestPatchI=iPart(indexK(1))+iMIN;
        TestPatchJ=jPart(indexK(1))+jMIN-1;
        checkval=0;
        stepval=2;
        while checkval==0;
            if abs(TestPatchI-iMIN)<2 | abs(TestPatchJ-jMAX)<2%;
                try
                    TestPatchI=iPart(indexK(stepval))+iMIN;
                    TestPatchJ=jPart(indexK(stepval))+jMIN-1;
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
        
        
        Bounded_tst_=CheckBound(iMIN,TestPatchI,TestPatchJ,jMAX,SubPixRes,IcrawlT);
        
        while (sum(sum(IpixelsTemp(iMIN+1:TestPatchI-1,TestPatchJ+1:jMAX-1)))~=0 | (Bounded_tst_==0)) && (ChkNxtCorner<length(indexK)) ...
                && YshortestDist(ChkNxtCorner)<FillInLimit+1
            ChkNxtCorner=ChkNxtCorner+1;
            TestPatchI=iPart(indexK(ChkNxtCorner))+iMIN;
            TestPatchJ=jPart(indexK(ChkNxtCorner))+jMIN-1;
            Bounded_tst_=CheckBound(iMIN,TestPatchI,TestPatchJ,jMAX,SubPixRes,IcrawlT);
        end
        
        
        
        
        
        if sum(sum(IpixelsTemp(iMIN+1:TestPatchI-1,TestPatchJ+1:jMAX-1)))==0 && (Bounded_tst_==1); 
            if sum(sum(IsparsePura{3}(iMIN:TestPatchI-1,TestPatchJ+1:jMAX)+IsparseNewWalls{3}(iMIN:TestPatchI-1,TestPatchJ+1:jMAX)))<=1
                
                
                
                IallRectsNew(rectsAllNew,:)=[iMIN TestPatchI TestPatchJ jMAX 1];
                
                rectsAllNew=rectsAllNew+1;
                
            end
        end
        
        
    end
    
    
end


for i=1:length(iSPNewWalls{4})
    
    
    searchlight=0;%checks if any corners have been found
    step_search=0; %increments area of search
    span_max=1; %checks if bounds of image have been reached
    
    
    
    N=FillInLimit+1; 
    iMIN=max(iSPNewWalls{4}(i),1);
    iMAX=min(iSPNewWalls{4}(i)+N,ImSize(1));
    jMIN=max(jSPNewWalls{4}(i),1);
    jMAX=min(jSPNewWalls{4}(i)+N,ImSize(2));
    
    
    
    
    
    if sum(sum(sign02(IsparsePura{2}(iMIN+1:iMAX,jMIN+1:jMAX)+IsparseNewWalls{2}(iMIN+1:iMAX,jMIN+1:jMAX))))~=0 
        [iPart,jPart,sPart]=find(IsparsePura{2}(iMIN+1:iMAX,jMIN+1:jMAX)+IsparseNewWalls{2}(iMIN+1:iMAX,jMIN+1:jMAX));
        KernSizeTempI=iMAX-iMIN+1;
        KernSizeTempJ=jMAX-jMIN+1;
        [YshortestDist,indexK]=sort(abs(iPart)+abs(jPart)); %Shortest Distance from (iMIN,jMIN)
        
        TestPatchI=iPart(indexK(1))+iMIN;
        TestPatchJ=jPart(indexK(1))+jMIN;
        
        checkval=0;
        stepval=2;
        while checkval==0;
            if abs(TestPatchI-iMIN)<2 | abs(TestPatchJ-jMIN)<2 ;
                try
                    TestPatchI=iPart(indexK(stepval))+iMIN;
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
        
        
        Bounded_tst_=CheckBound(iMIN,TestPatchI,jMIN,TestPatchJ,SubPixRes,IcrawlT);
        
        
        while (sum(sum(IpixelsTemp(iMIN+1:TestPatchI-1,jMIN+1:TestPatchJ-1)))~=0 | (Bounded_tst_==0)) && (ChkNxtCorner<length(indexK)) ...
                && YshortestDist(ChkNxtCorner)<FillInLimit+1
            ChkNxtCorner=ChkNxtCorner+1;
            TestPatchI=iPart(indexK(ChkNxtCorner))+iMIN;
            TestPatchJ=jPart(indexK(ChkNxtCorner))+jMIN;
            Bounded_tst_=CheckBound(iMIN,TestPatchI,jMIN,TestPatchJ,SubPixRes,IcrawlT);
        end 
        
        
        
        if sum(sum(IpixelsTemp(iMIN+1:TestPatchI-1,jMIN+1:TestPatchJ-1)))==0 && (Bounded_tst_==1); %&& (Inonsides(TestPatchI,TestPatchJ)+Inonsides(iMAX+1,jMAX+1)>0);%...
            
            if sum(sum(IsparsePura{4}(iMIN:TestPatchI-1,jMIN:TestPatchJ-1)+IsparseNewWalls{4}(iMIN:TestPatchI-1,jMIN:TestPatchJ-1)))<=1
                
                
                IallRectsNew(rectsAllNew,:)=[iMIN TestPatchI jMIN TestPatchJ 2];
                
                rectsAllNew=rectsAllNew+1;
                
            end
            
        end
        
    end
    
end

if exist('IallRectsNew')
IrectNewSize=IallRectsNew(:,2)+IallRectsNew(:,4)-IallRectsNew(:,1)-IallRectsNew(:,3);
[dumArray sortRect]=sort(IrectNewSize);
IallRectsNew=IallRectsNew(sortRect,:);
end

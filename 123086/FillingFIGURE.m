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
% as do the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('IallRects')
    IrectSize=IallRects(:,2)+IallRects(:,4)-IallRects(:,1)-IallRects(:,3);
    
    
    [dumArray sortRect]=sort(IrectSize);
    IallRects=IallRects(sortRect,:);
    
    NumTotalRects=size(IallRects,1);
    stpRects=1;
    
    while stpRects<=NumTotalRects
        tmpVect=IallRects(stpRects,:);
        
        sizRect=tmpVect(2)+tmpVect(4)-tmpVect(1)-tmpVect(3);
        EmptyRect = (isempty(tmpVect(1)+1:tmpVect(2)-1) || isempty(tmpVect(3)+1:tmpVect(4)-1));
        
      
        if (~EmptyRect) && sum(max(max(IfillIter(tmpVect(1)+1:tmpVect(2)-1,tmpVect(3)+1:tmpVect(4)-1))) == [1/dumArray(stpRects) 0])  && sum(sum(IpixelsFillA(tmpVect(1)+1:tmpVect(2)-1,tmpVect(3)+1:tmpVect(4)-1)))==0;
            tmp_dummy=0;
            
            [NborGround , tmp_dummy]=CheckOutBoundGround(tmpVect(1),tmpVect(2),tmpVect(3),tmpVect(4),SubPixRes,IpixelsFillA,FillNoFillVal,ImSize,tmp_dummy);
            
            NumFilled_Rule3(round((kk-1)/5)+1)=NumFilled_Rule3(round((kk-1)/5)+1)+tmp_dummy;
            
            filled_in_rule3_rects_store(rule3_step,:)=[round((kk-1)/5)+1 tmpVect(1) tmpVect(2) tmpVect(3) tmpVect(4)];
            rule3_step=rule3_step+1;
            
            gCorner=...
                (sum((IpixelsFillA(tmpVect(1)+1,tmpVect(3)-1)==[1  FIGR_VAL*(.15+((1-.15)*Iemergent_corners{4}(tmpVect(1),tmpVect(3))))])...
                +(IpixelsFillA(tmpVect(1)-1,tmpVect(3)+1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{4}(tmpVect(1),tmpVect(3))))]))>1)+...
                (sum((IpixelsFillA(tmpVect(2)+1,tmpVect(4)-1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{2}(tmpVect(2),tmpVect(4))))])...
                +(IpixelsFillA(tmpVect(2)-1,tmpVect(4)+1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{2}(tmpVect(2),tmpVect(4))))]))>1)+...
                (sum((IpixelsFillA(tmpVect(1)-1,tmpVect(4)-1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{3}(tmpVect(1),tmpVect(4))))])...
                +(IpixelsFillA(tmpVect(1)+1,tmpVect(4)+1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{3}(tmpVect(1),tmpVect(4))))]))>1)+...
                (sum((IpixelsFillA(tmpVect(2)+1,tmpVect(3)+1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{1}(tmpVect(2),tmpVect(3))))])...
                +(IpixelsFillA(tmpVect(2)-1,tmpVect(3)-1)==[1 FIGR_VAL*(.15+((1-.15)*Iemergent_corners{1}(tmpVect(2),tmpVect(3))))]))>1)+...
                NborGround;
            
            if  (gCorner==0) ;
                
                PosBreak=max(find(IallRects(:,5)==10));
                if isempty(PosBreak)
                    PosBreak=0;
                end
                
                if (stpRects>PosBreak) 
                    
                    IpixelsTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=max(.5*ones(size(I(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)))),I(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)));
                    
                    Iinterpol(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=InterpolateLines(Iinterpol(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)),tmpVect(5),SubPixRes);
                    NumFilled_FIGURE(round((kk-1)/5)+1)=NumFilled_FIGURE(round((kk-1)/5)+1)+1;
                    
                    
                    
                    filled_in_figr_rects_store(figr_step,:)=[round((kk-1)/5)+1 tmpVect(1) tmpVect(2) tmpVect(3) tmpVect(4)];
                    figr_step=figr_step+1;
                    
                    
                    %This is to create diagonals at a different
                    %resolution later on
                    
                    StoreDiagVals(stp_diag_,:)=tmpVect;
                    stp_diag_=stp_diag_+1;
                    
                    IpixelsFillTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=(IpixelsTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))==.5)*.5+(IpixelsTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))==FillNoFillVal)*FillNoFillVal;
                else
                    
                    %Rectangle designated as a FIGURE filling-in candidate was
                    %newly created. Thus, cannot be filled-in before other
                    %newly created rectangles are dealt with.
                    IallRects_top=[IallRects(1:stpRects-1,:) ; IallRects(stpRects+1:PosBreak-1,:)];
                    IallRects_bot=[IallRects(stpRects,:) ; IallRects(PosBreak+1:size(IallRects,1),:) ];
                    Irects_botSize=IallRects_bot(:,2)+IallRects_bot(:,4)-IallRects_bot(:,1)-IallRects_bot(:,3);
                    [dumArray sortRect]=sort(Irects_botSize);
                    IallRects_bot=IallRects_bot(sortRect,:);
                    IallRects=[IallRects_top; [13 13 13 13 10] ; IallRects_bot];
                    %The rectangle at stpRects is a new one because of the
                    %reshuffle. Cycle back to the previous one to check if
                    %the new rectangle is of a different size.
                    stpRects=stpRects-1;
                    resortArray=1;
                end
                
                
                
            else
               
                IpixelsTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=max(FillNoFillVal*ones(size(I(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)))),I(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)));
                IpixelsFillTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=(IpixelsTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))==.5)*.5+(IpixelsTempIter12(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))==FillNoFillVal)*FillNoFillVal;
                
                IpixelsFilledInTemp(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=1;
                
                NumFilled_GROUND(round((kk-1)/5)+1)=NumFilled_GROUND(round((kk-1)/5)+1)+1;
                
                filled_in_grnd_rects_store(grnd_step,:)=[round((kk-1)/5)+1 tmpVect(1) tmpVect(2) tmpVect(3) tmpVect(4)];
                grnd_step=grnd_step+1;
                
          
                
                
                if kk==2*5+1
                    IrectsShow(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=1;
                end
            end
            
            if (stpRects>0) && ((dumArray(stpRects)<dumArray(min(stpRects+1,size(IallRects,1)))) || (stpRects==size(IallRects,1)))
          
                
                IpixelsFilledIn=sign02(IpixelsFilledInTemp-sign02(IpixelsFillA));
                
                
                IcrawlTWall_temp{5}=sparse(IcrawlTWall{5}).*conv2(1-sign02(IpixelsFilledIn),[0 0 0; 0 0 0; 1 0 0],'same');
                IcrawlTWall_temp{6}=sparse(IcrawlTWall{6}).*conv2(1-sign02(IpixelsFilledIn),[0 0 0; 0 0 0; 0 0 1],'same');
                IcrawlTWall_temp{7}=sparse(IcrawlTWall{7}).*conv2(1-sign02(IpixelsFilledIn),[0 0 1; 0 0 0; 0 0 0],'same');
                IcrawlTWall_temp{8}=sparse(IcrawlTWall{8}).*conv2(1-sign02(IpixelsFilledIn),[1 0 0; 0 0 0; 0 0 0],'same');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %This is updated after all the walls eliminated at
                %this point have been accounted for
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                IpixelsFillA=IpixelsTempIter12;
                
                
                if sum(sum(IcrawlTWall_temp{5}+IcrawlTWall_temp{6}+IcrawlTWall_temp{7}+IcrawlTWall_temp{8}))...
                        <sum(sum(IcrawlTWall{5}+IcrawlTWall{6}+IcrawlTWall{7}+IcrawlTWall{8}))
                    
             
                    
                    IcrawlTNewWalls{5}=...
                        sparse(IcrawlTWall_temp{5}.*(1-sign02(IcrawlTWall_temp{6}+IcrawlTWall_temp{8})));
                    IcrawlTNewWalls{6}=...
                        sparse(IcrawlTWall_temp{6}.*(1-sign02(IcrawlTWall_temp{5}+IcrawlTWall_temp{7})));
                    IcrawlTNewWalls{7}=...
                        sparse(IcrawlTWall_temp{7}.*(1-sign02(IcrawlTWall_temp{6}+IcrawlTWall_temp{8})));
                    IcrawlTNewWalls{8}=...
                        sparse(IcrawlTWall_temp{8}.*(1-sign02(IcrawlTWall_temp{7}+IcrawlTWall_temp{5})));
                    
                    for i=1:4
                        IsparseNewWalls{i}=sparse(IcrawlTNewWalls{4+i});
                        [iSPNewWalls{i},jSPNewWalls{i},sSPNewWalls{i}]=find(IsparseNewWalls{i});
                    end
                    
                 
                    
                    AddNewRectangles
                    
                    
                    IcrawlTPura{5}=sparse(sign02(IcrawlTNewWalls{5}+IcrawlTPura{5}));
                    IcrawlTPura{6}=sparse(sign02(IcrawlTNewWalls{6}+IcrawlTPura{6}));
                    IcrawlTPura{7}=sparse(sign02(IcrawlTNewWalls{7}+IcrawlTPura{7}));
                    IcrawlTPura{8}=sparse(sign02(IcrawlTNewWalls{8}+IcrawlTPura{8}));
                    
                    
                   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %There are new corners, the process needs to be repeated all over again.
                    %Some efficiency may be obtained by setting
                    %stRects to first rectangle created by a new
                    %empty corner
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if exist('IallRectsNew')
                        IallRects = [IallRectsNew ; [13 13 13 13 10] ;IallRects];
                        rectSmallest=min(IallRectsNew(:,2)+IallRectsNew(:,4)-IallRectsNew(:,1)-IallRectsNew(:,3));
                        clear IallRectsNew
                        rectsAllNew=1;
                        
                        IrectSize=IallRects(:,2)+IallRects(:,4)-IallRects(:,1)-IallRects(:,3);
                   
                        
                        
                        [dumArray sortRect]=sort(IrectSize);
                        IallRects_cmp=IallRects(sortRect,:);
                        
                        NumTotalRects=size(IallRects,1);
                       
                        %Found new rectangles; start iterating from the
                        %start
                     
                        stpRects=0;
                    end
                    
                end
                
            end
            
            
            
            if (stpRects>0) && (resortArray==0)
                IfillIter(tmpVect(1)+1:tmpVect(2)-1,tmpVect(3)+1:tmpVect(4)-1)=1/dumArray(stpRects);
            end
            resortArray=0;
            %This prevents noting down on filling in when we are merely
            %resorting the array and not filling in
            
            
        end % Loop checking if the rectangle is empty, and defined by empty corners, not wall corners
        stpRects=stpRects+1;
    end  % End of While Loop 
    
    
    
    IfillIter=zeros(size(I));
    
    
    Ipixels=I+(IpixelsFillTempIter24==IpixelsFillTempIter12).*IpixelsFillTempIter24...
        +(IpixelsFillTempIter24~=IpixelsFillTempIter12).*(IpixelsFillTempIter24+IpixelsFillTempIter12);
    
    IpixelsTemp=Ipixels;
    IpixelsFill=Ipixels-I;
    IpixelsFillA=Ipixels;
end % End of Exist AllRects If Loop
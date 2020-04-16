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


        if (~EmptyRect) && sum(max(max(IfillIter(tmpVect(1)+1:tmpVect(2)-1,tmpVect(3)+1:tmpVect(4)-1))) == [1/dumArray(stpRects) 0]) &&  sum(sum(IpixelsFillA(tmpVect(1)+1:tmpVect(2)-1,tmpVect(3)+1:tmpVect(4)-1)))==0;
            tmp_dummy=0;

            [NborGround , tmp_dummy]=CheckOutBoundGround(tmpVect(1),tmpVect(2),tmpVect(3),tmpVect(4),SubPixRes,IpixelsFillA,FillNoFillVal,ImSize,tmp_dummy);

            NumFilled_Rule3(round((kk-1)/5)+1)=NumFilled_Rule3(round((kk-1)/5)+1)+tmp_dummy;

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
                
                %This rectangle is a candidate for filling-in as FIGURE.


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

            if (dumArray(stpRects)<dumArray(min(stpRects+1,size(IallRects,1)))) || (stpRects==size(IallRects,1))

                IpixelsFilledIn=sign02(IpixelsFilledInTemp-sign02(IpixelsFillA));


                IcrawlTWall_temp{5}=sparse(IcrawlTWall{5}).*conv2(1-sign02(IpixelsFilledIn),[0 0 0; 0 0 0; 1 0 0],'same');
                IcrawlTWall_temp{6}=sparse(IcrawlTWall{6}).*conv2(1-sign02(IpixelsFilledIn),[0 0 0; 0 0 0; 0 0 1],'same');
                IcrawlTWall_temp{7}=sparse(IcrawlTWall{7}).*conv2(1-sign02(IpixelsFilledIn),[0 0 1; 0 0 0; 0 0 0],'same');
                IcrawlTWall_temp{8}=sparse(IcrawlTWall{8}).*conv2(1-sign02(IpixelsFilledIn),[1 0 0; 0 0 0; 0 0 0],'same');



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
                        IallRects = [IallRectsNew ; IallRects];
                        rectSmallest=min(IallRectsNew(:,2)+IallRectsNew(:,4)-IallRectsNew(:,1)-IallRectsNew(:,3));
                        clear IallRectsNew
                        rectsAllNew=1;

                        IrectSize=IallRects(:,2)+IallRects(:,4)-IallRects(:,1)-IallRects(:,3);




                        [dumArray sortRect]=sort(IrectSize);
                        IallRects=IallRects(sortRect,:);

                        NumTotalRects=size(IallRects,1);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %This step sets the 'pointer' to the
                        %position of the rectangle with same size a
                        %as the smallest new rectangle
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        stpRects=max(min(find(dumArray==rectSmallest))-1,1);
                    end

                end


            end



            if gCorner==1
                IfillIter(tmpVect(1)+1:tmpVect(2)-1,tmpVect(3)+1:tmpVect(4)-1)=1/dumArray(stpRects);
            end
        end % Loop checking if the rectangle is empty, and defined by empty corners, not wall corners
        stpRects=stpRects+1;
    end  % End of While Loop



    IfillIter=zeros(size(I));


    Ipixels=I+(IpixelsFillTempIter24==IpixelsFillTempIter12).*IpixelsFillTempIter24...
        +(IpixelsFillTempIter24~=IpixelsFillTempIter12).*(IpixelsFillTempIter24+IpixelsFillTempIter12);

    IpixelsTemp=Ipixels;
    IpixelsFill=Ipixels-I;
    IpixelsFillA=Ipixels;
end
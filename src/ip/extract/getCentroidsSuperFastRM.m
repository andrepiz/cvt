function [cameraPoints]= getCentroidsSuperFastRM(I_rm,threshold,minROI,maxROI,minDC,height,width)

maxCentNum=400;
xCent_all=zeros(1,maxCentNum);
yCent_all=zeros(1,maxCentNum);
k=0;
idx = int32(2*width+4);
step=max(int32(1),idivide(minROI,int32(2)));
while idx < (height-1)*width-step
    idx=idx+step;
    I0=I_rm(idx);
    if I0<threshold
        continue
    end

    I_rm(idx)=uint8(0);

    y0=idivide(idx,width)+1;
    x0=mod(idx-1,width)+1;

    sidesClear=false;
    nUp=int32(1);
    nLeft=int32(1);
    nRight=int32(1);
    nDown=int32(1);
    DC_max=I0;
    DC2_tot=double(I0)^2;
    xCent_tmp=double(x0)*DC2_tot;
    yCent_tmp=double(y0)*DC2_tot;

    upLimit = y0-nUp < 0;
    if upLimit
        nUp = nUp-1;
    end

    downLimit = y0+nDown > height;
    if downLimit
        nDown = nDown-1;
    end

    leftLimit = x0-nLeft < 0;
    if leftLimit
        nLeft = nLeft-1;
    end

    rightLimit = x0+nRight > width;
    if rightLimit
        nRight = nRight-1;
    end


    while ~sidesClear && max(nLeft+nRight,nDown+nUp)+1<maxROI

        idx_TL=idx-nLeft-width*nUp;
        idx_TR=idx+nRight-width*nUp;
        idx_BL=idx-nLeft+width*nDown;
        idx_BR=idx+nRight+width*nDown;
        nUp_old=nUp;
        nDown_old=nDown;


        DC_TL=I_rm(idx_TL);
        DC_TR=I_rm(idx_TR);
        DC_BL=I_rm(idx_BL);
        DC_BR=I_rm(idx_BR);

        upClear = (DC_TL<threshold)&&(DC_TR<threshold);
        idxUp = idx_TL;
        for kUp = int32(1): nLeft+nRight-1
            %         vec_Up=zeros(1,nLeft+nRight-1,'int32');
            %         vec_Up(:)=idx_TL+1:idx_TR-1;
            idxUp = idxUp+1;
            DC_Up=I_rm(idxUp);
            I_rm(idxUp)=uint8(0);

            if ~(DC_Up <threshold)
                upClear = false;
                DC_max=max(DC_max,DC_Up);
                DC2_Up= double(DC_Up)^2;
                DC2_tot=DC2_tot+DC2_Up;
                xCent_tmp=xCent_tmp+double(x0-nLeft+kUp)*DC2_Up;
                yCent_tmp=yCent_tmp+double(y0-nUp)*DC2_Up;
            end
        end

        if ~upClear && y0-nUp >2
            nUp=nUp+1;
        end

        %%
        downClear = (DC_BL<threshold)&&(DC_BR<threshold);
        idxDown = idx_BL;
        for kDown = int32(1): nLeft+nRight-1
            idxDown = idxDown+1;
            DC_Down=I_rm(idxDown);
            I_rm(idxDown)=uint8(0);

            if ~(DC_Down <threshold)
                downClear = false;
                DC_max=max(DC_max,DC_Down);
                DC2_Down= double(DC_Down)^2;
                DC2_tot=DC2_tot+DC2_Down;
                xCent_tmp=xCent_tmp+double(x0-nLeft+kDown)*DC2_Down;
                yCent_tmp=yCent_tmp+double(y0+nDown)*DC2_Down;
            end
        end

        if ~downClear && y0+nDown < height-1
            nDown=nDown+1;
        end

        %%
        leftClear = true;
        if leftLimit
            I_rm(idx_TL) = uint8(0);
            I_rm(idx_BL) = uint8(0);
            DC2_BL= double(DC_BL)^2;
            DC2_TL= double(DC_TL)^2;
            DC2_tot = DC2_tot+ DC2_BL+DC2_TL;
            xCent_tmp=xCent_tmp+double((x0-nLeft))*DC2_TL+double((x0-nLeft))*DC2_BL;
            yCent_tmp=yCent_tmp+double((y0-nUp_old))*DC2_TL+double((y0+nDown_old))*DC2_BL;
        else

            for kLeft = int32(1): nDown_old+1+nUp_old
                idxLeft = idx_TL + (kLeft-1) * width;
                DC_Left=I_rm(idxLeft);
                I_rm(idxLeft)=uint8(0);

                if ~(DC_Left <threshold)
                    leftClear = false;
                    DC_max=max(DC_max,DC_Left);
                    DC2_Left= double(DC_Left)^2;
                    DC2_tot=DC2_tot+DC2_Left;
                    xCent_tmp=xCent_tmp+double(x0-nLeft)*DC2_Left;
                    yCent_tmp=yCent_tmp+double(y0-nUp_old-1+kLeft)*DC2_Left;
                end
            end

            if ~leftClear
                if  x0-nLeft > 2
                    nLeft=nLeft+1;
                else
                    leftLimit = true;
                end
            end
        end
        %%

        rightClear = true;
        if rightLimit
            I_rm(idx_TR) = uint8(0);
            I_rm(idx_BR) = uint8(0);
            DC2_BR= double(DC_BR)^2;
            DC2_TR= double(DC_TR)^2;
            DC2_tot = DC2_tot+ DC2_BR+DC2_TR;
            xCent_tmp=xCent_tmp+double((x0+nRight))*(DC2_TR+DC2_BR);
            yCent_tmp=yCent_tmp+double((y0-nUp_old))*DC2_TR+double((y0+nDown_old))*DC2_BR;
        else

            for kRight = int32(1): nDown_old+1+nUp_old
                idxRight = idx_TR + (kRight-1) * width;
                DC_Right=I_rm(idxRight);
                I_rm(idxRight)=uint8(0);

                if ~(DC_Right <threshold)
                    rightClear = false;
                    DC_max=max(DC_max,DC_Right);
                    DC2_Right= double(DC_Right)^2;
                    DC2_tot=DC2_tot+DC2_Right;
                    xCent_tmp=xCent_tmp+double(x0+nRight)*DC2_Right;
                    yCent_tmp=yCent_tmp+double(y0-nUp_old-1+kRight)*DC2_Right;
                end
            end

            if ~rightClear
                if x0+nRight < width-1
                    nRight=nRight+1;
                else
                    rightLimit = true;
                end
            end
        end


        sidesClear = (downClear||downLimit) && (leftClear||leftLimit)&&...
            (rightClear||rightLimit) && (upClear||upLimit);
    end

    if min(nLeft+nRight+1,nDown+nUp+1)+1>minROI && ...
            max(nLeft+nRight,nDown+nUp)+1<maxROI &&...
            DC_max>minDC

        k=k+1;

        xCent_all(k)=xCent_tmp/DC2_tot;
        yCent_all(k)=yCent_tmp/DC2_tot;
    end


end


%xCent=nonzeros(xCent_all)'-0.5;
%yCent=nonzeros(yCent_all)'-0.5;
 

cameraPoints = [nonzeros(xCent_all) , nonzeros(yCent_all)]'-0.5;

end
 
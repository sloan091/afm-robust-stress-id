function  [Astruct,H] = calcInfoPerformance(x1,x2,yo,ym,bins,outFlag)

 % Information diagnostics
        [H,I,~] = calcInfoDiagnostics(ym,ym,yo,bins);
        Ap = (H.Hy - I.Ix1_y)/H.Hy; 
        
        % Functional dependence of observed
         [~,Im,Partm] = calcInfoDiagnostics(x1,x2,ym,bins);
         [~,Io,Parto] = calcInfoDiagnostics(x1,x2,yo,bins);
        Afy = (Im.Ix1x2_y - Io.Ix1x2_y)./Io.Ix1x2_y;
        Afx1 = Partm.Ux1 - Parto.Ux1;
        Afx2 = Partm.Ux2 - Parto.Ux2;
        Afr = Partm.Rs - Parto.Rs;
        Afs = Partm.S - Parto.S;
        Afp = abs(Afx1) + abs(Afx2) + abs(Afr) + abs(Afs);
        
        switch outFlag
            case 1 % Structure
                Astruct = v2struct(Ap,Afy,Afx1,Afx2,Afr,Afs,Afp);
            case 2 % Array
                Astruct = [Ap,Afy,Afx1,Afx2,Afr,Afs,Afp];
        end
        
end
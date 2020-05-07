function Xout = advectV(Xu1, qtx, qty)

UUb1 = Xu1(:,:,1);
    [aa,bb] = meshgrid(1:length(qtx),1:length(qtx));

    qx = aa - qtx; qy = bb - qty;
bd = min(25,size(UUb1,1));
    UUb1 = Border(UUb1,bd); 
    Ub1 = interp2(UUb1,bd+qx,bd+qy,'bicubic');
  
    
    VVb1 = Xu1(:,:,2);

    bd = min(25,size(VVb1,1));
    VVb1 = Border(VVb1,bd); 
    Vb1 = interp2(VVb1,bd+qx,bd+qy,'bicubic');

    [qxx,qxy] = gradient(qx); [qyx,qyy] = gradient(qy);
  
   Vbt = qxx.*Vb1 -qyx .* Ub1;
   Ubt = qyy.*Ub1 -qxy .* Vb1;

   Xout(:,:,1) = Ubt; Xout(:,:,2) = Vbt;
   
   
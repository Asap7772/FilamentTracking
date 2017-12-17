Colx = 1;
Coly = 2;
for i = 1:size(I_orig,3)
    if i == 1
        excel_data(:,[Colx Coly]) = coord_reconstructed(:,[1 2],i);
    else
        Colx = Colx +2;
        Coly = Coly +2;
        excel_data(:,[(Colx) (Coly)]) = coord_reconstructed(:,[1 2],i);       
    end
    
    
end
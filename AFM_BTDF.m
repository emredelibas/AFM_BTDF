function [nwk,tree,sqrfrm] = AFM_MGMNF(filename,nameofSpecies)
clear all;
keys={'RR_D','YY_D','D_D','A','G','C','T'};
values={0,0,0,0,0,0,0};
vector=containers.Map(keys,values);
allvectors=vector.values(keys);
str = extractFileText(filename);
sequences = split(str,newline);
numofSpecies=length(sequences);
tic;
for i=1:numofSpecies
    RR_distance=0;
    YY_distance=0;
    D_distance=0;
    RR=0;YY=0;D=0;
    seqlen=strlength(sequences(i))-1;
    for j=1:seqlen-1
        if ((sequences{i}(j)=='A') && ((sequences{i}(j+1)=='A')||(sequences{i}(j+1)=='G')))   
            RR=RR+1;
            RR_distance=RR_distance+j;
        elseif ((sequences{i}(j)=='G') && ((sequences{i}(j+1)=='G')||(sequences{i}(j+1)=='A')))
            RR=RR+1;
            RR_distance=RR_distance+j;
        elseif ((sequences{i}(j)=='C') && ((sequences{i}(j+1)=='C')||(sequences{i}(j+1)=='T')))  
            YY=YY+1;
            YY_distance=YY_distance+j;
        elseif ((sequences{i}(j)=='T') && ((sequences{i}(j+1)=='T')||(sequences{i}(j+1)=='C')))  
            YY=YY+1;
            YY_distance=YY_distance+j;
        end
        if (sequences{i}(j)==sequences{i}(j+1))
            D=D+1;
            D_distance=D_distance+j;
        end
        vector(sequences{i}(j))=vector(sequences{i}(j))+(1);
    end
    vector('RR_D')=RR_distance/RR;
    vector('YY_D')=YY_distance/YY;
    vector('D_D')=D_distance/D;
    
    allvectors(i,:)=vector.values(keys);
    vector=containers.Map(keys,values);
end 
toc;
for j=1:numofSpecies
    for i=1:7
        A(j,i)=allvectors{j,i};
    end
end
hold off;
Z=pdist(A);
sqrfrm=squareform(Z);
tree = seqlinkage(Z,'average',nameofSpecies);
nwk=getnewickstr(tree,'Distances',false);
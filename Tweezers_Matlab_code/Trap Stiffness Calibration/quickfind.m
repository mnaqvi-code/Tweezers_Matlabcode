function ind=quickfind(ff,nbin)

%This function can only replace 'find' when the elements in ff are linearly  increasing from left to right
%finds indexes of ff when ff is in bin 1 to nbin and puts them in ind(i,:)
%
%Does the same as the built-in function 'find' eg. ind(i,:)=find(ff<c(i)+delta/2 &
%ff>c(i) -delta/2) where delta = c(i+1)-c(i) . But works faster.

edges=floor(length(ff)/nbin*[1:nbin]);
ind = zeros(nbin,length(1:edges(1)));
for i=1:nbin
    ind(i,:)=[edges(i)-edges(1)+1:edges(i)];
end

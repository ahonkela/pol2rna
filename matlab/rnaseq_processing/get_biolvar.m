function[b_v]=get_biolvar(y,m_groups,biol_var)

group_no=zeros(10,1);

for j=1:10
    aa=repmat(y(j),78,1);
    dist=abs(aa-m_groups);
    mindist=min(abs(aa-m_groups));
    group_no(j)=find(dist==mindist);
end

b_v=biol_var(group_no);

end


function edges=bear_boundary_masks(S,h)
% Put a segment between each pair of adjacent nodes with different values.
% Mark nodes within the specified Euclidean distance of these segments.
edges.widths=[.025,.05,.10];
edges.outer=cell(1,3); edges.inner=cell(1,3); edges.all=cell(1,3);
vertical_left=S(:,1:end-1); vertical_right=S(:,2:end);
horizontal_low=S(1:end-1,:); horizontal_high=S(2:end,:);
is_outer=@(a,b) xor(a==0,b==0);
is_inner=@(a,b) (a~=b) & (a~=0) & (b~=0);
for j=1:numel(edges.widths)
    radius=edges.widths(j)/h;
    edges.outer{j}=near_segments(is_outer(vertical_left,vertical_right), ...
        is_outer(horizontal_low,horizontal_high),size(S),radius);
    edges.inner{j}=near_segments(is_inner(vertical_left,vertical_right), ...
        is_inner(horizontal_low,horizontal_high),size(S),radius);
    edges.all{j}=edges.outer{j} | edges.inner{j};
end
% Check a straight interface halfway between two grid columns.
Sline=zeros(5,6); Sline(:,4:end)=1;
mask=near_segments(diff(Sline,1,2)~=0,diff(Sline,1,1)~=0,size(Sline),.5);
assert(isequal(mask,repmat([false false true true false false],5,1)));
end

function mask=near_segments(vertical,horizontal,grid_size,radius)
mask=false(grid_size);
[vr,vc]=find(vertical);
[hr,hc]=find(horizontal);
extent=ceil(radius+.5);
for row_offset=-extent:extent
    for col_offset=-extent:extent
        if hypot(max(abs(row_offset)-.5,0),col_offset-.5)<=radius+1e-10
            mask=mark(mask,vr+row_offset,vc+col_offset);
        end
        if hypot(row_offset-.5,max(abs(col_offset)-.5,0))<=radius+1e-10
            mask=mark(mask,hr+row_offset,hc+col_offset);
        end
    end
end
end

function mask=mark(mask,row,col)
keep=row>=1 & row<=size(mask,1) & col>=1 & col<=size(mask,2);
mask(sub2ind(size(mask),row(keep),col(keep)))=true;
end

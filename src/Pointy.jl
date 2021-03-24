module Pointy
    
using MinkowskiReduction, LinearAlgebra

#u = [2,1,0]; v = [3,4,6]; w = [0,-1,-2]

export pointgroup_basic

function pointgroup_basic(a1,a2,a3,debug=false)
    """
    Generate the symmetry operations of a lattice, defined by three 3-vectors
    (This function is _very_ simple but not efficient. Used to test more efficient algorithms.)

    ```juliadoctest
    julia> u = [1,0,0]; v = [.5,√3/2,0]; w = [0,0,√(8/3)];
    julia> pointgroup_basic(u,v,w)
    24-element Array{Array{Float64,2},1}:
    [-1.0 0.0 0.0; -1.0 1.0 0.0; 0.0 0.0 -0.9999999999999999]
    ...
    ```
    """
u,v,w = minkReduce(a1,a2,a3)
A = [u v w] # Put the lattice vectors as columns in matrix A
B = inv(A)*transpose(inv(A)) # Use this for checking for orthogonality
# A list of all possible lattice vectors in a rotated basis
c = [A*[i;j;k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# A list of all possible bases, (i.e., all combinations of c vectors)
R = [[i j k] for i ∈ c for j ∈ c for k ∈ c]
RT = [transpose(i) for i ∈ R]
# This is the U^T*U, where U transforms original basis to candidate basis
T = [R[i]*B*RT[i] for i ∈ 1:length(R)]
if debug return T end
# If T==identity then the U was a symmetry of the lattice
idx = findall([t≈I(3) for t ∈ T].==true)
Ai = inv(A)
ops = [Ai*R[i] for i in idx]
return ops
end


end


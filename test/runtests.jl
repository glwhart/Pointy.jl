using Pointy
using Test

@testset "Pointy.jl" begin
    u = [1,0,0]; v = [.5,√3/2,0]; w = [0,0,√(8/3)];
    @test length(pointgroup_basic(u,v,w))==24
    a = -u + v; b = v; c = w
    @test length(pointgroup_basic(a,b,c))==24

end

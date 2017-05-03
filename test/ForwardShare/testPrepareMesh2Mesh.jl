import jInv.ForwardShare
import jInv.Mesh
import Base.Test
import jInv.Utils

@everywhere begin
  using jInv.ForwardShare
  using jInv.Mesh
  using Base.Test
  using jInv.Utils
  type TestProb <: ForwardProbType
  	Mesh::AbstractMesh
  end
end

Minv = getRegularMesh([0 1 0 2 0 3],[12 14 16])
Mfor = getRegularMesh([0 1 0 2 0 3],[6 7 8])
pFor = TestProb(Mfor)
pFors = Array{RemoteChannel}(2)
workerList  = workers()
nw          = nworkers()
pFors[1]    = initRemoteChannel(identity,workerList[1%nw+1],pFor)
pFors[2]    = initRemoteChannel(identity,workerList[2%nw+1],pFor)


# prepare single mesh2mesh
M2Mc = prepareMesh2Mesh(pFor,Minv,true) # compact
M2Mf = prepareMesh2Mesh(pFor,Minv,false) # not compact
M2Mrc = prepareMesh2Mesh(pFors,Minv,true) # compact
M2Mrf = prepareMesh2Mesh(pFors,Minv,false) # compact

println("compare remote versions with local ones")
@test all(fetch(M2Mrc[1]).==M2Mc)
@test all(fetch(M2Mrf[1]).==M2Mf)

println("test interpGlobal2Local")
x   = randn(Minv.nc)
D1  = speye(Minv.nc)
D2  = UniformScaling(1.0)
xf1 = interpGlobalToLocal(x,M2Mc)
xf2 = interpGlobalToLocal(x,M2Mf)
xf3 = interpGlobalToLocal(D1,x,M2Mf)
xf4 = interpGlobalToLocal(D2,x,M2Mf)
@test norm(xf1-xf2)/norm(xf2) < 1e-12
@test norm(xf2-xf3)/norm(xf2) < 1e-12
@test norm(xf2-xf4)/norm(xf2) < 1e-12

println("test interpLocal2Global")
x = randn(Mfor.nc)
xf1 = interpLocalToGlobal(x,M2Mc)
xf2 = interpLocalToGlobal(x,M2Mf)
xf3 = interpLocalToGlobal(D1,x,M2Mf)
xf4 = interpLocalToGlobal(D2,x,M2Mf)
@test norm(xf1-xf2)/norm(xf2) < 1e-12
@test norm(xf2-xf3)/norm(xf2) < 1e-12
@test norm(xf2-xf4)/norm(xf2) < 1e-12

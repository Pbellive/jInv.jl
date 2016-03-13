export getFaceAverageMatrix, getEdgeAverageMatrix, getNodalAverageMatrix
export ndgrid, meshgrid
export getDivergenceMatrix, getNodalGradientMatrix, getCurlMatrix, getNodalLaplacianMatrix
export getMassMatrix, getdMassMatrix

# 1D operator
function av(n)
# A = av(n), 1D average operator
	return spdiagm((fill(.5,n),fill(.5,n)),(0,1),n,n+1)
end

function ddx(n)
# D = ddx(n), 1D derivative operator
	return spdiagm((fill(-1.0,n),fill(1.0,n)),(0,1),n,n+1)
end

# --- linear operator
"""
Af = getFaceAverageMatrix(Mesh::AbstractTensorMesh)

"""
function getFaceAverageMatrix(Mesh::AbstractTensorMesh)
# Mesh.Af = getFaceAverageMatrix(Mesh::AbstractTensorMesh) builds face-to-cc average operator
	if isempty(Mesh.Af)
		if Mesh.dim==2
			A1 = kron(speye(Mesh.n[2]),av(Mesh.n[1]))
			A2 = kron(av(Mesh.n[2]),speye(Mesh.n[1]))
			Mesh.Af = [A1 A2]
		elseif Mesh.dim==3
			A1 = kron(speye(Mesh.n[3]),kron(speye(Mesh.n[2]),av(Mesh.n[1]))) 
			A2 = kron(speye(Mesh.n[3]),kron(av(Mesh.n[2]),speye(Mesh.n[1]))) 
			A3 = kron(av(Mesh.n[3]),kron(speye(Mesh.n[2]),speye(Mesh.n[1])))
			Mesh.Af = [A1 A2 A3]
		end
	end
	return Mesh.Af
end

function getEdgeAverageMatrix(Mesh::AbstractTensorMesh)
# Mesh.Ae = getEdgeAverageMatrix(Mesh::AbstractTensorMesh) builds edge-to-cc average operator
	if isempty(Mesh.Ae)
		A1 = kron(av(Mesh.n[3]),kron(av(Mesh.n[2]),speye(Mesh.n[1]))) 
		A2 = kron(av(Mesh.n[3]),kron(speye(Mesh.n[2]),av(Mesh.n[1]))) 
		A3 = kron(speye(Mesh.n[3]),kron(av(Mesh.n[2]),av(Mesh.n[1])))
		Mesh.Ae = (1/3)*[A1 A2 A3]
	end
	return Mesh.Ae
end

function getNodalAverageMatrix(Mesh::AbstractTensorMesh)
# Mesh.An = getNodalAverageMatrix(Mesh::TensorMesh3D) builds nodal-to-cc average operator
	if isempty(Mesh.An)
		if Mesh.dim==2
			Mesh.An = kron(av(Mesh.n[2]),av(Mesh.n[1]))
		elseif Mesh.dim==3
			Mesh.An = kron(av(Mesh.n[3]),kron(av(Mesh.n[2]),av(Mesh.n[1])))
		else
			error("getNodalAverageMatrix: Dimension must be 2 or 3")
		end
	end
	return Mesh.An
end


# --- ndgrid
ndgrid(v::AbstractVector) = copy(v)
# function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
# 	m, n = length(v1), length(v2)
# 	v1 = reshape(v1, m, 1)
# 	v2 = reshape(v2, 1, n)
# 	(repmat(v1, 1, n), repmat(v2, m, 1))
# end

function ndgrid_fill(a, v, s, snext)
	for j = 1:length(a)
		a[j] = v[div(rem(j-1, snext), s)+1]
	end
end

function ndgrid{T}(vs::AbstractVector{T}...)
	n = length(vs)
	sz = map(length, vs)
	out = ntuple(i->Array(T, sz), n)
	s = 1
	for i=1:n
		a = out[i]::Array
		v = vs[i]
		snext = s*size(a,i)
		ndgrid_fill(a, v, s, snext)
		s = snext
	end
	out
end

# --- meshgrid
meshgrid(v::AbstractVector) = meshgrid(v, v)
function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
	m, n = length(vy), length(vx)
	vx = reshape(vx, 1, n)
	vy = reshape(vy, m, 1)
	(repmat(vx, m, 1), repmat(vy, 1, n))
end

function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T},
	                 vz::AbstractVector{T})
	m, n, o = length(vy), length(vx), length(vz)
	vx = reshape(vx, 1, n, 1)
	vy = reshape(vy, m, 1, 1)
	vz = reshape(vz, 1, 1, o)
	om = ones(Int, m)
	on = ones(Int, n)
	oo = ones(Int, o)
	(vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

# --- differential operators
function getDivergenceMatrix(Mesh::AbstractTensorMesh)
# Mesh.Div = getDivergenceMatrix(Mesh::AbstractTensorMesh) builds face-to-cc divergence operator
	if isempty(Mesh.Div)
		if Mesh.dim==2
			D1 = kron(speye(Mesh.n[2]),ddx(Mesh.n[1]))
			D2 = kron(ddx(Mesh.n[2]),speye(Mesh.n[1]))
			Div = [D1 D2]
		elseif Mesh.dim==3
			D1 = kron(speye(Mesh.n[3]),kron(speye(Mesh.n[2]),ddx(Mesh.n[1])))
			D2 = kron(speye(Mesh.n[3]),kron(ddx(Mesh.n[2]),speye(Mesh.n[1])))
			D3 = kron(ddx(Mesh.n[3]),kron(speye(Mesh.n[2]),speye(Mesh.n[1])))
			Div = [D1 D2 D3]
		end
		Vi = getVolumeInv(Mesh)
		F  = getFaceArea(Mesh)
		Mesh.Div = Vi*(Div*F)
	end
	return Mesh.Div
end

function getNodalGradientMatrix(Mesh::AbstractTensorMesh)
# Mesh.Grad = getNodalGradientMatrix(Mesh::AbstractTensorMesh) builds nodal-to-edge gradient operator
	if isempty(Mesh.Grad)
		if Mesh.dim==2
			G1 = kron(speye(Mesh.n[2]+1),ddx(Mesh.n[1]))
			G2 = kron(ddx(Mesh.n[2]),speye(Mesh.n[1]+1))
			Grad =[G1; G2]
		elseif Mesh.dim==3
			G1 = kron(speye(Mesh.n[3]+1),kron(speye(Mesh.n[2]+1),ddx(Mesh.n[1])))
			G2 = kron(speye(Mesh.n[3]+1),kron(ddx(Mesh.n[2]),speye(Mesh.n[1]+1)))
			G3 = kron(ddx(Mesh.n[3]),kron(speye(Mesh.n[2]+1),speye(Mesh.n[1]+1)))
			Grad =[G1; G2; G3]
		end
		Li  = getLengthInv(Mesh)
		Mesh.Grad = Li*Grad
	end
	return Mesh.Grad
end

function getCurlMatrix(Mesh::AbstractTensorMesh)
# Mesh.Curl = getCurlMatrix(Mesh::AbstractTensorMesh) builds edge-to-face curl operator
	if isempty(Mesh.Curl)
		if Mesh.dim==3
			# The Curl from edges to faces
			Dyz = kron(ddx(Mesh.n[3]),kron(speye(Mesh.n[2]),speye(Mesh.n[1]+1)))
			Dzy = kron(speye(Mesh.n[3]),kron(ddx(Mesh.n[2]),speye(Mesh.n[1]+1)))
        	
			Dxz = kron(ddx(Mesh.n[3]),kron(speye(Mesh.n[2]+1),speye(Mesh.n[1])))
			Dzx = kron(speye(Mesh.n[3]),kron(speye(Mesh.n[2]+1),ddx(Mesh.n[1])))
        	
			Dxy = kron(speye(Mesh.n[3]+1),kron(ddx(Mesh.n[2]),speye(Mesh.n[1])))
			Dyx = kron(speye(Mesh.n[3]+1),kron(speye(Mesh.n[2]),ddx(Mesh.n[1])))
        	
			# curl on the edges
			Curl = [
				 spzeros(Mesh.nf[1],Mesh.ne[1])  -Dyz   Dzy
				 Dxz   spzeros(Mesh.nf[2],Mesh.ne[2])  -Dzx
				-Dxy   Dyx   spzeros(Mesh.nf[3],Mesh.ne[3])]
			
				Fi = getFaceAreaInv(Mesh)
				L  = getLength(Mesh)
				Mesh.Curl = Fi*(Curl*L)
		else
			error("CURL is only implemented for 3D meshes.")
		end	

	end
	return Mesh.Curl
end

function getNodalLaplacianMatrix(Mesh::AbstractTensorMesh)
	if isempty(Mesh.nLap)
		G = getNodalGradientMatrix(Mesh)
		Mesh.nLap = G'*G
	end
	return Mesh.nLap	
end
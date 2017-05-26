export GlobalToLocal, getGlobalToLocal
export interpLocalToGlobal, interpGlobalToLocal
export prepareGlobalToLocal

abstract AbstractModel

"""
type jInv.InverseSolve.GlobalToLocal

Maps global model to local model

sigLocal = PForInv*sigGlobal + sigmaBackground

Fields:
	PForInv::SparseMatrixCSC - linear operator between inverse mesh and 
	                           forward mesh, e.g., linear interpolation matrix
                               and projection on active set.
	sigmaBackground::Vector  - background model
	
Constructors:
    getGlobalToLocal(P::SparseMatrixCSC)
	getGlobalToLocal(P::SparseMatrixCSC,sigmaBack::Vector)
	
Example: 
	Mesh2Mesh = getInterpolationMatrix(Minv,Mfwd)   # Minv and Mfwd have different resolutions
	sigmaBack = 1.2*ones(Minv.nc)                   # put background conductivity
	gloc      = getGlobalToLocal(Mesh2Mesh,sigmaBack)
	
	sigLocal     = gloc.PForInv\' * sigGlobal + gloc.sigmaBack
	# this call is equivalent, but needed in case the Mesh2Mesh matrix is compressed
	sigLocalFast = interpGlobalToLocal(sigGlobal,gloc.PForInv,gloc.sigmaBack)
"""
type GlobalToLocal
	PForInv::Union{SparseMatrixCSC,Vector{SparseMatrixCSC},AbstractFloat} # interpolation matrix from fwd mesh to inv mesh
	sigmaBackground::Union{Vector{Float64},AbstractModel,AbstractFloat} #  (# of cells fwd mesh)
end # type GlobalToLocal

# Constructors
getGlobalToLocal(P) = GlobalToLocal(P,1e-8)
getGlobalToLocal(P,sigBack) = GlobalToLocal(P,sigBack)
getGlobalToLocal(P,sigBack,fname) = GlobalToLocal(P,sigBack)



function prepareGlobalToLocal(Mesh2Mesh,Iact,sigmaBackground,fname)
	return getGlobalToLocal(Iact'*Mesh2Mesh,interpGlobalToLocal(sigmaBackground,Mesh2Mesh),fname)
end

function prepareGlobalToLocal(Mesh2Mesh::Vector,Iact::Vector,sigmaBackground::AbstractModel,fname)
        if ~( length(Mesh2Mesh) == length(Iact))
            error("prepareGlobalToLocal: Inconsistent Mesh2Mesh and Iact vector lengths")
        end
        n = length(Mesh2Mesh)
        PforInv = Vector{SparseMatrixCSC}()
        for i = 1:n
            push!(PforInv,Iact[i]'*Mesh2Mesh[i])
        end
        return GlobalToLocal(PforInv,interpGlobalToLocal(sigmaBackground,Mesh2Mesh))
end

function prepareGlobalToLocal(Mesh2Mesh,Iact::Vector,sigmaBackground::AbstractModel,fname)
        n = length(Mesh2Mesh)
        PforInv = Vector{SparseMatrixCSC}()
        for i = 1:n
            push!(PforInv,Iact[i]'*Mesh2Mesh)
        end
        return GlobalToLocal(PforInv,interpGlobalToLocal(sigmaBackground,Mesh2Mesh))
end


#############################################################################################################################
##### ERAN: I think that with the the update of misfitParam and the associated functions, the code below is no longer needed.
#####       Experimentally I comment out this code.
#############################################################################################################################

# function prepareGlobalToLocal(Mesh2Mesh::RemoteRef{Channel{Any}},Iact,sigmaBackground,fname)

	# M2M= take!(Mesh2Mesh)
	
	# G2L = prepareGlobalToLocal(M2M,Iact,sigmaBackground,fname)
   
	# put!(Mesh2Mesh,false)
	# return G2L
	
# end

# function prepareGlobalToLocal(Mesh2Mesh::RemoteRef{Channel{Any}},IactRef::RemoteRef{Channel{Any}},
	                           # sigmaBackgroundRef::RemoteRef{Channel{Any}},fname)
	
	# Iact = fetch(IactRef)
	# sigmaBackground = fetch(sigmaBackgroundRef)
	# return prepareGlobalToLocal(Mesh2Mesh,Iact,sigmaBackground,fname)
	
# end

# function prepareGlobalToLocal(Mesh2Mesh::Array{RemoteRef{Channel{Any}},1},Iact,sigmaBackground,fname="")
	# pMod               = Array(RemoteRef{Channel{Any}},length(Mesh2Mesh))
	# IactRef            = Array(RemoteRef{Channel{Any}},maximum(workers()))
	# sigmaBackgroundRef = Array(RemoteRef{Channel{Any}},maximum(workers()))
	# @sync begin
		# for p=workers()
			# @async begin
				# # communicate Iact and sigmaBackground
				# IactRef[p]            = remotecall_wait(p,identity,Iact)            # send active cells projection matrix
				# sigmaBackgroundRef[p] = remotecall_wait(p,identity,sigmaBackground) # send background conductivity
				# # call prepareGlobalToLocal with remote refs
				# for idx=1:length(Mesh2Mesh)
					# if p==Mesh2Mesh[idx].where
						# pMod[idx] = remotecall_wait(p,prepareGlobalToLocal,Mesh2Mesh[idx],IactRef[p],sigmaBackgroundRef[p],fname)
					# end
				# end
			# end
		# end
	# end
	# return pMod	
# end
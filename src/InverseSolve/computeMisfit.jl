export computeMisfit


"""
	Dc,F,dF,d2F = computeMisfit(...)

	Computes misfit for PDE parameter estimation problem.

	computeMisfit has several options.
"""
function computeMisfit(sig,     # conductivity on inv mesh (active cells only)
                       pMis::MisfitParam,doDerivative::Bool=true, doClear::Bool=false)
	#
	#	computeMisfit for a single forward problem
	#
	#	Notes: Everything is stored in memory on the node executing this function.
	#          The gradient and sensitivities returned are with respect to the conductivity on the inv mesh.
	#

	times = zeros(4)

	sigma,dsigma = pMis.modelfun(sig)

	tic();
		sigmaloc = interpGlobalToLocal(sigma,pMis.gloc.PForInv,pMis.gloc.sigmaBackground);
	times[1]=toq()
	tic()
		Dc,pMis.pFor  = getData(sigmaloc,pMis.pFor)      # fwd model to get predicted data
	times[2]=toq()
	tic()
		F,dF,d2F = pMis.misfit(Dc,pMis.dobs,pMis.Wd)
	times[3]=toq()

	if doDerivative
		tic()
		dF = interpLocalToGlobal(dsigma,getSensTMatVec(dF,sigmaloc,pMis.pFor),pMis.gloc.PForInv)
		times[4]=toq()
	end

	# free memory
	if doClear
		clear!(pMis.pFor.Ainv)
	end

	return Dc,F,dF,d2F,pMis,times
end # function computeMisfit


function computeMisfit(sigmaRef::Future,
                        pMisRef::RemoteChannel,
				      dFRef::RemoteChannel,
                  doDerivative,doClear::Bool=false)
	#
	#	computeMisfit for single forward problem
	#
	#	Note: model (including interpolation matrix) and forward problems are RemoteRefs
	#
	rrlocs = [ pMisRef.where  dFRef.where]
	if !all(rrlocs .== myid())
		warn("computeMisfit: Problem on worker $(myid()) not all remote refs are stored here, but rrlocs=$rrlocs")
	end

	sigma = fetch(sigmaRef)
	pMis  = take!(pMisRef)

	Dc,F,dFi,d2F,pMis,times = computeMisfit(sigma,pMis,doDerivative,doClear)

	put!(pMisRef,pMis)
	# add to gradient
	if doDerivative
		dF = take!(dFRef)
		put!(dFRef,dF += dFi)
	end
	# put predicted data and d2F into remote refs (no need to communicate them)
	Dc  = remotecall(identity,myid(),Dc)
	d2F = remotecall(identity,myid(),d2F)

	return Dc,F,d2F,times
end


function computeMisfit(sigma,
	pMisRefs::Array{RemoteChannel,1},
	doDerivative::Bool=true,
	indCredit=collect(1:length(pMisRefs)))
	#
	#	computeMisfit for multiple forward problems
	#
	#	This method runs in parallel (iff nworkers()> 1 )
	#
	#	Note: ForwardProblems and Mesh-2-Mesh Interpolation are RemoteRefs
	#								(i.e. they are stored in memory of a particular worker).
	#

	F   = 0.0
	dF  = (doDerivative) ? zeros(length(sigma)) : []
	d2F = Array{Any}(length(pMisRefs));
	Dc  = Array(Future, size(pMisRefs))

	indDebit = []
	updateRes(Fi,idx) = (F+=Fi;push!(indDebit,idx))
	updateDF(x) = (dF+=x)

	# find out which workers are involved
	workerList = []
	for k=indCredit
		push!(workerList,pMisRefs[k].where)
	end
	workerList = unique(workerList)

	# send sigma to all workers
	sigRef = Array(Future,maximum(workers()))
	dFiRef = Array(RemoteChannel,maximum(workers()))

	times = zeros(4);
	updateTimes(tt) = (times+=tt)

	@sync begin
		for p=workerList
			@async begin
				# communicate model and allocate RemoteRef for gradient
				sigRef[p] = remotecall(identity,p,sigma)   # send conductivity to workers
				dFiRef[p] = initRemoteChannel(zeros,p,length(sigma)) # get remote Ref to part of gradient
				# solve forward problems
				for idx=1:length(pMisRefs)
					if pMisRefs[idx].where==p
						Dc[idx],Fi,d2F[idx],tt = remotecall_fetch(computeMisfit,p,sigRef[p],pMisRefs[idx],dFiRef[p],doDerivative)
						updateRes(Fi,idx)
						updateTimes(tt)
					end
				end

				# sum up gradients
				if doDerivative
					updateDF(fetch(dFiRef[p]))
				end
			end
		end
	end
	return Dc,F,dF,d2F,pMisRefs,times,indDebit
end


function computeMisfit(sigma,pMis::Array,doDerivative::Bool=true,indCredit=collect(1:length(pMis)))
	#
	#	computeMisfit for multiple forward problems
	#
	#	This method runs in parallel (iff nworkers()> 1 )
	#
	#	Note: ForwardProblems and Mesh-2-Mesh Interpolation are stored on the main processor
	#		  and then sent to a particular worker, which returns an updated pFor.
	#
	numFor   = length(pMis)
 	F        = 0.0
    dF       = (doDerivative) ? zeros(length(sigma)) : []
 	d2F      = Array{Any}(numFor)
 	Dc       = Array{Any}(numFor)
	indDebit = []

	# draw next problem to be solved
	nextidx() = (idx = (isempty(indCredit)) ? -1 : pop!(indCredit))

 	updateRes(Fi,dFi,idx) = (F+=Fi; dF= (doDerivative)? dF+dFi : []; push!(indDebit,idx))

	times = zeros(4);
	updateTimes(tt) = (times+=tt)

 	@sync begin
 		for p = workers()
 				@async begin
 					while true
 						idx = nextidx()
 						if idx == -1
 							break
 						end
 							Dc[idx],Fi,dFi,d2F[idx],pMis[idx],tt = remotecall_fetch(computeMisfit,p,sigma,pMis[idx],doDerivative)
 							updateRes(Fi,dFi,idx)
							updateTimes(tt)
 					end
 				end
 		end
 	end

 	return Dc,F,dF,d2F,pMis,times,indDebit
 end

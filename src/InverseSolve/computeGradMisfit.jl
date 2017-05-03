export computeGradMisfit

function computeGradMisfit(sig,Dc::Array,pMis::MisfitParam)
	#
	#	gc = computeGradMisfit(mc,model,Dc,Dobs,Wd,misfit,pFor)
	#
	#	computes gradient of misfit for single forward problem
	#
	#	Note: all variables have to be in memory of the worker executing this method
	#
	try
		sigma,dsigma = pMis.modelfun(sig)
		sigmaloc = interpGlobalToLocal(sigma,pMis.gloc.PForInv,pMis.gloc.sigmaBackground)
		F,dF,d2F = pMis.misfit(Dc,pMis.dobs,pMis.Wd)
		return interpLocalToGlobal(dsigma,getSensTMatVec(dF,sigmaloc,pMis.pFor),pMis.gloc.PForInv)
	catch err
		if isa(err,InterruptException)
			return -1
		else
			throw(err)
		end
	end
end

function computeGradMisfit(sigmaRef::Future,
                           DcRef::Future,
                           pMisRef::RemoteChannel,
                           dFiRef::RemoteChannel)
	#
	#	gc = computeGradMisfit(mc,model,Dc,Dobs,Wd,misfit,pFor)
	#
	#	compute gradient of misfit for single forward problem
	#
	#	Note: model and forward problem are represented as RemoteReferences.
	#
	#	!! make sure that everything is stored on this worker to avoid communication !!
	#
	rrlocs = [sigmaRef.where pMisRef.where DcRef.where]
	if !all(rrlocs .== myid())
		warn("computeGradMisfit: Problem on worker $(myid()) not all remote refs are stored here, but rrlocs=$rrlocs")
	end

	tic()
	pMis  = take!(pMisRef) # this is a no-op if pFor is stored on this worker
	sigma = fetch(sigmaRef)
	Dc    = fetch(DcRef)
	commTime = toq()

	tic()
	dFt  = computeGradMisfit(sigma,Dc,pMis)
	compTime = toq()

	tic()
	dFi  = take!(dFiRef)
	put!(dFiRef,dFi+dFt)
	put!(pMisRef,pMis)    # does not require communication if PF lives on this worker
	commTime += toq()

#	if commTime/compTime > 1.0
#		warn("computeGradMisfit: Communication time larger than computation time! commTime/compTime = $(commTime/compTime)")
#	end
	return true,commTime,compTime
end


function computeGradMisfit(sigma,
	                   DcRef::Array{Future,1},
	                   pMisRefs::Array{RemoteChannel,1},
	                   indFors=1:length(pMisRefs))
	#
	#	gc = computeGradMisfit(mc,model,Dc,Dobs,Wd,misfit,pFor)
	#
	#	compute gradient of misfit for multiple forward problems
	#
	#	Note: models and forward problems are represented as RemoteReferences.
	#
	#
	tic()
	# find out which workers are involved
	workerList = []
	for k=1:length(pMisRefs)
		push!(workerList,pMisRefs[k].where)
	end
	workerList = unique(workerList)
	# send sigma to all workers
	sigmaRef = Array(Future,maximum(workers()))
	dFiRef   = Array(RemoteChannel,maximum(workers()))
	dF = zeros(length(sigma))

	commTime = 0.0
	compTime = 0.0
	updateTimes(c1,c2) = (commTime+=c1; compTime+=c2)
	updateDF(x) = (dF+=x)


	tic()
	@sync begin
		for p = workers()
				@async begin
					# send model to worker and get a remote ref
					tic()
						sigmaRef[p] = remotecall(identity,p,sigma)   # send sigma to workers
						  dFiRef[p] = initRemoteChannel(identity,p,zeros(length(sigma))) # get remote Ref to part of gradient
					c1 = toq()
					updateTimes(c1,0.0)

					# the actual computation
					for idx=indFors
						if pMisRefs[idx].where==p
							isDone,c1,c2 = remotecall_fetch(computeGradMisfit,p, sigmaRef[p],DcRef[idx],pMisRefs[idx],dFiRef[p])
							updateTimes(c1,c2)
						end
					end

					# fetch result and add
					tic()
					updateDF(fetch(dFiRef[p]))
					c1 = toq()
					updateTimes(c1,0.0)
				end
		end
	end
	chkTime=toq()

	totalTime = toq()
	return dF
end

###############################################################################
#                                                                             #
# Method below is old, currently broken and possibly obsolete. However,       #
# it seems to include functionality for interrupting timed-out workers        #
# that is not present anywhere else in jInv. For that reason I keep it around #
# but commented out in case anyone wants to resurrect it in the future.       #
#                                                                             #
###############################################################################


# function computeGradMisfit(mc,Dc::Array,pMis::Array{Future},M2M::Array{Future},indFors=1:length(PF))
# 	#
# 	#	gc = computeGradMisfit(mc,model,Dc,Dobs,Wd,misfit,pFor)
# 	#
# 	#	compute gradient of misfit for multiple forward problems
# 	#
# 	#	Note: forward problems and interpolation matrices are represented as RemoteReferences.
# 	#
# 	#	!! there is only one(!) model. Interpolation matrix has to be entered for each problem !!
# 	#
# 	#	NEW: First worker that finishes takes the time and interrupt hanging workers.
# 	#
#
# 	dF     = zeros(length(mc))
# 	times  = fill(-1.0,length(PF))
#
# 	updateRes(dFi,tm,idx) = (dF+=dFi; times[idx]=tm)
#
# 	fwid = 0; tmStop = 0.0
# 	stopTime() = tmStop
# 	function isFirstWorker(p)
# 		if fwid ==0 # p is first worker that finished
# 			fwid = p # store p's id
# 			tmStop = time() + 3 * mean(times[times.>-1])  # set stop time
# 			return true
# 		else
# 			return false
# 		end
# 	end
#
# 	sigma,dsigmadm = mfun(mc)
#
# 	@sync begin
# 		for p = workers()
# 			@async begin
# 				for idx=indFors
# 					if PF[idx].where==p
# 						tic()
# 						dFi = remotecall_fetch(p, computeGradMisfit,sigma,gloc,Dc[idx],Dobs[idx],Wd[idx],misfit,PF[idx],M2M[idx])
# 						tm = toq()
# 						if dFi != -1; updateRes(dFi,tm,idx); else; break; end
# 					end
# 				end
#
# 				if isFirstWorker(p)
# 					while true
# 					  wt = remotecall_fetch(p,sleep,.01)
# 					  if minimum(times[indFors])>0.0 # all workers are done
# 						  break
# 					  elseif time()  > stopTime() # time to stop
# 							for pw=workers();
# 								if p!=pw; interrupt(pw) end;
# 							end
# 						  break
# 					  end
# 				  	end
# 				end
# 			end
# 		end
# 	end
# 	dF = (dsigmadm' * dF) / length(times[times.>-1])
# 	if length(times[times.>-1]) != length(indFors)
# 		println("computeGradMisfit was interrupted after computing ", length(times[times.>-1]) ," of ", length(indFors), " gradients")
# 	end
# 	return dF
# end

function computeGradMisfit(sigma,Dc::Array,pMis::Array{MisfitParam},indFors=1:length(pMis))
	#
	#	gc = computeGradMisfit(mc,model,Dc,Dobs,Wd,misfit,pFor)
	#
	#	compute gradient of misfit for multiple forward problems
	#
	#	Note: models, interpolations and forward problems are stored on main process and then sent to workers
	#
	#	!! this method may lead to more communication than the ones above !!
	#
	#
	numFor = length(pMis);

	# get process map
	i      = 1; nextidx = p -> i;
	procMap = zeros(Int64,numFor)
	if isa(pMis[1].pFor.Ainv,MUMPSsolver) &&  (pMis[1].pFor.Ainv.Ainv.ptr !=-1)
		for ii=1:numFor
			if any(ii.==indFors)
				procMap[ii] = pMis[ii].pFor.Ainv.Ainv.worker
			end
		end
		nextidx = p -> ( ind = find(procMap.==p);
			if !isempty(ind);
				ind = ind[1];
				procMap[ind] = -1;
			end;
			return ind)
	else
		nextidx = p -> (idx=i; i+=1; idx)
	end
	dF = zeros(length(sigma))
	updateRes(dFi) = (dF+=dFi)

	@sync begin
		for p = workers()
				@async begin
					while true
						idx = nextidx(p)
						if isempty(idx) || idx > numFor
							break
						end
						if any(idx.==indFors)
							dFi = remotecall_fetch(computeGradMisfit,p,sigma,Dc[idx],pMis[idx])
							updateRes(dFi)
						end
					end
				end
		end
	end
	return dF
end

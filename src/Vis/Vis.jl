"""
module jInv.Vis

Visualization tools for jInv

"""
module Vis
	hasPyPlot = false
	try
		if myid()==1
			using PyPlot
			hasPyPlot = true
		end
	catch
	end
	
	using jInv.Mesh
	
	if hasPyPlot
		include("plotGrid.jl")
		include("viewImage2D.jl")
		include("viewOrthoSlices2D.jl")
	end
end
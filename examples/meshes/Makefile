
# Provide a series of meshes in which the number of vertices doubles every
# time. Since the relationship of the number of vertices |V| and the number
# of tetrahedral elements |E| is linear,
#    6*|V| \approx |E|,
# we want the number of elemnts to double.
# In Gmsh, one can scale the all characteristic lengths by the -clscale
# command line argument. When the length scales are halved, the number of 
# tetrahedral elements increases approximately by a factor of 6, i.e.,
#
#     |E*| / |E| = (|l*| / |l|)^alpha
#
# with
#
#     alpha = -ln(6)/ln(2) \approx = -2.58.
#
# To get double the number of elements, reduce the characterstic length by a
# factor of
#
#     2^(1/alpha) \approx 0.76479616424594177.
#     
cheese:
	@echo "Swiss cheese."
	@for scaling in 1.0 0.76479616424594177 0.58491317284530553 0.44733935100901323 0.34212341976796229 0.26165467913724177 0.20011249496116515 0.15304526856398445 0.11704843435372532 0.089518193624722034 0.068463171114412941 0.052360370660416572; do \
	    echo "scaling $$scaling..."; \
	    gmsh swiss-cheese3d.geo -3 -optimize -clscale $$scaling -o cheese-$$scaling.vtk > cheese-$$scaling-out.log 2> cheese-$$scaling-err.log; \
	    tail cheese-$$scaling-out.log | grep vertices; \
	done
	@echo "done."

clean:
	rm -f cheese-*-out.log
	rm -f cheese-*-err.log
	rm -f cheese-*.vtk

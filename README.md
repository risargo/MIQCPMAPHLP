**Mixed-integer quadratic conic formulations for two competitive multiple allocation p-hub location problems**

In the small-instances and large-instances folders, one can find the instances and source codes used to
assess the mixed-integer quadratic conic formulations for two competitive multiple allocation p-hub location problems.
One problem allows multiple paths to flow the demands of origin-destination nodes; while the other problem
restricts these flows to a single path.  

The formulations were implemented in C++ using the ILOG API to access to CPLEX solver. After compiling, 
just call the executable followed by the discount factor alpha, a scalar value between 0 <= alpha <= 1, to represent
the scale economics on the inter-hub connections.

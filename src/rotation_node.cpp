#include "rotation_node.hpp"



bool operator < (const RotNode & n1, const RotNode & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
            return n1.l > n2.l;
			//return n1.ub > n2.ub;
	}
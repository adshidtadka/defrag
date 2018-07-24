/* Given parameter */

param N integer, >0; 		/*Nodes number*/
param M integer, >0;		/* Spectrum */

set V := 0..N-1;  			/*list of nodes*/
set E := {V,V}; 			/*list of s/d couples*/
set EM within E; 			/*list of links*/
		
/* Decision variables */

var f{E,EM} binary >= 0;
var MS >= 0;
var TS >= 0;

/* Objective function */

minimize IND: TS + MS;

/* Constraints */

s.t. CI1 : 
		TS >= sum{(s,d) in E, (i,j) in EM: s != d} f[s,d,i,j];
		
s.t. CI2 {(i,j) in EM}: 
		MS >= sum{(s,d) in E: s != d} f[s,d,i,j];

s.t. TDS {(s,d) in E : s != d}: 
		sum{(i,j) in EM : i = s} f[s,d,i,j] = 2;
		
s.t. TDD {(s,d) in E : s != d}: 
		sum{(i,j) in EM : j = d} f[s,d,i,j] = 2;
		
s.t. SpecCont {(s,d) in E, j in V: s != d}: sum {i in V : (i,j) in EM && j != d} (f[s,d,i,j]) = sum {i in V: (j,i) in EM && j!=s} (f[s,d,j,i]);

solve;

display {(s,d) in E, (i,j) in EM: s != d && f[s,d,i,j] == 1} f[s,d,i,j];

end;

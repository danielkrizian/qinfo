ncount:{count[x]-$[type x;sum null x;0i {x+null y}/ x]}
nsum:{$[type x;sum x;0i {x+0i^y}/ x]}
navg:{$[type x;avg x;nsum[x]%ncount x]}
nvar:{$[type x;var x;navg[x*x]-m*m:navg x]}
ndev:(')[sqrt;nvar]
nsvar:{$[type x;svar x;(n*nvar x)%-1+n:ncount x]}
nsdev:(')[sqrt;nsvar]

norm:{[d;a;x]$[0h>type first x; d[x;a x]; d[;a x]peach x]}
demean:norm[-;navg] / centered
zscore:(')[norm[%;nsdev];demean] / feature normalization (centered/unit variance)
rnorm:{-6f+sum x cut (12*x)?1f}
cvm:{(x+flip(not n=\:n)*x:(n#'0.0),'(x$/:'(n:til count x)_\:x)%count first x)-a*\:a:avg each x} / covariance matrix (8 times faster than x cov/:\:x)
crm:{cvm[x]%u*/:u:dev each x} / correlation matrix
cm:{ / correl coefs to matrix
 $[0>type x;rotate[-1]\[1,(count x)#x];x]}
copula:{[n;c]
 x:#[r,n;] rnorm n*r:count first c:cm c;
 U:{flip .qml.mchol x};
 U[c]$inv[U cvm x]$x}

c2o:{$[type x;where x;cross[til count x;count[x]+til count first x] where raze x]} / generate (o)bservations from (c)ounts

x:(200+200?.1),(5+200?.1),-6000+200?.1
x:"f"$-469,(200?.1),200
(.info.H.p .info.pd.f[.info.bins.tails[10];x])%log 10     / proportional entropy

x:copula[10000;.8] / two correlated variables

t1: c2o (7 1 3; 87 18 84; 12 3 4; 9 1 7)
c1:4 2 3 0 2 4 0 0 2 1 1
"1.968382"~.Q.f[6] .info.H.x c2o c1

c2:(80 100 120;20 100 80)
c3: (1 2 3 ; 6 5 4 )
"1.701464" ~ .Q.f[6] .info.H.p .info.odds group c2o c2
"0.02646764"~.Q.f[8] .info.mi.d c2o c2
"0.03450492"~.Q.f[8] .info.mi.d c2o c3

c4:(41 0;38 0;0 92)
all (1f~U c4;"0.6833207"~.Q.f[7] U flip c4)
mutinf c2o m2
c2o m3
theilu[c2o m3;(c2o m3)[;0]]

mmm:c2o m2
mmm~c2o m2

chi2 count each group c2o (24 20;36 5)
chi2 count each group c2o (24 20;5 0;0 10) / case with zero counts
\t:100000 entropy[mmm[;0]]+entropy[mmm[;1]]-entropy mmm
\t:100000 mutinf mmm

(`stat`df`p!("5.4885459";"6.0000000";"0.4828422"))~.Q.f[7] each chi2 count each group t1

x:{(neg[x]?/:2#x)} 100000
x:gencop[4000000;-.5f]
x:2#enlist til 1000
n:100000; k:50;r:n cut (k*n)?1f;sum1234:sum r s:-4?k;sum12:sum r -2?s;sume:sum r e:til[k] except s;rand1:r first 1?e;sig1:r first 1?s
(.info.mi.f (sum1234;rand1);.info.mi.f (sum1234;sume);.info.mi.f (sum1234;sig1);.info.mi.f (sum1234;sum12))
.[cor;(sum1234;sig1)]

split:{
 0N!"spl ",-3!x;
 if[(5=first x)|3>=c:count x;:enlist x];
 0N!"piv ",-3!piv:"j"$c%2;
 (0,piv) cut x}


part:{[p;x] piv(piv:count[x]*0f,p) bin "f"$x};
grid2x2:{[p]
 raze {x group flip {y+part[.5]x-y}'[flip x;y]}'[value p;key p]}
/ raze {group flip x}each{y+piv[.5]x-y}''[flip each value p;key p]


uniform:{.05<=chi2[x]`p}
small:{4>count x} 
split:{$[small x; x;
  not uniform count each r:rect[x]; b;
  uniform count each rect r; x;
  r]}

x:{flip (neg[x]?/:2#x)} 10000
.[cor;flip x]
x:2#/:til 13
x:p:enlist[0f]!enlist x
{count first 0N!x}
count each grid2x2 x
x2:grid2x2 grid2x2 x
(where small each x2)_x2


\t split p
\ts uniform count each part[p:(enlist 0 0)!enlist x]
split:{r:rect[x] $[small x; x;
  not uniform count each r:rect[x]; b;
  uniform count each rect r; x;
  r]}

//////
x:flip rank each gencop[400000;.8f]
pf:reciprocal {x%2}\[1f<;count x] / partitioning factor
uniform:{.05<=chi2[count each x]`p} 
small:{30>count x} 
part:{[pf;x] x group floor x*pf}
split:{[pf;x]
 nonsplitable:$[
  not uniform p:part[pf] x;p where small each p:value p;
  small x;enlist x;
  not uniform part[pf*2;x];p where small each p:value p;
  enlist x];
 splitable:$[nonsplitable~enlist x;();p except nonsplitable];
 / 0N!-3!count each splitable;
 / 0N!-3!count each nonsplitable;
 (splitable;nonsplitable)}

dolevel:{[pf;x]
 stack:split[pf l:x 0] each x 1;
 :(l+1;raze stack[;0];x[2],raze stack[;1])}

sum {KL[count[x]%y;prd %[;y]{1+max[x]-min x} @/: flip x]}[;count x] each last dolevel[reciprocal {x%2}\[1f<;count x] ]/[{count x 1};(1;enlist x;())]
.[cor;flip x]

////////////
x:flip rank each gencop[4000000;.8f]
x:{flip (neg[x]?/:2#x)} 30
x:2#/:til 29
pf:reciprocal {x%2}\[1f<;count x] / partitioning factor
uniform:{.05<=chi2[count each x]`p} 
small:{30>count x} 
split:{[pf;x] x group floor x*pf}

part:{[pf;x]
 p:split[2*pf] x;
 nonsplitable:$[not u:uniform p; p where small each p;
  small x;enlist x;
  not uniform split[4*pf] x;p where small each p;
  enlist x];
 splitable:$[nonsplitable~enlist x;();p except nonsplitable];
 nonsplitable,raze part[pf*2] each splitable}

jmc:{count[x],{1+max[x]-min[x]}@/:flip x}; / (joint,marginal) counts of partition x
jminmax:{c1ount[x],(min;max)@/:\:flip x}
mi.f:{
 x:flip rank each flip x;
 jmd:%[;tc] jmc each part[reciprocal tc:count x] x; / (joint,marginal) probability of each partition
 .[KL;] {(x 0;prd 1_x)} flip jmd / Kullback-Leibler divergence between joint and product of marginals
 }
mi.f x
count each
split[2*reciprocal count x] x
.[KL;]prd
1_flip %[;tc]
jminmax
jmc each part[reciprocal tc:count x] x
.517 * log .517%.2675
mi.f gencop[400000;.0]
MIf flip x
{KL[x 0;prd 1_x]} flip %[;count x] jmc each part[reciprocal count x] x
reciprocal count x
xx:jmc each part[reciprocal tc:count x] x

nbyn:{[n;x;y]  / e.g. nbyn[4] splits x into 4-by-4 rectangles while keeping all ties in the same rectangle
 n:(2 4!1 2)n;
 pivots:{`s#1_raze flip ({x x binr floor {(x+y)%2}':[y]}[x;y]; y)};
 pivotsn:{x pivots[y]/z};
 bounds:{(x,'next x) y};
 split:{flip (bounds'[x;] each key[g]; flip .[y;] (::;) value g:group flip bin'[x;y])};
 split[pivotsn'[n;x;y 0];y 1]}

uniform:{
 c:{count first x} each x[;0]!x[;1];
 .05<=chi2[c]`p}
small:{30>count first x 1}
part:{[nt;x]
 p:nbyn[2;nt] x;
 nonsplitable:$[not u:uniform p; p where small each p;
  small x;enlist x;
  not uniform nbyn[4;nt] x;p where small each p;
  enlist x];
 splitable:$[nonsplitable~enlist x;();p except nonsplitable];
 nonsplitable,raze part[nt] each splitable}

jmc:{(count first x[1]), {x[1]-x 0}each x 0} / (joint;marginal) counts of the partition
mi.f:{
 nobs:count first x;
 i:iasc each x; / indices of sorted values
 nt:{[x;i] asc i where differ x i}'[x;i]; / first index of each distinct value, ignores ties
 x:iasc each i; / use integer ranks
 p0:(count[x]#enlist 0,nobs:count first x; x); / initial partition
 p:part[nt] p0;
 c: flip %[;nobs] jmc each p;
 KL[c 0;prd 1_c]}  

x:{(neg[x]?/:2#x)} 100000
x:gencop[4000000;-.5f]
x:2#enlist til 1000
r:1000000 cut 4000000?1f;
sum12:sum r[0 1]
sum34:sum r[2 3]
sum1234:sum r
mi.f x
.info.mi.f (sum1234;sum r[)
.[cor;x]
.[cor;(sum34;r 3)]

nb[9?9;4]
bins.equal[4;9?9]

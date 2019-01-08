\d .info
norm:{[d;a;x]$[0h>type first x; d[x;a x]; d[;a x]peach x]}
prb:norm[%;sum]                 / convert densities into probabilities
odds:{[g]prb count each g}
bins.equal:{[n;x] -0w,lo+1_til[n]*(max[x]-lo:min x)%n} / split range into n# equal width bins
bins.tails:{[n;x] / split into n-2 equal sized interior bins + 2 tail bins that minimize the interior range
 i:iasc x;istop:floor (c:count x)*1-2%n;
 xir:x i ir[;im:{first where x=min x} {y-x} . x i ir:(;s+istop) s:til c-istop];
 -0w,first[xir],(1_bins.equal[n-2;`s#xir]),last xir}
nbin.hgr:{"j"$(ee%6)+(1%3)+2%3*ee:xexp[;1%3] 8+(324*x)+12*sqrt (36*x) + 729*x xexp 2} / Hacine-Gharbi,Ravier number of bins

pd.f:{[b;x] / probability distribution of continuous variable discretized using (b)inning function
 `s#odds b#group b (b:b[x]) bin x}

pdf.chi2:{[df;x] reciprocal[xexp[2;df%2]*gamma df%2] * xexp[x;-1+df%2]*exp neg x%2}
cdf.chi2:{[df;x] cdf.gamma[df%2;2.;x]}
chi2:{[c]       / chi-squared test;
 if[0>type first key c;
  :chi2p[c; 1%count c]];
 if[0=df:prd {-1+count x}@/:dims:{enlist each distinct x}@'flip key c;
  :`stat`df`p!(0n;0;0n)];
 c:0^((cross) . dims)#c;
 e:%[;sum c] prd ((sum;c:value c) fby)@/:flip key c;
 yates:neg $[1=df;min .5,abs c-e;0];
 stat:sum %[;e] xexp[yates+abs c-e;2];
 `stat`df`p!(stat;df;1f-.qml.c2cdf[df;stat])}

chi2p:{[c;p]    / goodness-of-fit; univariate (c)ounts with expected (p)robabilities
 df:-1+count c;
 stat:sum %[;e] xexp[abs c-e:p*sum c;2];
 `stat`df`p!(stat;df;1f-.qml.c2cdf[df;stat])}

KL:{sum x*log x%y} / Kullback-Leibler divergence between joint distribution and the product of marginal distributions

E:{[p;x]sum p*x} / expected value of outcomes x occurring with probabilities p
I:{[p] neg log p} / information content/self-information/surprisal/signal of event with probability p
H.p:{[p] E[p] I p} / entropy
H.d:{H.p odds group x} / discrete variable entropy
H.f:{{H.p pd.f[x;y]}[bins[`equal;nbin.hgr count x ];x]} / continuous variable entropy
H.x:{H[`d`f`d 7 9 11h?type x] x}

mut.p:{[p]                             / mutual information given (margins combinations!joint (p)robabilities)
 m:((sum;value p) fby)@/: flip key p; / marginal probabilities (p(x);p(y); ...) as sums along margins/dimensions
 KL[p;prd m]} / Kullback-Leibler divergence btwn joint and product of marginals

mut.d:{[x] / mutual information between discrete (x;y; ...) variables
 j:odds group x; / joint probability p(x and y and ...)
 mut.p j}

\d .dv

uniform:{x<=.info.chi2p[y;1%count y]`p}
split:{(0,floor x*count y) cut y}
members:{(inter/) x}
counts:{count members x}'
grid:{raze .[(;)/:\:;] split'[x;y]}
grid2x2:grid[.5]
part:{[pf;pval;big;x]
 p:$[count x 2; x[2;0]; pf x 0];
 n:$[count x 2; x[2;1]; counts p];
 stack:flip $[not uniform[pval] n;(p;n;count[n]#());
  $[not uniform[pval] raze n2:counts each p2:pf each p;(p;n;flip (p2;n2));
   ()]];
 done:$[count stack;();enlist x];
 if[count stack;done:stack where small:n<big;
  stack@:where not small];
 (stack;done)}

ap:{[pf;pval;big;x]  / adaptive partitions
 unstack:{[f;stack;done] (-1_stack;done),'f last stack}; / process last on stack
 repeat:(unstack part[pf;pval;big;]).;      / this will be repeated until empty stack
 x:iasc each x;                             
 x:(enlist (x;n:count first x;());());      / initial stack
 @[;1] repeat/[{count x 0};x]}

mut:{n:sum x[;1];
 .info.KL . (x[;1]%n;(prd')%[;n](count'')x[;0])}

\d .

.info.mut.f:{.dv.mut .dv.ap[.dv.grid2x2;.05;30;x]}

/ TODO: ties
/ i:iasc each x; / indices of sorted values
/ nt:{[x;i] asc i where differ x i}'[x;i]; / first index of each distinct value, ignores ties
/ x:iasc each i; / use integer ranks

/ entropy:{neg sum x*log x:odds( group x})@
/ theilu:{mi.i[x]%entropy y}
/ f:{x%sum sum each x}                        / relative (f)requency
/ KLd:{sum sum each x*log x%y}                 / Kullback-Leibler divergence from random var X to Y with corresponding prob mass functions x and y
/ H:{ neg sum sum each {x*log x} f x}         / univariate or joint entropy, depending whether x is frequency vector or contingency matrix 
/ Hr:{ H[x]%log count x }                     / entropy relative
/ mi:{                                        / mutual information of bivariate matri(x) of frequencies
/  frq:f x;                                   / relative frequencies
/  outer:(sum each frq)*\: sum each flip frq; / outer product of row totals with column totals; result is rows*cols matrix
/  KLd[frq;outer] }                            / Kullback-Leibler divergence
/ U:{ mi[x] % H sum each flip x}              / Theil's U; aka (U)ncertainty coefficient. See wikipedia

/ NOTE:
/ If empirical frequencies are used for x (bin frequencies) then this is this is maximum likelihood estimator.
/ LIMITATION: If there are many zero counts and the sample size is small, the MLE is very inefficient and strongly biased.
/ Alternative but slower calculation of mutual information:
/ mi:{ H[sum peach x] + H[sum x] - H x}       / mutual information; H(X)+H(Y)-H(X,Y) where X/Y are row/column marginals


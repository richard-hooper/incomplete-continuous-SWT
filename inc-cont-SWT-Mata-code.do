/*  A Mata function in Stata to search for designs for incomplete     */
/*  stepped wedge trials in continuous time that achieve given power  */
/*  with the fewest participants.                                     */

/*  Richard Hooper [1], Jessica Kasza [2], Andrew Forbes [2]          */
/*  [1] Queen Mary University of London, London UK                    */
/*  [2] Monash University, Melbourne Australia                        */

/*  October 2019                                                      */


mata: mata set matastrict on

mata:
real matrix opti_ctsw(real scalar halfnseq,  ///
		real scalar nper, real scalar rho,  ///
		real scalar tau, real scalar effect,  ///
		real scalar order)

/*  halfnseq is half the number of clusters i.e. possible             */
/*  "sequences"; nper is the number of "periods" i.e. participants    */
/*  per cluster; rho is the correlation between outcomes of two       */
/*  participants from the same cluster sampled at the same time; tau  */
/*  is the parameter that determines the decay in this correlation    */
/*  over time; effect is the effect size; order is the order or       */
/*  degree of the polynomial for the time effect.                     */

{

	real matrix bestnmiss, bestnmiss1, trynmiss,  ///
			nmiss
	real colvector bestxover, bestxover1, tryxover,  ///
			xover, xul, nll, nul
	real scalar nseq, bestse, tryse, se, targse,  ///
			i, j, iopp, xll, n, n1, inc, tailseq, minn

/*  The algorithm works on the top half of the sequences (halfnseq    */
/*  sequences in total) and then copies these to the bottom half,     */
/*  but reversing time and swapping control and intervention.         */
	
	nseq = 2*halfnseq

/*  The target for the standard error of the treatment effect is      */
/*  calculated to achieve 90% power at the 5% significance level.     */

	targse = effect/(1.9600+1.2816)

/*  A design is defined by the matrix nmiss ("non-missing") which is  */
/*  1 if a participant is to be recruited and . if not, and the       */
/*  vector xover ("cross-over") which determines when each sequence   */
/*  crosses over.                                                     */

/*  STEP 1: For the backward search, initialise the design to a       */
/*  complete design with classic stepped wedge form:                  */

	inc = max((floor(nper/(nseq-1) + 0.5),1))
	tailseq = floor((nseq-1)/(2*nper)+1)
	minn = 2*inc*(nseq-tailseq)
	nmiss = J(nseq, nper, 1)
	xover = floor((0 :: nseq-1) :* (nper/(nseq-1))  ///
			:+ 0.5)
	xover[(halfnseq+1 .. nseq)] =  nper :-  ///
			xover[(halfnseq .. 1)]

/*  The program calls the function gls_polyt_se which calculates the  */
/*  standard error of the treatment effect for the current design.    */

	se = gls_polyt_se(xover, nmiss, rho, tau, order)

/*  STEP 2: Look for modifications which reduce the standard error    */
/*  without changing the sample size:                                 */
	
	do {

		xul = (xover[2 :: halfnseq] \ nper/2)
		xll = 0
		bestse = se

		for (i=1; i<=halfnseq; i++) {
		
			iopp = nseq+1-i
			bestxover = xover

			if (bestxover[i]-1 >= xll) {

/*  Try shifting the cross-over in cluster i to the left:             */

				tryxover = bestxover
				tryxover[i] = bestxover[i]-1
				tryxover[iopp] = bestxover[iopp]+1
				tryse = gls_polyt_se(tryxover,  ///
						nmiss, rho, tau, order)
				if (tryse < se) {
					xover = tryxover
					se = tryse
				}
			}	

			if (bestxover[i]+1 <= xul[i]) {

/*  Try shifting the cross-over in cluster i to the right:            */

				tryxover = bestxover
				tryxover[i] = bestxover[i]+1
				tryxover[iopp] = bestxover[iopp]-1
				tryse = gls_polyt_se(tryxover,  ///
						nmiss, rho, tau, order)
				if (tryse <= se) {
					xover = tryxover
					se = tryse
				}
			}
			
			xll=xover[i]
			
		}		

/*  Repeat until no further improvements can be found                 */
		
	} while (se < bestse)
	
	n = nper*nseq
	
	while ((bestse < targse) & (n > minn)) {

/*  The following command - "n" - just displays the current sample    */
/*  size so that the user can keep track of progress (and see that    */
/*  something is happening). It can be deleted if preferred.          */
	
		n
		
		bestnmiss = nmiss

		se = .

/*  STEP 3: Remove the participant with the least impact on standard  */
/*  error:                                                            */
	
		for (i=1; i<=halfnseq; i++) {
			for (j=1; j<=nper; j++) {

				if (bestnmiss[i,j] != .) {
					trynmiss = bestnmiss
					trynmiss[i,j] = .
					trynmiss[nseq+1-i,nper+1-j] = .
					tryse = gls_polyt_se(xover,  //
							trynmiss, rho, tau, order)
					if (tryse < se) {
						nmiss = trynmiss
						se = tryse
					}
				}
			
			}
		}
	
		n = n-2

/*  STEP 4: Look for modifications which reduce the standard error    */
/*  without changing the sample size:                                 */
	
		do {

			xul = (xover[2 :: halfnseq] \ nper/2)
			xll = 0
			bestse = se

			for (i=1; i<=halfnseq; i++) {
		
				iopp = nseq+1-i
				bestxover = xover
				bestnmiss = nmiss

/*  Try shifting the recruitment schedule in cluster i to the left    */
/*  (with wrap-around):                                               */

				trynmiss = bestnmiss
				trynmiss[i,.] =  ///
						(bestnmiss[i,(2 .. nper)],  ///
						bestnmiss[i,1])
				trynmiss[iopp,.] =  ///
						(bestnmiss[iopp, nper],  ///
						bestnmiss[iopp,(1 .. nper-1)])
				tryse = gls_polyt_se(bestxover,  ///
						trynmiss, rho, tau, order)
				if (tryse < se) {
					nmiss = trynmiss
					se = tryse
				}
				if (bestxover[i]-1 >= xll) {

/*  Try shifting the cross-over in cluster i to the left:             */

					tryxover = bestxover
					tryxover[i] = bestxover[i]-1
					tryxover[iopp] = bestxover[iopp]+1
					tryse = gls_polyt_se(tryxover,  ///
							bestnmiss, rho, tau, order)
					if (tryse < se) {
						xover = tryxover
						se = tryse
					}

/*  Try shifting the cross-over and the recruitment schedule in       */
/*  cluster i to the left:                                            */

					tryse = gls_polyt_se(tryxover,  ///
							trynmiss, rho, tau, order)
					if (tryse < se) {
						xover = tryxover
						nmiss = trynmiss
						se = tryse
					}
				}
	
/*  Try shifting the recruitment schedule in cluster i to the right   */
/*  (with wrap-around):                                               */

				trynmiss = bestnmiss
				trynmiss[i,.] =  ///
						(bestnmiss[i, nper],  ///
						bestnmiss[i,(1 .. nper-1)])
				trynmiss[iopp,.] =  ///
						(bestnmiss[iopp,  ///
						(2 .. nper)],  ///
						bestnmiss[iopp,1])
				tryse = gls_polyt_se(bestxover,  ///
						trynmiss, rho, tau, order)
				if (tryse < se) {
					nmiss = trynmiss
					se = tryse
				}
				if (bestxover[i]+1 <= xul[i]) {

/*  Try shifting the cross-over in cluster i to the right:            */

					tryxover = bestxover
					tryxover[i] = bestxover[i]+1
					tryxover[iopp] = bestxover[iopp]-1
					tryse = gls_polyt_se(tryxover,  ///
							bestnmiss, rho, tau, order)
					if (tryse <= se) {
						xover = tryxover
						se = tryse
					}

/*  Try shifting the cross-over and the recruitment schedule in       */
/*  cluster i to the right:                                           */

					tryse = gls_polyt_se(tryxover,  ///
							trynmiss, rho, tau, order)
					if (tryse < se) {
						xover = tryxover
						nmiss = trynmiss
						se = tryse
					}
				}
			
				xll=xover[i]
			
			}		
		
		} while (se < bestse)

/*  STEP 5: Return to STEP 3 (while the standard error is less than  */
/*  the target).                                                     */
	
	}
	
	bestnmiss = nmiss
	bestxover = xover

/*  STEP 1: For the forward search, initialise the design to a        */
/*  staircase design:                                                 */

	xover = floor((0 :: nseq-1) :* (nper/(nseq-1))  ///
			:+ 0.5)
	xover[(halfnseq+1 .. nseq)] =  nper :-  ///
			xover[(halfnseq .. 1)]
	nll = xover :- inc
	nul = xover :+ inc
	nmiss = 1 :/ (((J(nseq, 1, (1 .. nper)) :>  ///
			J(1, nper, nll)) :&  ///
			(J(nseq, 1, (1 .. nper)) :<=  ///
			J(1, nper, nul))))

	se = gls_polyt_se(xover, nmiss, rho, tau, order)

/*  STEP 2: Look for modifications which reduce the standard error    */
/*  without changing the sample size (as above):                      */
	
	do {

		xul = (xover[2 :: halfnseq] \ nper/2)
		xll = 0
		bestse = se

		for (i=1; i<=halfnseq; i++) {
		
			iopp = nseq+1-i
			bestxover1 = xover
			bestnmiss1 = nmiss

			trynmiss = bestnmiss1
			trynmiss[i,.] =  ///
					(bestnmiss1[i,(2 .. nper)],  ///
					bestnmiss1[i,1])
			trynmiss[iopp,.] =  ///
					(bestnmiss1[iopp, nper],  ///
					bestnmiss1[iopp,(1 .. nper-1)])
			tryse = gls_polyt_se(bestxover1,  ///
					trynmiss, rho, tau, order)
			if (tryse < se) {
				nmiss = trynmiss
				se = tryse
			}
			if (bestxover1[i]-1 >= xll) {
				tryxover = bestxover1
				tryxover[i] = bestxover1[i]-1
				tryxover[iopp] = bestxover1[iopp]+1
				tryse = gls_polyt_se(tryxover,  ///
						bestnmiss1, rho, tau, order)
				if (tryse < se) {
					xover = tryxover
					se = tryse
				}
				tryse = gls_polyt_se(tryxover,  ///
						trynmiss, rho, tau, order)
				if (tryse < se) {
					xover = tryxover
					nmiss = trynmiss
					se = tryse
				}
			}	

			trynmiss = bestnmiss1
			trynmiss[i,.] =  ///
					(bestnmiss1[i, nper],  ///
					bestnmiss1[i,(1 .. nper-1)])
			trynmiss[iopp,.] =  ///
					(bestnmiss1[iopp,(2 .. nper)],  ///
					bestnmiss1[iopp,1])
			tryse = gls_polyt_se(bestxover1,  ///
					trynmiss, rho, tau, order)
			if (tryse < se) {
				nmiss = trynmiss
				se = tryse
			}
			if (bestxover1[i]+1 <= xul[i]) {
				tryxover = bestxover1
				tryxover[i] = bestxover1[i]+1
				tryxover[iopp] = bestxover1[iopp]-1
				tryse = gls_polyt_se(tryxover,  ///
						bestnmiss1, rho, tau, order)
				if (tryse <= se) {
					xover = tryxover
					se = tryse
				}
				tryse = gls_polyt_se(tryxover,  ///
						trynmiss, rho, tau, order)
				if (tryse < se) {
					xover = tryxover
					nmiss = trynmiss
					se = tryse
				}
			}
			
			xll=xover[i]
			
		}		
		
	} while (se < bestse)
	
	n1 = minn
	
	while ((bestse > targse) & (n1 < nper*nseq)) {

/*  The following command - "n1" - just displays the current sample   */
/*  size so that the user can keep track of progress (and see that    */
/*  something is happening). It can be deleted if preferred.          */
		
		n1
		
		bestnmiss1 = nmiss

		se = .

/*  STEP 3: Add the participant with the greatest impact on standard  */
/*  error:                                                            */
	
		for (i=1; i<=halfnseq; i++) {
			for (j=1; j<=nper; j++) {
		
				if (bestnmiss1[i,j] == .) {
					trynmiss = bestnmiss1
					trynmiss[i,j] = 1
					trynmiss[nseq+1-i,nper+1-j] = 1
					tryse = gls_polyt_se(xover,  ///
							trynmiss, rho, tau, order)
					if (tryse < se) {
						nmiss = trynmiss
						se = tryse
					}
				}
			
			}
		}
		
		n1 = n1+2

/*  STEP 4: Look for modifications which reduce the standard error    */
/*  without changing the sample size (as above):                      */

		do {

			xul = (xover[2 :: halfnseq] \ nper/2)
			xll = 0
			bestse = se

			for (i=1; i<=halfnseq; i++) {
		
				iopp = nseq+1-i
				bestxover1 = xover
				bestnmiss1 = nmiss

				trynmiss = bestnmiss1
				trynmiss[i,.] =  ///
						(bestnmiss1[i,(2 .. nper)],  ///
						bestnmiss1[i,1])
				trynmiss[iopp,.] =  ///
						(bestnmiss1[iopp, nper],  ///
						bestnmiss1[iopp,(1 .. nper-1)])
				tryse = gls_polyt_se(bestxover1,  ///
						trynmiss, rho, tau, order)
				if (tryse < se) {
					nmiss = trynmiss
					se = tryse
				}
				if (bestxover1[i]-1 >= xll) {
					tryxover = bestxover1
					tryxover[i] = bestxover1[i]-1
					tryxover[iopp] = bestxover1[iopp]+1
					tryse = gls_polyt_se(tryxover,  ///
							bestnmiss1, rho, tau, order)
					if (tryse < se) {
						xover = tryxover
						se = tryse
					}
					tryse = gls_polyt_se(tryxover,  ///
							trynmiss, rho, tau, order)
					if (tryse < se) {
						xover = tryxover
						nmiss = trynmiss
						se = tryse
					}
				}	

				trynmiss = bestnmiss1
				trynmiss[i,.] =  ///
						(bestnmiss1[i, nper],  ///
						bestnmiss1[i,(1 .. nper-1)])
				trynmiss[iopp,.] =  ///
						(bestnmiss1[iopp,  ///
						(2 .. nper)],  ///
						bestnmiss1[iopp,1])
				tryse = gls_polyt_se(bestxover1,  ///
						trynmiss, rho, tau, order)
				if (tryse < se) {
					nmiss = trynmiss
					se = tryse
				}
				if (bestxover1[i]+1 <= xul[i]) {
					tryxover = bestxover1
					tryxover[i] = bestxover1[i]+1
					tryxover[iopp] = bestxover1[iopp]-1
					tryse = gls_polyt_se(tryxover,  ///
							bestnmiss1, rho, tau, order)
					if (tryse <= se) {
						xover = tryxover
						se = tryse
					}
					tryse = gls_polyt_se(tryxover,  ///
							trynmiss, rho, tau, order)
					if (tryse < se) {
						xover = tryxover
						nmiss = trynmiss
						se = tryse
					}
				}
			
				xll=xover[i]
			
			}		
		
		} while (se < bestse)

/*  STEP 5: Return to STEP 3 (while the standard error is greater    */
/*  than the target).                                                */
	
	}

/*  Choose the best of the two designs from the forward and          */
/*  backward searches.                                               */	
	
	if (n1 < n) {
		bestnmiss = nmiss
		bestxover = xover
	}

/*  Return a matrix which describes the design (based on bestnmiss   */
/*  and bestxover): 0 = recruitment under control; 1 = recruitment   */
/*  under intervention; . = no recruitment.                          */
	
	return ((J(1, nper, bestxover) :<  ///
			J(nseq, 1, (1 .. nper))) :* bestnmiss)
	
}
end


mata:
real scalar gls_polyt_se(real colvector xover,  ///
		real matrix nmiss, real scalar rho,  ///
		real scalar tau, real scalar order)

{
	real matrix seq, x, y, x1, y1, ccom, cincom,  ///
			tmat
	real colvector t
	real scalar nper, nseq, i
	
	nper = cols(nmiss)
	nseq = rows(nmiss)
	seq = (J(1, nper, xover) :<  ///
			J(nseq, 1, (1 .. nper))) :* nmiss
	
	y = J(order+2, 0, .)
	x = J(0, order+2, .)

	t = (0 :: (nper-1)) :/ (nper-1) :- 0.5
	tmat = J(1,order,t) :^ J(nper,1,(1 .. order))
	ccom = rho :* (J(nper, nper, tau) :^ abs  ///
		(J(nper, 1, t') - J(1, nper, t)))  ///
		+ (1-rho) :* I(nper)
		
	for (i=1; i<=nseq; i++) {		
		x1 = select((seq[i,.]' , J(nper, 1, 1),  ///
				tmat), seq[i,.]' :!= .)
		cincom = select(  ///
			select(ccom, seq[i,.]' :!= .),  ///
			seq[i,.] :!= .)
		y1 = x1' * invsym(cincom)
		x = x \ x1
		y = y , y1
	}
	
	return(sqrt(invsym(y*x)[1,1]))
}
end


/*  Here's how to run a toy example with 10 clusters (halfnseq=5),   */
/*  20 participants per cluster, rho = 0.05 and tau = 0.2, effect    */
/*  size 0.6, with time modelled as a sixth order polynomial. Note   */
/*  that it may take considerably more time to run examples with     */
/*  larger numbers of clusters or participants per cluster. The      */
/*  design is displayed in transposed form, with clusters in         */
/*  columns and participants in rows.                                */ 

mata:
a = opti_ctsw(5,20,0.05,0.2,0.6,6)
a'
end

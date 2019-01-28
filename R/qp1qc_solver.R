# To do:
# Fixed simDiagBasis function, but need to go through and remove language like "basis" since it's not really appropriate, the matrix is not orthonomal. 
# !! = important note to keep track of, or note for future changes.



#' binary search with arbitrary resolution
#'
#' The \code{\link{binsearch}} function in the \code{gtools} package searches only over integers. This is a wrapper for \code{\link{binsearch}} that also allows searching over a grid of points with distance \code{tol} between them. A target must also be entered (see \code{\link{binsearch}}).
#' @param fun a monotonic function over which we search (passed to \code{\link{binsearch}})
#' @param  tol resolution of the grid over which to search for \code{target} (see \code{\link{binsearch}})
#' @param  range a range over which to search for the input to \code{fun}
#' @param  ... passed to \code{\link{binsearch}}
#' @export
#' @import gtools
#' @seealso \code{\link{binsearch}}
#' @examples
#'  # best solution here is at x0,
#'  # which we can find with increasing precision
#'  x0 <- 10.241
#'  binsearchtol( function(x) x-x0, target=0, range=c(0,2*x0) , tol=1.00)
#'  binsearchtol( function(x) x-x0, target=0, range=c(0,2*x0) , tol=0.10)
#'  binsearchtol( function(x) x-x0, target=0, range=c(0,2*x0) , tol=0.05)
#'
binsearchtol <- function(fun, tol, range, ...){
  funTol <- function(x) fun(x*tol)

  soln <- binsearch(fun=funTol, range=range/tol, ...)
  soln$where <- soln$where*tol

  soln
}



#' Simultaneously diagonalize two symmetric matrices (1 being positive definite)
#' 
#' @param M1 a positive definite symmetric matrix
#' @param M2 a symmetric matrix
#' @param tol used for error checks
#' @param return_diags should Q^(-1); Q'M1Q; and Q'BQ be returned.
#' Resulting matrix Q gives  Q'M1Q=I (identity); Q'BQ = a diagonal matrix.
#' Q is *not* orthonormal or even symmetric, but it is invertible.
#' #########
# #TEST
# p <- 5
# M1 <- diag(p) + crossprod(matrix(rnorm(p^2),p,p))
# M1[1,] <- M1[,1] <- 0
# M2 <- diag(p)+2
# sdb <- simDiagBasis(M1=M1,M2=M2,tol=10^-10,return_diags=TRUE)
# Q <- sdb$diagonalizer
# QM1Q <- t(Q) %*% M1 %*% Q
# QM2Q <- t(Q) %*% M2 %*% Q
# range(QM1Q -diag(sdb$M1diag))
# range(QM2Q -diag(sdb$M2diag))
# range(Q %*% sdb$inverse_diagonalizer - diag(p))
# range(sdb$inverse_diagonalizer %*% Q - diag(p))
# #########
simDiagBasis <- function(M1,M2, eigenM1=NULL, eigenM2=NULL,tol, return_diags=FALSE){
	
	if(!isSymmetric(M1)) stop('M1 must be symmetric')
	if(!isSymmetric(M2)) stop('M2 must be symmetric')


	############
	#Central computation - part 1
	if(is.null(eigenM1)) eigenM1 <- eigen(M1)

	############
	#Error checks on inputs
	if(any(eigenM1$value<=0)){ #M1 not pos def.
		if(is.null(eigenM2)) eigenM2 <- eigen(M2)
		if(any(eigenM2$value<=0)){ #M2 not pos def.
			stop('Neither M1 nor M2 is positive definite')
		}
		if(all(eigenM2$value>0)){ #M2 not pos def.
			warning('M1 is not positive definite, but M2 is. Switching roles of M1 & M2')
			out_pre <- simDiagBasis(M1=M2,M2=M1,eigenM1=eigenM2,tol=tol,return_diags=return_diags)
			if(!return_diags) return(out_pre)
			out <- out_pre
			out$M1diag <- out_pre$M2diag
			out$M2diag <- out_pre$M1diag
			return(out)
		}
	}


	############
	#Central computation - part 2
	sqrtInvM1 <- eigenM1$vectors %*% (diag(1/sqrt(eigenM1$values)))
	Z <- t(sqrtInvM1) %*% M2 %*% sqrtInvM1
	# if(any( abs(t(Z)-Z) > tol)) stop('Symmetry error') #!! Computational instability!! This shouldn't be needed. This is flagging errors where it shouldn't (for smaller kernel size, errors are more common?)
	Z <- (Z + t(Z) )/ 2 #force matrices to by symmetric
	Q <- sqrtInvM1 %*% eigen(Z)$vectors # recall, real symmetric matrices are diagonalizable by orthogonal matrices (wikipedia)

	# prove that this is invertible
	# let V = eigen(Z)$vectors, then V'V=VV'=I. Let M^-1/2 = WD^-1/2,
	# Then Q_inverse =V' M1^(1/2)
	# Q [V' M1^(1/2)] = Q [V' D^(1/2)W'] = WD^(-1/2) V [V' D^(1/2)W'] = I
	# [V' M1^(1/2)] Q = [V' D^(1/2)W'] Q = [V'D^(1/2 )W'] WD^(-1/2) V = I
	Q_inv <- t(eigen(Z)$vectors) %*% diag(sqrt(eigenM1$values)) %*% t(eigenM1$vectors)
	inv_err <- max(abs(Q_inv %*% Q - diag(dim(M1)[1])))
	if( inv_err > tol * min( sqrt(abs(eigenM1$values))) ){
		warning(paste('possible inverting or machine error of order', signif(inv_err,4)))	
	} 

	if(!return_diags) return(Q)
	return(list(
		diagonalizer = Q,
		inverse_diagonalizer = Q_inv,
		M1diag = rep(1,dim(M1)[1]),
		M2diag = diag(t(Q) %*% M2 %*% Q)
	))
		# https://math.stackexchange.com/questions/1079627/simultaneously-diagonalization-of-two-matrices
			# Z = UDU'
			# Q'AQ   = U'[M1^(-1/2)' M1 M1^(-1/2)] U = U'[I]U = I
			# Q'BQ   = U'[M1^(-1/2)' M2 M1^(-1/2)] U = U'[Z]U = U'[UDU']U = D

	############
	#Workchecks
	# round(t(sqrtInvM1) %*% M1 %*% sqrtInvM1,12)
	# round(t(Q) %*% M1 %*% Q, 10)
	# round(t(Q) %*% M2 %*% Q, 10)
}

pseudo_invert_diagonal_matrix<- function(M){
	diag_M<-diag(as.matrix(M))
	M_pseudo_inv_diagonal <- 1/diag_M
	M_pseudo_inv_diagonal[diag_M==0] <- 0
	M_pseudo_inv <- diag(M_pseudo_inv_diagonal)
	M_pseudo_inv
}

is.diagonal <- function(M,tol=10^-8){
	if(!is.matrix(M)) stop('M must be a matrix')
	if(any(abs(M-t(M))>tol)) warning(paste('M must be symmetric, possible machine error of order', max(abs(M-t(M)) )))

	all( abs( M - drop_off_diagonals(M) ) < tol)
}

to_diag_mat <- function(vec){
	#input is a vector to put on the diagonal of a matrix
	if(length(vec)==1) return(as.matrix(vec))
	return(diag(vec))
}
drop_off_diagonals <- function(mat){
	to_diag_mat(diag(as.matrix(mat)))
}




#   ____  _____    ____            _
#  / __ \|  __ \  |  _ \          (_)
# | |  | | |__) | | |_) | __ _ ___ _  ___ ___
# | |  | |  ___/  |  _ < / _` / __| |/ __/ __|
# | |__| | |      | |_) | (_| \__ \ | (__\__ \
#  \___\_\_|      |____/ \__,_|___/_|\___|___/



# TEST # set_QP_unconstrained_eq_0(M=diag(0:3),v=0:3,k=-2,tol=10^-9)
set_QP_unconstrained_eq_0 <- function(M,v,k,tol){
	# find x s.t. x'Mx + v'x + k = 0,
	# or show that it is not possible
	# M should be diagonal
	# We find the saddle point
	# Then determine the direction we need to go (up or down) to get to zero
	# solve the resulting polynomial
	# (See notes in pdf)

	feasible <- FALSE
	soln <- value <- NA

	############################
	#Initialize a reference point where the derivative  = 0.
	diag_M <- diag(as.matrix(M))
	M_pseudo_inv <- pseudo_invert_diagonal_matrix(M)
	x0 <- -(1/2)*M_pseudo_inv %*% v 
	eval_x <- function(x){
		sum(x^2 * diag_M) + sum(x*v) + k
	}
	eval_x0 <- eval_x(x0)
	feasible <- eval_x0==0
	if(eval_x0==0){
		return(list(feasible=TRUE, soln=x0, value=eval_x0))
	}


	############################
	# Find an index of x0 along which we can move such that a new x0 solution evaluates to zero (with eval_x)
	move_ind <- NA

	v_candidates <- (diag_M==0)&(v!=0)
	if(sum(v_candidates)>0){
		move_ind <- which(v==max(abs(v[v_candidates])))[1]
	}else if(eval_x0 < 0 & max(diag_M) > 0){
		move_ind <- which(diag_M==max(diag_M))[1] 
	}else if(eval_x0 > 0 & min(diag_M) < 0){
		move_ind <- which(diag_M==min(diag_M))[1] 
	}

	if(!is.na(move_ind)){
		soln <- 
		soln_base <- x0
		soln_base[move_ind] <- 0
		coeffs <- c(eval_x(soln_base),
			v[move_ind],
			diag_M[move_ind])
		move_complex <- polyroot(coeffs)[1]
		if(Im(move_complex) > tol) stop('error in root finding')
		soln[move_ind] <- Re(move_complex)
		if(abs(eval_x(soln))>tol) stop('calculation error')
		feasible <- TRUE
	}	

	return(list(
		feasible=feasible,
		soln=soln,
		value=eval_x(soln)
	))
}






min_QP_unconstrained <- function(M,v,tol){
	# minimize x'Mx+v'x
	# if M can be invertible (after zero elements are discarded), then solution is when
	# 2Mx + v=0; (-1/2)M^-1 v = x
	BIG <- (1/tol)^4 #To avoid Inf*0 issues

	if(!all(is.finite(c(M,v)))) stop('M & v must be finite')
	if(!is.diagonal(M)) stop('M should be diagonal')
	if(length(v)!=dim(M)[1]) stop('dimension of M & v must match')

	diag_M <- diag(M)
	
	zero_directions <- 		 (diag_M==0) & c(v==0)
	moves_up_to_pos_Inf<- 	((diag_M==0) & c(v>0)) | (diag_M>0) #as we increase these indeces, do we approach +Inf?
	moves_up_to_neg_Inf <- 	((diag_M==0) & c(v<0)) | (diag_M<0) #as we increase these indeces, do we approach -Inf?
	moves_down_to_pos_Inf<- ((diag_M==0) & c(v<0)) | (diag_M>0)
	moves_down_to_neg_Inf<- ((diag_M==0) & c(v>0)) | (diag_M<0)
	if( any(
		moves_down_to_pos_Inf + moves_down_to_neg_Inf + zero_directions !=1 |
		moves_up_to_pos_Inf   + moves_up_to_neg_Inf   + zero_directions !=1 
		)){
		stop('direction error')#work check
	}
	finite_soln <- !any(moves_up_to_neg_Inf|moves_down_to_neg_Inf)

	if(finite_soln){ #only true if all zero elements of M correspond to zero elements of v
		M_pseudo_inv <- pseudo_invert_diagonal_matrix(M)
		soln <- -(1/2)*M_pseudo_inv %*% v 
	}else{
		soln <- rep(0,length(v))
		#need to use "big number" BIG to avoid multiplying by zero
		soln[moves_up_to_neg_Inf] <- BIG 
		soln[moves_down_to_neg_Inf] <- - BIG
	}
	value <- t(soln) %*% M %*% soln + crossprod(v,soln)

	return(list(
		zero_directions=zero_directions,
		moves_up_to_pos_Inf=moves_up_to_pos_Inf,
		moves_up_to_neg_Inf=moves_up_to_neg_Inf,
		moves_down_to_pos_Inf=moves_down_to_pos_Inf,
		moves_down_to_neg_Inf=moves_down_to_neg_Inf,
		finite_soln=finite_soln,
		value=value,
		soln=soln
		))
}










#   ____    _____    __    ____     _____
#  / __ \  |  __ \  /_ |  / __ \   / ____|
# | |  | | | |__) |  | | | |  | | | |
# | |  | | |  ___/   | | | |  | | | |
# | |__| | | |       | | | |__| | | |____
#  \___\_\ |_|       |_|  \___\_\  \_____|


#' Solve (non-convex) quadratic program with 1 quadratic constraint
#' 
#' Solves a possibly non-convex quadratic program with 1 quadratic constraint. Either \code{A_mat} or \code{B_mat} must be positive definite, but not necessarily both (see Details, below).
#' 
#' Solves a minimization problem of the form:
#' 
#' \deqn{ 	min_{x} x^T A_mat x + a_vec^T x }
#' \deqn{ s.t. x^T B_mat x + b_vec^T x + k \leq 0,}
#' 
#' where either \code{A_mat} or \code{B_mat} must be positive definite, but not necessarily both.
#' 
#' @param A_mat see details below
#' @param a_vec see details below
#' @param B_mat see details below
#' @param b_vec see details below
#' @param k see details below
#' @param tol a calculation tolerance variable used at several points in the algorithm.
#' @param eigen_A_mat (optional) the precalculated result \code{eigen(A_mat)}, where A_mat is defined in the optimization problem below.
#' @param eigen_B_mat (optional) the precalculated result \code{eigen(B_mat)}, where B_mat is defined in the optimization problem below.
#' @param verbose show progress from calculation
#' @import quadprog
#' @export
solve_QP1QC <- function(A_mat, a_vec, B_mat, b_vec, k, tol, eigen_B_mat=NULL, eigen_A_mat=NULL, verbose= TRUE){

	if(tol<0){tol <- abs(tol); warning('tol must be positive; switching sign')}
	if(tol>.01){tol <- 0.01; warning('tol must be <0.01; changing value of tol')}
	
	if(is.null(eigen_B_mat)) eigen_B_mat <- eigen(B_mat)
	if(any(eigen_B_mat$value<=0)){
		#B_mat not pos def.

		if(is.null(eigen_A_mat)) eigen_A_mat <- eigen(A_mat)
		if(any(eigen_A_mat$value<=0)){
			#A_mat not pos def.
			stop('Neither B_mat nor A_mat is positive definite')
		}

	}

	####### Diagonalize
	suppressWarnings({
		sdb <-  simDiagBasis(
		M1=B_mat, M2=A_mat, eigenM1=eigen_B_mat, eigenM2= eigen_A_mat, tol=tol, return_diags=TRUE)
	})
	Q  <- sdb$diagonalizer
	B <- crossprod(Q,B_mat)%*%Q
	b <- crossprod(Q,b_vec)
	A <- crossprod(Q,A_mat)%*%Q
	a <- crossprod(Q,a_vec)

	# Finish diagonalizing by dropping machine error	
	# if(!is.diagonal(A)) stop('A should be diagonal') #!! throwing bugs due to machine error, possibly for low rank matrices?
	# if(!is.diagonal(B)) stop('B should be diagonal')#!! throwing bugs due to machine error, possibly for low rank matrices?

	A <- drop_off_diagonals(A)
	B <- drop_off_diagonals(B)
	Bd <- diag(as.matrix(B))
	Ad <- diag(as.matrix(A))

	# Define helper functions on diagonal space
	calc_diag_Lagr <- function(x,nu){
		#output as vector
		sum(x^2 * diag(as.matrix(A + nu * B))) + c(crossprod(x,(a + nu* b)))
	}
	calc_diag_obj <- function(x){
		#output as vector
		sum(x^2 * diag(as.matrix(A))) + c(crossprod(x,a))
	}
	calc_diag_constraint <- function(x){
		#output as vector
		sum(x^2 * diag(as.matrix(B))) + c(crossprod(x,b)) + k
	}
	return_on_original_space <- function(x){
		list(
			soln = Q %*% x,
			constraint = calc_diag_constraint(x),
			objective = calc_diag_obj(x)
			)
	}


	#  _            _   _
	# | |          | | (_)
	# | |_ ___  ___| |_ _ _ __   __ _   _ __  _   _
	# | __/ _ \/ __| __| | '_ \ / _` | | '_ \| | | |
	# | ||  __/\__ \ |_| | | | | (_| | | | | | |_| |
	#  \__\___||___/\__|_|_| |_|\__, | |_| |_|\__,_|
	#                            __/ |
	#                           |___/

	test_nu <- function(nu, tol){
		#return: is nu too low, too high, or optimal
		if(any(Ad+nu*Bd < -tol)) warning(paste('invalid nu, should not have been submitted. Possible machine error of magnitude',abs(min(Ad+nu*Bd))))
		if(all(Ad+nu*Bd > 0)){
			return(test_nu_pd(nu))
		}else{
			return(test_nu_psd(nu,tol))
		}
	}

	# Returns high, low, optimal or non-optimal. Beacuse it's used for a binary search, it can't just say "non-optimal."
	test_nu_pd <- function(nu){

		soln <- -(1/2)*((Ad+nu*Bd)^-1)*(a+nu*b)
		constraint_value <- calc_diag_constraint(soln)
		
		if(constraint_value>0) out_type <- 'low' #contraint value is monotone decreasing in nu
		if(constraint_value<0) out_type <- 'high'
		if(constraint_value==0){
			out_type <- 'optimal'
		}

		return(list(
			type=out_type,
			soln=soln
		))
	}


	test_nu_psd <- function(nu,tol){

		diag_A_nu_B <- Ad + nu * Bd
		if(any(diag_A_nu_B< -tol)) warning(paste('possible error in pd, or machine error of order',abs(min(diag_A_nu_B))))
		diag_A_nu_B[diag_A_nu_B<0] <- 0
		mat_A_nu_B <- to_diag_mat(diag_A_nu_B)
		if(all(diag_A_nu_B>0)) stop('error in psd')

		I_nu <- diag_A_nu_B > 0
		N_nu <- diag_A_nu_B == 0

		##### Infer context in which function was called:
			#if this value of nu turns out to be non-optimal, infer whether we were testing nu_min or nu_max. If nu_min is non-optimal, then that means nu_min is lower than the optimum nu. Likewise, if we are testing nu_max and it proves non-optimal, than nu_max is too high.
		non_opt_value <- NA
		if(any(Bd[N_nu]<0) & any(Bd[N_nu]>0)) return('optimal') 
			#min=max!! need to say what the actual returned value is though. (!!) Maybe handle this separately
		if(any(Bd[N_nu]<0)) non_opt_value <- 'high' #we've inferred that we were testing nu_max.
		if(any(Bd[N_nu]>0)) non_opt_value <- 'low' #we've inferred that we were testing nu_min.
		if(is.na(non_opt_value)) stop('PSD test should not have been called') # !! Relevant for using this for initial checks of feasibility?
		#####


		#### Optimality check 1 (necessary) 
		if(max(abs(a + nu * b)[N_nu])>0){
			return(list(
				type=non_opt_value,
				soln=NA))
		}

		#### Optimality check 2 (necessity implied by 1st check not holding) 
		# First solve PD problem over I_nu
		A_nu_B_pseudo_inv <- pseudo_invert_diagonal_matrix(mat_A_nu_B)
		x_I <- -(1/2)*A_nu_B_pseudo_inv %*% (a + nu * b)
		if(any(x_I[N_nu]!=0)) stop('indexing error')

		free_constr_opt <- set_QP_unconstrained_eq_0(
			M = as.matrix(B[N_nu, N_nu]),
			v = b[N_nu],
			k = calc_diag_constraint(x_I),
			tol=tol
			)
		if(free_constr_opt$feasible){
			soln <- x_I 
			soln[N_nu] <- free_constr_opt$soln
			return(list(
				type = 'optimal',
				soln = soln
			))
		}else{
			return(list(
				type=non_opt_value,
				soln=NA
				))
		}
	}


	# Test constraint feasibility
	constr_prob <- min_QP_unconstrained(M=B, v=b, tol=tol)
	if(constr_prob$finite_soln){
		if(constr_prob$value > -k) stop('Constraint is not feasible')
	}
	if(constr_prob$value == -k) warning('Constraint may be too strong -- no solutions exist that strictly satisfy the constraint.')





	######  Step 1 ######
	# Check if unconstrained solution is feasible

	# Notes
		# Could we simplify by just doing test_nu(0)? No, later functions assume that constraint is met with equality.

	u_prob <- min_QP_unconstrained(M=A, v= a, tol=tol) # unconstrained problem
	min_constr_over_restricted_directions <- function(directions){
		
		if(length(directions)==0) return(u_soln)

		search_over_free_elements <- min_QP_unconstrained(
			M = to_diag_mat(Bd[directions]),
			v = b[directions],
			tol= tol
		)
		u_soln <- u_prob$soln
		u_soln[directions] <- search_over_free_elements$soln
		u_soln
	}

	u_soln <- u_prob$soln
	if(u_prob$finite_soln){
		if(any(u_prob$zero_directions)){ # we have a finite nonunique solution
			u_soln <- min_constr_over_restricted_directions(u_prob$zero_directions)
		}
		#otherwise we have a UNQIUE finite solution (u_soln), assigned above
	}
	# (Code commented out below) Don't bother to check the case when solution is inifite. Since we need A or B to PD in order to simultaneously diagonalize them (right now), the solution to the unconstrained problem has to either be finite, or lead to a non-feasible constraint value. Simply using u_soln should achieve this, unless BIG is too small, which is highly unlikely.

				# if(!u_prob$finite_soln){
				# 	#We have an infinite solution, so 1 dimension of soln must be Inf.
				# 	#Check each of these dimensions to see if this is possible while
				# 	#Meeting constraint
				# 	u_soln <- u_prob$soln
				# 	search_set <- 	u_prob$moves_up_to_neg_Inf | 
				# 					u_prob$moves_down_to_neg_Inf | 
				# 					u_prob$zero_directions

				# 	if(length(search_set)>1){ for(i in 1:p){ #If solution is not unique, then for each element in the search set...
				# 		BIG <- (1/tol)^4 #To avoid Inf*0 issues

				# 		if(u_prob$moves_up_to_neg_Inf[i]){
				# 			#Fix index i at +infinity, see if we can satisfy constraint.
				# 			if(constr_prob$moves_up_to_pos[i]) next #can't satisfy
				# 			search_set_up_i <- search_set 
				# 			search_set_up_i[i] <- FALSE #don't search over index i
				# 			u_soln_i <- min_constr_over_restricted_directions(search_set_up_i)
				# 			u_soln_i[i] <- BIG #fix at +Inf
				# 			if(c(calc_diag_constraint(u_soln_i) <= 0)){
				# 				u_soln <- u_soln_i
				# 				break
				# 			}
				# 		}

				# 		if(u_prob$moves_down_to_neg_Inf[i]){
				# 			#Fix index i at -Inf, see if we can satisfy constraint
				# 			if(constr_prob$moves_down_to_pos[i]) next #can't satisfy
				# 			search_set_down_i <- search_set
				# 			search_set_down_i[i] <- FALSE
				# 			u_soln_i <- min_constr_over_restricted_directions(search_set_down_i)
				# 			u_soln_i[i] <- -BIG
				# 			if(c(calc_diag_constraint(u_soln_i) <= 0)){
				# 				u_soln <- u_soln_i
				# 				break
				# 			}
				# 		}

				# 	}}

				# }


	if(c(calc_diag_constraint(u_soln) <= 0)){
		if(verbose) cat('\nUnconstrained solution also satisfies the constraint\n')
		return(return_on_original_space(u_soln))
	}


	###### Step 2 #####
	#Get upper/lower bounds on nu

	nu_opt <- x_opt <- NA # (yet) unknown optimal nu value

	nu_to_check <- c()
	nu_max <- Inf
	nu_min <- -Inf
	if(any(Bd>0)){
		nu_min <- max( (-Ad/Bd) [Bd>0] )
		nu_to_check <- c(nu_to_check,nu_min)
	}
	if(any(Bd<0)){ #non infinite check!! only relevant if not infinite Bd
		nu_max <- min( (-Ad/Bd) [Bd<0] )
		nu_to_check <- c(nu_to_check,nu_max)
	}
	if(length(nu_to_check)==0){
		if(any(Bd!=0)) stop('Error in Bd check')
		if(any(B_mat!=0)) stop('Error in Bd check')
		if(any(Ad<=0)) stop('Error in PD check')

		warning('Quadratic constraint not active')# (All diagonal elements of B are zero; a linear constraint is sufficient)
		
		qp_soln <- solve.QP(Dmat = 2*A_mat, dvec= -a_vec,
			Amat = matrix(-b_vec,ncol=1), bvec = k)$solution
			# solve.QP formulates their problem with different constants than us.
			# Requires PD 

		return( return_on_original_space( solve(Q)%*%qp_soln ) ) 
			# Need to cancel out multiplication by Q in return_on_original_space function
	}


	# Check nu_min, nu_max
	for(i in 1:length(nu_to_check)){
		test_nu_check <- test_nu(nu_to_check[i],tol=tol)
		if(test_nu_check$type=='optimal'){
			nu_opt <- nu_to_check[i]
		}
	}


	# Find endpoints for binary search.
	# print(c(nu_min,nu_max))
	if(is.infinite(nu_max)){
		nu_max <- abs(nu_min) + 10 #arbitrary number just to make it >1
		test_nu_max <- test_nu(nu_max, tol=tol)
		counter <- 0
		while(abs(nu_max) < 1/tol & test_nu_max$type=='low'){
			nu_max <- abs(nu_max) * 10
			test_nu_max <- test_nu(nu_max, tol=tol)

			# Extra safety/error check
			counter <- counter + 1
			if(counter > 1000) stop('While loop broken')
		}
		if(test_nu_max$type=='low'){nu_opt <- nu_max; warning('outer limit reached')}
	}
	if(is.infinite(nu_min)){
		nu_min <- -abs(nu_max) - 10 #arbitrary number just to make it < -1
		test_nu_min <- test_nu(nu_min, tol=tol)
		counter <- 0
		while(abs(nu_min) < 1/tol & test_nu_min$type=='high'){
			nu_min <- -abs(nu_min) * 10
			test_nu_min <- test_nu(nu_min, tol=tol)

			# Extra safety/error check
			counter <- counter + 1
			if(counter > 1000) stop('While loop broken')
		}	
		if(test_nu_min$type=='high'){nu_opt <- nu_min; warning('outer limit reached')}
	}
	# print(c(nu_min,nu_max))


	##### Step 3 #####
	# Binary Search (if nu_min or nu_max are not optimal)

	if(is.na(nu_opt)){
		bin_serach_fun <- function(nu){
			tested_type <- test_nu(nu, tol=tol)$type
			if(tested_type=='high') return(1)
			if(tested_type=='low') return(-1)
			if(tested_type=='optimal') return(0)
		}
		nu_opt <- binsearchtol(fun=bin_serach_fun, tol=tol, range=c(nu_min, nu_max), target=0)$where[1]
	}


	x_opt <- test_nu(nu_opt,tol=tol)$soln #either from binary search, or from nu_min or nu_max
	return(return_on_original_space(x_opt))

}


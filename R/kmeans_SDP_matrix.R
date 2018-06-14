

#' K-means SDP Matrices.
#' 
#' Constructs the SDP constraint matrices so they can be passed to a black-box solver.
#' 
#' @param D \eqn{d \times d} numeric matrix.
#' @param K number of clusters.
#' 
#' @export
gforce.kmeans_SDP_matrix <- function(D,K){
    
    # Initialization
    d <- ncol(D)
    # block_desc is a list with two elements - one of which is a vector of
    # variable types and the other a list of block sizes
    block_desc <- NULL
    block_desc$type <- c('s','l')
    block_desc$size <- c(d,d*d)

    # list of constraints and RHS
    A <- list()
    b <- c()


    # Objective Function
    C <- list()
    C[[1]] <- D
    C[[2]] <- rep(0,d*d)


    # Constraint: X1 = 1
    for(a in 1:d) {
        R_a_blocks <- list()
        R_a_blocks[[1]] <- R_a(d,a)
        R_a_blocks[[2]] <- rep(0,d*d)

        b[a] <- 2
        A[[a]] <- R_a_blocks
    }


    # Constraint: X >= 0, entrywise
    constraint_idx <- d
    for(a in 1:d){
        for(b in 1:a) {
            I_ab_blocks <- list()
            I_ab_blocks[[1]] <- I_ab(d,a,b)

            slack_mat <- rep(0,d*d)
            slack_mat[(a-1)*d + b] <- 0.5
            slack_mat[(b-1)*d + a] <- 0.5
            I_ab_blocks[[2]] <- slack_mat

            constraint_idx <- constraint_idx + 1 
            A[[constraint_idx]] <- I_ab_blocks
            b[constraint_idx] <- 0
        }
    }


    # Constraint: Tr(X) = K
    constraint_idx <- constraint_idx + 1
    I_d_blocks <- list()
    I_d_blocks[[1]] <- diag(d)
    I_d_blocks[[2]] <- rep(0,d*d)

    b[constraint_idx] <- K
    A[[constraint_idx]] <- I_d_blocks

    # Return Value
    res <- NULL
    res$block_desc <- block_desc
    res$A <- A
    res$b <- b
    res$C <- C

    return(res)
}
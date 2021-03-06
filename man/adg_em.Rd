% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screen.R
\name{adg_em}
\alias{adg_em}
\title{Estimate Q, screening latent attributes (Alternating Gibbs)}
\usage{
adg_em(
  X,
  Z_ini,
  Q_ini,
  max_iter,
  err_prob,
  must_maxiter = 0,
  D_mat = NULL,
  X1 = NULL,
  model = "DINA"
)
}
\arguments{
\item{X}{N by J2 binary data matrix - level 2}

\item{Z_ini}{N by K initial latent attributes}

\item{Q_ini}{J by K inital Q matrix}

\item{max_iter}{maximum iterations (e.g., 50)}

\item{err_prob}{noise level}

\item{must_maxiter}{1 to force maxiter; default is \code{0}}

\item{D_mat}{J1 by J2 binary matrix to indicate children in two-level trees.
\code{D_mat} is the \code{J1 * J2} binary adjacency matrix specifying how the
        trees are grown in the second layer, i.e., which second-level
        responses are "children" of each first-level response. Default is \code{NULL}}

\item{X1}{N by J1 binary data matrix - level 1; default is \code{NULL}}

\item{model}{"DINA" (default) or "DINO"}
}
\value{
\itemize{
\item Z_est Estimated latent attributes for all people
\item Z_candi candidate latent attribute patterns (unique)
\item Q_arr a list of Q matrices obtained from the algorithm
\item c J dimensional 1-slipping parameter
\item g J dimensional guessing parameter
}
}
\description{
This function implements the Alternating Direction Gibbs EM (ADG-EM)
algorithm in the scenario of responses observed over many taxonomies
(trees)
}

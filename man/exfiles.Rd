\name{ghap.exfiles}
\alias{ghap.exfiles}

\title{
 Example files
}
\description{
 This function retrieves the list of example files available.
}
\usage{
 ghap.exfiles()
}

\details{
This function requires internet connection. It returns a data table containing the list of example files in our github repository (\url{https://github.com/ytutsunomiya/GHap}). To get any of those files, please use \code{\link{ghap.makefile}}.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{
# # See list of example files
# exlist <- ghap.exfiles()
# View(exlist)
}

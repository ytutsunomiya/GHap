\name{ghap.manhattan}
\alias{ghap.manhattan}

\title{
 Manhattan plot
}
\description{
  Generate a Manhattan plot from a dataframe.
}
\usage{
 ghap.manhattan(data, chr, bp, y, colors = NULL, type = "p", pch = 20,
                cex = 1, lwd = 1, ylim = NULL, ylab = "", xlab = "", main = "",
                backcolor = "#F5EFE780", chr.ang = 0, hlines = NULL,
                hcolors = NULL, hlty = 1, hlwd = 1)
}

\arguments{
  \item{data}{
   A data.frame containing the data to be plotted.
  }
  \item{chr}{
   A character value with the name of the column containing chromosome labels.
  }
  \item{bp}{
   A character value with the name of the column containing base pair positions.
  }
  \item{y}{
   A character value with the name of the column containing the variable to be plotted in the y axis.
  }
  \item{colors}{
   A character value containing colors to be used for chromosomes.
  }
  \item{type}{
   What type of plot should be drawn (default = "p"). See \code{\link{plot}}.
  }
  \item{pch}{
   Either an integer specifying a symbol or a single character to be used as the default in plotting points (default = 20). See \code{\link{points}} for possible values and their interpretation. 
  }
  \item{cex}{
   A numeric value for the relative point size (default = 1).
  }
  \item{lwd}{
   A numeric value for the line width (default = 1). Only meaningful for type = "l".
  }
  \item{ylim}{
   A numeric vector of size 2 containing the lower and upper limits of the y-axis.
  }
  \item{ylab}{
   A chracter value for the y-axis label.
  }
  \item{xlab}{
   A chracter value for the x-axis label.
  }
  \item{main}{
   A chracter value for the plot title.
  }
  \item{backcolor}{
   The background color.
  }
  \item{chr.ang}{
   A numeric value representing the rotation of chromosome labels in degrees (default = 0).
  }
  \item{hlines}{
   A numeric vector containing y-axis positions for horizontal lines.
  }
  \item{hcolors}{
   A character vector containing colors for the horizontal lines.
  }
  \item{hlty}{
   A numeric vector containing types for horizontal lines.
  }
  \item{hlwd}{
   A numeric vector for the relative width of vertical lines.
  }
}

\details{
 This function takes a dataframe of genomic positions and generates a Manhattan plot. The chromosome column must be a vector of factors (the order of the chromosomes will be displayed according to the order of the factor levels).
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# ### RUN ###
# 
# # Generate some data
# set.seed(1988)
# genome <- c(1:22,"X")
# chr <- rep(genome, each = 1000)
# bp <- rep(1:1000, times = length(genome))
# y <- abs(rnorm(n = length(bp)))
# y <- 1 - (pnorm(y) - pnorm(-y))
# y <- -log10(y)
# y[10300:10350] <- sort(runif(n = 51, min = 0, max = 10))
# mydata <- data.frame(chr,bp,y)
# mydata$chr <- factor(x = mydata$chr, levels = genome, labels = genome)
# 
# # Minimal plot
# ghap.manhattan(data = mydata, chr = "chr", bp = "bp", y = "y")
# 
# # Customization example
# ghap.manhattan(data = mydata, chr = "chr", bp = "bp", y = "y",
#                main = "Genome-wide association analysis", ylab = "-log10(p)",
#                xlab = "Chromosome", chr.ang = 90, backcolor = "white", type = "l",
#                hlines = c(6,8), hlty = c(2,3), hcolors = c(1,2), hlwd = c(1,2))

}

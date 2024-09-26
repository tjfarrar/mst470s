library(berryFunctions)

linepointgraph <- function(fX, xset) {
  par(mar = c(4, 5, 0, 0.1))
  plot(xset, fX(xset), type = "p", xlab = "",
       ylab = "", pch = 20, ylim = c(0, 1), axes = F, cex = 1.5)
  axis(side = 1, at = xset)
  mtext(expression(x), side = 1, las = 1, line = 2.5)
  mtext(expression(f[X](x)), side = 2, las = 1, line = 3)
  axis(side = 2, at = seq(0, 1, 0.2), las = 1)
  invisible(lapply(xset, function(x) lines(c(x, x), c(0, fX(x)))))
}

probdensityplot <- function(fX, xsupport, val, ineq) {
  xset <- seq(xsupport[1], xsupport[2], 0.001)
  par(mar = c(4, 5, 0, 0.1))
  plot(xset, fX(xset), type = "l", xlab = "",
       ylab = "", ylim = c(0, ceiling(max(fX(xset)) * 10) / 10) * 1.05, 
       axes = F, cex = 1.5, xaxs = "i", yaxs = "i")
  axis(side = 1, at = seq(min(xset), max(xset), 1))
  mtext(expression(x), side = 1, las = 1, line = 2.5)
  mtext(expression(f[X](x)), side = 2, las = 1, line = 3)
  axis(side = 2, at = seq(0, 1, 0.1), las = 1)
  xsubset <- xset[which(do.call(ineq, args = list(xset, val)))]

  polygon(c(xsubset, rev(xsubset)), c(rep(0, length(xsubset)), fX(rev(xsubset))), border = NA, col = "grey")
  lines(c(val, val), c(0, fX(val)), lty = "dotted")
  lines(c(xsupport[1], xsupport[1]), c(0, fX(xsupport[1])), lty = "dotted")
  lines(c(xsupport[2], xsupport[2]), c(0, fX(xsupport[2])), lty = "dotted")
  lines(xset, fX(xset))

}

# draws region for double-integral for probability calculation
# of the form Pr(X < val) or Pr(X > val). or, if rv = "y", Pr(Y < val) or Pr(Y > val)
# if rv = "x":
# xminmax must be numeric
# yminmax is a list of functions of x; if one element is a numeric it will be coerced to a function
# if rv = "y":
# yminmax must be numeric
# xminmax is a list of functions of x; if any element is a numeric it will be coerced to a function

xyregionplot <- function(xminmax, yminmax, val, rv = "x", 
                         ineq) {
  if (val == as.integer(val)) {
    axsinc <- 1
  } else {
    axsinc <- 0.5
  }
  par(mar = c(2, 3, 0, 0.1))

  if (rv == "x") {
    xseq <- seq(xminmax[1], xminmax[2], 1e-3)
    yarefunctions <- vapply(yminmax, function(y) is.function(y), TRUE)
    yminmax <- lapply(yminmax, function(y) {
      if (is.function(y)) {
        y
      } else {
        function(x) rep(y, length(x))
      }
    })
    
    ymin <- vapply(1:2, function(i) 
      min(yminmax[[i]](xseq)), 0) %>% min
    ymax <- vapply(1:2, function(i) 
      max(yminmax[[i]](xseq)), 0) %>% max  
    plot(0, type = "n", axes = FALSE, xlim = c(xminmax[1] - 0.25, xminmax[2] * 1.05), 
         ylim = c(ymin - 0.25, ymax * 1.05),
         xlab = "", ylab = "",
         main = "", xaxs = "i", yaxs = "i")
    grid()
    xaxseq <- seq(xminmax[1], xminmax[2], axsinc)
    yaxseq <- seq(ymin, ymax, axsinc)
    
    xvertex <- NULL -> yvertex
    if (between(yminmax[[1]](val), ymin, ymax)) {
      # it is a vertex
      xvertex <- c(xvertex, val)
      yvertex <- c(yvertex, yminmax[[1]](val))
    }
    if (between(yminmax[[2]](val), ymin, ymax)) {
      # it is a vertex
      xvertex <- c(xvertex, val)
      yvertex <- c(yvertex, yminmax[[2]](val))
    }
    if (ineq %in% c("<", "<=")) {
      # check intersections of yminmax functions with xminmax[1]
      if (between(yminmax[[2]](xminmax[1]), ymin, ymax)) {
        xvertex <- c(xvertex, xminmax[1])
        yvertex <- c(yvertex, yminmax[[2]](xminmax[1]))
      }
      if (between(yminmax[[1]](xminmax[1]), ymin, ymax)) {
        xvertex <- c(xvertex, xminmax[1])
        yvertex <- c(yvertex, yminmax[[1]](xminmax[1]))
      }
    } else if (ineq %in% c(">", ">=")) {
      # check intersections of yminmax functions with xminmax[2]
      if (between(yminmax[[2]](xminmax[2]), ymin, ymax)) {
        xvertex <- c(xvertex, xminmax[2])
        yvertex <- c(yvertex, yminmax[[2]](xminmax[2]))
      }
      if (between(yminmax[[1]](xminmax[2]), ymin, ymax)) {
        xvertex <- c(xvertex, xminmax[2])
        yvertex <- c(yvertex, yminmax[[1]](xminmax[2]))
      }
    }
    polygon(x = xvertex, y = yvertex, border = NA, col = "grey")
    lines(x = c(val, val), y = c(ymin, ymax), lty = "dotted")
    invisible(lapply(yminmax, function(f) lines(xseq, f(xseq))))
    # lines(x = c(xminmax[1], xminmax[1]), y = c(ymin, ymax))
    # lines(x = c(xminmax[2], xminmax[2]), y = c(ymin, ymax))
  } else if (rv == "y") {
    yseq <- seq(yminmax[1], yminmax[2], 1e-3)
    xarefunctions <- vapply(xminmax, function(x) is.function(x), TRUE)
    xminmax <- lapply(xminmax, function(x) {
      if (is.function(x)) {
        x
      } else {
        function(y) rep(x, length(y))
      }
    })
    
    xmin <- vapply(1:2, function(i) 
      min(xminmax[[i]](yseq)), 0) %>% min
    xmax <- vapply(1:2, function(i) 
      max(xminmax[[i]](yseq)), 0) %>% max  
    plot(0, type = "n", axes = FALSE, ylim = c(yminmax[1] - 0.25, yminmax[2] * 1.05), 
         xlim = c(xmin - 0.25, xmax * 1.05),
         xlab = "", ylab = "",
         main = "", xaxs = "i", yaxs = "i")
    grid()
    xaxseq <- seq(xmin, xmax, axsinc)
    yaxseq <- seq(yminmax[1], yminmax[2], axsinc)
    
    xvertex <- NULL -> yvertex
    if (between(xminmax[[1]](val), xmin, xmax)) {
      # it is a vertex
      yvertex <- c(yvertex, val)
      xvertex <- c(xvertex, xminmax[[1]](val))
    }
    if (between(xminmax[[2]](val), xmin, xmax)) {
      # it is a vertex
      yvertex <- c(yvertex, val)
      xvertex <- c(xvertex, xminmax[[2]](val))
    }
    if (ineq %in% c("<", "<=")) {
      # check intersections of yminmax functions with xminmax[1]
      if (between(xminmax[[2]](yminmax[1]), xmin, xmax)) {
        yvertex <- c(yvertex, yminmax[1])
        xvertex <- c(xvertex, xminmax[[2]](yminmax[1]))
      }
      if (between(xminmax[[1]](yminmax[1]), xmin, xmax)) {
        yvertex <- c(yvertex, yminmax[1])
        xvertex <- c(xvertex, xminmax[[1]](yminmax[1]))
      }
    } else if (ineq %in% c(">", ">=")) {
      # check intersections of yminmax functions with xminmax[2]
      if (between(xminmax[[2]](yminmax[2]), xmin, xmax)) {
        yvertex <- c(yvertex, yminmax[2])
        xvertex <- c(xvertex, xminmax[[2]](yminmax[2]))
      }
      if (between(xminmax[[1]](yminmax[2]), xmin, xmax)) {
        yvertex <- c(yvertex, yminmax[2])
        xvertex <- c(xvertex, xminmax[[1]](yminmax[2]))
      }
    }
    polygon(x = xvertex, y = yvertex, border = NA, col = "grey")
    lines(x = c(xmin, xmax), y = c(val, val), lty = "dotted")
    invisible(lapply(xminmax, function(f) lines(f(yseq), yseq)))
    # lines(x = c(xmin, xmax), y = c(yminmax[1], yminmax[1]))
    # lines(x = c(xmin, xmax), y = c(yminmax[2], yminmax[2]))
  }
  axis(side = 1, at = xaxseq, labels = fractions(xaxseq), pos = 0)
  axis(side = 2, at = yaxseq, las = 1, labels = fractions(yaxseq), pos = 0)
  mtext(expression(x), side = 1, las = 1, line = 0.5)
  mtext(expression(y), side = 2, las = 1, line = 1.5)
}

# probdensityplot(fX = function(x) exp(-x), xsupport = c(0, 10), val = 2, ineq = "<=")
# linepointgraph(fX = function(x) 1/130 * (x^2 + 1), xset = 4:7)

f2latex <- function(f) {
  if (!is.function(f)) {
    f
  } else {
    if (exists("maxx") && exists("maxy")) {
      deparse(body(f)) %>% gsub(pattern = "\\*", replacement = "") %>% 
        gsub(pattern = "maxx", replacement = maxx) %>% 
        gsub(pattern = "maxy", replacement = maxy)
    } else {
      myf <- deparse(body(f)) %>% gsub(pattern = "\\*", replacement = "") %>% 
        gsub(pattern = "[[:space:]]", replacement = "")
      if (regexpr("/", myf) == -1) return(myf)
      lastnondigalphabefore <- suppressWarnings(max(gregexpr("[^[:digit:]|^[:alpha:]]", 
        myf)[[1]][gregexpr("[^[:digit:]|^[:alpha:]]", myf)[[1]] < regexpr("/", myf) - 1]))
      digalphabefore <- suppressWarnings(max(gregexpr("[[:digit:]|[:alpha:]]", 
            myf)[[1]][which(gregexpr("[[:digit:]|[:alpha:]]", myf)[[1]] < regexpr("/", myf) - 1 & 
            gregexpr("[[:digit:]|[:alpha:]]", myf)[[1]] > lastnondigalphabefore)]))
      if (is.infinite(digalphabefore)) digalphabefore <- regexpr("/", myf) - 1
      
      numer <- substr(myf, start = digalphabefore, stop = regexpr("/", myf) - 1)
      
      firstnondigalphaafter <- suppressWarnings(min(gregexpr("[^[:digit:]|^[:alpha:]]", 
                                            myf)[[1]][gregexpr("[^[:digit:]|^[:alpha:]]", myf)[[1]] > regexpr("/", myf) + 1]))
      digalphafter <- suppressWarnings(min(gregexpr("[[:digit:]|[:alpha:]]", 
            myf)[[1]][which(gregexpr("[[:digit:]|[:alpha:]]", myf)[[1]] > regexpr("/", myf) + 1 & 
            gregexpr("[[:digit:]|[:alpha:]]", myf)[[1]] < firstnondigalphaafter)]))
      if (is.infinite(digalphafter)) digalphafter <- regexpr("/", myf) + 1
      
      denom <- substr(myf, start = regexpr("/", myf) + 1, stop = digalphafter)
      fracsubs <- paste0("\\\\dfrac{", numer, "}{", denom, "}")
      gsub(x = myf, pattern = paste0(numer, "/", denom), replacement = fracsubs)
    }
  }
}

sympy2latex <- function(s, bracfrac = FALSE, double = TRUE) {
  s2 <- gsub(pattern = "\\*\\*", replacement = "^", x = s) %>% gsub(pattern = "\\*", replacement = " ")
  plusminusloc <- gregexpr(pattern = "\\+|\\-", s2) %>% unlist
  plusminus <- vapply(1:length(plusminusloc), function(p) substr(s2, start = plusminusloc[p], stop = plusminusloc[p]), "")
  
  if (bracfrac) {
    numdenom <- strsplit(s2, split = "/")[[1]]
    numdenom <- vapply(numdenom, function(n) gsub(pattern = "\\(|\\)", replacement = "", x = n), "")
    plusminusloc <- gregexpr(pattern = "\\+|\\-", numdenom)
    plusminus <- lapply(1:2, function(j) {
      vapply(1:length(plusminusloc[[j]]), function(p) 
        substr(numdenom[j], start = plusminusloc[[j]][p], stop = plusminusloc[[j]][p]), "")
    }) 
    bothparts <- rep(NA_character_, 2L)
    s3 <- strsplit(numdenom, split = "\\+|\\-")
    segs <- lapply(1:2, function(j) rep(NA_character_, length(s3[[j]])))
    for (j in 1:2) {
      for (i in 1:length(s3[[j]])) {
        if (regexpr("/", s3[[j]][i]) == -1) {
          segs[[j]][i] <- s3[[j]][i]
        } else {
          denom <- substr(s3[[j]][i], start = gregexpr("/", s3[[j]][i])[[1]] + 1, 
                          stop = nchar(s3[[j]][i])) %>% trimws
          numer <- substr(s3[[j]][i], start = 1, 
                          stop = regexpr("[^[:digit:]|^[:space:]]", text = s3[[j]][i]) - 1) %>% trimws
          numer[numer == ""] <- 1
          if (regexpr("[[:alpha:]]", s3[[j]][i]) != -1) {
            varpart <- substr(s3[[j]][i], 
                              start = regexpr("[[:alpha:]]", s3[[j]][i]), 
                              stop = regexpr("/", s3[[j]][i]) - 1)
            if (double) {
              segs[[j]][i] <- paste0("\\\\dfrac{", numer, "}{", denom, "}", varpart)
            } else {
              segs[[j]][i] <- paste0("\\dfrac{", numer, "}{", denom, "}", varpart)
            }
            
          } else {
            if (double) {
              segs[[j]][i] <- paste0("\\\\dfrac{", numer, "}{", denom, "}")
            } else {
              segs[[j]][i] <- paste0("\\dfrac{", numer, "}{", denom, "}")
            }
          }
        }
      }
      if (length(s3[[j]]) == 1) {
        bothparts[j] <- segs[[j]][1]
      } else {
        bothparts[j] <- paste(paste0(segs[[j]], c(plusminus[[j]], "")), collapse = "")
      }
    }
    if (double) {
      paste0("\\\\dfrac{", bothparts[1], "}{", bothparts[2], "}")
    } else {
      paste0("\\dfrac{", bothparts[1], "}{", bothparts[2], "}")
    }
  } else {
    s3 <- strsplit(s2, split = "\\+|\\-")
    segs <- rep(NA_character_, length(s3[[1]]))
    for (i in 1:length(s3[[1]])) {
      if (regexpr("/", s3[[1]][i]) == -1) {
        segs[i] <- s3[[1]][i]
      } else {
        denom <- substr(s3[[1]][i], start = gregexpr("/", s3[[1]][i])[[1]] + 1, stop = nchar(s3[[1]][i])) %>% trimws
        numer <- substr(s3[[1]][i], start = 1, stop = regexpr("[^[:digit:]|^[:space:]]", text = s3[[1]][i]) - 1) %>% trimws
        numer[numer == ""] <- 1
        if (regexpr("[[:alpha:]]", s3[[1]][i]) != -1) {
          varpart <- substr(s3[[1]][i], 
                            start = regexpr("[[:alpha:]]", s3[[1]][i]), 
                            stop = regexpr("/", s3[[1]][i]) - 1)
          if (double) {
            segs[i] <- paste0("\\\\dfrac{", numer, "}{", denom, "}", varpart)
          } else {
            segs[i] <- paste0("\\dfrac{", numer, "}{", denom, "}", varpart)
          }
        } else {
          if (double) {
            segs[i] <- paste0("\\\\dfrac{", numer, "}{", denom, "}")
          } else {
            segs[i] <- paste0("\\dfrac{", numer, "}{", denom, "}")
          }
        }
      }
    }
    if (length(s3[[1]]) == 1) {
      segs[1]
    } else {
      paste(paste0(segs, c(plusminus, "")), collapse = "")
    }
  }
}

frac2dfrac <- function(x, double = TRUE) {
  x <- as.character(x)
  if (regexpr("/", x[1]) == -1) {
    x
  } else {
    numerator <- substr(x, start = 1, stop = regexpr("/", x) - 1)
    denominator <- substr(x, start = regexpr("/", x) + 1, stop = nchar(x))
    if (double) {
      paste0("\\\\dfrac{", numerator, "}{", denominator, "}")
    } else {
      paste0("\\dfrac{", numerator, "}{", denominator, "}")
    }
  }
}

frac2latex <- function(x, double = TRUE) {
  x <- as.character(x)
  if (regexpr("/", x[1]) == -1) {
    x
  } else {
    numerator <- substr(x, start = 1, stop = regexpr("/", x) - 1)
    denominator <- substr(x, start = regexpr("/", x) + 1, stop = nchar(x))
    if (double) {
      paste0("\\\\frac{", numerator, "}{", denominator, "}")
    } else {
      paste0("\\frac{", numerator, "}{", denominator, "}")
    }
  }
}

frac2rational <- function(expr) {
  expr <- as.character(expr)
  slashes <- gregexpr("/", expr) %>% unlist
  if (slashes[1] == -1) {
    return(expr)
  } else {
    exponentpos <- gregexpr("\\*\\*", expr)[[1]]
    if (exponentpos[1] != -1) {
      charb4exponent <- vapply(exponentpos, function(e) 
        substr(expr, start = e - 1, stop = e - 1), "")
      patrep <- lapply(1:length(exponentpos), function(i) {
        if (!is.na(suppressWarnings(as.numeric(charb4exponent[i])))) {
          lastnondigitb4base <- max(gregexpr("[^[:digit:]]", expr)[[1]][gregexpr("[^[:digit:]]", expr)[[1]] < exponentpos[i]])
          pattern <- paste0(substr(expr, start = lastnondigitb4base + 1, stop = exponentpos[i] - 1),
                            "\\*\\*", substr(expr, start = exponentpos[i] + 2, stop = exponentpos[i] + 2))
          pattern2 <- substr(expr, start = lastnondigitb4base + 1, stop = exponentpos[i] + 2)
          replacement <- sympy(as.character.Sym(paste0("simplify(", pattern2, ")"))) %>% as.character
          list("pattern" = pattern, "replacement" = replacement)
        } else {
          list()
        }
      })
      for (i in length(patrep)) {
        if (length(patrep[[i]]) > 0) {
          expr <- gsub(pattern = patrep[[i]]$pattern, replacement = patrep[[i]]$replacement, 
                       x = expr)
        }
      }
    }
    slashes <- gregexpr("/", expr) %>% unlist
    startpoint <- c(1, slashes[-length(slashes)] - 1)
    stoppoint <- c(slashes[-1] - 1, nchar(expr))
    numb4slash <- vapply(slashes, function(s) 
      regexpr("[[:digit:]]", substr(expr, start = s - 1, stop = s - 1)) != -1, TRUE)
    lastnonnumb4slash <- lapply(1:length(slashes), function(i) {
      gregexpr("[^[:digit:]]", substr(expr, start = startpoint[i], 
          stop = slashes[i] - 1))[[1]] + startpoint[i] - 1})
    exponentb4slash <- vapply(1:length(slashes), function(i) {
      stuff <- lastnonnumb4slash[[i]]
      if (stuff[1] == -1 || length(stuff) == 1) {
        FALSE
      } else {
        substr(expr, start = stuff[length(stuff) - 1], stop = stuff[length(stuff)]) == "**"
      }
    }, TRUE)
    isnumericfraction <- which(numb4slash & !exponentb4slash)
    
    if (length(isnumericfraction) == 0) {
      expr
    } else {
      numerator <- vapply(isnumericfraction, function(j) {
        stuff <- lastnonnumb4slash[[j]]
        substr(expr, start = max(1, stuff[length(stuff)] + 1), 
               stop = slashes[j] - 1) %>% as.numeric
      }, 0)
      denominator <- vapply(isnumericfraction, function(j) {
        stringafterslash <- substr(expr, start = slashes[j] + 1, stop = stoppoint[j])
        firstnonnumberafterslash <- regexpr("[^[:digit:]]", stringafterslash)
        if (firstnonnumberafterslash == -1) {
          as.numeric(stringafterslash)
        } else {
          firstnonnumberafterslash <- firstnonnumberafterslash + slashes[j]
          substr(expr, start = slashes[j] + 1, 
                 stop = firstnonnumberafterslash - 1) %>% as.numeric
        }
      }, 0)
      fractionexpr <- paste0(numerator, "/", denominator)
      
      newexpr <- expr
      for (k in seq_along(isnumericfraction)) {
        newexpr <- gsub(pattern = fractionexpr[k], 
                        replacement = paste0("Rational(", numerator[k], ", ", denominator[k], ")"), 
                        x = newexpr)
      }
      newexpr
    } 
  }
}

num2rational <- function(x) {
  frac2rational(fractions(x))
}

ineq2latex <- function(i, double = TRUE) {
  if (i %in% c("<", ">")) {
    i
  } else if (i %in% c("<=", ">=")) {
    if (double) {
      paste0("\\\\", ifelse(i == ">=", "g", "l"), "e")
    } else {
      paste0("\\", ifelse(i == ">=", "g", "l"), "e")
    }
  } else if (i == "!=") {
    paste0(ifelse(double, "\\\\", "\\"), "ne")
  }
}

inequality2latex <- function(i, double = TRUE) {
  if (is.function(i)) i <- deparse(body(i))
  if (regexpr("<=|>=|!=", i) != -1) {
    if (double) {
      gsub(i, pattern = ">=", replacement = "\\\\ge") %>% 
        gsub(pattern = "<=", replacement = "\\\\le") %>% 
        gsub(pattern = "!=", replacement = "\\\\ne") 
    } else {
      gsub(i, pattern = ">=", replacement = "\\ge") %>% 
        gsub(pattern = "<=", replacement = "\\le") %>% 
        gsub(pattern = "!=", replacement = "\\ne")
    }
  } else if (regexpr("<|>", i) != -1) {
    i    
  } else if (regexpr("==", i) != -1) {
    gsub(i, pattern = "==", replacement = "=")
  }
}


pdfgraph <- function(fX, xlimits, xval, ineq) {
  par(mar = c(4, 5, 0, 1))
  xset <- seq(xlimits[1], xlimits[2], 1e-3)
  ymax <- max(fX(xset)) * 1.1
  
  plot(0, type = "n", axes = FALSE, xlim = c(min(xset) / 1.05, max(xset) * 1.05), 
       ylim = c(0, ymax),
       xlab = "", ylab = "",
       main = "", xaxs = "i", yaxs = "i")
  grid()
  axis(side = 1, at = seq(min(xset, 0), max(xset), 1))
  mtext(expression(x), side = 1, las = 1, line = 2.5)
  mtext(expression(f[X](x)), side = 2, las = 1, line = 3)
  axis(side = 2, at = seq(0, ceiling(ymax * 10) / 10, 0.1), las = 1)
  xsubset <- xset[which(do.call(ineq, list(xset, xval)))]
  polygon(c(xsubset, rev(xsubset)), c(rep(0, length(xsubset)), fX(rev(xsubset))), border = NA, col = "grey")
  lines(c(xval, xval), c(0, fX(xval)), lty = "dotted")
  lines(xset, fX(xset))
}

f2latex <- function(fsymp, double = TRUE, dfrac = FALSE) {
  d <- paste0("latex(simplify(", fsymp, "))") %>% sympy %>% 
    as.character %>% gsub(pattern = "\\$", replacement = "")
  
  if (dfrac) {
    d <- d %>% 
      gsub(pattern = "frac", replacement = "dfrac")
  }
  
  if (double) {
    gsub(d, pattern = "\\\\", replacement = "\\\\\\\\")
  } else {
    d
  }
}

afterlatex <- function(expr, double = FALSE, dfrac = FALSE) {
  expr <- gsub(x = expr, pattern = "\\$", replacement = "")
  if (dfrac) {
    expr <- gsub(x = expr, pattern = "frac", replacement = "dfrac")
  }
  if (double) {
    gsub(expr, pattern = "\\\\", replacement = "\\\\\\\\")
  } else {
    expr
  }
}

# Old version
# sympynlatex <- function(expr, double = FALSE) {
#   as.character.Sym(expr) %>% sympy %>% afterlatex(double = double)
# }

sympynlatex <- function(expr, double = TRUE) {
  paste0("latex(", expr, ")") %>% 
    as.character.Sym(expr) %>% 
    sympy %>% 
    afterlatex(double = double)
}


justsympy <- function(expr) {
  as.character.Sym(expr) %>% sympy
}

split_str_by_index <- function(target, index) {
  index <- sort(index)
  substr(rep(target, length(index) + 1),
         start = c(1, index),
         stop = c(index -1, nchar(target)))
}

#Taken from https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
interleave <- function(v1,v2)
{
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}

insert_str <- function(target, insert, index) {
  insert <- insert[order(index)]
  index <- sort(index)
  paste(interleave(split_str_by_index(target, index), insert), collapse="")
}
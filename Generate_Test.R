# Load packages
library(exams)
library(tidyverse)
library(knitr)
library(MASS)
library(pracma)
library(rSymPy)

source("mathstatsfunctions.R")

mytest <- list("example_essay_question.Rnw",
               "example_numerical_question.Rnw")

maxn <- 10L
set.seed(1234)
rm(cases2)
exams2html(mytest, n = maxn, mathjax = TRUE, converter = "pandoc-mathjax", 
           name = "Test1_", dir = "Test1")

set.seed(1234)
rm(cases2)
exams2blackboard(mytest, n = maxn, mathjax = FALSE, 
                 name = "Test1_",
                 points = 6, dir = "Test1")
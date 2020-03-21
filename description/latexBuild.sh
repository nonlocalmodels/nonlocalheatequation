#!/bin/bash
# 
# name of tex file
#
mainTex="problem_description"

# #
# # run pdflatex once
# #
# pdflatex -synctex=1 -interaction=nonstopmode "$mainTex"".tex" 

# #
# # run bibtex
# #
# bibtex "$mainTex"

#
# run pdflatex twice to make sure the bibliography works correctly
#
pdflatex -synctex=1 -interaction=nonstopmode "$mainTex"".tex"
pdflatex -synctex=1 -interaction=nonstopmode "$mainTex"".tex"
pdflatex -synctex=1 -interaction=nonstopmode "$mainTex"".tex"

# clean up the debug output files
rm *~  &> /dev/null
rm *#  &> /dev/null
rm *.gz &> /dev/null
rm *.out &> /dev/null
rm *.aux &> /dev/null
# rm *.bbl &> /dev/null
rm *.blg &> /dev/null
rm *.spl &> /dev/null
rm *.log &> /dev/null
rm *.thm &> /dev/null
rm *.nav &> /dev/null
rm *.snm &> /dev/null
rm *.toc &> /dev/null
rm *.bak &> /dev/null
rm *.brf &> /dev/null
	 

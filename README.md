# SuroModel
Engineering Design via Surrogate Modelling: A Practical Guide

# What
Downloaded from https://optimizationcodes.wordpress.com/ (2017)

# Original Readme

This file is part of the Engineering Design via Surrogate Modelling: A Practical Guide Maltab code.

    This code is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any 
    later version.

    This code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this toolset.  If not, see <http://www.gnu.org/licenses/>.


The code is split into the folders:

	\Sampling Plans,
	\Constructing a Surrogate, and
	\Exploring and Exploiting a Surrogate,

which follow the chapters in part I of the book,

	\Advanced Concepts,

which includes code related to part II,

	\Example Problems,

which contains the problems detailed in the Appendix, and

	\Example Scripts,

which contains a number of scripts to demonstrate the use of the code.

Start by copying these folders, for example into your MATLAB directory. 
Then try going to the 'Example Scripts' folder where you can study and run the various examples.


Reported errors in the code have been resolved as of March 2010 

(please see notes below and comments in code). 
Please report any further erros to Alexander.Forrester@soton.ac.uk 


********************************************************
*** NOTE ERROR multiei.m WHICH HAS NOW BEEN RESOLVED ***
********************************************************


* Main changes September 2009 *

multiei.m:
This code would not work properyl unless called via constrainedmultiei.m. 
This has now been resolved - see comments in function.

All kriging predictor functions:
The "small number" used to avoid floating point underflow / log of zero has been changed to "realmin".
This gives more accuracy, but may occasionally yeild complex numbers.

screeningplot.m:
There was an error which made it want another (non-existent) random orientation for certain inputs.
This error has been resolved.

ga.m
This has been "tidied up" and gademo.m added to \Example Scripts,
which calls a new example problem: rastrigin.m in \Example Problems.

* Main changes March 2010 *

multiei.m:
This had an error in the calculation of Ybar, which also appears in the book on page 185 (see errata for book). 
This has now been resolved.

* Main changes November 2010 *

buildcokriging.m
eps term removed from diagonal of PsicXeXc matrix

ga.m
Major overhaul, with more control of evolutionary parameters now possible

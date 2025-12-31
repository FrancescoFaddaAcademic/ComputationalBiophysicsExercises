# Computational Biophysics Assignments
This repository contains the required assignments, each assignment has the structure of a julia package to achieve a good level of portability and version control.

Each assigment contains:
- a `src/` directory containing most of the implementation required for the completion of the assignment.
- a `test/` directory containing scripts that should reproduce the required behaviour for each assigments.
- a `README.md` file that should provide the documentation of the code.
- `Project.toml` and `Manifest.toml`, two files that store useful information for dependencies and version control.

## How to run the code 
To get started I provide some instruction for the setup of the environment required for the correct functioning of the scripts.
### Install Julia
You can install julia from [julialang.org](https://julialang.org/)
### Setup the environment
NOTE: The following instruction are tested on linux, the proces may be slightly different for other OS.
As I was saying, each assignments has the structure of a julia Package, to run them one must first properly setup the julia environment.
#### Download the Packages
First download the packages from this directory, you can keep them together or separate, they are completely independent from one another.
You can do this by cloning this repository (using git), or by manually downloading the Zip file.
#### Position yourself
Locate the files and, from the terminal, using `cd` position yourself in the Package directory (`SantaLucia` or `HPModel` depending on which you want to run) 
#### Launch Julia
Julia scripts run in an environment called REPL, to start it, it should be sufficient to simply type `julia` in the terminal (to exit run `exit()`)
#### Enter Package mode and setup the environment
By typing `]` you will enter the julia package manager 
Run `activate .` to activate the local environment for the package manager.
Run `instatiate` to install all the required dependencies.
Then go back to the standard REPL mode (Ctrl + C).
You are now ready to run any script in the directory!
### Run the Packages
To run a julia file run `include("[relativepath]")` (For example to run the test reproducing the results from the Chan-Dill1989 paper you can run from HPModel `include("test/Chan_Dill_1989.jl")`)

## Why Julia?
I really enjoy a very bare metal approach to programming, I understand however that for most applications a code written in C may create many more problems than those it manages to solve. With this out of the way, I still really dislike the very high level of Python, in particular its non statically typed and non compiled nature, together with the forced indentation. I find it very difficult to read and mantain, moreover, the speed for general non vectorizable operation is often underwhelming. I've found that for me, Julia strikes a good balance between being easy to write and memeory safe, while still achieving impressive performances (if handled correclty) because of its JIT compilation. Moreover I find particularly useful that Julia lets you decide how low level you want to go. For rapid prototyping and scientific implementations it is great! This said the support and some of the features that come naturally with a higher level language are still a little annoying sometimes as what happens under the hood may not be very clear.

## Some additional (and personal) notes
The provided scripts aren't by any means perfect, in fact some of the functions defined in the main package files are not even accessed. The reason behind this is that I've tried to explore different approaches to solving the proposed problems together with some of the inner mechanisms of this language (with varying success) and now I am quite honestly running out of time. The result, while properly working, is a little incoherent, I thought however to keep everything as I would like to discuss some of the ideas that didn't make it into the final implementation due to technical and temporal problems.

Thank you for your attention :)

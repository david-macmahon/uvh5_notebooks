# Pluto notebooks for learning about UVH5

[Pluto.jl](https://plutojl.org) is a lightweight framework for creating reactive notebooks in [Julia](https://julialang.org).  The Pluto notebooks that accompany this `README` file are tutorials for interacting with the UVH5 file format used to store radio interferometer visibility data.  They also provide an introduction to some Julia programming concepts and radio astronomy concepts in general.

## First time setup

Before you can run these Pluto notebooks you need to have Julia and Pluto installed on your system.  The website for the *Introduction to Computational Thinking* course at MIT has great [step-by-step instructions](https://computationalthinking.mit.edu/Spring21/installation/) for getting Julia and Pluto up and running on your system.

In addition to the Firefox and Chrome browsers that you may see referenced in the installation instructions, the Safari browser also works well with these notebooks, so no need to panic if that's your favorite browser.

## Running the notebooks

Once you have Julia and Pluto installed, you can run the notebooks straight off the web without even cloning the repository containing this `README` and the notebooks!

Simply follow these steps:

1. Run Pluto as described in the instructions linked to above
2. Open the indicated URL in your browser (if it didn't opwn automatically)
3. Copy/paste the URL of the notebook "raw" format (see below) into the "Open from file:" text box
4. Click the "Open" button

The first time you load a notebook it will download any additional Julia packages needed by the notebook (e.g. for plotting, HDF5 file access, etc).  This may take a little while, but the next time you open the notebook it will be faster.  The downloaded packages are shared across notebooks, so the initialization time will be longest for the first notebook you run.  The notebooks also download an example UVH5 data file that is similarly shared between the notebooks so it only gets downloaded and stored on your system one time.  The example data file is a 2.3 MB download that expands to a 4.4 MB file.

All of the notebooks here have a Table of Contents that will appear in the upper right corner of the notebook when it has finished initializing.

If you run a notebook from the web, a copy of the notebook will be saved locally, but Pluto will give it an arbitrary name like `Random_notebook.jl`, so you may want to save it with a more meaningful name using the "Save notebook..." field at the top of the notebook.  Links to locally saved notebooks will appear in the "Recent sessions" section of Pluto's welcome page.  If you hover your cursor over these links, most browsers will show the path name of the linked notebook file.

Of course you can also clone this repository to get local copies of these notebooks all at once.  This is the recommended way to get and run the notebooks if you want to create pull requests to improve the notebooks.  Contributions welcome!

## Notebook Overview

The notebooks are presented here in the recommended reading order.  After the first notebook, **UVH5 explorations**, the others can be read in any order, but later notebooks tend to build on earlier notebooks.  If something in a later notebook seems unclear, you might want to go back to an earlier notebooks to see whether they clarify things for you.  Feel free to open a GitHub issue if you think something could use a better explanation.

The links here are suitable for pasting into the "open" text box on Pluto's welcome page.

- [`1/uvh5_explorations.jl`](1/uvh5_explorations.jl) is recommended as the starting notebook.  It covers the basics of the UVH5 file format, radio astronomy in general, and even some Julia programming techniques.
- [`2/uvh5_and_spatial_awareness.jl`](2/uvh5_and_spatial_awareness.jl) goes into detail about some of the various reference frames used with radio telescope arrays, how to convert between them, and how to read various positional data from a UVH5 data file.  It is far from exhaustive so you can look forward to future notebooks that delve into additional details on this topic.

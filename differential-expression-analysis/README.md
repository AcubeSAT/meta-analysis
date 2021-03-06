<div align="center">
<p>
    <a href="https://benchling.com/organizations/acubesat/">Benchling ๐๐งฌ</a> &bull;
    <a href="https://gitlab.com/acubesat/documentation/cdr-public/-/blob/master/DDJF/DDJF_PL.pdf?expanded=true&viewer=rich">DDJF_PL ๐๐งช</a> &bull;
    <a href="https://spacedot.gr/">SpaceDot ๐๐ช</a> &bull;
    <a href="https://acubesat.spacedot.gr/">AcubeSAT ๐ฐ๏ธ๐</a>
</p>
</div>

## Description

Coming soon :tm:

## Table of Contents

<details>
<summary>Click to expand</summary>

[[_TOC_]]

</details>


## Getting Started

### TL;DR

1. `git clone --depth 1 https://gitlab.com/acubesat/su/bioinformatics/meta-analysis.git`
2. `docker pull xlxs4/meta-analysis:latest`
3. `docker run -ti --entrypoint /bin/bash -v $PWD/meta-analysis:/meta-analysis xlxs4/meta-analysis`
4. `cd ../meta-analysis`
5. `Rscript differential-expression-analysis/src/flocculation.R --qc --plots --feather`

### Docker

First of all, you'll need [Docker](https://www.docker.com/) :whale:.
You can download and install [Docker for Desktop](https://www.docker.com/products/docker-desktop) for starters.

### Grab the Docker image

* Open your favorite terminal. If you think you don't have one, you're on Windows.
In that case, get [Windows Terminal](https://aka.ms/terminal)
* You can pull the image straight from the cloud: `docker pull xlxs4/meta-analysis:latest`
* ...  or you can build it yourself (something like `docker build -t meta-analysis .`)


### Clone the repository

* `Explorer` (`Ctrl + Shift + E`) -> `Clone Repository` [:camera:](https://prnt.sc/1bu52pf)
* Clone the `meta-analysis` repository [:camera:](https://prnt.sc/1bu5ttp)  [:camera:](https://prnt.sc/1bu602g)
* Open the cloned repository [:camera:](https://prnt.sc/1bu63we)  [:camera:](https://prnt.sc/1bu6990)
* :bulb: Hint: click on a :camera: when you find one

### Create a container

* `docker run -ti -v $PWD/meta-analysis:/meta-analysis xlxs4/meta-analysis` 
* `Ctrl + d` to exit. [:camera:](https://prnt.sc/1bu4lcr) [:camera:](https://prnt.sc/1bu4p8s)

### Start it up

* Everything below will help you get an IDE running inside the container.
  The cross-platform [Visual Studio Code](https://code.visualstudio.com/) will be used.
  If you want to use something else, you probably know what you're doing, so go to the [TL;DR](#tldr) section
* Get the VSCode [`remote-containers`](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) [:camera:](https://prnt.sc/1bueb0d)
  While not necessary, you could grab the [Docker extension](https://marketplace.visualstudio.com/items?itemName=ms-azuretools.vscode-docker) as well [:camera:](https://prnt.sc/1btzhi3)
* `Docker Containers: Start` [:camera:](https://prnt.sc/1bty5z9)  -> `Individual Containers` [:camera:](https://prnt.sc/1btyb4t) -> <container_hash>  [:camera:](https://prnt.sc/1btyimw)

### Get in

* `Open a Remote Window` [:camera:](https://prnt.sc/1btzqbd) -> `Attach to Container` [:camera:](https://prnt.sc/1btzy1h) -> <container_hash> [:camera:](https://prnt.sc/1bu0ffu)

### Development tips

* Installing `python-pip3` and using it to install [`radian`](https://github.com/randy3k/radian) inside the container is highly recommended
* If you want to update the image, e.g. for CI/CD purposes, but have made the changes to the `renv` lockfile inside an attached running Docker container, you can `git clone` the repository in your main working environment, and then grab the updated `renv.lock` file, like this: `docker cp <container_name>:"<path_to_renv.lock>" .`
  Then, just do `docker build -t meta-analysis .`, `docker image tag meta-analysis:latest <repo_name>/meta-analysis:latest` and `docker push <repo_name>/meta-analysis:latest`
* [`conflicted`](https://www.tidyverse.org/blog/2018/06/conflicted/) is a great `tidyverse` package to check for conflicts arising from ambiguous function names. From within the container, open the R terminal (e.g. `radian`), and `install.packages("conflicted")`. Then, you can just load it (`library(conflicted)`) in the running session. To re-check for conflicts, simply run `conflicted::conflict_scout()`
* Since you mounted the volume, any changes you do in `meta-analysis` while inside the container will persist in the host directory even after closing the container. This is bidirectional: you can keep updating the `meta-analysis` directory from outside, and the changes will be immediately reflected inside the container. You can also work inside the container, and sign your commits with `git commit -S` out of the box!

## Results

Can be found in `results/flocculation`.

### DE Genes

|PROBEID      | adj.P.Val|  P.Value|         t|        B|      logFC|   AveExpr|ENSEMBL   |GENENAME |
|:------------|---------:|--------:|---------:|--------:|----------:|---------:|:---------|:--------|
|1771832_at   | 0.0001419| 0.00e+00|  9.205550| 9.055121|  0.8654733| 11.989411|YDR418W   |RPL12A   |
|1776010_at   | 0.0001571| 1.00e-07|  8.762955| 8.357195|  1.0428640| 11.629032|YJL190C   |RPS22A   |
|1777989_at   | 0.0006174| 3.00e-07|  7.763274| 6.773180|  0.9994524| 11.189723|YLR406C   |RPL31B   |
|1773591_at   | 0.0007119| 5.00e-07|  7.531448| 6.383697|  1.0416106|  8.927194|YMR199W   |CLN1     |
|1772366_at   | 0.0015277| 1.40e-06|  7.013690| 5.483140|  1.1977834| 11.651799|YML063W   |RPS1B    |
|1776927_at   | 0.0019711| 2.40e-06|  6.716756| 4.947542|  0.8278819| 10.880286|YGL147C   |RPL9A    |
|1773292_at   | 0.0019711| 2.50e-06|  6.689588| 4.923252|  0.6378308| 12.136133|YGR034W   |RPL26B   |
|1777846_at   | 0.0019711| 3.40e-06|  6.540306| 4.638924|  0.6715883| 11.060921|YBL002W   |HTB2     |
|1773827_s_at | 0.0019711| 3.40e-06|  6.525859| 4.620072|  0.5766213| 12.068804|YHL033C   |RPL8A    |
|1776763_at   | 0.0019711| 3.80e-06|  6.488877| 4.530887|  0.7016451|  9.108503|YPL256C   |CLN2     |
|1772978_at   | 0.0019711| 3.80e-06| -6.472194| 4.519797| -0.5752284|  8.262176|YDR430C   |CYM1     |
|1769403_at   | 0.0041044| 8.70e-06| -6.083449| 3.759518| -0.7739060| 10.196870|YKL195W   |MIA40    |
|1772208_at   | 0.0041594| 9.60e-06|  6.026704| 3.670536|  0.5767523| 11.099817|YBR048W   |RPS11B   |
|1774263_at   | 0.0041965| 1.04e-05|  5.998287| 3.595161|  0.7473007|  9.645624|YNL327W   |EGT2     |
|1773091_at   | 0.0042034| 1.13e-05|  5.959616| 3.520182|  0.7765164| 11.678413|YFR031C-A |RPL2B    |
|1770023_at   | 0.0042034| 1.25e-05|  5.898664| 3.418212|  0.5908748| 13.024148|YDL055C   |PSA1     |
|1776370_at   | 0.0042034| 1.26e-05|  5.904253| 3.412461|  0.7537768| 12.057036|YFR031C-A |RPL2A    |
|1774667_at   | 0.0047992| 1.54e-05| -5.799085| 3.225383| -0.5450599|  9.621721|YGL049C   |TIF4632  |
|1776049_at   | 0.0047992| 1.61e-05|  5.777935| 3.182707|  0.5683761| 12.634236|YNL302C   |RPS19B   |
|1775495_at   | 0.0055168| 1.96e-05|  5.691768| 2.999526|  0.6100716|  7.731038|YDL003W   |MCD1     |
|1770653_at   | 0.0055168| 2.05e-05|  5.674419| 2.960671|  0.7209517| 11.434438|YGR214W   |RPS0A    |
|1778823_at   | 0.0056583| 2.20e-05| -5.640505| 2.893392| -0.6360340|  8.052246|YIL072W   |HOP1     |
|1775193_at   | 0.0057160| 2.35e-05| -5.607095| 2.829481| -0.6077096| 10.195407|YFL036W   |RPO41    |
|1776223_at   | 0.0057160| 2.42e-05| -5.585044| 2.800234| -0.3950174| 12.317287|YOR362C   |PRE10    |
|1775167_at   | 0.0058560| 2.59e-05|  5.564186| 2.741431|  0.9589474| 11.829023|YER131W   |RPS26B   |
|1776390_at   | 0.0062352| 2.95e-05| -5.493755| 2.617067| -0.4972437| 10.673441|YBL022C   |PIM1     |
|1779099_at   | 0.0062352| 2.98e-05| -5.489394| 2.608290| -0.4162036|  8.279802|YLR095C   |IOC2     |
|1775536_at   | 0.0064018| 3.17e-05|  5.469372| 2.551589|  1.0407435| 11.613532|YLR061W   |RPL22A   |
|1773427_at   | 0.0074819| 3.88e-05|  5.375099| 2.361702|  0.7031388| 10.371738|YLR367W   |RPS22B   |
|1779506_at   | 0.0074819| 3.97e-05|  5.364601| 2.340491|  1.0925925|  9.943237|YEL040W   |UTR2     |

### DE Genes (treat)

|PROBEID      | adj.P.Val| P.Value|        t|    logFC|  AveExpr|ENSEMBL   |GENENAME |
|:------------|---------:|-------:|--------:|--------:|--------:|:---------|:--------|
|1772279_at   |         0|       0| 401.1422| 15.07725| 15.07915|YMR251W-A |HOR7     |
|1777569_at   |         0|       0| 383.2136| 14.86353| 14.85811|YBR072W   |HSP26    |
|1774896_at   |         0|       0| 382.6178| 14.74570| 14.71289|YJR104C   |SOD1     |
|1776636_s_at |         0|       0| 378.0824| 14.61553| 14.61095|YBR118W   |TEF2     |
|1771705_at   |         0|       0| 376.5488| 14.73980| 14.71536|YER150W   |SPI1     |
|1775388_at   |         0|       0| 369.9227| 14.46565| 14.44869|YJL052W   |TDH1     |
|1774354_at   |         0|       0| 366.8373| 14.42229| 14.41623|YOL052C-A |DDR2     |
|1771908_at   |         0|       0| 366.1760| 14.68892| 14.67079|YNL160W   |YGP1     |
|1777549_at   |         0|       0| 355.4468| 14.61756| 14.60979|YPR036W-A |SPO24    |
|1776517_at   |         0|       0| 354.6041| 14.32870| 14.30089|YJL151C   |SNA3     |



## Technical

### Computational Reproducibility

What is the point in conducting research behind open doors? Research must be accessible and reproducible. But it must also be **easily** reproducible. Turns out the latter is way more difficult to achieve than the former. There's a couple steps taken to make sure this pipeline is 100% reproducible, but first:

* [Here](https://www.frontiersin.org/articles/10.3389/fninf.2017.00069/full)'s a quick rundown on **replicability**
* And a little bit [more](http://web.archive.org/web/20210908093208/https://ropensci.github.io/reproducibility-guide/sections/introduction/) from the great folks over at [rOpenSci](https://ropensci.org/)

#### renv

> The renv package is a new effort to bring project-local R dependency management to your projects. [...] Underlying the philosophy of renv is that any of your existing workflows should just work as they did before โ renv helps manage library paths (and other project-specific state) to help isolate your projectโs R dependencies, and the existing tools youโve used for managing R packages (e.g. install.packages(), remove.packages()) should work as they did before.

You code some analyses. To run these analyses, you use several packages (aka prewritten bundled code from other kind people). You upload the code. Everyone can download the source and play around tinkering with your creation! Right? Most likely no... Here's some of the reasons why:

* You use a package that depends on (requires) other packages that depend on other packages that depend... and everything is a tedious mess to go through in order to get every prerequisite installed
* You use a package that has undergone several modifications (or only one) since the time that you used it, and now its functionality has changed and your code no longer works and people might not even know why
* You use a package that has undergone several modifications (sound familiar?) since the time that you used it, and while the code seems to work, something is slightly off-kilter. Maybe no one spots there are mistakes in the first place
* You use a lot of different packages to develop your pipeline whilst not working on an isolated environment; you end up using pieces of code (mainly functions) from packages that you already had installed in your system but forgot to mention as dependencies for your pipeline to run... etc.

`renv` to the rescue! By employing renv we can hopefully alleviate some issues that can quickly become painstakingly difficult to resolve if not dealt with snappily enough. From [this](https://www.youtube.com/watch?v=yjlEbIDevOs) great little presentation by [Kevin](https://kevinushey.github.io/):

> You can use `renv` to make your projects more:
> * **Isolated**: Each project gets its own library of R packages, so you can feel free to upgrade and change package versions in one project without worrying about breaking your other projects. (read: and vice-versa)
> * **Portable**: Because `renv` captures the state of your R packages within a _lockfile_, you can more easily share and collaborate on projects with others, and ensure that everyone is working from a common base.
> * **Reproducible**: Use `renv::snapshot()` to save the state of your r  _library_ to the _lockfile_ `renv.lock`. You can later use `renv::restore()` to restore your R library exactly as specified in the lockfile. 

So, in brief:

* `renv` isolates the project; there's no external libraries bleeding in during development, meaning there's no way for you to miss explicitly stating a dependency due to an oversight
* There's a _lockfile_, called `renv.lock`, that captures all project packages that are directly used, including their dependencies. For every package, the version is also captured. Anyone trying to fiddle with the project can be sure to work with the exact same packages as you, or at least close to that
* Mainly due to the lockfile, all packages necessary can be installed automatically. This facilitates automation (e.g. for testing), makes it easier for the user, _and_ prevents them from any mistakes they could make during a manual installation
* and many more :telescope:

Without getting in too much detail, you can simply use `renv` like so:

1. Install `renv` in your local machine, by running something along the lines of `R -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv')"`
2. Just grab the `renv.lock` lockfile, and run `R -e "renv::restore()`
3. That's it!

#### Docker

From the `renv` webpage:

> It is important to emphasize that renv is not a panacea for reproducibility. Rather, it is a tool that can help make projects reproducible by solving one small part of the problem: it records the version of R + R packages being used in a project, and provides tools for reinstalling the declared versions of those packages in a project. Ultimately, making a project reproducible requires some thoughtfulness from the user: what does it mean for a particular project to be reproducible, and how can renv (and other tools) be used to accomplish that particular goal of reproducibility?

> There are a still a number of factors that can affect whether this project could truly be reproducible in the future โ for example,

> The results produced by a particular project might depend on other components of the system itโs being run on โ for example, the operating system itself, the versions of system libraries in use, the compiler(s) used to compile R and the R packages used, and so on. Keeping a โstableโ machine image is a separate challenge, but Docker is one popular solution. [...]

Using `renv` alone can not and will not guarantee your analysis is reproducible.
Let's try to investigate what's going on here. Imagine these few from the many possible scenarios:

* You use a package that depends on a non-R library that the user has not installed in their machine
* You use a package that depends on a non-R library that the user can not install in their machine
* You use a package that depends on a non-R library that the user has installed at a different version, and changing this could break other parts of their system
* A non-R dependency works different on the user's machine, even if using the exact same version
* The user's machine uses different computational routines (BLAS, OpenBLAS, LAPACK...) that produce slightly different results (e.g. see [this](vhttps://scicomp.stackexchange.com/questions/26137/are-blas-implementations-guaranteed-to-give-the-exact-same-result))
* Things on the user's machine were compiled in a different manner
* The user's machine runs on different hardware
* ...

So, as you can see, going from mostly reproducible to reproducible involves circumventing a lot of subtle and not-so-subtle pitfalls, not to mention the effort the user must go through to correct something that is awry. It is always better to put in the extra effort yourself, if it is to save others from committing 10x your effort.

In other words, it seems that to get rid of most of our problems, we would need to somehow ensure that the environment state in which we develop must remain identical for the user, not only at the R dependency scope, but many, many abstraction layers lower beyond that. Thankfully, [Docker](https://www.docker.com/) is waiting right around the corner and is what all the cool kids already do.

What Docker allows you to do is essentially bundle up your pipeline, together with all dependencies (libraries and/or binaries), inside a box. Then, you can fetch this box in the machine of your choosing, and look inside it to be transferred in an identical environment/state hustle-free. You can check [this](https://www.youtube.com/watch?v=TvnZTi_gaNc) presentation for a little bit more on the topic.

This box is called a Docker **image** and the blueprint to create it is called a **Dockerfile**. The isolated environment where the pipeline runs is called a Docker **container**.

To get a better look into how this works in our usecase, and why it facilitates replicability, let's take a look into our [`Dockerfile`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/Dockerfile):

```docker
FROM ubuntu:20.04
```

Here we choose to base our pipeline environment in Ubuntu Linux 20.04

```docker
# Dependencies for the R packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
libbz2-dev \
libcairo2-dev \
libcurl4-openssl-dev \
libfreetype6-dev \
libfribidi-dev \
libharfbuzz-dev \
libjpeg-dev \
libpng-dev \
libproj-dev \
libssl-dev \
libtiff5-dev \
libxml2-dev \
libxt-dev \
zlib1g-dev
```

Remember all the external library dependencies? We make sure they're installed.

```docker
# Compile dependencies for R and OpenBLAS.
RUN apt-get install -y --no-install-recommends build-essential \
cmake \
g++ \
gfortran \
make \
tk
```

And remember how what you use to compile things with is important? We take care of that, too.

```docker
# Install R.
RUN apt-get install -y --no-install-recommends r-base liblapack-dev

# Install OpenBLAS.
RUN apt-get install -y --no-install-recommends git
RUN git clone https://github.com/xianyi/OpenBLAS.git --branch v0.3.15 --depth=1
WORKDIR /OpenBLAS
RUN make install
```

Here we take care of the number-crunching software.

```docker
RUN R -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv')"

WORKDIR /project
COPY ./renv.lock .
# Install all dependencies based on the lock file.

RUN R -e "options(renv.consent = TRUE); options(timeout = 300); renv::restore()"
```

And, at long last, `renv` inside Docker for the win!
#### CI/CD

This project uses the [GitLab CI/CD](https://docs.gitlab.com/ee/ci/). Everything happens in the [`.gitlab-ci.yml`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/.gitlab-ci.yml) file:

```yaml
after_script:
  - R -e "sessionInfo()"

run:
  image: xlxs4/meta-analysis:latest
  script:
      - Rscript differential-expression-analysis/src/flocculation.R --qc --plots --feather --time

```

This simple set of directives ensures that each time code is pushed upstream, the script can run without crashing.
The script can run in the environment/container spawned by the Docker image, simplyfying things and also cutting down on CI/CD time.
There's a lot more that can be done here, but more on that later. Notice how almost everything has been abstracted away using the image.

### Path handling

We want to input and output a lot of files in our pipeline. Therefore, we need to handle paths in our code.
We want to mainly be able to:

1. Use relative instead of absolute, hard-wired paths (reproducibility)
2. Be able to set and change paths programmatically (reproducibility)
3. Have some sort of sanity checking when fiddling with paths

If only there was something like Python's [`pathlib`](https://docs.python.org/3/library/pathlib.html) for R...

But there is! And it's called [`here`]()! And [here](https://github.com/jennybc/here_here#tldr)'s (pun intended) some of the reasons why you should use it.

```r
library(here)
results_dir <- here(
    "differential-expression-analysis",
    "results",
    "flocculation"
)
plots_dir <- here(results_dir, "plots")
```

Boom.

### Logging

As mentioned [below](#input), [`logger`](https://daroczig.github.io/logger/) is used to log messages to the user. There are [numerous](https://daroczig.github.io/logger/#why-yet-another-logging-r-package) alternative R packages that can be used.

We can log to `stderr`:
```r
library(logger)

log_appender(appender_console)
```

Change the threshold level based on the levels from Apache [`Log4j`](https://logging.apache.org/log4j/2.x/), namely [these](https://daroczig.github.io/logger/reference/log_levels.html) ones.
```r
if (arguments$q) {
    log_threshold(SUCCESS)
}
```

Log messages using different log levels.
```r
log_info("Setting up paths...")

log_warn(capture.output(suppressMessages(arrayQualityMetrics(
            expressionset = cel_affybatch,
            outdir = here(qc_data_dir, "qc-report-affybatch"),
            force = TRUE,
            do.logtransform = TRUE
))))
```

Log messages in pretty print.
```r
if (!arguments$no_color) {
    log_layout(layout_glue_colors)
}
```

Control our messages programmatically.
```r
t <- tic.log()
log_success("The script finished running successfully! Time: {t}")
```

### Tibbles

All resulting dataframes meant to be part of the script output are first converted to [`tibble`](https://tibble.tidyverse.org/)s:

> A tibble, or tbl_df, is a modern reimagining of the data.frame, keeping what time has proven to be effective, and throwing out what is not. Tibbles are data.frames that are lazy and surly: they do less (i.e. they donโt change variable names or types, and donโt do partial matching) and complain more (e.g. when a variable does not exist). This forces you to confront problems earlier, typically leading to cleaner, more expressive code. Tibbles also have an enhanced print() method which makes them easier to use with large datasets containing complex objects.

```r
de_genes <- as_tibble(de_genes)
```

### Feather Files

The tibbles are then written on-disk as feather files. [Feather](https://github.com/wesm/feather) is a very fast, lightweight, and easy-to-use binary file format to read/write dataframes (tibbles, too!):

> Feather provides binary columnar serialization for data frames. It is designed to make reading and writing data frames efficient, and to make sharing data across data analysis languages easy.

> Feather uses the Apache Arrow columnar memory specification to represent binary data on disk.

```r
arrow::write_feather(de_genes, here::here(tibbles_dir, tibble_path))
```

I found `knitr::kable()` from `knitr` very helpful to convert the tibbles to markdown tables: `kable(df, "pipe")`.

### NASA GeneLab

[NASA Genelab](https://genelab.nasa.gov/) is an open-source, comprehensive space-related omics platform. To learn more about GeneLab you can begin from their [about](https://genelab.nasa.gov/about) page. Our analysis works on raw data from this platform. Information here is tidy, with most of the uploaded files having been processed by the GeneLab Team. There's some additional cool features, like the [visualization page](https://visualization.genelab.nasa.gov/data/GLDS-62).

### Original Publication

[1]: Goossens, K. V., Ielasi, F. S., Nookaew, I., Stals, I., Alonso-Sarduy, L., Daenen, L., ... & Willaert, R. G. (2015). Molecular mechanism of flocculation self-recognition in yeast and its role in mating and survival. MBio, 6(2), e00427-15.

TODO: @elsandal, give a brief description and mention why we can do our own DEA on these data maybe?

### GLDS-62

[`GLDS-62`](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-62/) is the GeneLab accession number for all the data related to [1]. From other potentially useful microarray data we came across, these were selected to analyze, for several reasons:

* The experiments were conducted using Affymetrix microarrays, specifically the Affymetrix Yeast Genome 2.0. We found data from these received greater support and were easier to work with than, say, the Hitachisoft Hypergene Yeast Chip Ver. 2.0
* How the data was processed and analyzed in the original publication is relatively clear and straightforward:
    > Raw .CEL files were read in and normalized using the r script 'affyNormQC.R' which utilizes the RMA algorithm through the oligo R package [rma() with default parameters]. Quality control reports were generated via the r script 'affyNormQC.R', with parameter 'do.logtransform' set to TRUE for the generating the raw report. This microarray experiment was annotated with the r script 'annotateProbes.R' which utilized Annotation-Db class probe annotations specific to each chip from the Bioconductor repository. In cases where multiple probes mapped to the same gene ID, representative probes were selected with the highest mean normalized intensity across all samples. Differential gene expression analysis was performed using the r script limmaDiffExp.R which utilizes the limma R package to perform pair-wise comparisons for all groups. For each probe set, the variance of mean signal intensities was estimated, improved by an empirical Bayes method for combining variances of probes showing similar variability, and the significance of the difference between the means was evaluated with a t-test to obtain p-values. P-values were adjusted for multiple hypothesis testing using the Benjamini and Hochberg method to control the false discovery rate.
* While the authors didn't provide the source code of their analysis, specific [GeneLab scripts](https://github.com/jdrubin91/GeneLab-Microarray/tree/master/GeneLab-Microarray/R_scripts) were used for all the pipeline stages
* Additional to the raw data and the DE genes list, the authors uploaded the data across the various analysis steps (preprocessing, normalization...), including their QC
* The experimental design was sound, and the data of adequate quality, with no significant outliers and good QC

### CEL Files

The `.CEL` files are produced by Affymetrix microarray platforms and contain experiment-specific information, such as intensity value and standard deviation, potential outliers as indicated by some preprocessing algorithms, etc. The microarray used in this study produces [Command Console](https://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html#calvin) CEL files.

### Bioconductor

[Bioconductor](https://bioconductor.org/about/) is a mostly R package repository centered on analyzing biological data (DNA microarray, sequence, SNP...). It is widely used in our pipeline to use various packages such as `affy`, `limma` and more.

## Script Rundown

### Input

[`flocculation.R`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/differential-expression-analysis/src/flocculation.R) uses [`docopt`](https://github.com/docopt/docopt.R) :books: to interface with the CLI and parse user input. You can read more on `docopt` [here](http://docopt.org/).

The complete interface for `flocculation.R` can be found inside the script:

```r
"DEA script for GLDS-62 GeneLab entry (raw data).

Usage:
  flocculation.R [-q | --no-color] [--time]
  flocculation.R (-h | --help)
  flocculation.R --version
  flocculation.R --qc [-r] [-n] [-t] [--plots] [-q | --no-color] [--feather] [--time]
  flocculation.R --plots [-q | --no-color] [--feather] [--time]
  flocculation.R --feather [--time]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --qc          Produce QC reports.
  --plots       Produce DEA plots.
  --feather     Save result tibbles as (arrow) feather files.
  --time        Output script execution time.

" -> doc
```

First, some `docopt` magic: :crystal_ball:

* `[ ]` square brackets indicate **optional elements**
* `( )` the parentheses mark **required elements**
* ` | ` the pipe operator is used to separate **mutually-exclusive elements**
* `-o` short options can be stacked: `-xyz โก -x -y -z`

So, you have a couple of different options:

* Call the script without any additional elements (`flocculation.R`). This will run the script beginning-to-end, excluding all code generating output, be it QC reports, various plots or `arrow:feather` files. Useful mainly to shorten execution time when developing, and to make profiling easier
* If you want to generate QC reports, you can set the `--qc` option. Note that each time the QC-generating code is called, new reports will be created regardless of whether there are already existing ones. To change this setting, open the script and remove the `force = TRUE` argument inside `arrayQualityMetrics()`. To speed up execution time, you can opt for running only a subset of the QC analyses, namely you can:
    * Additionally pass `-r` (`--qc -r`) to generate QC for the raw data (pre-RMA)
    * pass `-n` to generate QC for the normalized data (directly after RMA)
    * pass `-t` to generate QC for the statistical test results (after the `eBayes()` call on the linear model fit)
    * Example: `flocculation.R --qc -rt`
* To generate plots (as `.pdf`s), you can set the `--plots` option. Said plots include:
    * Various GSEA graphs
    * [Volcanoplot](https://rdrr.io/bioc/limma/man/volcanoplot.html)
    * [MD plot](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/plotMD) (log2 fold-change vs mean log2 expression)
    * [Venn diagram](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/venn.html)
    * [Heatmap](https://rdrr.io/bioc/limma/man/coolmap.html) (`coolmap`)
    * [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
* `flocculation.R` uses `tibble`s from [tidyverse](https://tibble.tidyverse.org/) to hold various dataframe-like results. You can have them be saved on-disk as `.feather` [files](https://arrow.apache.org/docs/python/feather.html) by passing `--feather`
* If you want to time the script call without using some external wrapper program, set the `--time` option
* [`logger`](https://daroczig.github.io/logger/) is used to log various events to the user. If you want to disable all `INFO`/`WARN` log messages, set the `-q` option
* By default, messages logged are colorful, using `glue` to format the message and ANSI escape codes, depending on the [`crayon`](https://github.com/r-lib/crayon) package. If that's not to your taste, feel free to use `--no-color`

### Reading in Files

First, the experimental data from the Affymetrix microarray are loaded by reading the [`.CEL` files](#cel-files). [`full.names`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/list.files.html) is required.

```r
cel_affybatch <- ReadAffy(filenames = list.celfiles(
    raw_data_dir,
    full.names = TRUE
))
```

All data from the `.txt`s are now loaded into a single object, of class [`AffyBatch`](https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/affy/html/AffyBatch-class.html). `AffyBatch` extends the [`eSet`](https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/Biobase/html/class.eSet.html) class. 

```
r$> cel_affybatch

AffyBatch object
size of arrays=496x496 features (25 kb)
cdf=Yeast_2 (10928 affyids)
number of samples=17
number of genes=10928
annotation=yeast2
notes=
```

An `eSet` and, by extension, an `AffyBatch` and an `ExpressionSet` (that we'll see [later on](#rma)), contains both the intensity readings and various metadata, as well as various methods to access and set the various `Slots`.

```
r$>str (cel_affybatch)                                                
Formal class 'AffyBatch' [package "affy"] with 10 slots
  ..@ cdfName          : chr "Yeast_2"
  ..@ nrow             : Named int 496
  ..@ ncol             : Named int 496
  ..@ assayData        :<environment: 0x4185ee98> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' 
  ..@ featureData      :Formal class 'AnnotatedDataFrame' 
  ..@ experimentData   :Formal class 'MIAME'
```

The Yeast Genome 2.0 Array [contains probe sets](https://www.affymetrix.com/support/technical/datasheets/yeast2_datasheet.pdf) to detect transcripts from both *Saccharomyces cerevisiae* and *Schizosaccharomyces pombe*. However, we are only interested in the *S. cerevisiae* data for our analysis. To filter out the *pombe* data, we can use a maskfile provided by [Affymetrix](http://www.affymetrix.com/Auth/support/downloads/maskfiles/scerevisiae.zip).

```r
s_cerevisiae_mask <- read.table(
    mask_data_dir,
    skip = 2,
    stringsAsFactors = FALSE
)
probe_filter <- s_cerevisiae_mask[[1]]
```

### Helper functions

Besides the [`flocculation.R`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/differential-expression-analysis/src/flocculation.R) main script, there's also [`helpers.R`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/differential-expression-analysis/src/helpers.R). Here, the two functions of note are `extract_ids()` and `remove_probes()`.

```r
extract_ids <- function(annotation_data_dir, probe_filter) {
    # Get both S. pombe & S. cerevisiae ids.
    genenames <- as.list(yeast2.db::yeast2GENENAME)
    probes <- names(genenames)

    # Get all transcript ids from yeast2annotation.csv.
    annotations <- read.csv(
        file = annotation_data_dir, header = TRUE,
        stringsAsFactors = FALSE
    )
    transcript_id <- annotations[, 3]
    probeset_id <- annotations[, 1]

    # Reorder the transcript_id to match probes.
    transcript_id <- transcript_id[match(probes, probeset_id)]

    # Retrieve the probeset and transcript ids for S. cerevisiae.
    c_probe_id <- probes[-match(probe_filter, probes)]
    c_transcript_id <- transcript_id[-match(probe_filter, probes)]

    # We need the TranscriptID if the gene name is NA.
    yeast_genenames <- transcript_id
    for (i in seq(along = probeset_id)) {
        gname <- genenames[i][[1]]
        if (!is.na(gname)) {
            yeast_genenames[i] <- gname
        }
    }

    # Set the gene name.
    c_genename <- yeast_genenames[-match(probe_filter, probes)]
    df <- tibble(
        probe = c_probe_id, transcript = c_transcript_id,
        genename = c_genename
    )
    return(df)
}
```

`remove_probes()` is more interesting:

```r
remove_probes <- function(list_out_probe_sets,
                          cdfpackagename,
                          probepackagename) {
    # list_out_probe_sets: Probe sets that are removed.

    # cdfpackagename: The cdf package name.
    # probepackagename: The probe package name.
    probe_env_org <- get(probepackagename)

    # Remove probesets from the CDF environment.
    rm(list = list_out_probe_sets, envir = get(cdfpackagename))

    # Set the PROBE env accordingly.
    tmp <- get("xy2indices", paste("package:", cdfpackagename, sep = ""))

    # cleancdf is from cdfpackagename... I think.
    new_affybatch <- new("AffyBatch", cdfName = cleancdf)
    pm_index <- unlist(affy::indexProbes(new_affybatch, "pm"))
    sub_index <- match(tmp(probe_env_org$x, probe_env_org$y,
        cdf = cdfpackagename
    ), pm_index)

    i_na <- which(is.na(sub_index))

    # Need to unlock the environment binding to alter the probes.
    ipos <- grep(probepackagename, search())
    env <- as.environment(search()[ipos])

    unlockBinding(probepackagename, env)
    assign(probepackagename, probe_env_org[-i_na, ], env = env)
    lockBinding(probepackagename, env)
}
```

We use `extract_ids()` to grab a *S. Cerevisiae*-only tibble and save it on-disk:

```r
s_cerevisiae_df <- extract_ids(annotation_data_dir, probe_filter)

if (arguments$feather) {
    log_info("Saving S.cerevisiae-only Yeast 2.0 probe sets tibble...")
    tibble_path <- "cerevisiae-only-probesets.feather"
    arrow::write_feather(s_cerevisiae_df, here::here(tibbles_dir, tibble_path))
}
```

Then, again using the mask, the `AffyBatch` object is edited in-place, in our current environment.

```r
remove_probes(probe_filter, cleancdf, probe_package_name)
```

### QC

**QC** stands for Quality Control and is a necessary step in every Differential Expression (and not only!) Analysis. It is the combination of all diferent methods used throughout the pipeline to gauge your data and ascertain that they are of adequate quality. More on what exactly that entails, and how it's integrated in our analysis as a pipeline step will be below (TODO). Here, let's take a look at the very simple code behind this:

We set up a path thanks to `here`:

```r
qc_data_dir <- here(
    "differential-expression-analysis",
    "qc",
    "flocculation"
)
```

And almost all of the heavy lifting is done under the scene by [`arrayQualityMetrics`](https://www.bioconductor.org/packages/release/bioc/html/arrayQualityMetrics.html). As is usual with R bioinformatics-related packages, everything is bundled in an all-seeing, aptly named wrapper function, `arrayQualityMetrics()`:

<p float="left">
  <img src="/differential-expression-analysis/qc/flocculation/qc-report-affybatch/msd.png" width="30%" />
  <img src="/differential-expression-analysis/qc/flocculation/qc-report-affybatch/ma.png" width="66%" /> 
</p>

The function calls generates a fully-fledged HTML QC report with various quality metrics about the dataset. Most of these metrics can be used to judge the quality of an _array_ across the whole sample spectrum. Also, some can be used to detect outliers/batch effects, giving you a better idea of the overall quality of the whole dataset:

```r
arrayQualityMetrics(
    expressionset = cel_affybatch,
    outdir = here(qc_data_dir, "qc-report-affybatch"),
    force = TRUE,
    do.logtransform = TRUE
)
```

The `force` argument is to ensure the report will be recreated, regardless of whether there already exists one on the output directory specified. `do.logtransform` makes sure to log transform the intensity values first. Note that when we call the function a second time, after [normalization](#rma), `do.logtransform` is set to `FALSE`, since thanks to the RMA algorithm we are now dealing with log-scale values. To browse the interactive report you can load the [`index.html`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/differential-expression-analysis/qc/flocculation/qc-report-affybatch/index.html) file.

Beside the report and the individual plots generated by `ArrayQualityMetrics`, there's some more QC plots generated:

1) A p-value distribution plot
2) An adjusted p-value distribution plot
3) And a Q-Q moderated t-statistics plot

### RMA

Again, not much is being done here. `rma()` is used from `affy` to normalize our data. The Robust Multichip Average expression measure is computed [2]. Note that the measure is in log2 scale, so `do.logtransform` is no longer needed in `arrayQualityMetrics()`. `rma()` converts our `AffyBatch` object into an [`ExpressionSet`](https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet), another extension of `eSet`. This `ExpressionSet` is what we're going to work with from now on.

```r
eset_rma <- rma(cel_affybatch, verbose = FALSE)
```

```
r$> eset_rma

ExpressionSet (storageMode: lockedEnvironment)
assayData: 5900 features, 17 samples 
  element names: exprs 
protocolData
  sampleNames: GSM1571870_hyb8616.CEL GSM1571871_hyb8617.CEL ... GSM1571886_hyb8637.CEL (17 total)
  varLabels: ScanDate
  varMetadata: labelDescription
phenoData
  sampleNames: GSM1571870_hyb8616.CEL GSM1571871_hyb8617.CEL ... GSM1571886_hyb8637.CEL (17 total)
  varLabels: sample
  varMetadata: labelDescription
featureData: none
experimentData: use 'experimentData(object)'
Annotation: yeast2
```

```
r$> str(eset_rma)

Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  ..@ assayData        :<environment: 0x3cf74d70> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  ..@ annotation       : chr "yeast2"
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
```

While we'll look into RMA (and other normalization methods) more in-depth later on, here is a quick visual (from `arrayQualityMetrics()`) of the effects the normalization has on our dataset:

<p float="left">
  <img src="/differential-expression-analysis/qc/flocculation/qc-report-affybatch/msd.png" width="49%" />
  <img src="/differential-expression-analysis/qc/flocculation/qc-report-rma/msd.png" width="49%" /> 
</p>

[2]: Irizarry, R. A., Hobbs, B., Collin, F., BeazerโBarclay, Y. D., Antonellis, K. J., Scherf, U., & Speed, T. P. (2003). Exploration, normalization, and summaries of high density oligonucleotide array probe level data. Biostatistics, 4(2), 249-264.

### Annotation

From [`affycoretools:annotateEset`](https://www.rdocumentation.org/packages/affycoretools/versions/1.44.2/topics/annotateEset):

> Annotating results is tedious, and can be surprisingly difficult to get right.

To annotate our `ExpressionSet`, we have to:

1) Find annotation data (ENTREZ ID, ENSEMBL ID, gene name, etc.) for our probes
2) Load the data in R
3) Merge the data with the `ExpressionSet`
4) Verify that the annotation was correctly carried out

The annotation data can be in SQL table(s) form, a CSV file, it can be fetched by an online database using their corresponding API, and more. Also, depending on the probes used, you might need to combine data from more than one sources to gather every thing needed. According to the form the data is in, how it is loaded and handled in R also has to change accordingly. To merge the data in R, it's usually one of (see [this](https://stackoverflow.com/questions/1299871/how-to-join-merge-data-frames-inner-outer-left-right)):

1) Outer join `merge(x = df1, y = df2, by = "CustomerId", all = TRUE)`
2) Left outer `merge(x = df1, y = df2, by = "CustomerId", all.x = TRUE)`
3) Right outer `merge(x = df1, y = df2, by = "CustomerId", all.y = TRUE)`
4) Cross join `merge(x = df1, y = df2, by = NULL)`

This is following standard [Relational Algebra](https://www.wikiwand.com/en/Relational_algebra).

Previously we [mentioned](#bioconductor) Bioconductor, which also provides "extensive annotation resources", that are either _gene centric_ or _genome centric_. For an intro into the various Bioconductor annotation packages you can check [here](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf). There are different kinds of packages, and each comes with a different kind of objects, for example:

1) the `AnnotationDb` objects
2) the `ChipDb` objects
3) the `OrgDb` objects
4) the `GODb` objects
5) the `TxDb` objects
6) the `EnsDb` objects

Another Bioconductor solution to this annotation problem comes in the form of a package: [`biomaRt`](https://bioconductor.org/packages/release/bioc/html/biomaRt.html):

> [...] biomaRt provides an interface to a growing collection of databases implementing the BioMart software suite (). The package enables retrieval of large amounts of data in a uniform way without the need to know the underlying database schemas or write complex SQL queries.

To annotate our _S. Cerevisiae_ probesets (note: probe _sets_, **not** probes; this will be important later on), we use:

1) [`yeast2.db`](https://www.bioconductor.org/packages/release/data/annotation/html/yeast2.db.html)
2) [`Yeast_2.na24.annot.csv`](https://gitlab.com/acubesat/su/bioinformatics/meta-analysis/-/blob/master/differential-expression-analysis/data/flocculation/annotation/Yeast_2.na24.annot.csv), as [provided](http://www.affymetrix.com/Auth/analysis/downloads/na24/ivt/Yeast2.na24.annot.csv.zip) by Affymetrix (requires user login/registration)
3) `biomaRt`
4) `AnnotationDbi`

However, still it is a tedious, prone to error process. Helpfully, this is where `annotateEset` comes into play:

```r
annotateEset(eset_rma, yeast2.db, columns = c("PROBEID", "ENSEMBL", "GENENAME"))
```

### Removing probesets

We must take care of some lingering probesets. The Affymetrix platforms include two types of probesets that can be removed post-normalization; the `AFFX` and the `RPTR`-labeled probesets, respectively. The `AFFX` probesets are control probesets, while the `RPTR` have to do with reporter genes:

```r
control_affymetrix <- grep("AFFX", featureNames(eset_final))
eset_final <- eset_final[-control_affymetrix, ]

control_reporter_genes <- grep("RPTR", featureNames(eset_final))
eset_final <- eset_final[-control_reporter_genes, ]
```

Since this analysis is conducted to help select the genes the expression of which will be probed on-board our nanosatellite, all probesets that do not map to an ENSEMBL/ENTREZ ID (and therefore do not map to a distinct `GENENAME`) are also removed from our study:

```r
no_ensembl_ids <- is.na(fData(eset_final)$ENSEMBL)
eset_final <- eset_final[!no_ensembl_ids, ]
```

### Selecting representative probes

Another somewhat common occurence when dealing with microarray data is coming across probesets that get annotated to multiple gene identifiers.
This might be confusing to some, and I find that there is an helpful distinction to be made (also highlighted above, remember?).
At least when working with Affymetrix arrays, each ID corresponds to a probe _set_, **not** an individual probe. Thus, each ID is a bundle of individual probes. See [here](https://www.reddit.com/r/bioinformatics/comments/544zqi/multiple_probes_for_one_gene/d7zpcxj) for an explanatory post, and [here](https://www.affymetrix.com/support/help/faqs/mouse_430/faq_8.affx) to better understand the different Affymetrix probe sets and suffixes in the probe IDs.
Before conducting a DEA, it is best to select some representative probes, if possible. See [this](https://www.biostars.org/p/47421/) Biostars post for more information on the matter and alternative approaches.

### Creating the design matrix

We'll go over Linear Regression Models and DEA below.
How these are implemented in R is very straightforward, beginning with creating the design matrix.

First, we need to label our microarray samples with the respective group they belong in, as per the experimental design. In our case, this can be found in the [respective GEO (Gene Expression Omnibus) entry](https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE64468). Since we are only interested in comparing between the on-ground control samples and the ones exposed to the simulated low-shear microgravity environmental conditions, we don't distinguish between wildtype, FLO1 and FLO8 samples. You can find more on the `factor` object [here](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/factor).

```r
group_membership_ground <- "00011100011000111"
sml <- base::strsplit(group_membership_ground, split = "")[[1]]

gs <- factor(sml)
```

```
r$> gs
[1] 0 0 0 1 1 1 0 0 0 1 1 0 0 0 1 1 1
Levels: 0 1
```

We label for each sample (total of 17):

```r
groups <- make.names(c("onground", "micro"))

levels(gs) <- groups
eset_final$group <- gs
```

```
r$> gs  
 [1] onground onground onground micro    micro    micro    onground onground onground micro    micro    onground onground onground micro    micro    micro   
Levels: onground micro
```

After this preperation, we can create our independent t-test design matrix by calling `model.matrix()`:

```r
design_matrix <- model.matrix(~ group + 0, eset_final)
colnames(design_matrix) <- levels(gs)
```

We get this design matrix:

```
r$> design_matrix    
                       onground micro
GSM1571870_hyb8616.CEL        1     0
GSM1571871_hyb8617.CEL        1     0
GSM1571872_hyb8620.CEL        1     0
GSM1571873_hyb8618.CEL        0     1
GSM1571874_hyb8621.CEL        0     1
GSM1571875_hyb8622.CEL        0     1
GSM1571876_hyb8624.CEL        1     0
...
```

Which corresponds to the following linear equation:

`y = mean(on ground) + mean(micro)`

The design matrix is essentially a incident matrix with `1` if the sample was subjected in simulated microgravity conditions and `0` otherwise. By default, `model.matrix()` includes a column of all 1's representing $`\mu`$ in the ANOVA model $`Y_{ij} = \mu + \alpha_{i} + \mathrm{error}`$ (see [this](https://www.bioconductor.org/packages/release/data/experiment/vignettes/ChimpHumanBrainData/inst/doc/DiffExpressVignette.pdf)).

### Fitting the linear models

`limma::lmFit()` is used to fit the multiple linear models, according to the design matrix we generated.
`lmFit()` fits linear models with generalized or weighted least squares. The `method=robust` argument can be passed to have robust regression instead (see the [man](https://bioconductor.org/packages/devel/bioc/manuals/limma/man/limma.pdf) for more info).
As stated in [this BioC thread](https://support.bioconductor.org/p/25749/), we opt for least squares, since there is a low number of replicates (triplicates and a duplicate) and we want to avoid removing real variation.

```r
fit <- lmFit(eset_final, design_matrix)
```

### Creating the contrast matrix

See this [StackExchange thread](https://stats.stackexchange.com/questions/78354/what-is-a-contrast-matrix) for some excellent answers on what a contrast matrix is. In our case, all we need to do is create our contrast matrix (control vs microgravity), recaluclate the model coefficients and re-orientate the previously fitted model to the set of contrasts of the original coefficients:

```r
contrast <- paste(groups[1], groups[2], sep = "-")
contrast_matrix <- makeContrasts(contrasts = contrast, levels = design_matrix)

fit2 <- contrasts.fit(fit, contrast_matrix)
```

### Empirical Bayes

We can employ Empirical Bayes (EB) statistical tests that use moderated genewise variances through [`limma::eBayes()`](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/ebayes.html), which ranks genes in order of evidence for differential expression, by:
1) using an empirical Bayes method to shrink the probe-wise sample variances towards a common value and
2) augmenting the degrees of freedom for the individual variances

See [3] for more. Additionally, the `robust = TRUE` argument is passed, to robustify the hyperparameter estimation as seen in [4]. This deals with outlier genes (genes with very large or small variances) greatly affecting the Empirical Bayes statistical tests.

```r
fit_eb <- eBayes(fit2, robust = TRUE)
```

[3]: Smyth, G. K. (2004). Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Statistical applications in genetics and molecular biology, 3(1).

[4]: Phipson, B., Lee, S., Majewski, I. J., Alexander, W. S., & Smyth, G. K. (2016). Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. The annals of applied statistics, 10(2), 946.

### Selecting DE genes

`fit_eb` now contains everything our `fit2` linear model fit comes with, plus some added components, such as `t`, the matrix of moderated t-statistics, and `p.value`, the p-values corresponding to the t-statistics.
`limma::topTable()` uses these statistics to extract the top-ranked genes from our linear model fit:

```r
de_genes <- topTable(fit_eb,
    number = 30,
    adjust.method = "BH",
    sort.by = "B",
    p.value = pval_cutoff
)
```

* `adjust.method` sets the method used by [stats::p.adjust()](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust) to adjust the p-values from `eBayes()`. The commonly used Benjamini & Hochberg [5] method is selected here. The Benjamini & Hochberg method controls the FDR, a less stringent condition than the family-wise error rate, so it's more powerful than the other methods available except Benjamini & Yekutieli [6]
For an intro to `BH` and `BY` see [these slides](http://www.stat.cmu.edu/~genovese/talks/hannover1-04.pdf) by CMU's Christopher R. Genovese.
* The criterion used to select the top genes is the B-statistic, the log-odds that the gene is differentially expressed (see `limma`'s [Users Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf))
* Lastly, we filter our DE genes according to an adjusted p-value threshold. An absolute log2 fold-change cutoff was not used, since if the fold changes and the p-values are not highly correlated, the use of a fold-change cutoff on top of a p-value cutoff can increase the FDR above the nominal level. See the [limma reference manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).

We end up with a `tibble` with a row for each gene considered differentially expressed and with the following columns:

```r
de_genes <- subset(de_genes,
    select = c(
        "PROBEID",
        "adj.P.Val",
        "P.Value",
        "t",
        "B",
        "logFC",
        "ENSEMBL",
        "GENENAME"
    )
)

de_genes <- as_tibble(de_genes)
```

```
r$> de_genes     
# A tibble: 242 ร 8
   PROBEID      adj.P.Val      P.Value     t     B logFC ENSEMBL GENENAME
   <chr>            <dbl>        <dbl> <dbl> <dbl> <dbl> <chr>   <chr>   
 1 1771832_at    0.000142 0.0000000251  9.21  9.06 0.865 YDR418W RPL12A  
 2 1776010_at    0.000157 0.0000000555  8.76  8.36 1.04  YJL190C RPS22A  
 3 1777989_at    0.000617 0.000000327   7.76  6.77 0.999 YLR406C RPL31B  
 4 1773591_at    0.000712 0.000000503   7.53  6.38 1.04  YMR199W CLN1    
 5 1772366_at    0.00153  0.00000135    7.01  5.48 1.20  YML063W RPS1B   
 6 1776927_at    0.00197  0.00000242    6.72  4.95 0.828 YGL147C RPL9A   
 7 1773292_at    0.00197  0.00000248    6.69  4.92 0.638 YGR034W RPL26B  
 8 1777846_at    0.00197  0.00000337    6.54  4.64 0.672 YBL002W HTB2    
 9 1773827_s_at  0.00197  0.00000344    6.53  4.62 0.577 YHL033C RPL8A   
10 1776763_at    0.00197  0.00000379    6.49  4.53 0.702 YPL256C CLN2    
# โฆ with 232 more rows
```

[5]: Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal statistical society: series B (Methodological), 57(1), 289-300.

[6]: Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of statistics, 1165-1188.

#### decideTests

`limma` also provides [`limma::decideTests()`](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/decideTests). `decideTests` also takes our `fit_eb` `eBayes()` fit and identifies significantly differentially expressed genes:

```r
results <- decideTests(fit_eb,
    adjust.method = "BH",
    p.value = pval_cutoff,
    method = "nestedF"
)
```

A multiple testing procedure and a significance level cutoff is applied to the statistics from `fit_eb` to determine whether each statistic should be considered significantly different from zero.

The Benjamini & Hochberg [5] is also used here to control the FDR.
The `nestedF` method is selected, which seems to achieve better consistensy between contrasts than the other available methods (first used by [7]). `nestedF` provides formal error rate control at the gene level but not for individual contrasts.

The `nestedF` multiple testing option tests across contrasts for each gene individually, and, when called through `decideTests()`, a multiple testing adjustment across genes is also incorporated.
A nested F-test approach is used with a focus on correctly classifying genes that have two or more significant t-statistics (i.e. are differentially expressed in two or more conditions). See [`limma:classifyTestsF`](https://rdrr.io/bioc/limma/man/classifytestsF.html) for more.

`decideTests` returns an object that belongs to the `TestResults` class; essentially what is a numeric matrix with elements `-1`, `0`, or `1`, depending on whether each t-statistic is classified as significantly negative, not significant, or significantly positive. We can use `summary()`:

```
r$> summary(results)

       onground-micro
Down              122
NotSig           5416
Up                120
```

Again, as is the case with `topTable`, an absolute log2 fold-change cutoff was not used, since if the fold changes and the p-values are not highly correlated, the use of a fold-change cutoff on top of a p-value cutoff can increase the FDR or family-wise error rate above the nominal level.

[7]: Michaud, J., Simpson, K. M., Escher, R., Buchet-Poyau, K., Beissbarth, T., Carmichael, C., ... & Scott, H. S. (2008). Integrative analysis of RUNX1 downstream pathways and target genes. BMC genomics, 9(1), 1-17.

#### topTreat

Because we would like to also run the analysis using fold change thresholding, we also use `treat` instead of `eBayes`:
```r
fit_treat <- treat(fit, fc = 1.1, robust = TRUE)
```

From the `limma` reference guide:
> treat computes empirical Bayes moderated-t p-values relative to a minimum fold-change threshold. Instead of testing for genes that have true log-fold-changes different from zero, it tests whether the true log2-fold-change is greater than lfc in absolute value [8]. In other words, it uses an interval null hypothesis, where the interval is [-lfc,lfc]. When the number of DE genes is large, treat is often useful for giving preference to larger fold-changes and for prioritizing genes that are biologically important. treat is concerned with p-values rather than posterior odds, so it does not compute the B-statistic lods. The idea of thresholding doesnโt apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed. When fc=1 and lfc=0, treat is identical to eBayes, except that F-statistics and B-statistics are not computed. The fc threshold is usually chosen relatively small, because genes need to have fold changes substantially greater than the testing threshold in order to be considered statistically significant. Typical values for fc are 1.1, 1.2 or 1.5.

Now all we have to do is use `topTreat` instead of `topTable` to summarize the linear fit produced by `lmfit` and processed by `treat`:
```r
de_genes_treat <- topTreat(fit_treat)
```

[8]: McCarthy, D. J., & Smyth, G. K. (2009). Testing significance relative to a fold-change threshold is a TREAT. Bioinformatics, 25(6), 765-771.

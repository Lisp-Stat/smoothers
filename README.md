

<!-- PROJECT SHIELDS -->

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MS-PL License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/lisp-stat/lisp-stat">
    <img src="https://lisp-stat.dev/images/stats-image.svg" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">Smoothers in Common Lisp</h3>

  <p align="center">
	Data smoothers for Common Lisp
	<br />
    <a href="https://lisp-stat.dev/docs/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/lisp-stat/smoothers/issues">Report Bug</a>
    ·
    <a href="https://github.com/lisp-stat/smoothers/issues">Request Feature</a>
    ·
    <a href="https://lisp-stat.github.io/smoothers/">Reference Manual</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About the Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
	<li><a href="#resources">Resources</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About the Project

Nonparametric regression offers a flexible alternative to classic (parametric) methods for regression.  Unlike parametric methods, which assume that the regression relationship has a known form that depends on a finite number of unknown parameters, nonparametric regression models attempt to learn the form of the regression relationship from a sample of data.

All nonparametric regression models involve finding some balance between fitting the observed sample of data (model fit) and “smoothing” the function estimate (model parsimony).  Typically, this balance is determined using some form of cross-validation, which attempts to find a function estimate that does well for predicting new data.  As a result, nonparametric regression models can be useful for discovering relationships between variables, as well as for developing generalizable prediction rules.


### Built With

* [lla](https://github.com/Lisp-Stat/lla)
* [numerical-utilities](https://github.com/Lisp-Stat/numerical-utilities)
* [alexandria](https://gitlab.common-lisp.net/alexandria/alexandria)
* [alexandria+](https://github.com/Symbolics/alexandria-plus)
* [array-operations](https://github.com/Lisp-Stat/array-operations)
* [let-plus](https://github.com/sharplispers/let-plus)


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these steps:

### Prerequisites

An ANSI Common Lisp implementation. Developed and tested with [SBCL](https://www.sbcl.org/).

### Installation

#### Getting the source

To make the system accessible to [ASDF](https://common-lisp.net/project/asdf/) (a build facility, similar to `make` in the C world), clone the repository in a directory ASDF knows about.  By default the `common-lisp` directory in your home directory is known.  Create this if it doesn't already exist and then:

1. Clone the repositories
```sh
cd ~/common-lisp && \
git clone https://github.com/Lisp-Stat/smoothers
```
2. Reset the ASDF source-registry to find the new system (from the REPL)
   ```lisp
   (asdf:clear-source-registry)
   ```
3. Load the system
   ```lisp
   (asdf:load-system :smoothers)
   ```

If you have installed the slime ASDF extensions, you can invoke this with a comma (',') from the slime REPL.

#### Getting dependencies

To get the third party systems that `smoothers` depends on, you can use a dependency manager, such as [Quicklisp](https://www.quicklisp.org/beta/) or [CLPM](https://www.clpm.dev/). Once installed, get the dependencies with either of:

```lisp
(clpm-client:sync :sources "smoothers") ;sources may vary
```

```lisp
(ql:quickload :smoothers)
```

You need do this only once. After obtaining the dependencies, you can
load the system with `ASDF` as described above without first syncing
sources.


<!-- USAGE EXAMPLES -->
<!--
## Usage

Create a data frame from a file named `sg-weather.csv` on the local disk:

```lisp
(defparameter *df*
	(read-csv #P"LS:DATA;sg-weather.csv"))

```

For more examples, please refer to the
[Documentation](https://lisp-stat.dev/docs/).
-->

## References

[Nonparametric Regression (Smoothers) in R](http://users.stat.umn.edu/~helwig/notes/smooth-notes.html)

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/lisp-stat/smoothers/issues) for a list of proposed features (and known issues).

## Resources

This system is part of the [Lisp-Stat](https://lisp-stat.dev/) project; that should be your first stop for information. Also see the [resources](https://lisp-stat.dev/docs/resources) and [community](https://lisp-stat.dev/community) page for more information.

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create.  Any contributions you make are **greatly appreciated**.  Please see [CONTRIBUTING](CONTRIBUTING.md) for details on the code of conduct, and the process for submitting pull requests.  Good first issues include porting the [xlisp-stat super-smoother](https://github.com/Lisp-Stat/xls-archive/tree/master/statistics/smoothers/supersmoother) to Common Lisp and bringing in Gary King's [smoothing.lisp](https://github.com/hraban/cl-mathstats/blob/master/dev/smoothing.lisp).

<!-- LICENSE -->
## License

Distributed under the MS-PL License. See [LICENSE](LICENSE) for more information.



<!-- CONTACT -->
## Contact

Project Link: [https://github.com/lisp-stat/smootherfs](https://github.com/lisp-stat/smoothers)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/lisp-stat/smoothers.svg?style=for-the-badge
[contributors-url]: https://github.com/lisp-stat/smoothers/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/lisp-stat/smoothers.svg?style=for-the-badge
[forks-url]: https://github.com/lisp-stat/smoothers/network/members
[stars-shield]: https://img.shields.io/github/stars/lisp-stat/smoothers.svg?style=for-the-badge
[stars-url]: https://github.com/lisp-stat/smoothers/stargazers
[issues-shield]: https://img.shields.io/github/issues/lisp-stat/smoothers.svg?style=for-the-badge
[issues-url]: https://github.com/lisp-stat/smoothers/issues
[license-shield]: https://img.shields.io/github/license/lisp-stat/smoothers.svg?style=for-the-badge
[license-url]: https://github.com/lisp-stat/smoothers/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/company/symbolics/

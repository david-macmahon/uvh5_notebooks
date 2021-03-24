### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 782d9788-77da-11eb-357f-e3857ed0f0ec
begin
	ARTIFACT_DIR = Dict{String,String}()
	artifact(aname,fname="") = joinpath(ARTIFACT_DIR[aname], fname)

	let
		deps = [
			(name="PlutoUI",  version="0.7"),
			(name="HDF5",     version="0.15"),
			(name="Plots",    version="1.11"),
			(name="PlotlyJS", version="0.14"),
			(url="https://github.com/david-macmahon/RadioInterferometry.jl",
			 rev="v0.2.0")
		]
		# Activate temporary environment and add deps
		import Pkg, Base.SHA1
		Pkg.activate(mktempdir())
		Pkg.add(deps)
		# Define artifacts (inline version of Artifacts.toml)
		ARTIFACTS = Dict{String, Any}()
		ARTIFACTS["uvh5_examples"] = Dict(
			"git-tree-sha1" => SHA1("618542fbe9049cad22b653de4804e83cd64c9b64"),
			"download" => [
				Dict(
"url" => "http://astro.berkeley.edu/~davidm/meerkat_example.uvh5.tar.xz",
"sha256" => "27372c8f28c4fb98513ec02a77fa0a4cf4bba69ac34fc8283035a16c3e9c7c21"						)
			]
		)
		# Ensure artifacts are downloaded and valid
		for (name, meta) in ARTIFACTS
			if !Pkg.Artifacts.verify_artifact(meta["git-tree-sha1"])
				for d in meta["download"]
					Pkg.Artifacts.download_artifact(
						meta["git-tree-sha1"], d["url"], d["sha256"]
					) && break
				end
			end
			if !Pkg.Artifacts.artifact_exists(meta["git-tree-sha1"])
				throw(ErrorException("unable to resolve artifact $name"))
			end
			ARTIFACT_DIR[name] = Pkg.Artifacts.artifact_path(meta["git-tree-sha1"])
		end
	end

	using HDF5, Plots, PlutoUI, Dates, RadioInterferometry
	plotlyjs(legend=false)
	HDF5.SHOW_TREE_MAX_CHILDREN[] = 100

	md"""
	## UVH5 explorations
	"""
end

# ╔═╡ 440400b0-7be9-11eb-0efb-c723129eeef2
md"""
The UVH5 file format is used to store radio interferometer visibility data in HDF5 data files.  A UVH5 file is an HDF5 data file that can be read using standard HDF5 libraries, though higher level libraries can make some tasks easier.  This notebook uses Julia and its HDF5 package to explore the UVH5 file format with a test UVH5 data file from BLUSE@MeerKAT as well as describe some basic concepts related to radio interferometry in general.  The UVH5 format is described in the [UVH5 file format memo](https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/main/docs/references/uvh5_memo.pdf).

The UVH5 file format was developed by the *Radio Astronomy Software Group* comprised of radio astronomers and cosmologists working with radio telescope arrays.  While several existing file formats already exist for storing visibility data (e.g. *MIRIAD*, *UVFITS*, *Measurement Sets*, *FHD*), their underlying file formats are older and/or highly specialized, making them harder to work with than modern file formats such as HDF5.  The Radio Astronomy Software Group provides [software](https://github.com/RadioAstronomySoftwareGroup/pyuvdata) that can be used to interact with UVH5 data files as well as convert between UVH5 and these other file formats.
"""

# ╔═╡ 2c9addb8-7be9-11eb-0adf-452042d88a4c
md"""
!!! note "Disclaimer"
    Being an introductory level notebook, some of the more advanced (and even some basic) features of UVH5 will be glossed over, ignored outright, or deferred to future notebooks.  Some of the techniques used here take advantage of *a priori* knowledge of the conventions used when creating the example UVH5 file being used (e.g. there is only one common integration time for all baselines) and are not generally applicable to more complicated but equally legitimate UVH5 file layouts.
"""

# ╔═╡ 99f47c64-78df-11eb-0ee3-af29fbc46edf
md"""
### Getting started
The `H5` in `UVH5` stands for `HDF5` and the `H ` in `HDF5` stands for *hierarchical*.  After opening a UVH5 file via the standard `h5open()` call, Julia (by default) shows the structure of the HDF5 file contents in a familiar "file explorer" hierarchical display.  In the next notebook cell we open the example UVH5 file named `meerkat_example.uvh5` and assign the returned object to the `uvh5` variable for later use.

!!! note
    This notebook is one of a series of notebooks that use the same example UVH5 file.  These notebooks automatically downloaded the example file during notebook initializtion using Julia's *artifacts* functionality.  This is the same mechanism used to download various components of Julia packages.  The example UVH5 file will only be downloaded the first time you run any of the notebooks in the series.  After that they will all use the same local copy.  The `artifact()` function defined in this notebook can used to obtain the local path to the example UVH5 file.  The exact location on your system can be seen on the `HDF5.File` line of the UVH5 file structure shown below.
"""

# ╔═╡ 83191f96-77da-11eb-16bb-8bf06f6b60ed

uvh5 = h5open(artifact("uvh5_examples", "meerkat_example.uvh5"))

# ╔═╡ b05be61e-77da-11eb-3a82-edc1be04b007
begin
	# Set some convenience variables
	visdata = uvh5["Data/visdata"]::HDF5.Dataset
	freqs_mhz = uvh5["Header/freq_array"][:,1] ./ 1e6
	B = uvh5["Header/channel_width"][]
	τ = uvh5["Header/integration_time"][1]
	Bτ = B*τ

	md"""
	As you can see, the UVH5 file contains two HDF5 *groups*: `Data` and `Header`.  The `Data` group contains three HDF5 *datasets* holding the visibility data (`visdata`), validity flags (`flags`), and fractional sample weights (`nsamples`). The `Header` group contains not a small number of HDF5 datasets that store various metadata required to properly interpret the visibility data.
	"""
end

# ╔═╡ d893adaa-7afb-11eb-143d-1fa4d8054566
let
	hdr = uvh5["Header"]
	telescope = hdr["telescope_name"][]
	instrument = hdr["instrument"][]
	latitude = d2dmsstr(hdr["latitude"][])
	longitude = d2dmsstr(hdr["longitude"][])
	altitude = hdr["altitude"][]
	srcname = hdr["object_name"][]
	ra  = hdr["phase_center_ra"][]/15  |> rad2deg |> h2hmsstr
	dec = hdr["phase_center_dec"][] |> rad2deg |> d2dmsstr
	obsdate = julian2datetime(hdr["time_array"][1])
	lst = h2hmsstr(rad2deg(hdr["lst_array"][1])/15)
	nants_data = hdr["Nants_data"][]
	nants_telescope = hdr["Nants_telescope"][]
	fc_mhz = (freqs_mhz[1] + freqs_mhz[end]) / 2
	nchan = hdr["Nfreqs"][]
	chanbw_mhz = hdr["channel_width"][] / 1e6
	bw_mhz = chanbw_mhz * nchan
	inttime = hdr["integration_time"][1]
	md"""
	### Observational metadata (UVH5 Header)
	Here is a simple table summarizing some of the details of the observational metadata contained in the test dataset's `Header` datasets:
	
	|Property          |Value                            |
	|:-----------------|:--------------------------------|
	| Telescope        | $telescope                      |
	| Instrument       | $instrument                     |
	| Latitude         | $latitude                       |
	| Longitude        | $longitude                      |
	| Altitude         | $altitude m                     |
	| Source           | $srcname                        |
	| RA               | $ra                             |
	| Dec              | $dec                            |
	| Date             | $obsdate                        |
	| LST              | $lst                            |
	| Antennas         | $nants_data of $nants_telescope |
	| Center frequency | $fc_mhz MHz                     |
	| Bandwidth        | $bw_mhz MHz                     |
	| Nchan            | $nchan                          |
	| Channel width    | $chanbw_mhz MHz                 |
	| Integration time | $inttime s                      |
	"""
end

# ╔═╡ 42167ea4-7b26-11eb-180b-138a95419a71
md"""
### Visibility data in a nutshell (UVH5 Data)
The bulk of the data in a UVH5 file is in the `Data` datasets, especially the `Data/visdata` dataset, which contains the so called *interferometric visibility data*, often referred to as just *visibilities*.  This term originates from the relationship these values have with the spatial distribution of RF energy that is incident on the radio telescope array.  We will not be going into that type of analysis too deeply in this notebook.  We will focus primarily on the computational genesis and statistical properties of these data, but we will cover some some basic interpretive concepts once we get to plotting data (aka the fun stuff).
"""

# ╔═╡ 6f2807e6-7b2b-11eb-0844-39b3e53b116c
md"""
###### Where do visibilities come from?
Visibilities come from a device known as a *correlator*.  Correlators have tended to use high performance digital electronics customized for a specific radio telescope array to match optimally the capabilities of the array, but increasingly correlators are starting to resemble HPC compute clusters.  Glossing over decades of correlator design history, correlators are essentially computers that compute the visibilities as the integrated products of pairs of input signals, one of which is conjugated when the signals are complex.  For complex input signals $Z$ from an array with ``A`` antennas, this can be written as:

```math
\mathrm{visibility}_{ij} = \sum_{n=1}^{N}(Z_i^{}(n) \cdot Z_j^*(n)),\ \mathrm{where}\ i,j \in [1, A]
```

Despite mathematically resembling a dot product more than a correlation function, a visibility is often referred to as a *correlation*.  When ``i = j``, it is called an *auto-correlation*; when ``i \neq j`` it is called a *cross-correlation*.  The summation over ``N`` terms is more commonly referred to in radio astronomy as an *integration*.  The choice of ``N`` is somewhat arbitrary.  Generally, ``N`` is chosen to be as large as possible, but the changing phase difference between the two inputs, when ``i \neq j``, tends to provide an upper bound.  We won't be getting into the cause or interpretation of the changing phase difference in this notebook, but it is worth mentioning that the cross-correlation phase is the difference of the phases of the two complex inputs.  This can be seen by expressing the complex inputs in polar form and recalling that complex conjugation negates the phase of the input.

```math
\begin{align}
Z_i &= V_i \enclose{phasorangle}{\theta_i} \\
Z_j &= V_j \enclose{phasorangle}{\theta_j} \\
Z_i^{} \cdot Z_j^* &= V_i \enclose{phasorangle}{\theta_i} \cdot V_j \enclose{phasorangle}{-\theta_j} \\
&= V_i V_j \enclose{phasorangle}{\theta_i - \theta_j}
\end{align}
```

!!! note "Tip"
    The phases of the cross-correlations carry much of the information in visibility data.

###### Correlator architectures
The visibility formula shown above is a *wideband* correlation function.  It provides the correlation of the full bandwidth of input ``Z``.  Two different but equivalent approaches can be used to gain insights into the frequency structure of the correlation.  One approach is to compute correlation of the two inputs at various *lags* (i.e. with signal delayed relative to the other) and then compute the FFT of the resultant *lag spectrum* to produce the *frequency spectrum* of the correlation.  This is known as an *XF* architecture because the correlation operation (*X*) happens before the frequency spectrum (*F*) is computed.  The other approach is to compute the frequency spectrum of the input signal first and then compute the wideband correlations for each frequency channel.  This in known as an *FX* architecture because the frequency spectrum is computed before the correlation operation happens.  Although the implementations of these architectures are quite different, the output visibilities are equivalent.

###### Baselines
Each pair of signals corresponds to a pair of input antennas (including "self pairs" where the two signals correspond to the same antenna).  The antennas of each pair form a *baseline vector* (often just *baseline*) from one antenna position to the other antenna position.  Auto-correlation baselines always have zero length.

The total number of baselines would appear to be ``A^2``, but half the cross-correlation baseline vectors are simple opposite signed to others (and the corresponding visibilities are conjugate redundant since ``\mathrm{visibility}_{ij}^{} = \mathrm{visibility}_{ji}^*``), so the total number of unique baselines is ``\frac{A(A+1)}{2} = A`` auto-correlation baselines plus ``\frac{A(A-1)}{2}`` cross-correlation baselines.
"""

# ╔═╡ 906c6b04-7b35-11eb-1ba6-110b111b728e
let
	nbls = uvh5["Header/Nbls"][]
	nants_data = uvh5["Header/Nants_data"][]
	Markdown.parse("""
	In a UVH5 file, the quantity ``A`` is available in the `Header/Nants_data` dataset ($(nants_data))
	and the number of unique baselines is available in `Header/Nbls` ($(nbls))
	or can be computed as `Nants_data*(Nants_data+1)÷2` ($(nants_data*(nants_data+1)÷2)).
	""")
end

# ╔═╡ d3a8122c-7b3c-11eb-1bc3-51a1da2515a3
md"""
###### Zero mean signals
We assume that the signals ``Z`` have have *zero mean*.  In other words, their expected value ``E(Z) = 0``.  This is generally a valid assumption when dealing with radio astronomy data, especially if the input signals are channelized by a *Fast Fourier Transform* (FFT) or *Polyphase Filter Bank* (PFB) prior to correlation.  This property of the input signals simplifies the calculation of various statistical properties such as the correlation coefficients described in a later section.

Zero mean signals have another very important property.  The correlation of uncorrelated zero mean signals grows by ``\sqrt{N}`` whereas the correlation of correlated zero mean signals grows by ``N``.  In other words, the correlation of correlated signals (e.g. a radio wave from space arriving at different radio telescopes) grows faster than the correlation of uncorrelated signals (e.g. independent noise from different radio telescope receivers) by a factor of ``\sqrt{N}``.  This is what makes correlators so very sensitive.
"""

# ╔═╡ b79230c6-7f93-11eb-046d-330e79e39964
md"""
### Dimensionality of the "Data/visdata" dataset

In a UVH5 data file, the `Data/visdata` dataset contains the actual visibility data.  This dataset has four dimensions.  Before describing what the four dimensions are,
we first need to cover a small detail regarding the order of dimensions.

!!! note "Note about order of dimensions"
    Julia uses *column major ordering* like Fortran and Matlab rather than *row major ordering* like C and Python.  Julia also indexes into multidimensional arrays with the first index changing the fastest like Fortran and Matlab rather the with the last index changing fastest like C and Python.  The effect of this is that the order of dimensions in Julia/Fortran/Matlab will be reversed from the order of dimensions used in C/Python (and the HDF5 command line tools like `h5ls` and `h5dump`).  The layout of data in memory and on disk remains the same, just the dimensions are flipped end to end.  Because this is a Julia notebook, the order of dimensions presented and used here will follow the column major ordering convention.

The dimensions of the `Data/visdata` dataset are ordered in Julia (Fortran/Matlab) as:

```math
(N_{pols}, N_{freqs}, N_{spws}, N_{blts})
```

with ``N_{pols}`` changing the fastest.  Each of these dimensions is described in one of the following sections.  The elements of this dataset are complex numbers, but the real and imaginary components can be an integer or floating point types of various sizes.  The real and imaginary components must be of the same type.  Because HDF5 datasets are self describing, the Julia HDF5 package is able to return the appropriate Julia type for a given dataset automatically.  Julia types for commonly used for visibility data types are `Complex{Float32}` and `Complex{Int32}`.  Care must be used when performing computations on integer data types to avoid overflow which will silently "wrap around" giving modulo results that may not be desired.
"""

# ╔═╡ 53e0fc08-7f8c-11eb-3f3f-e7c01518e7d7
md"""
###### Polarization products

Radio waves from space are sometimes polarized.  The polarization properties of an object's radio emission provides additional information about the object's intrinsic properties and surrounding environment.  Many radio telescopes contain *dual polarization feeds*.  Such radio telescopes are essentially two telescopes in one since each feed is sensitive to one of two orthogonal polarizations.  Lower frequency feeds tend to be linearly (X and Y) polarized and higher frequency feeds tend to be circularly (left-handed and right-handed) polarized, but using Jones matrices allows for conversion between polarization types.  Visibilities can be computed for inputs from the same polarization or from different polarizations.  For a given pair of dual polarization antennas there are four possible *cross-polarization* visibility products.  For example, if antennas A and B have linearly polarized feeds, the four cross-polarization visibility products are ``(A_X^{} B_{X}^*)``, ``(A_X^{} B_{Y}^*)``, ``(A_Y^{} B_{X}^*)``, and ``(A_Y^{} B_{Y}^*)``.  The cross polarization products can be referred to more generally without specific antennas, for example XX, XY, YX, or YY for linearly polarized feeds or LL, LR, RL, RR for circularly polarized feeds.  Polarization cross products can also be linearly combined to produce *Stokes parameters* I, Q, U, and V.

A convention used by multiple radio interferometry software packages is to assign integer constants to these various cross polarization products.  UVH5 has adopted/continued this convention.  The number of cross-polarizations present in a UVH5 file is stored in the `Header/Npols` field.  This value, ``N_{pols}``, is one of the dimensions of the four dimensional `Data/visdata` and related datasets.  The order in which the polarizations are presented in the data file is given by the `Header/polarization_array` dataset, which contains ``N_{pols}`` integer values indicating the cross-polarization products in the file and the order in which they are stored.  The values follow the convention described in [AIPS Memo 117](ftp://ftp.aoc.nrao.edu/pub/software/aips/TEXT/PUBL/AIPSMEM117.PDF) and summarized here:

| Cross-polarization | Integer code |
|-------------------:|:-------------|
| Stokes I           | +1           |
| Stokes Q           | +2           |
| Stokes U           | +3           |
| Stokes V           | +4           |
| Circular RR        | -1           |
| Circular LL        | -2           |
| Circular RL        | -3           |
| Circular LR        | -4           |
| Linear XX          | -5           |
| Linear YY          | -6           |
| Linear XY          | -7           |
| Linear YX          | -8           |

Cross referencing the values in the `Header/polarization_array` dataset of our UVH5 data file with this table, we can see that the cross-polarization products of each baseline are ordered as XX, XY, YX, YY:
"""

# ╔═╡ 810b605c-7f92-11eb-1825-9dbd89dcc5a4
pol_array = uvh5["Header/polarization_array"][]

# ╔═╡ df4cc7b2-7f79-11eb-1255-ab2dacc1053c
md"""
###### Frequency, bandwidth, and spectral windows

Correlators can produce visibilities covering multiple regions of frequency space.  These regions, often referred to as *spectral windows*, each have some number of contiguous frequency channels, but the windows are not contiguous with each other.  The UVH5 file format supports multiple spectral windows and every UVH5 file has at least one spectral window.  The number of spectral windows is stored in the `Header/Nspws` dataset.  The UVH5 file format currently constrains all spectral spectral windows in a data file to all have the same number of channels.  This value is stored in the `Header/Nfreqs` dataset.  These numbers, ``N_{freqs}`` and ``N_{spws}``, define the dimensions of the two dimensional `Header/freq_array` which stores the center frequency of each frequency channel for each window.  ``N_{freqs}`` and ``N_{spws}`` also define two dimensions of the four dimensional `Data/visdata` and related datasets.

The UVH5 file format also constrains all frequency channels of all spectral windows in a data file to have a common bandwidth.  This value is available in the `Header/channel_width` dataset.

!!! note "Future evolutions"
    The UVH5 conventions regarding frequency channels and spectral windows are an area where the specification is evolving.  Future versions of the UVH5 specification are likely to flatten the ``N_{freqs}`` and ``N_{spws}`` dimensions into a single dimension to allow the number of channels per spectral window to vary.  The bandwidth related parameters are also being enhanced to allow for more flexibility.  This future version of the UVH5 specification will also include new metadata to distinguish between various conventions.

Reading the `Header/freq_array` dataset from our UVH5 data file shows that it contains 1 spectral window of 256 frequency channels ranging from 989.75 MHz to 1043.04 MHz:
"""

# ╔═╡ 0be7e252-7f98-11eb-04f7-f5efb30cb004
freq_array = uvh5["Header/freq_array"][]

# ╔═╡ d862e16a-7f98-11eb-3422-31fa60d487b1
md"""
###### The baseline-time dimension

The final dimension of the `Data/visdata` dataset, *baseline-time*, is a flattening of what would otherwise be two separate dimensions of *baseline* and *time*.  This allows for more flexibility and file space efficiency.  Specifically, flattening these two dimensions into a single baseline-time dimension allows for time samples to have any number of baselines.  To put it another way, the UVH5 file format does not require that the dataset have space allocated for every baseline at every time sample.  This allows data files with different subsets of antennas to be combined more easily.

!!! note "Historical note"
    Early radio interferometers had as few as two physical antennas.  The antennas would be placed to form one baseline for some period of time (days or weeks), then one (or both) would be moved to another position to form a different baseline.  Data from these different baseline configurations would be combined to emulate an array of multiple telescopes with multiple baselines, just not with all baselines present at the same time.  For measuring non-transient objects this time difference was insignificant (other than perhaps some calibration challenges).

Many other datasets in a UVH5 file share this common baseline-time dimension, specifically those datasets which provide context for the data at each baseline-time position.  This includes:
- `Header/ant_1_array`
- `Header/ant_2_array`
- `Header/time_array`
- `Header/uvw_array`
- `Header/lst_array` (optional)
- All the datasets in the `Data` group

The size of this baseline-time dimension is available, as a scalar, in the `Header/Nblts` dataset.
"""

# ╔═╡ 27bd5ec8-7fd2-11eb-063b-d776aa9423a6
nblts = uvh5["Header/Nblts"][]

# ╔═╡ 72af3ff4-7a6f-11eb-14a6-eb0a67e3a741
md"""
### Correlation coefficients
For input signals ``Z_i`` and ``Z_j`` with expected value of zero, normalizing their cross-correlation by the geometric mean of their auto-correlations produces the *correlation coefficient*, often represented as ``\rho_{ij}``, for that pair of input signals.

```math
\rho_{ij} = \frac{\displaystyle \sum^NZ_i^{}Z_j^*}{\sqrt{\left(\displaystyle \sum^NZ_i^{}Z_i^*\right)\left(\displaystyle \sum^NZ_j^{}Z_j^*\right)}},\ \mathrm{where}\ i,j \in [1, A]
```

For complex inputs, the correlation coefficient is a complex value with magnitude between $0$ and $1$ and the same phase as the cross-correlation.  For 100% correlated signals, the correlation coefficient is ``1``.  For uncorrelated signals (e.g. independent noise), the correlation coefficient approaches ``0`` as ``N=B\tau`` approaches infinity, where ``B`` is the bandwidth, ``\tau`` is the integration time, and their product, ``N``, is the number of critically sampled values integrated.  The expected value of the magnitude of the correlation coefficient for uncorrelated inputs is ``\frac{1}{\sqrt{N}}=\frac{1}{\sqrt{B\tau}}``.
"""

# ╔═╡ e308c260-7ade-11eb-11a3-2fa5afe66a22
md"""
###### Expressing correlation coefficients in decibels (dB)
It is often useful to express the magnitude of the correlation coefficient in terms of decibels (dB).  When expressed as dB, it is implicitly understood to be of the magnitude.  Decibels are defined as ``10 \log_{10}{\frac{P_a}{P_b}}``, where ``P_a`` and ``P_b`` are powers in a common unit.  It is worth noting that decibels are unitless.  Often times it is more convenient and efficient to use properties of the logarithm to perform the  computation using an alternate but equivalent calculation. For example, instead of computing the expected value of the correlation coefficient of uncorrelated signals in dB as ``10\log_{10}{\frac{1}{\sqrt{B\tau}}}``, the square root and inversion can be avoided by equivalently computing ``-5\log_{10}B\tau``.  Likewise, calculating the correlation coefficient in dB can be done efficiently as:

```math
\rho_{ij} = 5 \log_{10}\left(\frac{\left\lvert\displaystyle \sum^NZ_i^{}Z_j^*\right\rvert^2}{\left(\displaystyle \sum^NZ_i^{}Z_i^*\right)\left(\displaystyle \sum^NZ_j^{}Z_j^*\right)}\right)\ \mathrm{dB}
```

The square of the absolute value in the numerator can be performed with the `abs2` function which is more efficient for complex inputs than squaring the output of the `abs` function because `abs(z)` for complex `z` inputs requires an extra `sqrt` call.  Likewise, the denominator has no square roots so this calculation avoids square roots altogether.
"""

# ╔═╡ 3b6083bc-7af8-11eb-1c69-8df28a0ac02f
md"""
### Plotting correlation coefficients in dB
Here you can explore the correlation coefficients for various baselines in the test UVH5 file.  This file includes 16 of 64 antennas antennas, so there are ``\frac{16\cdot15}{2} = 120`` possible cross correlation baselines.  The correlation coefficients across the frequency range of the data file are plotted in both magnitude (dB) and phase (degrees).  The expected value of the  correlation coefficient (in dB) for uncorrelated inputs with the baseline's channel width and integration time is also shown as a horizontal line.
"""

# ╔═╡ 59ed297e-7a66-11eb-2a4e-ff07618aa387
begin
	a1s = uvh5["Header/ant_1_array"][:]
	a2s = uvh5["Header/ant_2_array"][:]
	blts = collect(zip(a1s, a2s))
	ants = unique([a1s; a2s])
	md"""
	Select two antennas of a baseline:
	$(@bind a1str Select(string.(ants))) $(@bind a2str Select(string.(ants), default=string(ants[2])))
	"""
end

# ╔═╡ b8c340ec-7a67-11eb-2710-091fd2edabfd
begin
	a1 = parse(Int, a1str)
	a2 = parse(Int, a2str)
	# UVH5 allows baselines to have a1 > a2,
	# but this dataset does not have any such baselines
	a1 <= a2 || ((a1, a2) = (a2, a1))
	a11blts = findall(==((a1,a1)), blts)
	a22blts = findall(==((a2,a2)), blts)
	a12blts = findall(==((a1,a2)), blts)
	time_select = @bind bltidxstr Select(
		map(i->string(i)=>string(julian2datetime(uvh5["Header/time_array"][a12blts[i]])), 1:length(a12blts)))
	md"""
	Select an integration of that baseline:
	$time_select
	"""
end

# ╔═╡ 77ce98d8-7fd6-11eb-1c4d-cf2fa6fc50a0
begin
	polmap = Dict{Any,Any}(
		+1=>"I",  +2=>"Q",  +3=>"U",  +4=>"V",
		-1=>"RR", -2=>"LL", -3=>"RL", -4=>"LR",
		-5=>"XX", -6=>"YY", -7=>"XY", -8=>"YX"
	)
	pols = [polmap[p] for p in pol_array]
	md"""
	Select a polarization product:
	$(@bind polstr Select(pols))
	"""
end

# ╔═╡ 1404ed9e-7a6c-11eb-2883-13056a01a8f0
begin
	pol12idx = findfirst(==(polstr), pols)
	pol11idx = findfirst(==(polstr[1]^2), pols)
	pol22idx = findfirst(==(polstr[2]^2), pols)
	bltidx = parse(Int, bltidxstr)
	a11data = Float32.(abs.(visdata[pol11idx,:,1,a11blts[bltidx]]))
	a22data = Float32.(abs.(visdata[pol22idx,:,1,a22blts[bltidx]]))
	a12data = ComplexF32.(visdata[pol12idx,:,1,a12blts[bltidx]])
	corrcoeffdb = 5*log10.(abs2.(a12data)./(a11data.*a22data))
	uncorrdb = -5*log10(Bτ)
	p1 = plot(freqs_mhz, corrcoeffdb,
		xlims=extrema(freqs_mhz), xticks=:native,
		ylims=(-20.2,0.2),
		yticks=-20:5:0,
		title="Baseline $a1-$a2 ($polstr)",
		ylabel="Corr Coeff (dB)")
	plot!(freqs_mhz[[begin, end]], [uncorrdb, uncorrdb])
	
	p2 = scatter(freqs_mhz, rad2deg.(angle.(a12data)),
		xlims=extrema(freqs_mhz), xticks=:native,
		ylims=(-183,183), yticks=-180:90:180,
		xlabel="Frequency (MHz)",
		ylabel="Phase (deg)"
	)
	plot(p1, p2, layout=(2,1), link=:x)
end

# ╔═╡ 57406bdc-7b3f-11eb-15d9-b99031aa956a
md"""
###### Things to notice
As you explore the plots for the different baselines and integrations, you may notice certain recurring features, such as:

1. The correlation coefficient is generally higher than the expected value for uncorrelated noise when XX or YY cross-polarizations polarizations are selected, but not when XY or YX polarizations are selected.

2. For XX and YY cross-polarizations, the phases tend to cluster along a horizontal or slightly diagonal line, but for XY and YX cross polarizations the phases are generally uniformly scattered.

Having a correlation coefficient higher than the expected value and having structure in the phases (i.e. not uniformly scattered) are two qualities that indicate correlation, also sometimes called *coherence*.  The object of this observation, $(uvh5["Header/object_name"][]), is a quasar.  Quasars tend to have high and/or well defined flux (radio brightness) and relatively compact structure.  These two properties make quasars well suited for use as calibration sources for radio interferometers.  The fact that XX and YY cross-polarizations have coherence, but XY and YX do not suggests that this quasar is unpolarized.

In addition to detecting the desired correlation due to radio waves arriving from space, it should also be noted that correlators are very good at detecting unwanted correlation due to common mode signals that leak into or between the electronic signal paths at the observatory.  This unwanted correlation reduces the overall sensitivity of the interferometer.

The phase angle of a cross-correlation is the phase difference (or *relative phase*) between the two inputs.  The clustering of phases in a straight line is due to the linear (with frequency) relationship between delay and phase.  Phase (``\phi``) is essentially frequency (``\omega``) times time (``t``).  It follows that ``\Delta \phi = \omega \Delta t``.  If the correlated signal arrives at two antennas at different times, the time offset (aka delay) will result in a *phase slope* across the band.  The longer the delay, the steeper the phase slope.  Instrumental factors also influence the relative phases between pairs of input.  Phase calibration is the process of measuring and compensating for these various phase effects.
"""

# ╔═╡ ea0b9abe-7b98-11eb-33c5-0fd12d130aed
md"""
### Calculating input RMS values
Understanding the power levels of the correlator's input signals is important to ensure/verify that the system is operating in a linear regime.  A correlator's input signals (voltages) are generally noise dominated with zero-mean Gaussian statistics.  The meaningful quantity for understanding input levels is the power averaged over some time interval.  Power is proportional to the square of the voltage, so taking the mean of the square of the voltages produces a value that is proportional to the power.  Taking the square root of the mean of the squared voltages gives a value in the same scaled units of the input voltage.  This value is called the *RMS voltage*, where RMS stands for *root mean square* (named after the exact steps used to create it).  The RMS voltage conveys the same information as the average power (i.e. the mean of squares of voltages), but being in the same scaled units as the input samples tends to make the RMS voltage more useful in a practical sense.

Fortunately, the auto-correlations have done the hard part of summing up the squares of voltages.  In the case of complex inputs, the auto-correlations contain the sum of the squares of the real and imaginary components of the input voltages: ``\sum^N (Z_{re}^2+Z_{im}^2)``.  If we assume that the real and imaginary components of an input signal have equal average power (a valid assumption for voltages from an FFT or PFB), then the input RMS voltage, ``V_{RMS}``, can be calculated as:

```math
V_{RMS} = \sqrt{\frac{1}{2N} \sum^N (Z_{re}^2 + Z_{im}^2)}
```
"""

# ╔═╡ c8a2cb68-7b9d-11eb-2f6e-b9981c5282d1
md"""
### Plotting input RMS values
The auto-correlations from the UVH5 file's `Data/visdata` dataset provide the values for the summation factor in the input RMS formula shown above.  All we need to do is divide by ``2N`` and take the square root.  As mentioned elsewhere, ``N`` is not stored directly in the UVH5 file, but it can be calculated as ``N=B\tau``, where bandwidth ``B`` is available in `Header/channel_width` and integration time ``\tau`` is available in `Header/integration_time`.
"""

# ╔═╡ 1c3f6b7e-7bb6-11eb-2c1c-f34e81d4805e
begin
	md"""
	Select an antenna to see its input RMS:
	$(@bind autoantstr Select(string.(ants)))
	"""
end

# ╔═╡ 1f008a4a-7bb7-11eb-12a3-d3e7963003b6
begin
	autoant = parse(Int, autoantstr)
	autoblts = findall(==((autoant,autoant)), blts)
	auto_select = @bind autoidxstr Select(
		map(i->string(i)=>string(julian2datetime(uvh5["Header/time_array"][i])), autoblts))
	md"""
	Select an integration of that baseline:
	$time_select
	"""
end

# ╔═╡ 3d882432-7fe0-11eb-0ddf-ed7307e60383
md"""
Select a polarization product:
$(@bind autopolstr Select(filter(pp->pp[1]==pp[2], pols)))
"""

# ╔═╡ 3059d072-7bb8-11eb-190a-ab79c0e1068d
begin
	polidx = findfirst(==(autopolstr), pols)
	autoidx = parse(Int, autoidxstr)
	autodata = Float32.(real.(visdata[polidx,:,1,autoidx]))
	autorms = sqrt.(autodata ./ 2Bτ)
	plot(freqs_mhz, autorms,
		xlims=extrema(freqs_mhz), xticks=:native,
		title="Antenna $autoant ($autopolstr)",
		xlabel="Frequency (MHz)",
		ylabel="Input RMS (counts)")
end

# ╔═╡ 0e67cd02-7bbb-11eb-00aa-bfc7e8ef4455
md"""
###### Things to notice
Except for one channel around 1,019 MHz, the input RMS for these inputs is generally between 10.5 and 11.5.  Is that good or bad?  If the input samples were 5 bits each, that would be far too high, but if the input samples were 8 bits this is a very comfortable level.  The number of bits per input sample is not present in the UVH5 file, so this information must be obtained elsewhere.  For the test data file being used here, the real and imaginary components of the input samples were 8 bits each.

With this extra information, we can express the input levels in terms of *dBFS* (dB Full Scale) which is defined as:
```math
\mathrm{dBFS} = 20 \log_{10} \left(\frac{V_{RMS}}{2^{nbits-1}-1}\right)
```

Using 11 for ``V_{RMS}`` and 8 for ``nbits``, we can express the input level as approximately $(round(20*log10(11/(2^7-1)), digits=3)) dBFS.  Systems that quantize to even fewer bits need to pay even more attention to input levels to ensure that the quantization does not result in a non-linearity.  Most calibration and measurement done with an interferometer rely upon the system being linear.

!!! note
    While this discussion about input RMS levels applies to both XF and FX correlators, it is most germane to FX correlators which typically include a digital re-quantization of the signals between the F and X stages of the correlator.
"""

# ╔═╡ 5ba19006-7bc0-11eb-08b2-838f828988dd
md"""
### Next steps
I hope this notebook has given you a good introduction to the UVH5 data format and visibility data in general.  There are numerous topics that were glossed over or were not even covered at all, but this was intended as an overview and not an in depth treatise.  Some possible next steps could be exploring:

1. More details about the structure of a UVH5 file and how to navigate it
2. The geometric relationship between baselines, sky positions, delays, and phases
3. Instrumental effects and how to calibrate them or verify that they are already properly calibrated
4. Imaging or other localization techniques involving the visibility data

Stay phased!
"""

# ╔═╡ 32b5ab62-89c0-11eb-08ef-571a2c7c0e1a
PlutoUI.TableOfContents()

# ╔═╡ Cell order:
# ╟─782d9788-77da-11eb-357f-e3857ed0f0ec
# ╟─440400b0-7be9-11eb-0efb-c723129eeef2
# ╟─2c9addb8-7be9-11eb-0adf-452042d88a4c
# ╟─99f47c64-78df-11eb-0ee3-af29fbc46edf
# ╠═83191f96-77da-11eb-16bb-8bf06f6b60ed
# ╟─b05be61e-77da-11eb-3a82-edc1be04b007
# ╟─d893adaa-7afb-11eb-143d-1fa4d8054566
# ╟─42167ea4-7b26-11eb-180b-138a95419a71
# ╟─6f2807e6-7b2b-11eb-0844-39b3e53b116c
# ╟─906c6b04-7b35-11eb-1ba6-110b111b728e
# ╟─d3a8122c-7b3c-11eb-1bc3-51a1da2515a3
# ╟─b79230c6-7f93-11eb-046d-330e79e39964
# ╟─53e0fc08-7f8c-11eb-3f3f-e7c01518e7d7
# ╠═810b605c-7f92-11eb-1825-9dbd89dcc5a4
# ╟─df4cc7b2-7f79-11eb-1255-ab2dacc1053c
# ╠═0be7e252-7f98-11eb-04f7-f5efb30cb004
# ╟─d862e16a-7f98-11eb-3422-31fa60d487b1
# ╠═27bd5ec8-7fd2-11eb-063b-d776aa9423a6
# ╟─72af3ff4-7a6f-11eb-14a6-eb0a67e3a741
# ╟─e308c260-7ade-11eb-11a3-2fa5afe66a22
# ╟─3b6083bc-7af8-11eb-1c69-8df28a0ac02f
# ╟─59ed297e-7a66-11eb-2a4e-ff07618aa387
# ╟─b8c340ec-7a67-11eb-2710-091fd2edabfd
# ╟─77ce98d8-7fd6-11eb-1c4d-cf2fa6fc50a0
# ╟─1404ed9e-7a6c-11eb-2883-13056a01a8f0
# ╟─57406bdc-7b3f-11eb-15d9-b99031aa956a
# ╟─ea0b9abe-7b98-11eb-33c5-0fd12d130aed
# ╟─c8a2cb68-7b9d-11eb-2f6e-b9981c5282d1
# ╟─1c3f6b7e-7bb6-11eb-2c1c-f34e81d4805e
# ╟─1f008a4a-7bb7-11eb-12a3-d3e7963003b6
# ╟─3d882432-7fe0-11eb-0ddf-ed7307e60383
# ╟─3059d072-7bb8-11eb-190a-ab79c0e1068d
# ╟─0e67cd02-7bbb-11eb-00aa-bfc7e8ef4455
# ╟─5ba19006-7bc0-11eb-08b2-838f828988dd
# ╟─32b5ab62-89c0-11eb-08ef-571a2c7c0e1a

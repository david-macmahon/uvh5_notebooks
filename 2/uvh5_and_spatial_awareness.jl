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
			(name="PlutoUI",    version="0.7"),
			(name="HDF5",       version="0.15"),
			(name="Plots",      version="1.11"),
			(name="PlotlyJS",   version="0.14"),
			(name="DataFrames", version="0.22"),
			(name="Geodesy",    version="1.0"),
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

	using HDF5, Plots, PlutoUI, Dates, DataFrames, Geodesy, RadioInterferometry
	plotlyjs(legend=false)
	HDF5.SHOW_TREE_MAX_CHILDREN[] = 100
	md"""
	## UVH5 and spatial awareness
	"""
end

# ╔═╡ 440400b0-7be9-11eb-0efb-c723129eeef2
md"""
The [UVH5 explorations](open?path=uvh5_explorations.jl) notebook provided an introduction to the UVH5 file format and some basic interferometry concepts related to correlation.  In the discussion of auto-correlations and cross-correlations it made brief mention of a *baseline vector* as the vector from one antenna position to another antenna position (or, for auto-correlation baselines, the same antenna position).  This notebook will go into a bit more depth about baseline vectors and discuss various frames of reference that they can be expressed in.  Knowledge of the antennas' physical positions at the observatory and the observatory's physical position on planet Earth (or elsewhere) is critical to being able to interpret interferometric data properly.  Not surprisingly, the metadata in the UVH5 `Header` contains several pieces of information on this topic.  Also not surprisingly, there is a lot to talk about here, so let's get spatial!
"""

# ╔═╡ 2c9addb8-7be9-11eb-0adf-452042d88a4c
md"""
!!! note "Disclaimer"
    Being an introductory level notebook, some of the more advanced (and even some basic) features of UVH5 will be glossed over, ignored outright, or deferred to future notebooks.  Some of the techniques used here take advantage of *a priori* knowledge of the conventions used when creating the test data file being used (e.g. there is only one common integration time for all baselines) and are not generally applicable to more complicated but equally legitimate UVH5 file layouts.
"""

# ╔═╡ 99f47c64-78df-11eb-0ee3-af29fbc46edf
md"""
### Getting started
We are going to use the same UVH5 file as the **UVH5 explorations** notebook.  After opening, we are going to go straight to the `Header/antenna_positions` dataset.  This dataset is dimensioned ``(3, N_{ants\_telescope})``.

!!! note "Note about order of dimensions"
    Julia uses *column major ordering* like Fortran and Matlab rather than *row major ordering* like C and Python.  The order of dimensions in Julia/Fortran/Matlab will be reversed from the order of dimensions used in C/Python (and the HDF5 command line tools like `h5ls` and `h5dump`).  The layout of data in memory and on disk remains the same, just the dimensions are flipped end to end.  Because this is a Julia notebook, the order of dimensions presented and used here will follow the column major ordering convention.

!!! note "Note about indexing"
    As a reminder, Julia indexes Arrays from 1 like Fortran and Matlab rather than from 0 like C and Python.  It can be a bit of an adjustment coming from a language with 0-based indexing, but you're smart and will pick it up in no time!
"""

# ╔═╡ 0416bb34-7d7c-11eb-3060-b1c01edbf110
md"""
### Antenna positions
"""

# ╔═╡ 83191f96-77da-11eb-16bb-8bf06f6b60ed
antpos = begin
	uvh5 = h5open(artifact("uvh5_examples", "meerkat_example.uvh5"))
	uvh5["Header/antenna_positions"][]
end

# ╔═╡ b05be61e-77da-11eb-3a82-edc1be04b007
md"""
The `Header/antenna_positions` dataset of our UVH5 file contains a 3x64 Array of Float64 values, which we have stored in the `antpos` variable.  As you may recall from the **UVH5 explorations** notebook, this UVH5 file is from the MeerKAT radio telescope array which has 64 antennas.  That explains why we have 64 columns.  By default, Julia only shows the first and last few columns to fit the display width.  Each column contains three values that represent the three dimensional antenna position, but to pin down what these coordinates really mean, we need to know their units and reference frame.  Another detail we need to know is which column corresponds to which antenna!
	
###### The ECEF frame
	
The UVH5 memo specifies that the antenna position units are meters and that their reference frame is the *Earth-centered, Earth-fixed* ([ECEF](https://en.wikipedia.org/wiki/ECEF)) reference frame of the International Terrestrial Reference System (ITRS).  This is why the hypotenuse of each point is approximately six million meters, i.e. the approximate radius of Earth.  This frame is often referred to as a *geocentric XYZ frame* or just an *XYZ frame* .  These latter terms are not so canonical in that *geocentric XYZ frame* has some ambiguity about where the X axis points and *XYZ frame* also ambiguity about where it is centered.  For example, an *XYZ frame* may have axes that are parallel to the XYZ axes of the ITRS ECEF frame, but it may have been translated to have its origin at another location (e.g. the observatory's reference position).  In all reasonable uses, the Z axis of an XYZ frame is parallel to the polar Z axis of the ITRS ECEF frame.  In practice the context is usually clear, but it never hurts to be explicit in comments or accompanying text.
	
To know which column of the `antpos` Array corresponds to which antenna, we need to read the corresponding element from the `Header/antenna_numbers` dataset and/or the `Header/antenna_names`.  MeerKAT numbers the antennas from 0 to 63 and names them from `m000` to `m063`.  This particular data file happens to order the `Header/antenna_number` and `Header/antenna_names` in the same order, so the only thing we need to worry about is that we need to use column indices 1 to 64 to reference the position of antenna numbers 0 to 63, but generally speaking one should look up the index of the desired antenna number or name in the appropriate dataset to find the column index of the desired antenna position.  UVH5 allows antenna numbers to be stored in any order (out of order, non-sequential, etc) and the other antenna related Arrays must be indexed in the corresponding order.
"""

# ╔═╡ f3378276-7d7b-11eb-2495-030f08a01f0c
md"""
###### Plotting antenna positions

Let's make a 3D scatter plot of the antenna positions!
"""

# ╔═╡ 11499c58-7d7d-11eb-36d5-7f793f9d7d5d
scatter(antpos[1,:], antpos[2,:], antpos[3,:])

# ╔═╡ c5ebb850-7d7d-11eb-11e9-29094ba7160b
md"""
If you rotate the plot around, you will see that the antennas are on a plane that is somewhat inclined relative to the axes.  That's to be expected since the topocentric (local) plane of the observatory is itself inclined with respect to the XYZ axes of the ITRS ECEF frame.  As cool as this plot is, it is probably more useful to make a two dimensional bird's eye view looking down on the radio telescope array.

###### The ENU frame

To do that we need to transform the ITRS ECEF coordinates to a new topocentric coordinate system that is centered at the observatory's reference position and rotated to be aligned to local East, North, and Up vectors.  This new frame is known as a *topocentric ENU frame*.  One small wrinkle is that the observatory's reference position is specified in the UVH5 `Header` metadata as the more familiar geodetic `latitude` (degrees), `longitude` (degrees), and `altitude` (meters), which are not directly usable with coordinates in ITRS ECEF frame.
"""

# ╔═╡ eac37b0a-7d84-11eb-1430-95092833fe33
md"""
### Geodesy odyssey

To handle this conversion, we will use the `Geodesy.jl` Julia package.  This package provides transformations between ECEF and ENU coordinates for a given latitude, longitude, altitude (LLA) position.  It's going to be a great journey!

First we have to create an instance of Geodesy's `LLA` structure to store the observatory's LLA reference position.
"""

# ╔═╡ 91707832-7d7f-11eb-0ca8-fbd7a47adda4
begin
	latitude = uvh5["Header/latitude"][]
	longitude = uvh5["Header/longitude"][]
	altitude = uvh5["Header/altitude"][]
	reflla = LLA(latitude, longitude, altitude)
end

# ╔═╡ 86308de8-7d86-11eb-353e-d75a37cb7d1f
md"""
Now we can use Geodesy's `ENUfromECEF` function generator to create a transformation function that is customized to transform from the ECEF frame to to an ENU frame with the origin at the observatory's reference position.  We will name this customized transformation function `enufromecef`.  Notice that we need to specify not only the `reflla::LLA` object created above for the observatory's reference position but also the corresponding *reference ellipsoid*.  UVH5 files (or at least this one) use the WGS84 ellipsoid, which we specify with `wgs84` (defined in Geodesy).
"""

# ╔═╡ 13bb054c-7d89-11eb-34ca-0f84e57aa525
enufromecef = ENUfromECEF(reflla, wgs84)

# ╔═╡ 5c84f00a-7d8c-11eb-31cd-e3f4157f534d
md"""
Now we need to perform the actual conversion.  This involves a few extra steps because we like to work with Arrays (easier slicing and dicing), but Geodesy (and our customized transformation function) work with `ECEF` instances.  The line of code in the following notebook cell does these things:
1. Converts the 3x64 `antpos` Array to a Vector id `ECEF` instances
2. Transforms the Vector of `ECEF` instances to a Vector of `ENU` instances
3. Reduces the Vector of `ENU` instances to a 3x64 Array containing ENU values using the `hcat` function.
We store the results in the `antenu` variable.  The dots in the expression below invoke Julia's "broadcast" functionality, which is quite powerful.  We use it here to perform the function calls element-wise over the Vectors.
"""

# ╔═╡ c2d0b58e-7d86-11eb-1812-b96f2ce720b1
antenu = reduce(hcat, enufromecef.(ECEF.(eachcol(antpos))))

# ╔═╡ 9367d2fe-7d87-11eb-1eda-7f466558ff4d
md"""
###### Plotting antenna positions (take 2)

Now we can make the two dimensional bird's eye view of the antenna positions.  To make the plot a bit more engaging, we will show the antenna name when the cursor "hovers" over an antenna.  To provide a bit more context for this particular data file, we will use the color red for the antennas that are present in the UVH5 file and blue for antennas that were omitted from this observation.  As you may recall from the **UVH5 explorations** notebook, this UVH5 file contains data from 16 of MeerKAT's 64 antennas.
"""

# ╔═╡ 5ec2e7fc-7d8f-11eb-2f0a-11afe337cb4d
begin
	# This is probably overkill since all present antennas are likely
	# to be in just ant_1_array (or just ant_2_array), but UVH5 does
	# allow for files without auto-correlations so we use the union
	# operator (aka "\cup<TAB>") to ensure that we get all antennas
	# that are present in this data file.
	obsantnums = uvh5["Header/ant_1_array"][] ∪ uvh5["Header/ant_2_array"][]
	allantnums = uvh5["Header/antenna_numbers"][]
	scatter(antenu[1,:], antenu[2,:],
		xlims=(-4000,4000), ylims=(-4000,4000), aspect_ratio=:equal,
		xticks=:native, yticks=:native,
		title="""$(uvh5["Header/telescope_name"][]) Antenna Positions""",
		xlabel="East (m)",
		ylabel="North (m)",
		hover=uvh5["Header/antenna_names"][:],
		color=ifelse.(allantnums .∈ Ref(obsantnums), :red, :blue)
)
end

# ╔═╡ 7550d5ee-7ed7-11eb-3ce0-51c68ce5bf5a
md"""
### Baseline vectors

Baseline vectors can be represented in either of the reference frames we have already used for antenna positions, the XYZ reference frame or the ENU reference frame.  A baseline vector for two antennas can be calculated simply by subtracting one antenna's position from the other antenna's position (in the same frame, of course).  By convention, a baseline vector points from the first antenna to the second antenna, so the position of the first antenna is subtracted from the position of the second antenna to get the baseline vector.

Here is the baseline vector for baseline 0-1 of MeerKAT in the XYZ frame used for antenna positions in the UVH5 file and in the ENU frame that we transformed the ECEF antenna positions to.
"""

# ╔═╡ 2c52245a-7ed8-11eb-0f76-d1591196f9d1
bl_0_1_ecef = antpos[:,2] - antpos[:,1]

# ╔═╡ f9627760-7ed8-11eb-009f-45e3f92053b7
bl_0_1_enu = antenu[:,2] - antenu[:,1]

# ╔═╡ 8eba61ce-7ed9-11eb-08c5-4ffcb667f664
md"""
These vectors look quite different because they have different reference frames, but the length of the two vectors are the same as shown here.
"""

# ╔═╡ 1b46bd5a-7ed9-11eb-2661-75f2f5d4987f
hypot(bl_0_1_ecef...), hypot(bl_0_1_enu...)

# ╔═╡ e87c0060-7eec-11eb-3607-0beac189151e
md"""
Antenna positions and baseline vectors are not limited to these two reference frames.  We can use any three dimensional reference frame to describe antenna positions and baseline vectors.  When antennas are Earth-based, like MeerKAT's, the two Earth-based frames discussed thus far tend to be the most useful for antenna positions and therefore the most commonly used.  The XYZ frame tends to be used more for arrays with geographically distant antennas (e.g. VLBI) while the ENU frame tends to be used more for arrays on a single geographic site where a single concept of "up" can be used.
"""

# ╔═╡ ea4a03f6-7eda-11eb-0499-77b6062ea04e
md"""
### The UVW frame

The XYZ and ENU Earth-based frames are not particularly useful for baseline vectors.  Baseline vectors are more useful in relation to the astronomical objects or regions we wish to observe.  This leads to the introduction of yet another reference frame: the *UVW reference frame*.  Like the ENU frame, the UVW frame is really more like a meta-frame because there is no single canonical origin.  Unlike both the XYZ and ENU frames, the orientation of a UVW frame is fixed relative to a reference point on the *celestial sphere* rather than on (or in) Earth.  To an Earth-based observer, the UVW frame will appear to rotate due to Earth's rotation.  To a UVW-based observer, the Earth-based frames will appear to rotate along with the rotating Earth.  The UVW frame consists of the so-called *UV plane* and the *W axis* that completes the right handed three dimensional coordinate system.  The UV plane is parallel to the plane tangent to the celestial sphere at the reference point.  The UV plane is oriented such that the positive U axis points east (in *right ascension*) and the positive V axis points north (in *declination*).  The W axis is normal to (i.e. perpendicular to) the UV plane and is essentially the line of sight to the reference point (also called the *tangent point*).

Because the Earth-based antenna positions and baseline vectors are rotating with respect to a given UVW frame, the UVW coordinates of the antenna positions and baseline vectors change over time.  In other words, the visibility data for a given baseline sampled at different times will have different UVW coordinates for a given UVW frame.  UVH5 files store these UVW coordinates in the `Header/uvw_array` dataset.  The *(u,v)* coordinates of a baseline sample give the *UV spacing* of that baseline sample.  It is possible for different pairs of antennas to sample the same UV spacing at different times (or possibly the same time).  The *w* coordinate of a baseline gives the *line-of-sight distance* (also called the *delay distance*) that the EM wavefront traveled from one of the baseline's antennas to the other.  The *w* component is closely related to delay and delay is closely related to phase.

Coordinates in the UVW frame are commonly given in frequency independent units of meters or nanoseconds (based on speed of light) or the frequency dependent units of wavelengths, typically indicated by the Greek letter lambda ``\lambda``).  In UVH5 files, the UVW coordinates are in meters.
"""

# ╔═╡ d0c7bd6e-7f10-11eb-07e4-65b181fa87b2
md"""
### The UVW array
The UVH5 `Header/uvw_array` dataset contains the UVW coordinates of each occurrence of each baseline in the file.  Let's take a look at the UVW data in our sample data file.
"""

# ╔═╡ 0dda68f2-7f14-11eb-31eb-b5b3f0e0acef
uvw_array = uvh5["Header/uvw_array"][]

# ╔═╡ 4fc456e2-7f14-11eb-3b11-3159926f69fc
md"""
Notice that this is a two dimensional dataset.  Instead of using three dimensions of (3, *baseline*, *time*), the UVH5 file format flattens the *baseline* and *time* dimensions into a single *baseline-time* dimension.  The size of this baseline-time dimension is available, as a scalar, in the `Header/Nblts` dataset.  The dimensions of the `Header/uvw_array` dataset are ``(3, N_{blts})``.
"""

# ╔═╡ 9db1fed6-7f19-11eb-21a7-bb7fb8e19108
md"""
The antenna numbers of the two antennas corresponding to a particular baseline-time can be read from the corresponding position in the `Header/ant_1_array` and `Header/ant_2_array`, which are each `Nblts` long.  Likewise the time corresponding a particular baseline-time can be read from the corresponding position in the `Header/time_array`, which also `Nblts` long and contains Julian dates stored as `Float64`s values (which limits them to approximately millisecond precision).  We can combine all these datasets into a `DataFrame`:
"""

# ╔═╡ 6f4a42be-7f1a-11eb-0ee4-fb4f2108ba47
begin
	a1s = uvh5["Header/ant_1_array"][]
	a2s = uvh5["Header/ant_2_array"][]
	times = julian2datetime.(uvh5["Header/time_array"][])
	bltuvw_df = DataFrame("Ant 1"=>a1s, "Ant 2"=>a2s, "Time"=>times,
		"U"=>uvw_array[1,:], "V"=>uvw_array[2,:], "W"=>uvw_array[3,:]
	)
end

# ╔═╡ c91745c0-7f16-11eb-3b9c-6dd85327d520
md"""
You will see that the first 16 baseline-time positions correspond to auto-correlations for time 2020-08-13T16:16:55.01.  Auto-correlations have zero-length baseline vectors, so they are all (0,0,0) in the UVW frame.  The next 120 baseline-time positions correspond to cross-correlations baselines from the same time.

Closer examination will reveal that this data file contains 4 groups of 136 baseline-times, with each group corresponding to a unique time.  Furthermore, the baseline order within each group is consistent across groups.  In other words, this UVH5 data file happens to have a very regular structure that is not required by the UVH5 specification and generally should not be assumed when working with UVH5 data files.
"""

# ╔═╡ f7f47dea-7f1d-11eb-1d53-97a4821b4623
md"""
### Plotting UV spacings on the UV plane
We can create a scatter plot on the UV plane showing all the UV spacings that a dataset has sampled.  This is called the *UV coverage* of the observation.  Note that sampling UV spacing ``(u,v)`` means that we have also sampled UV spacing ``(-u,-v)`` because the visibility data we sampled for ``(u,v)`` is the conjugate of the visibility data for ``(-u,-v)``.  The plot includes both ``(u,v)`` points and ``(-u,-v)`` points in different colors.  Auto-correlations are not shown.  Only the UV points for cross-correlation baselines for the data file's first time value are plotted as the file is too short in duration to see any significant change in the UVW coordinates over time.  If an observation is long enough, each baseline's UV coverage will trace out a curve on the UV plane.
"""

# ╔═╡ 33992d3a-7f1e-11eb-1be4-9526e6e17743
begin
	scatter(uvw_array[1,17:136], uvw_array[2,17:136],
		title="UV Coverage",
		xlabel="U (m)", ylabel="V (m)",
		xticks=:native, yticks=:native,
		hover=string.(zip(a1s[17:136], a2s[17:136]))
	)
	scatter!(-uvw_array[1,17:136], -uvw_array[2,17:136],
		hover=string.(zip(a1s[17:136], a2s[17:136]))
	)
end

# ╔═╡ 3b6083bc-7af8-11eb-1c69-8df28a0ac02f
md"""
### Plotting amplitudes vs UV distance

One type of plot commonly used when working with visibility data is a scatter plot with *UV distance* on the horizontal axis and baseline amplitude on the vertical axis.  As the name suggests, UV distance is length of the baseline's UV spacing (i.e. the distance from the UV plane origin to the baseline's ``(u,v)`` coordinate.  The baseline amplitudes are typically the aggregate amplitude of the baseline's visibility data across all frequency channels.  This aggregate amplitude can be the amplitude of the sum of the complex visibility values (sometimes called the *vector sum*) or the sum of the amplitudes of complex visibility values (sometimes called the *scalar sum*).  The vector sum will always be less than or equal to the scalar sum.  The vector sum is generally preferable as it is a scaled version of the complex average, but it requires that the phases be relatively flat (i.e. constant) across the band.  If the baseline has coherence with significant phase slope across the band the vector sum can go to 0.  The scalar sum is unaffected by phase slope, but it is also unaffected by incoherence so it should be used with caution.  Comparing the scalar sum with the vector sum is one way to identify baselines with phase slope or lack of correlation.

For a point source at the phase center, such as the calibrator in this dataset, the cross-correlation amplitudes should be independent of UV distance.  For other types of observations, baseline amplitudes can be and often are dependent on UV distance, which can reveal information about the morphology of the objects being observed.

Here is the amplitude vs UV distance plot for the first integration (i.e. time sample) in the data file.  You can choose to plot either the vector sum or the  scalar sum of the frequency channels of each baseline.  Only the XX and YY visibilities are plotted (in different colors).
"""

# ╔═╡ b696f722-7fef-11eb-0be8-23af42681850
md"""
Choose the type of frequency channel summing for the visibility amplitudes:

$(@bind sumtype Radio(["Scalar", "Vector"], default="Scalar"))
"""

# ╔═╡ 1404ed9e-7a6c-11eb-2883-13056a01a8f0
begin
	visdata = uvh5["Data/visdata"]
	# We cheat and use a priori knowledge of the
	# polarization order and baseline-time order.
	xxdata, yydata = if sumtype == "Scalar"
		sum(abs.(visdata[1,:,1,17:136]), dims=1),
		sum(abs.(visdata[4,:,1,17:136]), dims=1)
	else
		abs.(sum(visdata[1,:,1,17:136], dims=1)),
		abs.(sum(visdata[4,:,1,17:136], dims=1))
	end
	ampmax = max(maximum(xxdata), maximum(yydata))
	uvws = uvw_array[:,17:136]
	uvdists = hypot.(uvws[1,:], uvws[2,:])
	scatter(uvdists, xxdata[1,:],
		ylims=(0, 1.1*ampmax),
		title="Amplitude vs UV Distance ($(lowercase(sumtype)) sum)",
		xlabel="UV Distance (m)",
		ylabel="Amplitude (counts)",
		legend=true,
		label="XX",
		hover=string.(zip(a1s[17:136], a2s[17:136]))
	)
	
	scatter!(uvdists, yydata[1,:],
		label="YY",
		hover=string.(zip(a1s[17:136], a2s[17:136]))
	)
end

# ╔═╡ 2e28c1a8-7ff0-11eb-1455-0be7a7848000
md"""
As we saw in the **UVH5 explorations** notebook, the baselines in this dataset have a range of different phase slopes across the band.  You can use the vector sum vs UV distance plot to identify those that have substantial phase slope.  You go back to the **UVH5 explorations** notebook to plot the correlation coefficients for various baselines and compare their phase slope across the band with where they appear on the vector sum vs UV distance plot.
"""

# ╔═╡ 5ba19006-7bc0-11eb-08b2-838f828988dd
md"""
### Next steps
I hope this notebook has helped you feel more spatially aware and more familiar with the UVH5 file format.  It probably won't surprise you to learn that even more spatial/geometric interferometric topics await your discovery!  Plus we haven't really talked yet about the different types of calibration that need to be done on interferometers and the interferometric data they produce.  Some possible next steps could be exploring:

1. More spatial awareness, including astrometric transformation from catalog position to observed direction
2. A closer look at the direction/delay/phase relationship (also very relevant to coherent beamforming)
3. Deriving antenna-based parameters, such as system temperature and instrumental phase, from baseline visibility data
4. Various types of *phase closure* and *amplitude closure*
5. Instrumental and observations effects and how to calibrate them or verify that they are already properly calibrated
6. Imaging or other localization techniques involving the visibility data

Stay spatial!
"""

# ╔═╡ 47af0e5e-8a25-11eb-2b6b-f14ef17fe0f7
PlutoUI.TableOfContents()

# ╔═╡ Cell order:
# ╟─782d9788-77da-11eb-357f-e3857ed0f0ec
# ╟─440400b0-7be9-11eb-0efb-c723129eeef2
# ╟─2c9addb8-7be9-11eb-0adf-452042d88a4c
# ╟─99f47c64-78df-11eb-0ee3-af29fbc46edf
# ╟─0416bb34-7d7c-11eb-3060-b1c01edbf110
# ╠═83191f96-77da-11eb-16bb-8bf06f6b60ed
# ╟─b05be61e-77da-11eb-3a82-edc1be04b007
# ╟─f3378276-7d7b-11eb-2495-030f08a01f0c
# ╠═11499c58-7d7d-11eb-36d5-7f793f9d7d5d
# ╟─c5ebb850-7d7d-11eb-11e9-29094ba7160b
# ╟─eac37b0a-7d84-11eb-1430-95092833fe33
# ╠═91707832-7d7f-11eb-0ca8-fbd7a47adda4
# ╟─86308de8-7d86-11eb-353e-d75a37cb7d1f
# ╠═13bb054c-7d89-11eb-34ca-0f84e57aa525
# ╟─5c84f00a-7d8c-11eb-31cd-e3f4157f534d
# ╠═c2d0b58e-7d86-11eb-1812-b96f2ce720b1
# ╟─9367d2fe-7d87-11eb-1eda-7f466558ff4d
# ╟─5ec2e7fc-7d8f-11eb-2f0a-11afe337cb4d
# ╟─7550d5ee-7ed7-11eb-3ce0-51c68ce5bf5a
# ╠═2c52245a-7ed8-11eb-0f76-d1591196f9d1
# ╠═f9627760-7ed8-11eb-009f-45e3f92053b7
# ╟─8eba61ce-7ed9-11eb-08c5-4ffcb667f664
# ╠═1b46bd5a-7ed9-11eb-2661-75f2f5d4987f
# ╟─e87c0060-7eec-11eb-3607-0beac189151e
# ╟─ea4a03f6-7eda-11eb-0499-77b6062ea04e
# ╟─d0c7bd6e-7f10-11eb-07e4-65b181fa87b2
# ╠═0dda68f2-7f14-11eb-31eb-b5b3f0e0acef
# ╟─4fc456e2-7f14-11eb-3b11-3159926f69fc
# ╟─9db1fed6-7f19-11eb-21a7-bb7fb8e19108
# ╟─6f4a42be-7f1a-11eb-0ee4-fb4f2108ba47
# ╟─c91745c0-7f16-11eb-3b9c-6dd85327d520
# ╟─f7f47dea-7f1d-11eb-1d53-97a4821b4623
# ╟─33992d3a-7f1e-11eb-1be4-9526e6e17743
# ╟─3b6083bc-7af8-11eb-1c69-8df28a0ac02f
# ╟─b696f722-7fef-11eb-0be8-23af42681850
# ╟─1404ed9e-7a6c-11eb-2883-13056a01a8f0
# ╟─2e28c1a8-7ff0-11eb-1455-0be7a7848000
# ╟─5ba19006-7bc0-11eb-08b2-838f828988dd
# ╟─47af0e5e-8a25-11eb-2b6b-f14ef17fe0f7

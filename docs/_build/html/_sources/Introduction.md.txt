(Introduction)=
Introduction
============

Most X-ray binaries (XRBs) exhibit spectral state transitions during their outburst phases. When these sources awaken from quiescence, their brightness increases dramatically, and they trace out a characteristic q-shaped pattern in the hardness-intensity diagram as they evolve from the hard to the soft state and eventually back to the hard state {cite:p}`2010LNP...794...53B`.

Typically, the transition begins in the *hard state*, where the spectrum is dominated by a power-law-like spectrum with a high-energy cutoff, produced by Compton up-scattering of seed photons in a hot corona {cite:p}`1975ApJ...195L.101T,1979Natur.279..506S`. In this state, reflection signatures also emerge as coronal photons irradiate and are reprocessed by the accretion disk {cite:p}`1991A&A...247...25M`. As the outburst progresses, the system enters the disk-dominated *soft state*, characterized by a multi-temperature blackbody spectrum from the accretion disk {cite:p}`1973blho.conf..343N`. During the transition from the hard state to the soft state, sources often pass through an *intermediate state*, in which quasi-periodic oscillations (QPOs) are frequently observed {cite:p}`2019NewAR..8501524I`. The presence of different types of QPOs further allows this intermediate state to be subdivided into the hard-intermediate state (HIMS) and the soft-intermediate state (SIMS).

The model `DAO` focuses on X-ray reflection, which is generally observed in the *hard state* {cite:p}`2023ApJ...951..145L`. It computes the reprocessed radiation from disk–corona systems, similar to well-known reflection models such as `reflionx` {cite:p}`2005MNRAS.358..211R` and `xillver` {cite:p}`2010ApJ...718..695G,2013ApJ...768..146G,2014ApJ...782...76G`. The main goal of `DAO` is to provide an accurate, flexible, and community-accessible reflection code that enables researchers to perform studies tailored to their specific scientific needs.

In this document, I aim to include all relevant information about the model—such as its usage, approximations, and underlying physics—and to continuously update this user guide. Please feel free to contact me if you encounter any questions or issues.

After you have installed, we recommend the following steps to get started with `DAO`:

1. Consult the [Quick Start](quickstart_guide) for a basic execution tutorial.

2. Review the [Model Input Parameters](InputParameters) and [Output Files](outputfile) sections for detailed specifications.

3. Utilize [Python scripts](runcodewithpython) to manage your workflow and see [Example](example) for more details practices.

4. For a detailed description of the physics and mathematics behind DAO, please refer to the [Theory](theory) section.
# The Impact of Satellite Constellations on Solar System Science by LSST
This project aims to assess the potential loss of Solar System Object (SSO) discoveries for the Vera C. Rubin Observatory's Legacy Survey of Space and Time (LSST) due to the presence of bright satellite trails from Low Earth Orbit (LEO) and Medium Earth Orbit (MEO) satellite constellations, such as Starlink and OneWeb.

## Overview
With over 400,000 satellites planned for launch in the near future, ground-based astronomical surveys like the LSST face challenges from satellite trails in their observations. This project simulates the impact of satellite constellations on LSST observations, quantifying the data loss and potential missed SSO discoveries.

## Objective
The primary objective is to quantify the data loss for the LSST survey caused by artificial satellite constellations, specifically SpaceX's Starlink V1 & V2, and OneWeb. The project involves the following steps:

Propagate Two-Line Element (TLE) data for satellite constellations.
Check each satellite against each field of the LSST survey.
Identify observations affected by satellite trails within a 3-arcsecond distance.
Account for Earth shadowing effects.
Determine the loss of SSO discoveries based on the minimum required observations.

## Findings
Preliminary results for the Starlink V2 constellation (29,988 satellites) propagated against the first 20,000 LSST fields (approximately 29 days) show:

A total of 836 observations (0.1% of all observations) are affected by satellite trails.
Assuming providers limit satellite brightness below the recommended 7th magnitude, about 6% of SSOs that should be discovered will be missed in the presence of the planned 400,000 satellites per month.
This translates to a detection deficit of about 0.15% per 10,000 satellites per month.

## Future Scope
The project aims to extend the analysis to cover the entirety of the LSST duration and create a census of all affected observations and missed objects. Additionally, the following aspects will be explored:

Analyze the impact of the entire Starlink constellation and the OneWeb constellation.
Simulate a worst-case scenario with 400,000 satellites and measure its impact on SSO data loss.
Inform mitigation strategies and provide assessments for specific SSO populations.
Improve the model's fidelity by better modeling satellite brightness and its effects on trail width.

## Acknowledgments
The authors acknowledge the support of the Heising-Simons foundation through the LSST Kickstarter Grant #KSI-9, collaboration with the LSST Solar System Science Collaboration, and the IAU Centre for the Protection of the Dark and Quiet Sky From Satellite Constellation Interference. This work utilized the Illinois Campus Cluster and various software tools like SGP4, Skyfield, Astropy, Matplotlib, Numpy, and Spiceypy.

## References
Tyson, J.A. et al. (2020) "Mitigation of LEO satellite brightness and trail effects on the Rubin Observatory LSST," The Astronomical Journal, 160(5), p. 226.
Hu, J.A. et al. (2022) "Satellite constellation avoidance with the Rubin Observatory Legacy Survey of Space and Time," The Astrophysical Journal Letters, 941(1).
Grav, T. et al. (2011) "The Pan-STARRS Synthetic Solar System Model: A tool for testing and efficiency determination of the Moving Object Processing System," Publications of the Astronomical Society of the Pacific, 123(902), pp. 423–447.

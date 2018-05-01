We thank the referee for a detailed and thorough reading of our paper.
We have corrected all the major deficiencies in the presentation of
the methods and the results.
We have also included an explicit estimate of the main results as a
function of satellite number. 
All the major conclusions quoted in the first version of the paper remain
unchanged. 
In the following we quote the referee with the symbol ">". 
Our reply to each comment follows immediately after.
We attach the file paper_diff.pdf highlighting the changes.

With best regards, 
Jaime and Veronica

> Major comments:
> section 2.2, page 4 and section 2.3, page 4, end of first column:
> The authors say Appendix A shows physical properties like e.g. stellar
> masses of the halo pairs. I can't find any stellar mass plots or
> listings. Furthermore, ELVIS is a DMO simulation suite so the only
> simulation with stellar masses is Illustris. Please, add the stellar
> mass information for solely Illustris or be consistent throughout the
> paper and remove the reference to stellar masses. 

We have removed the reference to stellar masses.

> * section 3.1:
> The authors select 11 up to 15 satellites. While for MW 11 satellites
> corresponds to the "classical MW satellites" nowadays we know many
> more satellites (>30). The number of 15 corresponds to "the minimum
> number of M31 satellites usually included in M31 studies". While this
> might be true for the number of plane members of studies of the plane
> of satellites around M31 by e.g. Ibata+2013 these authors indeed used
> ALL the satellites found in the Pandas footprint. However, there are
> many more satellites known around M31 and the 15 in the plane are not
> necessarily the most massive/luminous ones. Thus, it would be far more
> informative to calculate the asphericity of MW and M31 satellites for
> the classical/conventional number of satellites (11 or 15) and the
> whole (known) satellite population to see if asphericity is the result
> of low number statistics or if there is indeed a dependence on
> satellite mass. E.g. are the most massive satellites preferentially
> distributed in a planar fashion? 

We calculate now the asphericity at a fixed number of satellites. 
The results we examine in detail correspond to N=11, but we also quote
the results for N=12, 13, 14 and 15 to examine a possible dependence
on satellite number. 
We do not make the calculation for the whole population because the
resolution of the Illustris1 simulation does not allow us to make such
an estimate. 

For the MW we do not find any obvious trend. However, for the M31
there seems to be a trend for N>=12.  The new Figure 4 summarizes these results.

We also emphasize in the introduction that our study does not
address questions about M31 planes. It only addresses the
distribution of a given set of satellites with trivial selection
criteria. 


>* section 3.3:
>The authors randomize the position of satellites while keeping their
>radial distribution fixed to compare satellite distribution against
>their own randomized distribution. Different halo pairs have
>different masses. How do the authors account for the fact, that a
>more massive halo has on average a 'wider' satellite distribution
>with on average larger radii for the satellites? Would the results
>change if a rescaling, like e.g. rescaling/normalizing all satellite
>radii by R200 of the host would be done? Could this take out the halo
>mass dependence? 

Rescaling/normalizing all satellite radii by R200 produces the same
results after normalizing by the randomized distribution.  
Scaling positions by R200 keeps the asphericity of a halo unchanged.
The main motivation to renormalize each satellite distribution by its
own randomized distribution is precisely to factor out the effect of
the physical size. 

To put it in a different way, the normalized results
are sensitive only to the degree of asphericity and not to the
physical units.
This is more evident now in Figure 2 and Figure 3. 
Figure 3 shows the aesphericity in physical quantities, there is
visible a big difference between the MW and M31 results. 
Figure 3 shows the same in normalized quantities, there the
differences are more subtle to spot due to the rescaling. 

> * section 4.4 and 5, end of third paragraph:
> The authors say the results of the multivariate gaussian fitting are
> robust agianst changes in the simulation. However, as Figure 5 clearly
> shows there are large differences in the number of LG-like systems in
> the simulations (e.g. M31: 2700 Illustris, 5688 Illustris dark and
> 3700 ELVIS). Furthermore comparing the corner plots (Fig. 4, C2, C4)
> for the different simulations one can clearly identify major
> differences among the simulation setups. Correlations in the plane
> parameters are much weaker in the HYDRO case and for the b/a vs. c/a
> case correlations even change from slightly correlated to
> anti-correlated. This can be discussed in detail at the end of section
> 4.4 and 4.5 where already an attempt of discussion was made which can
> clearly be expanded.

What we mean is that the results of different simulations are always
consistent with gaussians, not that the parameters are similar.
 The results for the covariance matrix and mean vector differ widely
 among simulations, as the referee points out. We have changed the 
discussion section accodingly.


> Also the DMO simulations Illustris-dark and ELVIS show vastly
>  different correlations (e.g. b/a vs. w or b/a vs c/a). Elaborating
>  on these differences and working out/understanding where they are
>  coming from is needed. In the case of Illustris-dark
>  vs. Illustris-hydro cosmological parameters, simulation code setup,
>  etc. are the same so its solely hydro effects which might vastly
>  alter the correlations and their strength. However, the differences
>  between Illustris-dark and ELVIS must be of different
>  nature. Undertsanding and quantifying these effects is needed to
>  interpret results obtained from DMO simulations in order to compare
>  them with observations. An attempt is made in section 4.5 where
>  median halo masses are compared. This could be expanded a bit.
> In the case of Illustris vs. ELVIS the reason might lie in
>  resolution (simulation code, difference of WMAP 9 vs WMAP 7) or
>  selection effects introduced in the selection of ELVIS haloes. Here
>  it would be greatly appreciable if the vastly different
>  correlations for e.g. b/a vs. w could be discussed and
>  investigated. 

We agree that the results from ELVIS are puzzling given the vast
difference with respect to the Illustris simulations. 
We consider such a study far beyond the scope of this paper. 

This would require, at least, 1) having larger halo samples of halos in the
same resolution as Illustris1 to explore asphericity as a function of halo
mass. Right now we only have ~5 halos in the same mass range as
ELVIS. 2) Similar simulations as Illustris1 in volume with higher
resolution to explore the effect of numerical resolution. 

These two points could be addresed in detaile with the new IllustrisTNG
simulations, unfortunately the data is still not yet public. 
We highlight that our main objective is presenting a new tool and demonstrating 
its power to spot differences as the one we have just discussed between ELVIS and Illustris1. 

> * section 4.5, line 59:
> Here and in the abstract it says at most 2% of the pairs have
> asphericity comparable to the LG. This doesn't match the values stated
> in Figure caption 5 which says that only 1% of the LG analogues have
> properties similar to what is observed. Please clarify how the
> probability is derived and make sure the paper is coherent and not
> contradicting itself. 

We clarify now these probabilities as a function of satellite number.
We have checheked that we don't have errorrs in the number quoted in
the different sections.

Moderate comments:

>* section 3.1:
> The authors select satellites in the simulation ranked by maximum
> circular velocity while in observations the satellites are selected by
> luminosity. Could they please elaborate a bit on why this selection
> results in a fair comparison or at least state some references. Which
> assumptions go into this? 

We have expanded and added some references.


> * section 4, second paragraph, discussion Table 3,4:
> ".. confirms ... more spherical distribution for the M31." How does
> this go together with the claim that there is a plane of satellites
> around M31? How exactly are the 15 satellites for the M31 case
> chosen? The ones in the plane or the 15 most massive/luminous? If the
> 15 most massive/luminous are chosen, this explains why M31 is not
> planar. The satellites in the plane of M31 are not the most massive
> ones. This shows how important it is to clearly state the selection
> procedure of satellites and that a general characterization of
> asphericity (including all known satellites) is more informative. 

We stress now in the introduction that we do not find planes of satellites
and our results for M31 should not be compared against sophisticated
satellite selections, but only characterize asphericity in distributions ranked by
vmax, luminosity or mstar.



> * Figure 1, 2, 3:
> (a) What exactly do the 5 stars show? Satellite numbers from 11 to
> 15? If so, > please state that clearly. It is not obvious.

In the older version of the paper each star showed the results for 11,
12, 13, 14 or 15 satellites. We now present the results at fixed
N_satellite throughout the paper.

> (b) While it is a nice idea to show MW and M31 on x and y axis, it
  is not at all clear why this is informative. It implies a
  correlation/dependence between, e.g. the plane width of MW and
  M31. Furthermore, as the ELVIS paper (Garison-Kimmel+2014) states,
  there are no differences in satellite abundances and kinematics
  between paired and isolated haloes. Naively, from this one would
  expect that the width of planes in paired haloes should be
  independent. So, could the authors please spend some words on
  explaining the choice of visualization? I understand that this is a
  compact way of showing the results but a bit more discussion helps
  the reader to understand the authors's point. 

To avoid confusion we have changed the visualization and show separately
the results for the MW and M31. 
Indeed, we have checked that there isn't any correlation between the
two and mention that in the paper.

> * Figure 4 and Figure C2, C4:
> Axis limits are different for all corner plots. This makes a direct
> comparison difficult. I would advice to make the axis limits the
> same for all plots (or otherwise remind the reader to pay attention
> to this). 

Axis limits are now the same for all plots.

> * section 5, last two paragraphs:
>If true, this result is very nice and fits the other hints that MW
>might be special. But then the question is, do any of the
>particularities stated here correlate with asphericity? How do these
>findings influence the probability of finding aspherical satellite
>distributions? 

To address this question using simulations one would need a large
sample of galaxies to estimate these correlations. We have added an estimate in
the discussion about the volume required from a cosmological
simulation to study these correlations. 

> * section 4.5, last paragraph:
> values for plane width and c/a ratios do not match the values in the
> corner plots. Please be consistent!!! 

We have fixed this inconsistency. 

>* Figure A.1:
> The quality of the plots would highly be improved if the authors
> would add lin> es for the observed/estimated values and their
> uncertainties for the MW and M31. 

We have included now these uncertainties. 


Minor comments:

* Abstract:
> The authors say "not enough systems fully resembling the LG". Could
> you in the main text explain what you exactly mean by this? Number of
>satellites? Radial distribution of satellites? Galaxy pairs? Any clue
>why this is the case?  

Clarified.




* introduction, third paragraph:
> observational studies have found... There is another study by
> Collins+2015 finding no difference between satellites in plane and off
> plane for M31 satellites. 

We have included this reference. in our discussion about that plane.

* section 1, line 36:
> What exactly do the authors mean by "explicit probability
> distribution"? Can this term be elaborated>  in the main text?

We have changed "explicit" by "analytical". 



* section 2, second paragraph:
> Why are only the brightest 11 to 15 satellites used? And what does
> 11 to 15 mean? 11 satellites for MW and 15 for M31 or is the whole
> range from 11 to 15 been used for both galaxies (in observations and
> simulations)? 

We clarify in that section that is due both to historical reasons (the
11 classical satellites, the 15 satellites that are at minimum
considered for M31) and resolution of the simulations we use.

* section 2.2, second paragraph:
> baryonic mass resolution of Illustris is 1.2e6 Msun not! 8.0e7!!!

We corrected this error in the text.

* section 2.2, selection of LG counterparts:
>The authors discard haloes with a number of satellites smaller than
>15, which are the lowest mass haloes. Do the authors have any
>intuition how this might affect the results of the study? The authors
>say already in the abstract, there is a dependence on halo mass but is
>there any trends already seen which might hint at the effect of
>discarding lower mass LG analogues? 

In the current form of the test we do not discard those halos. We
include all halos with the same number of satellites as we
require. The main results do not change.

* Table 3:
> How exactly are the normalized values calculated. If I follow the
> description in the text I get the following result for e.g. M31: 
> (59-65)/12=-0.5 not -0.48 as the table says. Is the error on the
> normalized value just the standard error propagation or how is this
> calculated? 

THe error on the normalized value is the error on the physical
quantity normalized by the standard deviation in the randomized
sample.


* Table 4:
> It would be good if normalized results for the simulation (e.g. from
> the corner plots Fig.4, C2, C4) appear as well in the table. 

Now we only include plots for all quantities (physical, normalized,
gaussian model)

* section 4.1, last two sentences:
> missing word, general confusing, what do you want to say with these two sentences?

We corrected this in the text.

* section 4.1:
> Can the fact that MW (as the lower mass halo of the two in the pair)
> has a thinner plane be explained by a halo mass dependence? (see also
> point above of rescaling by R200) 

The problem is not only a thinner plane, but also a small c/a ratio,
i.e. the problem is the extreme asphericity. This is present in the
normalized quantities which do not depend on the physical units (see
our response of a rescaling by R200) and already take into account the
different sizes of the different halos.

* section 4.2, last sentence:
> close to between two and three... Is it close to or between? Please rewrite.

It was "between" and is already corrected in the text.

* section 4.5, first paragraph:
> When drawing samples from the multivariate gaussians you assume halo
> pairs can be drawn independently? Is this assumption valid? Do you
>find any dependence of the results on the companion, its mass,
>distance, relative velocity, etc.? E.g. does a more massive M31 result
>in a thinner plane? Please discuss a bit more the procedure and
>assumptions. 

We checked explicitly that the correlation between the halos can be
discarded, i.e. we ran everything using a six dimensional gaussian
distribution instead of two 3d gaussian distribution.
We don't find any clear trend of the asphericity with the mass or
companion distance.
This is probably due to the small number of pairs we have.
We mention that in the text. 

* section 4.5, second last paragraph:
> median halo masses are stated for the simulated analogues. Can the
> authors compare to observed/estimated values for MW and M31. It looks
> like the MW and M31 analogues are on the very lower end of estimat
> ed values for the LG (~220 km/s for MW and ~260 km/s for M31
>  e.g. Sofue+2016). 

We have included that comparison the figures of Appendix A.
Such a difference can be expected for 

* section 4.5, last paragraph:
> haloes get more spherical when baryonic effects are included, this was already found before (e.g. Bryan+2012).

It was found in the shape sampled by all the particles in the halo but
not in the satellite distribution. We have included that reference in
the paper. 

Spelling errors:

> * abstract, second sentence: to to
> * abstract, second sentence: atipicallity --> atypicality
> * abstract, sentence 6: at most of 2% the --> at most 2% of the
> * section 1, first sentence: a explicit --> an explicit
> * section 1, more recent numerical experiment --> more recent numerical experiments
> * section 1, second page: spherically symmetric distribution --> spherically symmetric distributions
> * section 1, second page: seemly --> seemingly
> * section 1, second page: incovenients --> inconveniences
> * section 3.1, first pragraph, last sentence: to to
> * section 3.1, second last sentence: included M31... --> included in M31...
> * section 3.2, equation 1: usage of k and i as indices, be consistent
> * section 3.3, dividing between --> dividing by
> * section 3.4, third paragraph: out of the the --> out of the
> * Figure caption 1: full stop missing after "pair".
> * Figure caption 4: spaned --> spanned
> * section 4.4, second paragraph: this result summarizes the result... doubled the word result
> * section 4.5, line 52: sampe --> sample
> * section 4.5, line 54: as aspherical as the observed in M31 --> as aspherical as the observed one in M31
> * section 5, page 9, line 35: distributions --> distribution
> * section 5, page 9, line 39: needing --> the need
> * section 5, page 9, line 45/46: of influence     doubled
> * section 5, page 9, line 46: incluedes --> includes
> * Figure caption A1: all lots --> all plots

The spelling erros were corrected in the text.

> Further references the authors might want to look at: *Gillet+2015,
> *Buck+2015, Buck+2016, Lipnicky+2017, Ahmed+2017, Fernando+2017  
> the ones marked with * are particularly interesting for table 2.  

All these reference focus on finding/charachterizing planes similar to
the one found in M31 which is a different direction from the work we
are presenting and don't belong to Table 2. However, we have included
them in the introduction.   


# Scientific basis for visualizing the PLA–ChAT exports

The available data can describe where HALO detected objects, how those objects differ in size and intensity, and how their reported overlap with ChAT varies across four named regions. The strongest presentation combines spatial maps, clearly defined overlap summaries, full distributions, and detection-quality checks. Biological interpretation requires the sample identities, image evidence, and experimental controls that connect these measurements to the experiment.

The evidence for this dataset is [abstract.txt](../abstract.txt), the four local CSV exports, and their documented inventory in [README.md](../README.md). The abstract describes a tissue/animal/condition structure, but the exports do not identify that structure explicitly. `22M`, `24M`, `Pfkp-TDP43 PLA`, and `Hk1-TDP43 PLA` therefore remain literal labels. One shared image reference does not determine how many specimens, sections, or animals contributed.

## Biological context and interpretation

In situ proximity ligation assay uses antibody-linked DNA probes to generate localized amplified signals when the probes bind sufficiently close together. The original method demonstrated visualization of endogenous protein complexes in cells and tissue. This explains why the exported PLA objects are scientifically interesting; it does not establish the exact protocol, specificity, or molecular counting accuracy of these particular exports. [Söderberg et al., 2006](https://pubmed.ncbi.nlm.nih.gov/17072308/).[^1]

**PLA target proximity and PLA–ChAT overlap answer different questions.** The target-pair label names the proteins assayed for proximity. ChAT overlap describes where the resulting segmented signal lies relative to a second marker. A positive overlap flag cannot by itself establish direct binding, a neuron identity, mitochondrial localization, or glycolytic activity. HALO's module documentation describes object measurements and colocalization, consistent with treating these as image-derived objects rather than automatically interpreting every row as a cell. [Indica Labs, Object Colocalization FL](https://indicalab.com/halo/halo-modules/object-colocalization-fl/).[^2]

Two primary studies explain why TDP-43 and glycolytic enzymes are plausible research targets. Barone and colleagues reported that cytoplasmic TDP-43 binds and sequesters HK1, disrupts its mitochondrial association, and reduces glycolytic capacity in their investigated ALS models. Their evidence included cell models, patient-derived motor neurons, mouse models and human postmortem tissue, with biochemical and functional assays alongside microscopy. This is useful motivation for follow-up experiments; it does not identify the present exports as belonging to that study. [Barone et al., 2026](https://link.springer.com/article/10.1007/s00401-026-02996-6).[^3]

Manzo and colleagues found evidence for compensatory glycolysis upregulation and protective effects of enhancing glucose metabolism in TDP-43 models, particularly Drosophila, with additional observations in human-derived material. Their PFK results do not establish direct PFKP–TDP-43 binding. The different models and measurements in these studies also make a universal “more PLA means more glycolysis” interpretation untenable. The present CSVs contain no metabolic flux, enzyme-activity, disease-status, species, or perturbation measurements. [Manzo et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6557627/).[^4]

## A spatial reconstruction that remains faithful to the evidence

Use the midpoint of each exported bounding box, `(XMin + XMax)/2` and `(YMin + YMax)/2`, as an **approximate object position**. These midpoints are neither measured centroids nor complete object masks. A grayscale raster of ChAT object positions can provide a useful background for a colored PLA count or overlap map, but it must be labeled an **object-derived spatial reconstruction**, not a microscopy image. True cell shapes, tissue texture, original fluorescence, and missing anatomy cannot be recovered from these fields.

Keep equal X/Y aspect and show coordinate axes. A downward Y direction is a display convention until the original coordinate frame is verified. Use “export coordinate units”; integer coordinates alone do not prove pixel units. OME metadata can specify physical pixel sizes, dimensions, channel identities, and spatial information, making the original OME-TIFF and its metadata the correct sources for calibration. Do not derive a micrometer scale bar from object area or bounding-box dimensions. [OME schema documentation](https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html).[^5]

For an eventual image overlay, establish image series, pyramid level, crop origin, axis order, channel identity, and any registration transform first. Review several exported boxes against the actual image at distinct positions. Preserve an unmodified image panel next to the overlay, and disclose display transformations. Publication guidance requires image processing to remain faithful to the underlying evidence. [Nature image-integrity policy](https://www.nature.com/nature/editorial-policies/image-integrity).[^6]

## The visual sequence and its questions

| View | Question it answers | Required explanation |
| --- | --- | --- |
| Object inventory | How many PLA and ChAT objects were exported? | Counts are detections, not tissue-normalized abundance or cell numbers. |
| Spatial object reconstruction | Where are detections located within each region? | Background and overlay both come from coordinates; anatomical identity is unknown. |
| PLA count heatmap | Where do many PLA objects occur within equal coordinate bins? | Color means objects per bin, not fluorescence or objects per mm². |
| ChAT-associated PLA fraction heatmap | Where is a larger share of local PLA associated with ChAT? | Denominator is local PLA objects; low-support bins are masked. |
| Overlap distribution | Are objects nonoverlapping, partially overlapping, or nearly fully overlapping? | Zero and 100% values remain visible. |
| Area and intensity distributions | Are regional differences broad shifts or driven by a tail? | Full distributions and medians complement totals. |
| Area–intensity hexbin | Do size and brightness covary, or form suspicious bands? | Color encodes object count; the pattern is descriptive. |
| Quality and sensitivity views | Could export quirks or display choices explain a pattern? | Show zero diameters, collapsed bounds, support and bin-size changes. |

A count heatmap and fraction heatmap should be read together. A bright count bin can arise from more detected material or greater sampling. A bright fraction bin can arise from only a handful of objects. Neither alone demonstrates molecular enrichment above an appropriate background.

## Spatial binning, support and smoothing

Use square bins with a stated width and origin, and retain the same coordinate-bin width when comparing counts across regions. Independently fitting a fixed number of bins to each region changes the sampled area per bin and compromises count-color comparisons. Keep the color normalization shared when a figure invites numerical comparison.

For bin `b`, define `N_b` as the number of PLA-type objects and `K_b` as the number of those objects with the ChAT-present flag. Display `100 × K_b/N_b` only where `N_b` meets a stated minimum. A threshold such as 10 or 20 objects is a display-support choice, not a biological or statistical significance cutoff. Show how the picture changes with this threshold and with finer/coarser bins; report the fraction of PLA objects retained in visible bins.

An empty bin has an undefined PLA-associated fraction. A bin with PLA objects and no ChAT association has a valid zero fraction. Give these cases different appearances. Observed ChAT occupancy is also **not a tissue mask**: tissue can exist without a detected ChAT object. The outer rectangle or convex hull of detections is not a validated analyzed-tissue boundary.

Prefer unsmoothed bins for quantitative reading. If smoothing is offered, retain an unsmoothed comparison and disclose kernel width in bins or verified physical units. For a fraction, smooth counts separately and then divide, `G(K)/G(N)`, rather than smoothing already calculated percentages; the latter gives sparse bins excessive influence. This is an analysis recommendation derived from the ratio's definition. Keep low-support regions masked and avoid interpreting apparent signal spread across empty space. SciPy documents the strong influence of bandwidth on density estimates and warns that smoothing can conceal multimodal structure. [SciPy KDE documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html).[^7]

Formal clustering or cross-type association tests should await calibrated positions and valid observation windows. Lagache and colleagues' SODA study shows why geometry, object density, and boundary correction matter: confined objects can be close by chance. Its molecular-localization assumptions do not automatically apply to irregular ChAT objects represented by box midpoints. A nearest-neighbor line is not an exported object match; a randomization over a bounding rectangle is not an adequate tissue-aware null model. [Lagache et al., 2018](https://www.nature.com/articles/s41467-018-03053-x).[^8]

## Distribution summaries and denominators

Use an empirical cumulative distribution function (ECDF) to answer “what fraction of objects is at or below this value?” Its full-data version avoids histogram-bin choices and kernel smoothing. Keep object types separate and use a logarithmic area axis when the positive area range is wide. Preserve zero intensity or diameter observations explicitly rather than silently dropping them for a logarithmic axis. Matplotlib distinguishes its exact ECDF from the binned approximation produced by a cumulative histogram. [Matplotlib cumulative distributions](https://matplotlib.org/stable/gallery/statistics/histogram_cumulative.html).[^9]

For a million rows, ordinary scatterplots hide points beneath one another. A hexbin aggregates nearby observations; logarithmic count coloring can reveal both dense and sparse areas. State axis transformations and color units, and expose the number of excluded or nonfinite observations. Use the same axes for comparable panels. [Matplotlib hexbin documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hexbin.html).[^10]

Three overlap quantities deserve separate labels:

- **Associated-object fraction:** number of PLA objects with ChAT present divided by all PLA objects.
- **Mean object overlap:** average of individual PLA overlap percentages; each object receives equal weight.
- **Area-weighted overlap:** sum of reported PLA overlap areas divided by summed PLA object areas; larger objects receive greater weight.

The third is a summary of exported object areas, not tissue coverage or necessarily the union of overlapping pixels. Do not add the PLA-side and ChAT-side overlap totals. Their denominators describe different object populations, and explicit object-to-object matches are absent.

The CSVs cannot recover pixelwise Pearson or Manders colocalization coefficients: those require the appropriate channel images and analysis masks. Correlating row-level area and average intensity answers a different measurement question. Lagache et al. distinguish fluorescence correlation, physical overlap, and distance-based association; these should not be presented as interchangeable validation of protein interaction. [Lagache et al., 2018](https://www.nature.com/articles/s41467-018-03053-x).[^8]

The documented associated-object fractions are 59.04% and 68.91% for the `22M` and `24M` Pfkp exports, and 42.44% and 55.46% for the corresponding Hk1 exports. These are observed regional differences. They do not quantify an age effect or establish that one target pair interacts more strongly, because region coverage, specimen identity, target accessibility and acquisition/analysis comparability are unresolved. [Local CSV inventory](../README.md#counting-pla-and-colocalization-correctly)

## Replication, controls and the next analysis stage

The million object measurements are not a million independent biological replicates. Present exact descriptive counts, fractions and distributions without object-level significance stars or confidence intervals implying animal-level precision. Once sample metadata is available, retain the hierarchy from animal or donor through section, region and object. Summarize at the independent experimental unit, or use a justified hierarchical model; display individual biological-replicate summaries alongside within-replicate distributions. This follows the SuperPlots approach to separating measurement variability from experimental reproducibility. [Lord et al., 2020](https://rupress.org/jcb/article/219/6/e202001064/151717/SuperPlots-Communicating-reproducibility-and).[^11]

Prioritize image review before imposing size or brightness filters. Review randomly chosen locations as well as sparse, dense, unusually bright and extreme-size examples, with the selection rule disclosed. Single-antibody omission controls, appropriate biological negative/positive controls, and matched acquisition settings help distinguish assay signal from background. Duolink's manufacturer notes that amplification, antibody concentration and overexposure can merge signals, illustrating why spot counts and areas require protocol-specific validation. This guidance does not establish use of Duolink in the repository. [Duolink troubleshooting guide](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/protein-biology/protein-and-nucleic-acid-interactions/duolink-troubleshooting-guide).[^12]

Then recover analyzed-tissue masks for physical density, validated cell masks for puncta per cell, and compartment markers for localization questions. Investigating glycolysis needs separate functional measurements; investigating binding needs orthogonal evidence. Plan biological replication before hypothesis testing, with any power calculation based on variability among independent units.

For the notebook, keep a clear reading order, plain-language captions, visible parameters, exported supporting tables, and saved figures. Use perceptually ordered sequential color maps for nonnegative quantities and reserve diverging colors for a meaningful central reference. Avoid rainbow and red–green encodings, which can create false visual boundaries or exclude readers with color-vision deficiencies. [Crameri et al., 2020](https://www.nature.com/articles/s41467-020-19160-7).[^13] Preserve provenance, display transformations, parameter choices and analysis versions so every figure can be interpreted and recreated; these are central themes of the microscopy community's publication checklists. [Schmied et al., 2024](https://www.nature.com/articles/s41592-023-01987-9).[^14]

## Sources

The repository files are the primary evidence for dataset-specific statements. External sources provide mechanism, methods and reporting context. Online sources were accessed on September 13, 2026.

[^1]: Söderberg, O. et al. [Direct observation of individual endogenous protein complexes in situ by proximity ligation](https://pubmed.ncbi.nlm.nih.gov/17072308/). *Nature Methods* 3, 995–1000 (2006). doi:10.1038/nmeth947. Original assay paper; the linked PubMed record provides its abstract.
[^2]: Indica Labs. [Object Colocalization FL module](https://indicalab.com/halo/halo-modules/object-colocalization-fl/). Official module overview, undated. It does not specify this export's thresholds or configuration.
[^3]: Barone, C. et al. [TDP-43 impairs glycolysis by sequestering hexokinase 1 in amyotrophic lateral sclerosis](https://link.springer.com/article/10.1007/s00401-026-02996-6). *Acta Neuropathologica* 151, 26 (March 16, 2026).
[^4]: Manzo, E. et al. [Glycolysis upregulation is neuroprotective as a compensatory mechanism in ALS](https://pmc.ncbi.nlm.nih.gov/articles/PMC6557627/). *eLife* 8:e45114 (June 10, 2019). doi:10.7554/eLife.45114.
[^5]: Open Microscopy Environment. [OME 2016-06 schema: Pixels and related metadata](https://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2016-06/ome_xsd.html). Official schema documentation.
[^6]: Nature Portfolio. [Image integrity and standards](https://www.nature.com/nature/editorial-policies/image-integrity). Official editorial policy, undated.
[^7]: SciPy developers. [scipy.stats.gaussian_kde](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html). Official API documentation, live version.
[^8]: Lagache, T. et al. [Mapping molecular assemblies with fluorescence microscopy and object-based spatial statistics](https://www.nature.com/articles/s41467-018-03053-x). *Nature Communications* 9, 698 (February 15, 2018).
[^9]: Matplotlib developers. [Cumulative distributions](https://matplotlib.org/stable/gallery/statistics/histogram_cumulative.html). Official example and explanation, live stable documentation.
[^10]: Matplotlib developers. [matplotlib.pyplot.hexbin](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hexbin.html). Official API documentation, live stable documentation.
[^11]: Lord, S. J., Velle, K. B., Mullins, R. D. and Fritz-Laylin, L. K. [SuperPlots: Communicating reproducibility and variability in cell biology](https://rupress.org/jcb/article/219/6/e202001064/151717/SuperPlots-Communicating-reproducibility-and). *Journal of Cell Biology* 219:e202001064 (2020).
[^12]: Merck / Sigma-Aldrich. [Duolink PLA Troubleshooting Guide](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/protein-biology/protein-and-nucleic-acid-interactions/duolink-troubleshooting-guide). Official manufacturer guidance, undated.
[^13]: Crameri, F., Shephard, G. E. and Heron, P. J. [The misuse of colour in science communication](https://www.nature.com/articles/s41467-020-19160-7). *Nature Communications* 11, 5444 (2020).
[^14]: Schmied, C. et al. [Community-developed checklists for publishing images and image analyses](https://www.nature.com/articles/s41592-023-01987-9). *Nature Methods* 21, 170–181 (2024 issue; published online September 14, 2023).

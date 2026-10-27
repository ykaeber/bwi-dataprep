# BWI data notes

Everything here was verified against the actual Thünen release
(`bwi20150320_*`, data status 2015) rather than taken from documentation.
Where a decision was made rather than discovered, the reasoning is given.

## 1. What the data is

The Bundeswaldinventur is a large-area sample inventory of German forest.
Four inventories exist:

| Inventory | Field years | Coverage |
|---|---|---|
| BWI 1 | 1986–88 | old federal states only, different design |
| BWI 2 | 2001–02 | all of Germany |
| BWI 3 | 2011–12 | all of Germany |
| BWI 4 | 2021–22 | all of Germany (not in this release) |

This project uses the **BWI 2 → BWI 3 transition**, which is the only pair in
the public release with a consistent design, permanent sample trees and
harmonised change tables.

The *erweiterter Grunddatensatz* ("extended raw dataset") is published
anonymised, without coordinates, and **needs no data-use agreement**:
<https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/>

Tract coordinates are released separately, snapped to the 1 km INSPIRE grid:
`ShapeFile_Tnr_INSPIRE_Poly1000.zip`.

Citation: see `data/raw/Zitiervorschlag-BWI-Basisdaten-ACCESS2003-(2015).txt`.

## 2. Sampling design

- **Trakt** — a square with 4 corners (*Ecken*, `Enr` 1–4 = A–D) 150 m apart,
  on a 4 × 4 km national grid, denser in states with a *Verdichtung* (`VBl`).
- **Ecke** — the sampling point. Where it falls in accessible forest, trees are
  tallied by **Winkelzählprobe with Zählfaktor 4 m²/ha** for all stems with
  BHD ≥ 7 cm.
- **Regeneration** is recorded separately on fixed-radius circles of r = 1 m and
  r = 2 m (`b3_jung`, `Bnr` 1xxx and 2xxx).

Angle-count sampling has no fixed plot area. A tree is tallied when its
horizontal distance from the corner is below 25 · d (d in m), so inclusion
probability is proportional to basal area. What the design yields is a density:

```
n_ha = K / g          stems per hectare represented by one tallied tree
```

with `g` the stem's basal area in m² and `K` the count factor.

**Verified:** `N_ha_Hb == 4 / G` exactly for 98.1 % of BWI3 stand trees. The
remaining 1.9 % carry a Grenzkreis (edge) correction, with a median corrected
`K` of 4.96 m²/ha. The pipeline therefore uses the published `N_ha` rather than
recomputing `4/G`, so Thünen's edge correction is preserved.

`N_ha` vs `N_ha_Hb`: `N_ha` is the **within-stand** density and is 0 for stems
flagged as not belonging to the stand (`Bz = 0`); `N_ha_Hb` is the across-stand
density on wooded land. The pipeline uses `N_ha`, because FINN models a
homogeneous patch and the within-stand normalisation is what "the stand this
corner samples" means. `N_ha` is available in every table we need, at both
dates; `N_ha_Hb` is not.

## 3. Tables used

From `bwi20150320_alle_daten2002-2012.mdb`:

| Table | Contents |
|---|---|
| `b0_ecke` | permanent corner info; `Use23` is Thünen's own flag for whether a corner may enter a BWI2→BWI3 change analysis (100 = usable) |
| `b23_baeume_a_s` | tree state at the **start** of the period, measured where the tree was a BWI2 sample tree, back-modelled (Sloboda) where it was not |
| `b23_baeume_m_s` | same key, carries `Pk` — the tree's fate over the period |
| `b23_ecke_w` | forest-cover change over the period |

From `bwi20150320_alle_daten2012.mdb`: `b3_baeume` (end-of-period state),
`b3_ecke_w`, `b3_jung`, `b3_bestock`, `b3_tab_w` (field dates), `b0_tab`
(per-tract elevation `HoeheNN`).

Key throughout: `(Tnr, Enr, Bnr)`. Verified that `b23_baeume_a_s` and
`b23_baeume_m_s` have **identical key sets** (484,796 trees), and that survivors
and ingrowth appear in `b3_baeume` while every removal and death does not.

`VegJahreBaum` in `b23_baeume_m_s` is entirely NULL in this release — period
lengths come from `Datum_feld` in `b2_tab_w` / `b3_tab_w` instead (realised
intervals 8–12 years, mode 10).

## 4. Fate codes (`Pk`, code table `x_pk`)

| `Pk` | Meaning | Treated as |
|---|---|---|
| 0 | neuer Probebaum | recruit **or** lateral ingrowth — see §5 |
| 1 | wiederholt aufgenommen | survivor |
| 2 | selektiv entnommen | **harvest** |
| 3 | bei Kahlschlag entnommen | **harvest** |
| 4 | am Ort verblieben (fallen) | natural mortality |
| 5 | abgestorben (standing dead) | natural mortality |
| 9 | nicht auffindbar | dropped, corner flagged |
| 6, 7, 8, 10, 11 | out of sample / invalid | dropped, corner flagged |

Separating harvest (2, 3) from natural mortality (4, 5) is essential: FINN has
no harvest process. Over 2002–2012, on usable corners, harvest removed basal
area at **1.56 %/yr** against **0.26 %/yr** for natural mortality — harvest is
6× the signal FINN is meant to learn.

Corners carrying unclassifiable trees are flagged rather than silently pruned:
dropping a live tree would read downstream as mortality that never happened.

## 5. The ingrowth trap

**This is the single most important correction in the pipeline.**

`Pk == 0` ("neuer Probebaum") does *not* mean the tree is a recruit. Because the
angle-count inclusion zone grows as d², a tree that was merely too far from the
corner in 2002 gets tallied in 2012 without anything having happened to it
biologically.

Measured on the release:

- 113,237 trees carry `Pk == 0`
- of these, only **8.9 %** had a back-modelled 2002 diameter below the 7 cm
  caliper threshold — genuine recruits
- the other **91.1 %** already existed as measurable trees, median back-modelled
  2002 dbh **23.6 cm**

The design predicts exactly this: ~15 % diameter growth over the period inflates
the inclusion area by ~32 %, so ~24 % of the final tally should be newly
included. Observed: 29.5 %.

### Direct geometric confirmation

The BWI records each tallied stem's position (`Azi` [gon], `Hori` [cm] in
`b2_wzp` / `b3_wzp`), so the mechanism can be tested rather than argued. For
every stem, compare its distance from the corner against its inclusion radius
`25 · d` at each date:

| Class | inside the **2002** inclusion zone |
|---|---|
| survivors (tallied at both dates) | **100.0 %** |
| lateral ingrowth (`Pk` = 0, back-modelled 2002 dbh ≥ 7 cm) | 29.9 % |

The survivors' hard cut-off at exactly 1.0 validates the geometry itself
(`Hori ≤ 25·d` holds for 99.98 % of WZP4 stand trees, max ratio 1.000).

For lateral ingrowth, **70 % lie outside** their 2002 inclusion zone — they
demonstrably grew into the sample. The other 30 % sit just inside it, and the
ratio distribution is tightly pinned to the boundary (10th–90th percentile
0.81–1.22, median 1.05). That concentration *is* the evidence: a class of trees
entering for unrelated reasons would scatter, whereas these pile up on the
inclusion boundary. The 30 % inside are best read as back-modelling error — the
2002 diameter is a Sloboda estimate, not a measurement, so marginal stems fall
on either side of a boundary they are sitting on.

Either way the conclusion for FINN is unchanged: these are established trees,
not recruits. Figure `outputs/figures/fig1_wzp_stemmap.png`, panel C.

Taking `Pk == 0` at face value would overstate FINN's regeneration response by
roughly an order of magnitude. The pipeline therefore splits on the
back-modelled start-of-period diameter — which is precisely what Thünen
produced it for:

- back-modelled dbh **< 7 cm** → genuine recruit; appears only at t1, `"new"`
- back-modelled dbh **≥ 7 cm** → lateral ingrowth; present at **both** dates,
  using the back-modelled state at t0

## 5b. Two selection biases in the demographic cohorts

Both were found by noticing that observed mortality had the wrong *shape*, not
the wrong magnitude — small trees appeared to be the safest stems in the forest.
Both come from the same source: laterally ingrown stems are not ordinary trees,
they are trees selected on their own future.

### Immortal time bias in the mortality denominator

A stem is classed as lateral ingrowth **because it was tallied at BWI3** — so it
is conditioned on having survived and cannot, by construction, die. Leaving it
in the at-risk cohort is textbook immortal time bias.

It is **26.6 %** of the naive at-risk cohort, and concentrated in the small size
classes where lateral ingrowth lives:

| dbh class | % of cohort that is ingrowth | biased mortality | corrected | inflation |
|---|---|---|---|---|
| 7–15 | 42.9 % | 0.89 %/yr | **1.60 %/yr** | 1.81× |
| 15–25 | 29.9 % | 0.38 | 0.55 | 1.44× |
| 25–35 | 24.5 % | 0.20 | 0.27 | 1.33× |
| 35–50 | 20.6 % | 0.15 | 0.18 | 1.26× |
| 50+ | 18.1 % | 0.18 | 0.22 | 1.22× |

Uncorrected, mortality *rose* with tree size. Corrected, it recovers the
expected shape. Pooled over the calibration cohort the observed rate moves from
0.77 %/yr to **1.21 %/yr** — which reverses the sign of the apparent model
error, from 27 % over-prediction to 19 % under-prediction.

### Model-derived growth

The same stems carry a second problem. Their t0 diameter is a **Sloboda
back-model**, not a measurement, so their apparent increment is Thünen's growth
model rather than an observation. Matched on size (25–35 cm at t0):

| class | n | mean growth | sd |
|---|---|---|---|
| lateral ingrowth | 24,239 | 1.427 %/yr | 0.0053 |
| survivor | 72,303 | 1.263 %/yr | 0.0075 |

13 % higher and visibly smoother — the reduced variance is the giveaway, since
a model trajectory has none of the measurement scatter a real remeasurement
does. They were 21.8 % of all growth observations.

**Both are fixed in `wzp_effective_n()` / `scripts/04`:** the mortality and
growth cohorts admit only survivors and natural deaths — the stems genuinely
measured twice. Harvested stems remain censored as a competing risk.

## 6. Angle-count → fixed-area patches

FINN needs tree counts on a patch of known area. Each tallied tree is expanded
into `round(n_ha × A)` identical virtual trees on a nominal patch of
`A = PATCH_SIZE_HA` (default 0.1 ha), allocated by **largest remainder** within
(patch × species) so group totals survive the rounding.

A tree only risks vanishing when `n_ha × A < 0.5`, i.e. `g > 8A`; at A = 0.1 ha
that is dbh > 101 cm and affects 0.14 % of BWI3 stems, all recovered at group
level by the remainder step.

Measured basal-area conservation on real corners:

| A (ha) | median abs. error | 90th pct | mean error |
|---|---|---|---|
| 0.05 | 2.3 % | 11.2 % | −0.19 % |
| **0.10** | **1.1 %** | **5.0 %** | **−0.01 %** |
| 0.25 | 0.4 % | 2.0 % | +0.01 % |

The expansion is essentially **unbiased**, and the residual scatter is an order
of magnitude below the angle-count's own sampling error — a corner tallies a
median of 7 trees, i.e. roughly 38 % standard error on BA/ha.

### Sample sizes are deliberately *not* expanded

Expansion is right for densities and wrong for information content: one measured
tree must not tell a binomial likelihood it is forty. So rates are computed from
the expanded, density-weighted table, while `n_at_risk`, `n_died` and `growth_n`
are overwritten with counts of **tallied** trees. In practice this changes the
median `n_at_risk` from 43 to 3 — a 14× difference in how confidently the loss
treats each observation.

Harvested stems are **censored** out of the at-risk denominator rather than
counted as survivors.

## 7. Time

FINN advances one timestep per row of the time axis; there is no sub-stepping.
Left as `{BWI2, BWI3}` the fitted rates would be per decade.

The pipeline therefore lays env over years 0…P and puts the single observation
at year P, carrying a `period_length` column. FINN then averages its simulated
annual growth and mortality over the window `(i-P+1):i`, and **sums**
regeneration over it, before comparing.

Two consequences, both handled in `scripts/04_prepare_finn_inputs.R`:

1. `period_length` must be **constant across sites** — FINN reads it as
   `period = unique(...)` and its own comment says "currently we assume that
   they are constant over sites". Sites are filtered to one realised interval
   (default 10 years, the mode).
2. Because FINN compares against a *mean of annual rates*, the observed `growth`
   (relative increment over the whole period) and `mort` (proportion dying over
   the whole period) are **annualised**:
   `growth → (1+g)^(1/P) − 1`, `mort → 1 − (1−m)^(1/P)`.
   Regeneration is a period total and already matches the summed prediction, so
   it is left alone.

Do **not** pad the observation table with NA-only years: a year whose
observations are entirely NA produces a loss with no gradient, and `backward()`
then fails with *"element 0 of tensors does not require grad"*.

## 8. Coordinates and elevation

Exact tract positions are not released. What is public is the 1 km INSPIRE
(ETRS89-LAEA, EPSG:3035) cell containing the tract centre, so the pipeline uses
the cell centroid and the positional error is bounded at ~700 m. WorldClim 30s
cells are ~1 km here, so climate extraction is at worst displaced by one cell.

Elevation comes from the inventory's own per-tract `HoeheNN`, not from a DEM
lookup at the anonymised coordinate, because terrain varies far faster than 1 km.

**Independent validation of the whole coordinate chain:** inventory elevation vs
WorldClim elevation at the derived coordinates gives r = 0.991 with a median
absolute difference of 12 m over 18,719 corners. If the INSPIRE join were wrong,
this would not hold.

## 9. Species

`x_ba` carries `Gattung` and `Art`, so proper binomials can be built. 66 species
codes occur in the BWI2/BWI3 tree tables; the legacy 9xx BWI1 codes do not
appear. Aggregate classes ("Sonstige Laubbäume") become `<Genus> spp.`.

`x_ba` is a forestry code list, not a nomenclatural authority, and carries
misspelt epithets — `Fraxinus exelsior` among them. These are corrected in
`bwi_species_table()`; left alone they would split a species across two names in
any join against an external trait or distribution database.

Top 15 species by corner occurrence cover **95.3 %** of national basal area:
Picea abies, Fagus sylvatica, Pinus sylvestris, Quercus petraea, Quercus robur,
Betula pendula, Abies alba, Fraxinus excelsior, Larix decidua, Carpinus betulus,
Acer pseudoplatanus, Pseudotsuga menziesii, Alnus glutinosa, Prunus avium,
Larix kaempferi.

## 10. Sample sizes

| Filter | Corners |
|---|---|
| in the change database | 52,178 |
| usable (`Use23` = 100, measured both dates, classifiable) | 47,232 |
| + tier M0 (no tree harvested) | 23,378 |
| + basal-area window | 22,568 |
| + top-15 species ≥ 99 % of corner BA | 18,759 |
| + complete environment | 18,745 |
| + realised period exactly 10 years | ~9,900 |

Tracts with all four corners usable: 5,887 — but only 733 with all four at M0,
which is why `patch_mode = "ecke"` (site = corner) is the default. FINN
simulates `patches` stochastic replicates per site and compares their mean to
the observation, so one observed patch per site is not a limitation.

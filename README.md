# p057p034

Code and data for Küey et al. manuscript `#p057p034`

**Recruitment of clathrin to intracellular membranes is sufficient for vesicle formation**

Cansu Küey, Méghane Sittewelle, Gabrielle Larocque, Miguel Hernández-González & Stephen J. Royle

*bioRxiv* [doi: 10.1101/2022.03.21.485220](https://doi.org/10.1101/2022.03.21.485220)

## Data

- CD8Dynole - text files (Fig 5 - supp 1)
- Clathrin - superplot data (Fig 3)
- Dynamin - superplot data (Fig 5)
- EM - IMOD models and text files (Fig 4)
- EpsinFCHO - superplot data (Figs 6 and 7)
- Hooks - alternate hooks (Fig 3 - supp 1)
- MitoPitLocations - raw counts (Fig 2)
- MitoPitProfiles - 3 channel line profiles (Fig 1)
- Vps4Actin - superplot data (Fig 5 - supp 2)

## R

- `MitoGraph` R project to process outputs from [MitoGraph](https://github.com/vianamp/MitoGraph). Dataset from the paper and associated outputs are included.

## Scripts

- `ImageJ_Mitopits_FreeSpotsQuantification.ijm` ImageJ macro to quantify number of free spots for EpsinKD and FchoKO experiments.
- `MitoPitEM.ipf` Igor code to quantify diameter of MitoPits from EM images (in `Data/EM/`).
- `MitoPits.ipf` Igor code to analyse line profiles through MitoPits (in `Data/MitoPitLocations/`).
- `SpotSize.ipf` Igor code to calculate spot size for hot-wired CME/dynole experiments (in `Data/CD8Dynole/`).

Superplots are generated using [SuperPlot](https://github.com/quantixed/SuperPlot) in IgorPro.
